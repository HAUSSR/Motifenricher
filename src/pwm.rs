//! Position matrices: PPM construction, log-odds PWM, exact score
//! distribution via dynamic programming (used to calibrate site thresholds
//! to a target per-position significance), multithreaded site scanning,
//! IUPAC consensus and matrix similarity for motif dedup.

use crate::fasta::SeqRec;
use crate::kmer::{for_each_window_pos, revcomp_enc, BASE_CHARS};
use rayon::prelude::*;

/// Discretization of log-odds scores: 100 units = 1 bit.
pub const GRANULARITY: i64 = 100;

#[derive(Clone, Debug)]
pub struct Site {
    pub motif: usize,
    /// sequence index
    pub seq: usize,
    /// 0-based start of the window on the forward strand
    pub pos_fwd: usize,
    /// b'+' or b'-'
    pub strand: u8,
    /// discretized log-odds score
    pub units: i64,
    /// per-position P-value of this score under the background
    pub p: f64,
}

impl Site {
    /// 1-based closed interval of the plus-strand footprint of this site
    /// (the same window for both strands; `strand` carries the orientation).
    pub fn start_end(&self, _seq_len: usize, w: usize) -> (usize, usize) {
        (self.pos_fwd + 1, self.pos_fwd + w)
    }

    /// Site sequence (plus-strand orientation).
    pub fn seq_string(&self, enc: &[u8], w: usize) -> String {
        let win: Vec<u8> = enc[self.pos_fwd..self.pos_fwd + w].to_vec();
        let out = if self.strand == b'+' { win } else { revcomp_enc(&win) };
        out.iter().map(|&e| BASE_CHARS[e as usize] as char).collect()
    }
}

/// A position probability matrix plus its discretized log-odds form and the
/// exact background score distribution.
pub struct Pwm {
    pub width: usize,
    pub ppm: Vec<[f64; 4]>,
    pub bg: [f64; 4],
    /// units[col][base] = round(GRANULARITY * log2(ppm/bg))
    pub units: Vec<[i64; 4]>,
    /// score-unit of masses[0]
    pub offset: i64,
    /// masses[i] = P(score = offset + i) under the 0-order background
    pub masses: Vec<f64>,
    /// tail[i] = P(score >= offset + i)
    pub tail: Vec<f64>,
    pub threshold_units: i64,
    pub threshold_p: f64,
}

impl Pwm {
    pub fn build(ppm: Vec<[f64; 4]>, bg: [f64; 4]) -> Self {
        let width = ppm.len();
        let units: Vec<[i64; 4]> = ppm
            .iter()
            .map(|col| {
                let mut u = [0i64; 4];
                for b in 0..4 {
                    let bits = (col[b] / bg[b]).log2();
                    u[b] = (bits * GRANULARITY as f64).round() as i64;
                }
                u
            })
            .collect();
        // exact score distribution by DP over columns
        let mut lo: i64 = 0;
        let mut hi: i64 = 0;
        let mut masses: Vec<f64> = vec![1.0];
        for col in &units {
            let cmin = *col.iter().min().unwrap();
            let cmax = *col.iter().max().unwrap();
            let new_lo = lo + cmin;
            let new_len = (hi + cmax - new_lo + 1) as usize;
            let mut next = vec![0.0f64; new_len];
            for (i, &m) in masses.iter().enumerate() {
                if m <= 0.0 {
                    continue;
                }
                for b in 0..4 {
                    let abs_unit = lo + i as i64 + col[b];
                    next[(abs_unit - new_lo) as usize] += m * bg[b];
                }
            }
            lo = new_lo;
            hi += cmax;
            masses = next;
        }
        // renormalize to absorb rounding loss
        let total: f64 = masses.iter().sum();
        if total > 0.0 {
            for m in masses.iter_mut() {
                *m /= total;
            }
        }
        let mut tail = vec![0.0f64; masses.len()];
        let mut acc = 0.0;
        for i in (0..masses.len()).rev() {
            acc += masses[i];
            tail[i] = acc.min(1.0);
        }
        Pwm { width, ppm, bg, units, offset: lo, masses, tail, threshold_units: i64::MIN, threshold_p: 1.0 }
    }

    /// Smallest score whose background tail probability is <= `target`,
    /// plus the actual tail at that threshold.
    pub fn set_threshold(&mut self, target: f64) {
        let mut chosen = self.tail.len() - 1;
        for i in 0..self.tail.len() {
            if self.tail[i] <= target {
                chosen = i;
                break;
            }
        }
        self.threshold_units = self.offset + chosen as i64;
        self.threshold_p = self.tail[chosen];
    }

    #[inline]
    pub fn tail_at(&self, units: i64) -> f64 {
        let idx = units - self.offset;
        if idx < 0 {
            1.0
        } else if idx as usize >= self.tail.len() {
            0.0
        } else {
            self.tail[idx as usize]
        }
    }

    /// Score a single window (encoded bases); returns None on N.
    #[inline]
    pub fn score_window(&self, enc: &[u8], start: usize) -> Option<i64> {
        let mut s = 0i64;
        for j in 0..self.width {
            let e = enc[start + j];
            if e >= 4 {
                return None;
            }
            s += self.units[j][e as usize];
        }
        Some(s)
    }

    pub fn max_bits(&self) -> f64 {
        (0..self.width)
            .map(|j| {
                let u = self.units[j].iter().max().unwrap();
                *u as f64 / GRANULARITY as f64
            })
            .sum()
    }
}

/// Normalize raw counts into a PPM with a background-mixed pseudocount.
pub fn ppm_from_counts(counts: &[[u32; 4]], bg: [f64; 4], alpha: f64) -> Vec<[f64; 4]> {
    counts
        .iter()
        .map(|col| {
            let n: f64 = col.iter().sum::<u32>() as f64;
            let mut out = [0.0f64; 4];
            for b in 0..4 {
                out[b] = (col[b] as f64 + alpha * bg[b]) / (n + alpha);
            }
            out
        })
        .collect()
}

/// Rebuild column counts from observed sites in motif orientation.
pub fn counts_from_sites(sites: &[Site], seqs: &[SeqRec], width: usize) -> Vec<[u32; 4]> {
    let mut counts = vec![[0u32; 4]; width];
    for s in sites {
        let enc = &seqs[s.seq].enc;
        for j in 0..width {
            let mut e = enc[s.pos_fwd + j];
            if s.strand == b'-' {
                e = crate::kmer::comp_base(e);
            }
            if e < 4 {
                counts[j][e as usize] += 1;
            }
        }
    }
    counts
}

/// Collect all exact occurrences of the word family `words` (canonical
/// codes in both-strand mode) as sites. Each site is scored (log-odds in
/// motif orientation) against `pwm`; per-site P-values come from the exact
/// background score distribution. Overlapping hits on the same strand are
/// suppressed greedily by score.
pub fn scan_word_sites(
    seqs: &[SeqRec],
    words: &std::collections::HashSet<u64>,
    k: u32,
    both_strands: bool,
    pwm: &Pwm,
    motif: usize,
) -> Vec<Site> {
    let w = k as usize;
    // matrix that scores a plus-strand window against the motif's reverse
    // complement: rc_units[j][b] = units[w-1-j][comp(b)]
    let rc_units: Vec<[i64; 4]> = (0..w)
        .map(|j| {
            let o = &pwm.units[w - 1 - j];
            [o[3], o[2], o[1], o[0]]
        })
        .collect();
    let per_seq: Vec<Vec<Site>> = seqs
        .par_iter()
        .enumerate()
        .map(|(si, rec)| {
            let mut cands: Vec<(usize, u8, i64)> = Vec::new();
            for_each_window_pos(&rec.enc, k, both_strands, &mut |idx, start| {
                if !words.contains(&idx) {
                    return;
                }
                let win = &rec.enc[start..start + w];
                let forward = if both_strands {
                    let code: u64 = win.iter().fold(0u64, |c, &e| (c << 2) | e as u64);
                    code == idx
                } else {
                    true
                };
                // score the window in motif orientation
                let units = if forward {
                    pwm.score_window(&rec.enc, start)
                } else {
                    let mut s = 0i64;
                    for j in 0..w {
                        s += rc_units[j][rec.enc[start + j] as usize];
                    }
                    Some(s)
                };
                if let Some(u) = units {
                    let strand = if forward { b'+' } else { b'-' };
                    cands.push((start, strand, u));
                }
            });
            cands.sort_by(|a, b| b.2.cmp(&a.2).then(a.0.cmp(&b.0)));
            let mut sites = Vec::new();
            let mut taken_p: Vec<(usize, usize)> = Vec::new();
            let mut taken_m: Vec<(usize, usize)> = Vec::new();
            for (pos, strand, u) in cands {
                let taken = if strand == b'+' { &mut taken_p } else { &mut taken_m };
                if taken.iter().all(|&(s, e)| pos + w <= s || pos >= e) {
                    taken.push((pos, pos + w));
                    let p = pwm.tail_at(u).max(1e-300);
                    sites.push(Site { motif, seq: si, pos_fwd: pos, strand, units: u, p });
                }
            }
            sites
        })
        .collect();
    let mut all: Vec<Site> = per_seq.into_iter().flatten().collect();
    all.sort_by(|a, b| a.seq.cmp(&b.seq).then(a.pos_fwd.cmp(&b.pos_fwd)));
    all
}

/// IUPAC consensus string from a PPM.
pub fn consensus(ppm: &[[f64; 4]]) -> String {
    const IUPAC: [char; 16] = [
        '-', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N',
    ];
    let mut s = String::new();
    for col in ppm {
        let mut bits = 0u8;
        let maxf = col.iter().cloned().fold(0.0f64, f64::max);
        if maxf <= 0.0 {
            s.push('N');
            continue;
        }
        for b in 0..4 {
            if col[b] >= 0.7 * maxf && col[b] >= 0.1 {
                bits |= 1 << b;
            }
        }
        if bits == 0 {
            bits = 15;
        }
        s.push(IUPAC[bits as usize]);
    }
    s
}

/// Reverse-complement of a PPM (columns reversed, bases complemented).
pub fn ppm_revcomp(ppm: &[[f64; 4]]) -> Vec<[f64; 4]> {
    ppm.iter()
        .rev()
        .map(|col| [col[3], col[2], col[1], col[0]])
        .collect()
}

fn col_corr(a: &[f64; 4], b: &[f64; 4]) -> f64 {
    let ma: f64 = a.iter().sum::<f64>() / 4.0;
    let mb: f64 = b.iter().sum::<f64>() / 4.0;
    let mut num = 0.0;
    let mut da = 0.0;
    let mut db = 0.0;
    for i in 0..4 {
        let x = a[i] - ma;
        let y = b[i] - mb;
        num += x * y;
        da += x * x;
        db += y * y;
    }
    if da <= 0.0 || db <= 0.0 {
        return 0.0;
    }
    num / (da * db).sqrt()
}

/// Best average column Pearson correlation between two PPMs over any
/// alignment (both orientations). Used to drop redundant motifs.
pub fn ppm_similarity(a: &[[f64; 4]], b: &[[f64; 4]]) -> f64 {
    let (short, long) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    if short.len() < 3 || long.len() < 3 {
        return 0.0;
    }
    let min_overlap = ((0.7 * short.len() as f64).ceil() as usize).clamp(3, short.len());
    let mut best = 0.0f64;
    // slide the short matrix along the long one, in both orientations
    for s_var in [short.to_vec(), ppm_revcomp(short)] {
        for off in 0..long.len() {
            let ov = (long.len() - off).min(short.len());
            if ov < min_overlap {
                break;
            }
            let mut sum = 0.0;
            for i in 0..ov {
                sum += col_corr(&s_var[i], &long[off + i]);
            }
            best = best.max(sum / ov as f64);
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seqrec(s: &str) -> SeqRec {
        SeqRec::new("x".into(), s.as_bytes().to_vec())
    }

    #[test]
    fn score_dist_matches_direct_enumeration() {
        let ppm = vec![[0.9, 0.03, 0.03, 0.04], [0.1, 0.3, 0.3, 0.3]];
        let bg = [0.25; 4];
        let pwm = Pwm::build(ppm, bg);
        // enumerate all 16 dinucleotides directly
        let mut dist = std::collections::BTreeMap::<i64, f64>::new();
        for a in 0..4u8 {
            for b in 0..4u8 {
                let u = pwm.units[0][a as usize] + pwm.units[1][b as usize];
                *dist.entry(u).or_insert(0.0) += 0.25 * 0.25;
            }
        }
        for (u, m) in &dist {
            let got = pwm.masses[(*u - pwm.offset) as usize];
            assert!((got - m).abs() < 1e-12, "unit {} got {} want {}", u, got, m);
        }
        assert!((pwm.tail.iter().sum::<f64>() - 0.0).abs() < 1e-9 || pwm.tail[0] > 0.99);
    }

    #[test]
    fn threshold_and_tail_consistent() {
        let ppm = vec![
            [0.9, 0.03, 0.03, 0.04],
            [0.05, 0.05, 0.85, 0.05],
            [0.3, 0.3, 0.2, 0.2],
            [0.9, 0.03, 0.03, 0.04],
        ];
        let mut pwm = Pwm::build(ppm, [0.25; 4]);
        // this matrix is so informative that even its best score has
        // P > 1e-3 under the background: the strictest possible threshold
        // must be selected then
        pwm.set_threshold(1e-3);
        let idx = (pwm.threshold_units - pwm.offset) as usize;
        assert_eq!(pwm.tail[idx], pwm.threshold_p);
        assert!(pwm.tail[idx] > 1e-3 || idx == pwm.tail.len() - 1 || pwm.tail[idx + 1] > 1e-3);
        // a reachable target is honoured exactly
        pwm.set_threshold(0.5);
        assert!(pwm.threshold_p <= 0.5);
        let idx = (pwm.threshold_units - pwm.offset) as usize;
        if idx > 0 {
            assert!(pwm.tail[idx - 1] > 0.5);
        }
        // tails are monotonically non-increasing
        for w in pwm.tail.windows(2) {
            assert!(w[0] >= w[1]);
        }
    }

    #[test]
    fn scan_finds_planted_site() {
        // motif = ACGT..AT positions; build a strong PWM for "ACGTAG"
        let counts: Vec<[u32; 4]> = "ACGTAG"
            .bytes()
            .map(|b| {
                let mut c = [0u32; 4];
                c[crate::kmer::enc_base(b) as usize] = 100;
                c
            })
            .collect();
        let ppm = ppm_from_counts(&counts, [0.25; 4], 0.25);
        let mut pwm = Pwm::build(ppm, [0.25; 4]);
        pwm.set_threshold(1e-4);
        let mut word = 0u64;
        for &b in b"ACGTAG" {
            word = (word << 2) | crate::kmer::enc_base(b) as u64;
        }
        let words: std::collections::HashSet<u64> =
            [crate::kmer::canonical(word, 6)].into_iter().collect();
        let seq = seqrec("TTTTTTACGTAGTTTTTTACGTAGTTTT");
        let sites = scan_word_sites(std::slice::from_ref(&seq), &words, 6, true, &pwm, 0);
        assert_eq!(sites.len(), 2, "expected 2 sites, got {:?}", sites.len());
        for s in &sites {
            let ss = s.seq_string(&seq.enc, 6);
            assert_eq!(ss, "ACGTAG");
            assert_eq!(s.strand, b'+');
        }
        // minus-strand occurrence: coordinates are the plus-strand footprint
        let seq2 = seqrec("GGGGCTACGTGGGG");
        let sites2 = scan_word_sites(std::slice::from_ref(&seq2), &words, 6, true, &pwm, 0);
        assert_eq!(sites2.len(), 1);
        let (st, en) = sites2[0].start_end(seq2.len(), 6);
        assert_eq!(st, 5);
        assert_eq!(en, 10);
        assert_eq!(sites2[0].strand, b'-');
        assert_eq!(sites2[0].seq_string(&seq2.enc, 6), "ACGTAG");
    }

    #[test]
    fn consensus_iupac() {
        let ppm = vec![[0.9, 0.05, 0.03, 0.02], [0.25, 0.25, 0.25, 0.25], [0.05, 0.45, 0.45, 0.05]];
        assert_eq!(consensus(&ppm), "ANS");
    }

    #[test]
    fn similarity_detects_rc() {
        let a: Vec<[f64; 4]> = "ACGTAG"
            .bytes()
            .map(|b| {
                let mut c = [0.0f64; 4];
                c[crate::kmer::enc_base(b) as usize] = 0.95;
                c
            })
            .collect();
        let b: Vec<[f64; 4]> = "CTACGT"
            .bytes()
            .map(|b| {
                let mut c = [0.0f64; 4];
                c[crate::kmer::enc_base(b) as usize] = 0.95;
                c
            })
            .collect();
        assert!(ppm_similarity(&a, &b) > 0.9);
        let c: Vec<[f64; 4]> = "GGGGGG"
            .bytes()
            .map(|b| {
                let mut c = [0.0f64; 4];
                c[crate::kmer::enc_base(b) as usize] = 0.95;
                c
            })
            .collect();
        assert!(ppm_similarity(&a, &c) < 0.5);
    }
}

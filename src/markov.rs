//! Markov background model estimated from sequence data. Used to derive the
//! expected probability of each k-mer word under a null of "no enrichment".

use crate::fasta::SeqRec;
use crate::kmer::for_each_window;

/// Order-`r` Markov model with counts of all words of length 1..=r+1.
pub struct MarkovBg {
    pub order: u32,
    /// counts[l][code] = count of l-mers (index l = word length)
    pub counts: Vec<Vec<u64>>,
    pub totals: Vec<u64>,
    pub letter: [f64; 4],
}

impl MarkovBg {
    pub fn build(seqs: &[SeqRec], order: u32, both_strands: bool) -> Self {
        let mut counts: Vec<Vec<u64>> = (1..=order + 1)
            .map(|l| vec![0u64; 1usize << (2 * l)])
            .collect();
        let mut totals = vec![0u64; order as usize + 1];
        for rec in seqs {
            for l in 1..=order + 1 {
                for_each_window(&rec.enc, l, false, &mut |code| {
                    counts[l as usize - 1][code as usize] += 1;
                    totals[l as usize - 1] += 1;
                });
            }
            if both_strands {
                let rc = crate::kmer::revcomp_enc(&rec.enc);
                for l in 1..=order + 1 {
                    for_each_window(&rc, l, false, &mut |code| {
                        counts[l as usize - 1][code as usize] += 1;
                        totals[l as usize - 1] += 1;
                    });
                }
            }
        }
        let letter = {
            let mut f = [0.0f64; 4];
            if totals[0] > 0 {
                for b in 0..4 {
                    f[b] = counts[0][b] as f64 / totals[0] as f64;
                }
            } else {
                f = [0.25; 4];
            }
            f
        };
        MarkovBg { order, counts, totals, letter }
    }

    fn freq_l(&self, code: u64, l: usize) -> f64 {
        let t = self.totals[l - 1];
        if t == 0 {
            return 0.25f64.powi(l as i32);
        }
        self.counts[l - 1][code as usize] as f64 / t as f64
    }

    /// Probability of the k-mer `code` (length k) under the Markov model.
    pub fn prob(&self, code: u64, k: u32) -> f64 {
        let r = self.order as usize;
        let kk = k as usize;
        if kk == 0 {
            return 1.0;
        }
        let mut p = if kk <= r {
            self.freq_l(code, kk)
        } else {
            let prefix = code >> (2 * (kk - r));
            let mut pp = self.freq_l(prefix, r.max(1));
            for j in r..kk {
                // base index j (0-based from the left)
                let b_j = ((code >> (2 * (kk - 1 - j))) & 3) as usize;
                let ctx = (code >> (2 * (kk - j))) & ((1u64 << (2 * r)) - 1);
                if r == 0 {
                    pp *= self.letter[b_j];
                } else {
                    // counts[l-1] holds l-mer counts: ctx is an r-mer,
                    // (ctx<<2)|b_j is an (r+1)-mer
                    let cnt0 = self.counts[r - 1][ctx as usize];
                    let cnt1 = self.counts[r][((ctx << 2) | b_j as u64) as usize];
                    pp *= if cnt0 > 0 {
                        cnt1 as f64 / cnt0 as f64
                    } else {
                        self.letter[b_j]
                    };
                }
            }
            pp
        };
        // floor: never allow exactly zero so P-values stay defined
        if p <= 0.0 {
            let t = self.totals[self.order as usize].max(1) as f64;
            p = 1.0 / (t + 4.0f64.powi(self.order as i32 + 1));
        }
        p
    }

    /// Strand-summed probability of the canonical word (P(fwd) + P(rc)),
    /// i.e. the per-position probability of observing the canonical k-mer
    /// when scanning both strands.
    pub fn prob_canonical(&self, code: u64, k: u32) -> f64 {
        let p = self.prob(code, k) + self.prob(crate::kmer::revcomp_code(code, k), k);
        p.min(0.999_999)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mk(s: &str) -> SeqRec {
        SeqRec::new("x".into(), s.as_bytes().to_vec())
    }

    #[test]
    fn uniform_sequence() {
        // alternating ACGT => roughly uniform letter frequencies
        let seq = mk("ACGTACGTACGTACGTACGT");
        let bg = MarkovBg::build(std::slice::from_ref(&seq), 2, false);
        let mut a = 0u64;
        for &b in b"ACGT" {
            a = (a << 2) | crate::kmer::enc_base(b) as u64;
        }
        let p = bg.prob(a, 4);
        assert!(p > 0.0 && p <= 1.0);
    }

    #[test]
    fn at_rich_sequence_favors_at() {
        let seq = mk("AAAAAATTTTTTAAAAAATTTTTT");
        let bg = MarkovBg::build(std::slice::from_ref(&seq), 1, true);
        let mut word = 0u64;
        for &b in b"AAAA" {
            word = (word << 2) | crate::kmer::enc_base(b) as u64;
        }
        let mut gc = 0u64;
        for &b in b"GGCC" {
            gc = (gc << 2) | crate::kmer::enc_base(b) as u64;
        }
        assert!(bg.prob(word, 4) > bg.prob(gc, 4) * 10.0);
    }
}

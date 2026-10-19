//! 2-bit k-mer encoding, canonical (strand-independent) representation,
//! sliding-window iteration and multithreaded per-sequence presence counting.

use crate::fasta::SeqRec;
use rayon::prelude::*;

pub const BASE_CHARS: [u8; 4] = [b'A', b'C', b'G', b'T'];
const COMP_TBL: [u8; 4] = [3, 2, 1, 0];

#[inline]
pub fn enc_base(b: u8) -> u8 {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

#[inline]
pub fn comp_base(e: u8) -> u8 {
    if e < 4 {
        COMP_TBL[e as usize]
    } else {
        4
    }
}

pub fn revcomp_enc(enc: &[u8]) -> Vec<u8> {
    enc.iter().rev().map(|&e| comp_base(e)).collect()
}

/// Reverse complement of a packed k-mer code.
pub fn revcomp_code(code: u64, k: u32) -> u64 {
    let mut rc = 0u64;
    let mut c = code;
    for _ in 0..k {
        rc = (rc << 2) | COMP_TBL[(c & 3) as usize] as u64;
        c >>= 2;
    }
    rc
}

/// Canonical code = min(code, revcomp(code)): the same value for a word and
/// its reverse complement, so both strands collapse into one index.
#[inline]
pub fn canonical(code: u64, k: u32) -> u64 {
    code.min(revcomp_code(code, k))
}

pub fn decode_code(code: u64, k: u32) -> String {
    let mut s = String::with_capacity(k as usize);
    for i in (0..k).rev() {
        s.push(BASE_CHARS[((code >> (2 * i)) & 3) as usize] as char);
    }
    s
}

/// Iterate over all valid k-length windows of `enc`, invoking `f` with the
/// canonical (both-strand mode) or plain (plus-strand mode) k-mer code and
/// the 0-based start position of the window.
pub fn for_each_window_pos<F: FnMut(u64, usize)>(enc: &[u8], k: u32, both_strands: bool, f: &mut F) {
    let kk = k as usize;
    if kk == 0 || enc.len() < kk {
        return;
    }
    let mask: u64 = (1u64 << (2 * k)) - 1;
    let mut code: u64 = 0;
    let mut n_bad: u32 = 0;
    for i in 0..enc.len() {
        let e = enc[i];
        if e < 4 {
            code = ((code << 2) | e as u64) & mask;
        } else {
            code = (code << 2) & mask;
            n_bad += 1;
        }
        if i >= kk {
            if enc[i - kk] >= 4 {
                n_bad -= 1;
            }
        }
        if i + 1 >= kk && n_bad == 0 {
            let idx = if both_strands { canonical(code, k) } else { code };
            f(idx, i + 1 - kk);
        }
    }
}

/// Position-free convenience wrapper.
pub fn for_each_window<F: FnMut(u64)>(enc: &[u8], k: u32, both_strands: bool, f: &mut F) {
    for_each_window_pos(enc, k, both_strands, &mut |code, _pos| f(code));
}

/// For each k-mer: number of sequences that contain it at least once
/// (on either strand when `both_strands`). Parallel over sequences.
/// Returns `(canonical_code, count)` entries sorted by code. Dense arrays
/// are used for k <= 10 and a sparse hash map above, where 4^k would
/// exhaust memory.
pub fn count_presence(seqs: &[SeqRec], k: u32, both_strands: bool) -> Vec<(u64, u32)> {
    let chunk = (seqs.len() / (rayon::current_num_threads() * 4)).max(1);
    if k <= 10 {
        let n = 1usize << (2 * k); // 4^k
        let partials: Vec<Vec<u32>> = seqs
            .par_chunks(chunk)
            .map(|chunk| {
                let mut counts = vec![0u32; n];
                let mut seen = vec![u32::MAX; n];
                for (local_i, rec) in chunk.iter().enumerate() {
                    let li = local_i as u32;
                    for_each_window(&rec.enc, k, both_strands, &mut |idx| {
                        if seen[idx as usize] != li {
                            seen[idx as usize] = li;
                            counts[idx as usize] += 1;
                        }
                    });
                }
                counts
            })
            .collect();
        let mut total = vec![0u32; n];
        for p in partials {
            for (t, v) in total.iter_mut().zip(p) {
                *t += v;
            }
        }
        total
            .into_iter()
            .enumerate()
            .filter(|&(_, c)| c > 0)
            .map(|(i, c)| (i as u64, c))
            .collect()
    } else {
        use std::collections::HashMap;
        let partials: Vec<HashMap<u64, (u32, u32)>> = seqs
            .par_chunks(chunk)
            .map(|chunk| {
                // (count, last sequence that contained the word)
                let mut map: HashMap<u64, (u32, u32)> = HashMap::new();
                for (local_i, rec) in chunk.iter().enumerate() {
                    let li = local_i as u32;
                    for_each_window(&rec.enc, k, both_strands, &mut |idx| {
                        let e = map.entry(idx).or_insert((0, u32::MAX));
                        if e.1 != li {
                            e.1 = li;
                            e.0 += 1;
                        }
                    });
                }
                map
            })
            .collect();
        let mut total: HashMap<u64, u32> = HashMap::new();
        for p in partials {
            for (code, (c, _)) in p {
                *total.entry(code).or_insert(0) += c;
            }
        }
        let mut v: Vec<(u64, u32)> = total.into_iter().collect();
        v.sort_unstable_by_key(|&(code, _)| code);
        v
    }
}

/// Presence count of `code` in a sorted entry list from `count_presence`.
pub fn count_of(entries: &[(u64, u32)], code: u64) -> u32 {
    entries
        .binary_search_by_key(&code, |&(c, _)| c)
        .map(|i| entries[i].1)
        .unwrap_or(0)
}

/// Hamming distance between two same-length k-mer codes.
pub fn hamming(a: u64, b: u64, k: u32) -> u32 {
    let mut d = 0;
    for i in 0..k {
        if ((a >> (2 * i)) & 3) != ((b >> (2 * i)) & 3) {
            d += 1;
        }
    }
    d
}

/// Orientation-independent distance: also compares against the reverse
/// complement so palindromic/near-palindromic variants cluster together.
pub fn word_dist(a: u64, b: u64, k: u32) -> u32 {
    hamming(a, b, k).min(hamming(a, revcomp_code(b, k), k))
}

/// Count PPM column counts over every exact occurrence of `words` (given as
/// canonical codes in both-strand mode, or plain codes in plus-only mode).
/// Minus-strand occurrences are counted in motif orientation (i.e. the
/// reverse complement of the plus-strand window), so the resulting PPM is
/// orientation-normalized.
pub fn ppm_counts_from_words(
    seqs: &[SeqRec],
    words: &std::collections::HashSet<u64>,
    k: u32,
    both_strands: bool,
) -> Vec<[u32; 4]> {
    let w = k as usize;
    let mut counts = vec![[0u32; 4]; w];
    let partials: Vec<Vec<[u32; 4]>> = seqs
        .par_iter()
        .map(|rec| {
            let mut local = vec![[0u32; 4]; w];
            for_each_window_pos(&rec.enc, k, both_strands, &mut |idx, start| {
                if !words.contains(&idx) {
                    return;
                }
                let win = &rec.enc[start..start + w];
                // in both-strand mode idx is canonical: if the raw window
                // code equals it the window is already in motif orientation,
                // otherwise the motif occurs on the minus strand
                let forward = if both_strands {
                    let code: u64 = win.iter().fold(0u64, |c, &e| (c << 2) | e as u64);
                    code == idx
                } else {
                    true
                };
                if forward {
                    for j in 0..w {
                        local[j][win[j] as usize] += 1;
                    }
                } else {
                    for j in 0..w {
                        local[j][comp_base(win[w - 1 - j]) as usize] += 1;
                    }
                }
            });
            local
        })
        .collect();
    for p in partials {
        for (c, l) in counts.iter_mut().zip(p) {
            for b in 0..4 {
                c[b] += l[b];
            }
        }
    }
    counts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_roundtrip() {
        let s = b"ACGTTA";
        let mut code = 0u64;
        let k = 6;
        for &b in s {
            code = (code << 2) | enc_base(b) as u64;
        }
        assert_eq!(decode_code(code, k), "ACGTTA");
        assert_eq!(decode_code(revcomp_code(code, k), k), "TAACGT");
    }

    #[test]
    fn canonical_min() {
        // ACGT is its own reverse complement
        let mut code = 0u64;
        for &b in b"ACGT" {
            code = (code << 2) | enc_base(b) as u64;
        }
        assert_eq!(canonical(code, 4), code);
        // AAAT vs ATTT: canonical is min = AAAT (0)
        let mut a = 0u64;
        for &b in b"AAAT" {
            a = (a << 2) | enc_base(b) as u64;
        }
        assert_eq!(canonical(a, 4), a);
    }

    #[test]
    fn word_dist_orientation() {
        let enc = |s: &[u8]| {
            let mut c = 0u64;
            for &b in s {
                c = (c << 2) | enc_base(b) as u64;
            }
            c
        };
        let a = enc(b"AAAAAAA");
        let b = enc(b"TTTTTTT");
        assert_eq!(word_dist(a, b, 7), 0); // rc of AAAAAAA is TTTTTTT
        let c = enc(b"AAAAACC");
        assert_eq!(word_dist(a, c, 7), 2);
    }

    #[test]
    fn window_iteration_skips_n() {
        let seq = SeqRec::new("t".into(), b"ACGTNACGT".to_vec());
        let mut seen = Vec::new();
        for_each_window(&seq.enc, 4, true, &mut |c| seen.push(decode_code(c, 4)));
        // windows: ACGT (canonical ACGT), then N blocks CGTN.., then ACGT
        assert_eq!(seen, vec!["ACGT", "ACGT"]);
        let mut n_all = 0;
        for_each_window(&seq.enc, 3, false, &mut |_| n_all += 1);
        // ACG CGT [GTN] [TN.] .. ACG CGT => 4 valid
        assert_eq!(n_all, 4);
    }

    #[test]
    fn presence_counting() {
        let mk = |s: &str| SeqRec::new("x".into(), s.as_bytes().to_vec());
        let seqs = vec![mk("AAAAACGT"), mk("ACGTAAAA"), mk("TTTT")];
        let counts = count_presence(&seqs, 4, true);
        let mut a = 0u64;
        for &b in b"AAAA" {
            a = (a << 2) | enc_base(b) as u64;
        }
        // AAAA and TTTT are reverse complements -> same canonical word,
        // so all three sequences contain it
        assert_eq!(count_of(&counts, a), 3);
        // ACGT appears in seq0 and seq1 (presence is once per sequence even
        // though both strands collapse to the same canonical code)
        let mut c = 0u64;
        for &b in b"ACGT" {
            c = (c << 2) | enc_base(b) as u64;
        }
        assert_eq!(count_of(&counts, c), 2);
        // sparse path (k > 10) agrees with the dense path on presence
        let seqs2: Vec<SeqRec> = (0..20)
            .map(|i| {
                let s: String = std::iter::repeat_with(|| match i % 4 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    _ => 'T',
                })
                .take(11)
                .collect::<String>()
                    + "ACGCGGATCCA";
                mk(&s)
            })
            .collect();
        let sparse_k = count_presence(&seqs2, 11, true);
        let mut w = 0u64;
        for &b in b"ACGCGGATCCA" {
            w = (w << 2) | enc_base(b) as u64;
        }
        // "ACGCGGATCCA" occurs once in every sequence (canonical form)
        assert_eq!(count_of(&sparse_k, canonical(w, 11)), 20);
    }
}

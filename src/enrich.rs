//! Word-level enrichment statistics and greedy star-clustering of
//! significant words into motif seeds.

use crate::kmer::word_dist;

#[derive(Clone, Debug)]
pub struct WordStat {
    pub code: u64,
    pub k: u32,
    /// number of input sequences containing this word
    pub n_with: u32,
    /// per-sequence expected occurrence fraction under the null
    pub exp_frac: f64,
    /// exact binomial upper-tail P-value
    pub p: f64,
    /// Benjamini-Hochberg q-value across all tested words
    pub q: f64,
}

/// Greedy star clustering: repeatedly take the most significant unclustered
/// word as a seed and attach all remaining significant words of the same
/// length within `nmis` mismatches (orientation-aware). Returns seeds sorted
/// by seed P-value.
pub fn cluster_words(words: &[WordStat], nmis: u32, max_clusters: usize) -> Vec<SeedCluster> {
    let mut order: Vec<usize> = (0..words.len()).collect();
    order.sort_by(|&a, &b| words[a].p.partial_cmp(&words[b].p).unwrap_or(std::cmp::Ordering::Equal));
    let mut used = vec![false; words.len()];
    let mut out = Vec::new();
    for &i in &order {
        if used[i] {
            continue;
        }
        let seed = &words[i];
        used[i] = true;
        let mut members = vec![seed.code];
        for &j in &order {
            if used[j] || words[j].k != seed.k {
                continue;
            }
            if word_dist(seed.code, words[j].code, seed.k) <= nmis {
                used[j] = true;
                members.push(words[j].code);
            }
        }
        out.push(SeedCluster {
            k: seed.k,
            seed: seed.code,
            words: members,
            seed_p: seed.p,
            seed_q: seed.q,
        });
        if out.len() >= max_clusters {
            break;
        }
    }
    out
}

#[derive(Clone, Debug)]
pub struct SeedCluster {
    pub k: u32,
    pub seed: u64,
    pub words: Vec<u64>,
    pub seed_p: f64,
    pub seed_q: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ws(code: u64, k: u32, p: f64) -> WordStat {
        WordStat { code, k, n_with: 10, exp_frac: 0.01, p, q: p }
    }

    fn enc(s: &str) -> u64 {
        let mut c = 0u64;
        for b in s.bytes() {
            c = (c << 2) | crate::kmer::enc_base(b) as u64;
        }
        c
    }

    #[test]
    fn star_clustering_groups_neighbors() {
        let words = vec![
            ws(enc("AAAAAAA"), 7, 1e-9),
            ws(enc("AAAAAAC"), 7, 1e-8),
            ws(enc("AAAAACC"), 7, 1e-7),
            ws(enc("GGGGGGG"), 7, 1e-6),
        ];
        // nmis=1: seed AAAAAAA grabs AAAAAAAC (dist 1); AAAAACC (dist 2) starts
        // its own cluster; GGGGGGG another.
        let cl = cluster_words(&words, 1, 10);
        assert_eq!(cl.len(), 3);
        assert_eq!(cl[0].words.len(), 2);
        assert_eq!(cl[0].seed, enc("AAAAAAA"));
    }

    #[test]
    fn clusters_respect_max() {
        let words: Vec<WordStat> = (0..20)
            .map(|i| {
                let code = enc("AAAAAAA") ^ (i as u64 + 1) << 8;
                ws(code, 7, 1e-5 / (i as f64 + 1.0))
            })
            .collect();
        let cl = cluster_words(&words, 0, 5);
        assert_eq!(cl.len(), 5);
    }
}

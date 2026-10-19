//! Motifenricher: multithreaded de novo motif enrichment for promoter FASTA.
//!
//! Pipeline: enumerate k-mers on both strands -> exact binomial enrichment
//! test against a Markov (or empirical background-FASTA) null -> BH-FDR ->
//! greedy clustering of significant words -> PPM/PWM construction ->
//! threshold calibration by exact score-distribution DP -> iterative
//! rescan/refine -> motif-level binomial P-value -> minimal MEME + TSV output.

pub mod enrich;
pub mod fasta;
pub mod kmer;
pub mod markov;
pub mod meme;
pub mod pwm;
pub mod stats;

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use std::collections::HashSet;
use std::path::Path;
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct Config {
    pub input: String,
    pub out_dir: String,
    pub bg_fasta: Option<String>,
    /// 0..=4, only used when no background FASTA is given
    pub markov_order: u32,
    pub kmin: u32,
    pub kmax: u32,
    /// mismatches allowed when clustering significant words
    pub nmis: u32,
    pub max_motifs: usize,
    pub fdr: f64,
    pub min_word_p: f64,
    /// minimum number of sequences that must contain a word
    pub min_occ: u32,
    /// target per-position P-value used to report the PWM score threshold
    pub scan_p: f64,
    pub pseudocount: f64,
    pub dedup_corr: f64,
    pub strand_both: bool,
}

pub struct MotifReport {
    pub id: usize,
    pub consensus: String,
    pub width: usize,
    pub nsites: usize,
    pub n_seq_with: usize,
    /// per-sequence expected fraction at the final threshold
    pub exp_frac: f64,
    /// motif-level binomial P-value
    pub p_value: f64,
    /// BH q-value inherited from the seed word
    pub q_value: f64,
    /// seed-word P-value x number of tested words
    pub e_value: f64,
    pub log2_fe: f64,
    pub seed_word: String,
    pub seed_p: f64,
    pub threshold_bits: f64,
    pub ppm: Vec<[f64; 4]>,
    pub sites: Option<Vec<pwm::Site>>,
}

pub struct RunResult {
    pub motifs: Vec<MotifReport>,
    pub n_seqs: usize,
    pub mean_len: f64,
    pub bg_desc: String,
    pub bg_freq: [f64; 4],
    pub n_words_tested: u64,
    pub out_dir: String,
}

/// Fraction of the smaller site list that coincides (same sequence, start
/// within +-2 bp) with the larger one. Both lists must be sorted by
/// (seq, pos). High containment marks a motif as a shifted rediscovery of
/// an already-kept one.
fn site_containment(a: &[pwm::Site], b: &[pwm::Site]) -> f64 {
    let (mut i, mut j, mut matched) = (0usize, 0usize, 0usize);
    while i < a.len() && j < b.len() {
        let (sa, pa) = (a[i].seq, a[i].pos_fwd as i64);
        let (sb, pb) = (b[j].seq, b[j].pos_fwd as i64);
        if sa == sb && (pa - pb).abs() <= 2 {
            matched += 1;
            i += 1;
            j += 1;
        } else if (sa, pa) < (sb, pb) {
            i += 1;
        } else {
            j += 1;
        }
    }
    let smaller = a.len().min(b.len());
    if smaller == 0 {
        return 0.0;
    }
    matched as f64 / smaller as f64
}

pub fn run(cfg: &Config) -> Result<RunResult> {
    if cfg.kmin < 4 || cfg.kmax > 16 || cfg.kmin > cfg.kmax {
        bail!("kmin/kmax 参数无效 (要求 4 <= kmin <= kmax <= 16)");
    }
    if cfg.nmis > 3 {
        bail!("--nmis 最大为 3");
    }
    if cfg.markov_order > 4 {
        bail!("--markov-order 最大为 4");
    }
    if cfg.kmax > 12 {
        eprintln!(
            "提示: k > 12 时 word 长度较长, 出现次数少, 结果对 --min-occ / --min-frac 较敏感"
        );
    }
    let t0 = Instant::now();

    // ---- 1. read sequences -------------------------------------------------
    let seqs = fasta::read_fasta(&cfg.input)
        .with_context(|| format!("读取输入 FASTA 失败: {}", cfg.input))?;
    if seqs.len() < 5 {
        bail!("输入序列过少 ({} 条); 至少需要 5 条", seqs.len());
    }
    let n_seq = seqs.len();
    let mean_len = seqs.iter().map(|s| s.len()).sum::<usize>() as f64 / n_seq as f64;
    eprintln!(
        "[1/6] 读入 {} 条序列, 平均长度 {:.0} bp (总 {:} bp), 双链模式: {}",
        n_seq,
        mean_len,
        seqs.iter().map(|s| s.len()).sum::<usize>(),
        if cfg.strand_both { "是" } else { "否" }
    );

    // ---- 2. background -----------------------------------------------------
    let bg_seqs = match &cfg.bg_fasta {
        Some(p) => {
            let b = fasta::read_fasta(p).with_context(|| format!("读取背景 FASTA 失败: {p}"))?;
            if b.len() < 5 {
                bail!("背景 FASTA 序列过少 ({} 条)", b.len());
            }
            Some(b)
        }
        None => None,
    };
    let bg_freq = fasta::letter_freqs(bg_seqs.as_deref().unwrap_or(&seqs));
    let bg_desc = match &bg_seqs {
        Some(b) => format!("经验背景 ({} 条序列)", b.len()),
        None => format!("{} 阶马尔可夫自背景 (由输入序列估计)", cfg.markov_order),
    };
    eprintln!("[2/6] 背景模型: {}", bg_desc);

    let markov = if bg_seqs.is_none() {
        Some(markov::MarkovBg::build(&seqs, cfg.markov_order, cfg.strand_both))
    } else {
        None
    };

    // ---- 3. word-level enrichment ------------------------------------------
    eprintln!("[3/6] 枚举 {}..{}-mer 并计算富集 P 值 ...", cfg.kmin, cfg.kmax);
    let ks: Vec<u32> = (cfg.kmin..=cfg.kmax).collect();
    // average number of scannable positions per sequence, per k
    let mean_windows: Vec<f64> = ks
        .iter()
        .map(|&k| {
            let mut m = 0.0f64;
            for s in &seqs {
                m += (s.len() as f64 - k as f64 + 1.0).max(0.0);
            }
            m / n_seq as f64
        })
        .collect();
    let fg_counts: Vec<Vec<(u64, u32)>> = ks
        .iter()
        .map(|&k| kmer::count_presence(&seqs, k, cfg.strand_both))
        .collect();
    let bg_counts: Option<Vec<Vec<(u64, u32)>>> = match &bg_seqs {
        Some(b) => Some(ks.iter().map(|&k| kmer::count_presence(b, k, cfg.strand_both)).collect()),
        None => None,
    };
    let n_bg = bg_seqs.as_ref().map(|b| b.len()).unwrap_or(0) as i64;

    let mut all_words: Vec<enrich::WordStat> = Vec::new();
    for (ki, &k) in ks.iter().enumerate() {
        let counts = &fg_counts[ki];
        let words_k: Vec<enrich::WordStat> = counts
            .par_iter()
            .filter(|e| e.1 >= cfg.min_occ)
            .map(|&(code, c)| {
                let exp_frac = match (&bg_counts, &markov) {
                    (Some(bgc), _) => {
                        ((kmer::count_of(&bgc[ki], code) as f64) + 0.5) / (n_bg as f64 + 1.0)
                    }
                    (None, Some(mk)) => {
                        let p_word = mk.prob_canonical(code, k);
                        1.0 - (1.0 - p_word).powf(mean_windows[ki])
                    }
                    (None, None) => unreachable!(),
                };
                let exp_frac = exp_frac.clamp(1e-9, 1.0);
                let p = stats::binom_sf(c as i64, n_seq as i64, exp_frac);
                enrich::WordStat {
                    code,
                    k,
                    n_with: c,
                    exp_frac,
                    p: p.max(1e-300),
                    q: 1.0,
                }
            })
            .collect();
        all_words.extend(words_k);
    }
    let n_words_tested = all_words.len() as u64;
    {
        let ps: Vec<f64> = all_words.iter().map(|w| w.p).collect();
        let qs = stats::bh_fdr(&ps);
        for (w, q) in all_words.iter_mut().zip(qs) {
            w.q = q;
        }
        all_words.sort_by(|a, b| a.p.partial_cmp(&b.p).unwrap_or(std::cmp::Ordering::Equal));
    }

    let sig: Vec<enrich::WordStat> = all_words
        .iter()
        .filter(|w| w.q <= cfg.fdr && w.p <= cfg.min_word_p)
        .cloned()
        .collect();
    eprintln!(
        "      测试 word 数: {}, 显著 (q <= {}, p <= {}): {} 个",
        n_words_tested, cfg.fdr, cfg.min_word_p, sig.len()
    );

    let out_dir = Path::new(&cfg.out_dir).to_path_buf();
    std::fs::create_dir_all(&out_dir).with_context(|| format!("创建输出目录失败: {}", out_dir.display()))?;

    if sig.is_empty() {
        eprintln!("      未发现显著富集的 motif。可尝试: 提高 --fdr / --min-word-pvalue, 降低 --min-occ, 调整 k 范围, 或提供 --bg 背景文件");
        meme::write_meme(out_dir.join("motifs.meme"), &[], bg_freq, cfg.strand_both)?;
        meme::write_motifs_tsv(out_dir.join("motifs.tsv"), &[], n_seq)?;
        let _ = meme::write_sites_tsv(out_dir.join("sites.tsv"), &[], &seqs);
        return Ok(RunResult {
            motifs: Vec::new(),
            n_seqs: n_seq,
            mean_len,
            bg_desc,
            bg_freq,
            n_words_tested,
            out_dir: out_dir.display().to_string(),
        });
    }

    // ---- 4. cluster significant words --------------------------------------
    eprintln!("[4/6] 聚类显著 word (错配半径 {}) -> motif 种子", cfg.nmis);
    let clusters = enrich::cluster_words(&sig, cfg.nmis, cfg.max_motifs * 3 + 5);

    // ---- 5. PWM construction from anchored word-family sites ----------------
    eprintln!("[5/6] 构建 {} 个 motif 的 PWM 与位点表 ...", clusters.len());
    let motif_results: Vec<Option<MotifReport>> = clusters
        .par_iter()
        .map(|cl| {
            let ki = (cl.k - cfg.kmin) as usize;
            let words: HashSet<u64> = cl.words.iter().cloned().collect();
            // exact occurrences of the word family, orientation-normalized
            let counts = kmer::ppm_counts_from_words(&seqs, &words, cl.k, cfg.strand_both);
            let total_sites: u32 = counts.first().map(|c| c.iter().sum()).unwrap_or(0);
            if total_sites == 0 {
                return None;
            }
            let ppm = pwm::ppm_from_counts(&counts, bg_freq, cfg.pseudocount);
            let mut pwm_obj = pwm::Pwm::build(ppm.clone(), bg_freq);
            pwm_obj.set_threshold(cfg.scan_p);
            let sites = pwm::scan_word_sites(&seqs, &words, cl.k, cfg.strand_both, &pwm_obj, 0);
            if sites.is_empty() {
                return None;
            }
            let n_seq_with = sites.iter().map(|s| s.seq).collect::<HashSet<_>>().len();
            // expected per-sequence fraction of the whole word family: at a
            // given position the exact-match events are disjoint, so their
            // probabilities sum exactly (Markov mode); for the empirical
            // background the sum is a conservative upper bound of the union
            let exp_frac = match (&bg_counts, &markov) {
                (Some(bgc), _) => {
                    let mut q = 0.0f64;
                    for wcode in &cl.words {
                        q += ((kmer::count_of(&bgc[ki], *wcode) as f64) + 0.5)
                            / (n_bg as f64 + 1.0);
                    }
                    q.min(1.0)
                }
                (None, Some(mk)) => {
                    let mut p_fam = 0.0f64;
                    for wcode in &cl.words {
                        p_fam += mk.prob_canonical(*wcode, cl.k);
                    }
                    let p_fam = p_fam.min(0.999_999);
                    1.0 - (1.0 - p_fam).powf(mean_windows[ki])
                }
                (None, None) => unreachable!(),
            };
            let p_value = stats::binom_sf(n_seq_with as i64, n_seq as i64, exp_frac).max(1e-300);
            let log2_fe = ((n_seq_with as f64 / n_seq as f64) / exp_frac).log2();
            Some(MotifReport {
                id: 0,
                consensus: pwm::consensus(&ppm),
                width: pwm_obj.width,
                nsites: sites.len(),
                n_seq_with,
                exp_frac,
                p_value,
                q_value: cl.seed_q,
                e_value: (cl.seed_p * n_words_tested as f64).max(1e-300),
                log2_fe,
                seed_word: kmer::decode_code(cl.seed, cl.k),
                seed_p: cl.seed_p,
                threshold_bits: pwm_obj.threshold_units as f64 / pwm::GRANULARITY as f64,
                ppm,
                sites: Some(sites),
            })
        })
        .collect();

    // ---- 6. dedup, sort, cap and write --------------------------------------
    eprintln!("[6/6] 去重与输出 ...");
    let mut motifs: Vec<MotifReport> = motif_results.into_iter().flatten().collect();
    motifs.sort_by(|a, b| {
        a.q_value
            .partial_cmp(&b.q_value)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(
                a.p_value
                    .partial_cmp(&b.p_value)
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
    });
    let mut kept: Vec<MotifReport> = Vec::new();
    'cand: for m in motifs {
        for k in &kept {
            if pwm::ppm_similarity(&k.ppm, &m.ppm) >= cfg.dedup_corr {
                continue 'cand;
            }
            // shifted rediscoveries of the same sites (same motif found at
            // an offset) share most of their site positions
            let empty: Vec<pwm::Site> = Vec::new();
            let sa = k.sites.as_ref().unwrap_or(&empty);
            let sb = m.sites.as_ref().unwrap_or(&empty);
            if site_containment(sa, sb) >= 0.5 {
                continue 'cand;
            }
        }
        kept.push(m);
        if kept.len() >= cfg.max_motifs {
            break;
        }
    }
    for (i, m) in kept.iter_mut().enumerate() {
        m.id = i + 1;
        if let Some(s) = m.sites.as_mut() {
            for site in s.iter_mut() {
                site.motif = i + 1;
            }
        }
    }

    let meme_path = out_dir.join("motifs.meme");
    let tsv_path = out_dir.join("motifs.tsv");
    let sites_path = out_dir.join("sites.tsv");
    meme::write_meme(&meme_path, &kept, bg_freq, cfg.strand_both)
        .with_context(|| "写入 motifs.meme 失败")?;
    meme::write_motifs_tsv(&tsv_path, &kept, n_seq).with_context(|| "写入 motifs.tsv 失败")?;
    meme::write_sites_tsv(&sites_path, &kept, &seqs).with_context(|| "写入 sites.tsv 失败")?;

    eprintln!(
        "完成: {} 个 motif, 耗时 {:.1}s",
        kept.len(),
        t0.elapsed().as_secs_f64()
    );
    Ok(RunResult {
        motifs: kept,
        n_seqs: n_seq,
        mean_len,
        bg_desc,
        bg_freq,
        n_words_tested,
        out_dir: out_dir.display().to_string(),
    })
}

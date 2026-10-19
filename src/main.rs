use clap::Parser;
use motifenricher::{run, Config};

#[derive(Parser)]
#[command(
    name = "Motifenricher",
    version,
    about = "多线程启动子 motif de novo 富集分析工具 (word 枚举 + PWM 精化)",
    after_help = "示例:\n  Motifenricher -i promoters_2kb.fasta -o result_dir -t 16\n  Motifenricher -i target_promoters.fasta --bg all_promoters.fasta --kmin 6 --kmax 10"
)]
struct Cli {
    /// 输入启动子序列 FASTA
    #[arg(short, long)]
    input: String,

    /// 输出目录
    #[arg(short, long, default_value = "motifenricher_out")]
    out_dir: String,

    /// 线程数 (默认: 全部逻辑核心)
    #[arg(short, long)]
    threads: Option<usize>,

    /// 背景 FASTA (推荐: 全基因组启动子集等); 缺省使用输入序列的马尔可夫自背景
    #[arg(short, long)]
    bg: Option<String>,

    /// 马尔可夫背景阶数 (0-4, 仅在未提供 --bg 时生效)
    #[arg(long, default_value_t = 2)]
    markov_order: u32,

    /// 最小 word 长度
    #[arg(long, default_value_t = 5)]
    kmin: u32,

    /// 最大 word 长度 (<= 16; k > 10 时自动切换稀疏计数)
    #[arg(long, default_value_t = 9)]
    kmax: u32,

    /// 聚类显著 word 时允许的错配数 (0-3)
    #[arg(long, default_value_t = 1)]
    nmis: u32,

    /// 最多报告的 motif 数
    #[arg(long, default_value_t = 25)]
    max_motifs: usize,

    /// word 水平 FDR 阈值 (BH 校正)
    #[arg(long, default_value_t = 0.05)]
    fdr: f64,

    /// word 水平显著性 P 值上限
    #[arg(long, default_value_t = 1e-4)]
    min_word_pvalue: f64,

    /// word 必须出现在至少多少条序列中 (绝对数)
    #[arg(long, default_value_t = 3)]
    min_occ: u32,

    /// word 必须出现在至少这个比例的序列中
    #[arg(long, default_value_t = 0.02)]
    min_frac: f64,

    /// PWM 位点阈值的单位置假阳性率目标
    #[arg(long, default_value_t = 1e-5)]
    scan_p: f64,

    /// PPM 伪计数 (与背景碱基组成混合)
    #[arg(long, default_value_t = 0.25)]
    pseudocount: f64,

    /// motif 去重的矩阵相似度阈值 (Pearson)
    #[arg(long, default_value_t = 0.8)]
    dedup_corr: f64,

    /// 仅分析正链 (默认双链)
    #[arg(long, default_value_t = false)]
    plus_strand_only: bool,
}

fn fmt_p(p: f64) -> String {
    if p == 0.0 {
        "0".to_string()
    } else {
        format!("{:.2e}", p)
    }
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    let threads = cli.threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    });
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build_global()
        .map_err(|e| anyhow::anyhow!("初始化线程池失败: {e}"))?;

    let n_seq = count_fasta_records(&cli.input)?;
    let min_occ_eff = cli
        .min_occ
        .max((cli.min_frac * n_seq as f64).ceil() as u32);

    let cfg = Config {
        input: cli.input.clone(),
        out_dir: cli.out_dir.clone(),
        bg_fasta: cli.bg.clone(),
        markov_order: cli.markov_order,
        kmin: cli.kmin,
        kmax: cli.kmax,
        nmis: cli.nmis,
        max_motifs: cli.max_motifs,
        fdr: cli.fdr,
        min_word_p: cli.min_word_pvalue,
        min_occ: min_occ_eff,
        scan_p: cli.scan_p,
        pseudocount: cli.pseudocount,
        dedup_corr: cli.dedup_corr,
        strand_both: !cli.plus_strand_only,
    };

    eprintln!(
        "Motifenricher v{} | 线程: {} | 输入: {}",
        env!("CARGO_PKG_VERSION"),
        threads,
        cli.input
    );
    let result = run(&cfg)?;

    // ---- stdout summary -----------------------------------------------------
    println!();
    println!("背景模型: {}  |  测试 word 数: {}", result.bg_desc, result.n_words_tested);
    println!(
        "序列: {} 条, 平均长度 {:.0} bp",
        result.n_seqs, result.mean_len
    );
    if result.motifs.is_empty() {
        println!("未发现显著富集的 motif。");
    } else {
        println!(
            "共 {} 个富集 motif (输出目录: {}):",
            result.motifs.len(),
            result.out_dir
        );
        println!(
            "{:<4} {:<14} {:>3} {:>7} {:>7} {:>8} {:>8} {:>11} {:>11}",
            "#", "consensus", "w", "nsites", "seqs", "%seqs", "log2FE", "p_value", "q_value"
        );
        for m in &result.motifs {
            println!(
                "{:<4} {:<14} {:>3} {:>7} {:>7} {:>8.1} {:>8.2} {:>11} {:>11}",
                m.id,
                m.consensus,
                m.width,
                m.nsites,
                m.n_seq_with,
                100.0 * m.n_seq_with as f64 / result.n_seqs as f64,
                m.log2_fe,
                fmt_p(m.p_value),
                fmt_p(m.q_value),
            );
        }
        println!();
        println!("输出文件:");
        println!("  {}/motifs.meme  (minimal MEME 格式, 可用于 TOMTOM/FIMO/ceqlogo)", result.out_dir);
        println!("  {}/motifs.tsv   (含 P 值 / q 值 / E 值 / 富集倍数)", result.out_dir);
        println!("  {}/sites.tsv    (每个 motif 的序列位点坐标)", result.out_dir);
    }
    Ok(())
}

fn count_fasta_records(path: &str) -> anyhow::Result<usize> {
    use std::io::{BufRead, BufReader};
    let f = std::fs::File::open(path)
        .map_err(|e| anyhow::anyhow!("无法打开输入 FASTA {path}: {e}"))?;
    let mut n = 0usize;
    let mut has_seq = false;
    for line in BufReader::new(f).lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            if has_seq {
                n += 1;
            }
            has_seq = false;
        } else {
            has_seq = true;
        }
    }
    if has_seq {
        n += 1;
    }
    Ok(n.max(1))
}

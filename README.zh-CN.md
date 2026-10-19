# Motifenricher — 启动子 motif 富集分析工具（Rust）

[English README](README.md)

一个多线程的 de novo motif 富集分析命令行工具：输入一批启动子（或任何 genomic region）的
FASTA 序列，自动发现在这批序列中显著富集的 DNA motif 元件，输出 **minimal MEME 格式矩阵**、
**精确的二项检验 P 值**、**BH-FDR q 值**、**E 值**以及每个 motif 在每条序列上的**位点坐标表**。

纯 Rust 实现（仅依赖 `rayon`/`clap`/`anyhow`），无外部库要求，Linux 服务器上 `cargo build --release` 即可编译。

```
Motifenricher 0.2.0
用法: Motifenricher -i promoters_2kb.fasta -o result_dir -t 16
      Motifenricher -i target_promoters.fasta --bg all_gene_promoters.fasta --kmin 6 --kmax 10
```

---

## 编译

```bash
# 服务器上（推荐 rustup 安装 Rust >= 1.70）
cargo build --release
# 产物: target/release/Motifenricher
```

Windows 下若使用 GNU 工具链编译，本项目自带 `.cargo/config.toml`（`link-self-contained`），无需额外配置。

## 输入

- `-i/--input`：目标启动子 FASTA（大小写均可，非 ACGT 字符按 N 处理，含 N 的窗口自动跳过）。
- `-b/--bg`（**强烈推荐**）：背景 FASTA，例如全基因组所有基因的 2kb 启动子。
  提供时使用**经验背景**（每个 word 在背景序列中的出现比例）做检验；
  缺省时使用由输入序列自身估计的 **r 阶马尔可夫背景**（`--markov-order`，默认 2）。

## 输出（写入 `-o/--out-dir`，默认 `motifenricher_out`）

| 文件 | 内容 |
|---|---|
| `motifs.meme` | minimal MEME 格式（version 4），可直接提交 TOMTOM 比对 JASPAR/植物 TF 库，或用 `ceqlogo`/`FIMO` 出图、扫描 |
| `motifs.tsv` | 每个 motif 的统计表（见下） |
| `sites.tsv` | 每个 motif 在每条序列上的位点：`motif_id, seq_id, start, end, strand, score_bits, p_value, sequence`（坐标为正链 1-based 闭区间） |

`motifs.tsv` 列含义：

| 列 | 含义 |
|---|---|
| `p_value` | motif 家族水平的精确二项尾概率 P(≥观察序列数)，motif 级检验 |
| `q_value_fdr` | 发现（word 级）检验的 Benjamini–Hochberg FDR q 值，取自种子 word |
| `e_value` | 种子 word P 值 × 测试 word 总数（MEME 风格 E 值，也写入 meme 文件的 `E=` 字段） |
| `frac_with_site` / `expected_frac` | 观察/期望的含 motif 序列比例；`log2_fold_enrichment` 为二者之比的 log2 |
| `seed_word` / `seed_p` | 构成该 motif 的最显著 word 及其 P 值 |
| `threshold_bits` | PWM 扫描阈值（bit），对应单位置假阳性率 `--scan-p`（默认 1e-5，由精确打分分布 DP 求得） |

## 算法

1. **word 枚举**：对 k = `--kmin..--kmax`（默认 5–9，最大 16；k>10 用稀疏哈希计数），双链合并（正反向互补视为同一 word）统计每个 word 出现在多少条序列中（多线程）。
2. **富集检验**：对每个出现次数 ≥ `--min-occ` 的 word，用**精确二项检验**（正则化不完全 Beta 函数）计算 P(含有该 word 的序列数 ≥ 观察值)。零假设概率：经验背景模式用背景序列出现比例；马尔可夫模式用 word 的链合并概率与每序列窗口数。全部测试 word 做 **BH-FDR** 校正。
3. **聚类成 motif**：显著 word（q ≤ `--fdr` 且 P ≤ `--min-word-pvalue`）按显著性贪心星型聚类（同长、错配 ≤ `--nmis`、双向互补感知），每个簇形成一个 motif 种子。
4. **位点与矩阵**：收集簇内所有 word 的精确匹配位点（**链取向归一化**后）构建 PPM（伪计数 `--pseudocount` 与背景组成混合），并用 log-odds PWM 为每个位点打分；单位置 P 值来自 PWM 打分分布的**精确动态规划**。
5. **motif 级 P 值**：word 家族联合概率（簇内 word 概率之和，位置级互斥事件）的二项检验。
6. **去重**：矩阵 Pearson 相关（任意偏移、含反向互补取向，≥ `--dedup-corr`）或位点包含率（±2bp 容差 ≥ 0.5）判定冗余，按 q 值排序输出前 `--max-motifs` 个。

**统计口径说明**：`q_value_fdr` 控制的是 word 枚举发现阶段的 FDR（多重检验家族 = 所有测试 word）；`p_value` 是 motif 家族水平的富集强度检验。经验背景模式下家族期望概率取簇内 word 比例之和（并集的上界），P 值偏保守。

## 主要参数

| 参数 | 默认 | 说明 |
|---|---|---|
| `-t/--threads` | 全部核心 | 线程数 |
| `--kmin/--kmax` | 5 / 9 | word 长度范围（≤16；k>10 自动切换稀疏计数。word 越长出现越稀少，结果受 `--min-occ/--min-frac` 影响越大） |
| `--nmis` | 1 | 聚类错配半径（0–3） |
| `--fdr` | 0.05 | word 级 BH-FDR 阈值 |
| `--min-word-pvalue` | 1e-4 | word 显著性 P 上限 |
| `--min-occ/--min-frac` | 3 / 0.02 | word 最少出现序列数（取二者较大值） |
| `--max-motifs` | 25 | 最多报告 motif 数 |
| `--markov-order` | 2 | 自背景马尔可夫阶数（0–4） |
| `--scan-p` | 1e-5 | PWM 阈值的单位置假阳性率 |
| `--pseudocount` | 0.25 | PPM 伪计数 |
| `--plus-strand-only` | 关 | 仅分析正链（默认双链） |

## 性能参考

1000 条 × 2 kb（2 Mbp），8 线程，默认 k=5–9：**约 0.3 s**；k=4–14 约 2 s；k=10–12 约 1.6 s。
十万条序列规模（200 Mbp）预计数十秒至数分钟。

## 局限与建议

- 基于 word 枚举（oligo-analysis 风格），适合 4–12 bp 的连续 motif；**不找有间隔（gapped）的 motif**，可变区较长的退化 motif 灵敏度不如 MEME 的 EM。建议重要结果用 MEME Suite 复核，meme 文件可直接给 TOMTOM 做 TF 注释。
- 自背景（马尔可夫）能校正序列组成偏好，但如果目标集与全基因组背景存在整体组成差异（如 GC 含量），优先提供 `--bg`。
- 强组成偏好数据（如 AT-only 区段）中，富集结果可能反映二核苷酸组成而非 TF 结合偏好——可提高 `--markov-order` 至 3–4 缓解。

## 测试

```bash
cargo test        # 18 个单元测试（统计函数、编码、打分分布 DP、扫描坐标）+ 3 个端到端测试
```

端到端测试在合成启动子中植入已知 motif（8–13 bp，随机位置/双链），验证工具能以 p < 1e-10 完整找回；纯随机阴性对照数据集（4 万余 word）在 q ≤ 0.05 下零假阳性。

# Motifenricher

**Motifenricher** is a multithreaded, dependency-light command-line tool for *de novo* motif
enrichment analysis in promoter (or any genomic) sequences. Given a FASTA file of target
sequences, it discovers DNA motif elements that are significantly enriched in the set and
outputs **minimal MEME-format matrices** with exact **binomial P-values**, **BH-FDR q-values**
and **E-values**, plus per-sequence **site coordinates**.

Written in pure Rust (only `rayon`/`clap`/`anyhow`), it builds anywhere with
`cargo build --release` — no external libraries required. A 1,000 × 2 kb promoter set with the
default settings completes in well under a second on a modern CPU.

[中文文档 (Chinese README)](README.zh-CN.md)

## Installation

```bash
git clone https://github.com/HAUSSR/Motifenricher.git
cd Motifenricher
cargo build --release
# binary: target/release/Motifenricher
```

or install directly into `~/.cargo/bin`:

```bash
cargo install --path /path/to/Motifenricher
```

Requires Rust >= 1.70 (via [rustup](https://rustup.rs)).

## Quick start

```bash
# self-background: r-th order Markov model estimated from the input set
Motifenricher -i promoters_2kb.fasta -o result_dir -t 16

# recommended: empirical enrichment against a background promoter set
Motifenricher -i target_promoters.fasta --bg all_promoters_2kb.fasta -o result_dir -t 16
```

By default both strands are analyzed: a word and its reverse complement are counted as the
same element, and motif sites are orientation-normalized. Use `--plus-strand-only` to
restrict to the forward strand.

## Input / Output

**Input** — `-i/--input`: target FASTA (case-insensitive; non-ACGT characters are treated as
N and windows containing them are skipped). `-b/--bg`: background FASTA (recommended, e.g.
2 kb promoters of all genes genome-wide). Without `--bg`, an r-th order Markov background
(`--markov-order`, default 2) is estimated from the input sequences themselves.

**Output** (written to `-o/--out-dir`, default `motifenricher_out/`):

| File | Content |
|---|---|
| `motifs.meme` | minimal MEME format (version 4) — ready for TOMTOM against JASPAR/plant TF databases, or `ceqlogo`/FIMO for logos and scanning |
| `motifs.tsv` | per-motif statistics table (see below) |
| `sites.tsv` | per-sequence site table: `motif_id, seq_id, start, end, strand, score_bits, p_value, sequence` (1-based closed coordinates on the plus strand) |

Key columns of `motifs.tsv`:

| Column | Meaning |
|---|---|
| `p_value` | exact binomial upper-tail P-value at the motif (word-family) level |
| `q_value_fdr` | Benjamini–Hochberg FDR q-value of the discovery (word-level) test |
| `e_value` | seed-word P-value × number of tested words (MEME-style; also written to the `E=` field of the meme file) |
| `frac_with_site` / `expected_frac` | observed vs. expected fraction of sequences containing the motif |
| `log2_fold_enrichment` | log2 of the observed/expected ratio |
| `seed_word` / `seed_p` | the most significant word of the motif cluster and its P-value |
| `threshold_bits` | PWM scan threshold (bits) calibrated by exact score-distribution DP to `--scan-p` per-position error |

## Algorithm

1. **Word enumeration** — for k = `--kmin..--kmax` (default 5–9, up to 16), count in how many
   sequences each word occurs (both strands merged via canonical encoding). k > 10 switches
   to sparse hash counting.
2. **Enrichment test** — exact binomial tail probability (regularized incomplete beta) for
   every word occurring in ≥ `--min-occ` sequences; null probability from the empirical
   background fraction or the Markov model. All tested words are jointly corrected with
   Benjamini–Hochberg FDR.
3. **Clustering** — significant words (q ≤ `--fdr` and P ≤ `--min-word-pvalue`) are greedily
   star-clustered (same length, ≤ `--nmis` mismatches, reverse-complement aware) into motif
   seeds.
4. **Sites and matrices** — all exact occurrences of a cluster's words are collected with
   orientation normalization to build the PPM; sites are scored by a log-odds PWM whose
   per-position significance comes from an exact dynamic program over the score distribution.
5. **Motif-level P-value** — binomial test on the summed probability of the word family.
6. **Deduplication** — redundant motifs are removed by matrix Pearson correlation over any
   alignment (both orientations) or by site containment (±2 bp tolerance).

See the [Chinese README](README.zh-CN.md) for a detailed statistical interpretation of the
P/q/E values.

## Key options

| Option | Default | Description |
|---|---|---|
| `-t/--threads` | all cores | number of threads |
| `--kmin/--kmax` | 5 / 9 | motif word length range (≤ 16; k > 10 uses sparse counting) |
| `--nmis` | 1 | mismatch radius when clustering words (0–3) |
| `--fdr` | 0.05 | word-level BH-FDR threshold |
| `--min-word-pvalue` | 1e-4 | word significance threshold |
| `--min-occ` / `--min-frac` | 3 / 0.02 | minimum sequences containing a word (max of both) |
| `--max-motifs` | 25 | maximum number of motifs reported |
| `--markov-order` | 2 | Markov order of the self-background (0–4) |
| `--scan-p` | 1e-5 | per-position false-positive rate of the PWM threshold |
| `--plus-strand-only` | off | analyze the forward strand only |

## Limitations

Motifenricher is a word-enumeration method (oligo-analysis style), best suited for continuous
motifs of 4–12 bp. It does not discover gapped motifs, and sensitivity for highly degenerate
motifs is lower than EM-based MEME. For important results we recommend cross-checking with
the MEME suite; the emitted `motifs.meme` file can be fed to TOMTOM directly for TF
annotation. When the target set differs globally in composition from the background
(e.g. GC content), prefer providing `--bg`; for strongly composition-biased data raise
`--markov-order` to 3–4.

## Testing

```bash
cargo test   # 18 unit tests (statistics, encoding, score-distribution DP, site coordinates)
             # + 3 end-to-end tests
```

End-to-end tests plant known motifs (8–13 bp, random positions, both strands) into synthetic
promoters and verify full recovery at p < 1e-10; a purely random negative-control dataset
(> 40k words) yields zero false positives at q ≤ 0.05.

## License

MIT — see [LICENSE](LICENSE).

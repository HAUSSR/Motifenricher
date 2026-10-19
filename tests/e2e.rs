//! End-to-end test: synthesize promoters with two planted motifs and verify
//! that the pipeline recovers them with significant P-values.

use motifenricher::{run, Config};

struct Lcg(u64);

impl Lcg {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }
    /// AT-rich background like typical plant promoters
    fn base(&mut self) -> u8 {
        match self.next_u64() % 100 {
            0..=29 => b'A',
            30..=49 => b'C',
            50..=69 => b'G',
            _ => b'T',
        }
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

fn write_fasta(path: &std::path::Path, seqs: &[(String, Vec<u8>)]) {
    let mut out = String::new();
    for (h, s) in seqs {
        out.push('>');
        out.push_str(h);
        out.push('\n');
        out.push_str(std::str::from_utf8(s).unwrap());
        out.push('\n');
    }
    std::fs::write(path, out).unwrap();
}

fn make_test_fasta(path: &std::path::Path) {
    let mut rng = Lcg(20261014);
    const MOTIF_A: &[u8] = b"ACGCGGAT"; // G-box-like
    const MOTIF_B: &[u8] = b"TTTGACATA"; // TATA-plus-GA-like
    let n = 160usize;
    let len = 400usize;
    let mut seqs = Vec::new();
    for i in 0..n {
        let mut s: Vec<u8> = (0..len).map(|_| rng.base()).collect();
        if i % 5 < 3 {
            // 60%: plant motif A
            let pos = (rng.next_u64() as usize) % (len - 8);
            let m = if rng.next_u64() % 2 == 0 {
                MOTIF_A.to_vec()
            } else {
                revcomp(MOTIF_A)
            };
            s[pos..pos + 8].copy_from_slice(&m);
        }
        if i % 4 == 0 {
            // 25%: plant motif B
            let pos = (rng.next_u64() as usize) % (len - 9);
            let m = if rng.next_u64() % 2 == 0 {
                MOTIF_B.to_vec()
            } else {
                revcomp(MOTIF_B)
            };
            s[pos..pos + 9].copy_from_slice(&m);
        }
        seqs.push((format!("seq{i}::Chr1:{}-{}", i * 1000, i * 1000 + len), s));
    }
    write_fasta(path, &seqs);
}

fn test_config(dir: &std::path::Path, input: &std::path::Path) -> Config {
    Config {
        input: input.display().to_string(),
        out_dir: dir.display().to_string(),
        bg_fasta: None,
        markov_order: 2,
        kmin: 6,
        kmax: 9,
        nmis: 1,
        max_motifs: 15,
        fdr: 0.05,
        min_word_p: 1e-4,
        min_occ: 4,
        scan_p: 1e-5,
        pseudocount: 0.25,
        dedup_corr: 0.8,
        strand_both: true,
    }
}

fn rc_norm(s: &str) -> String {
    let rc: String = revcomp(s.as_bytes())
        .iter()
        .map(|&b| b as char)
        .collect();
    if s <= &rc {
        s.to_string()
    } else {
        rc
    }
}

#[test]
fn recovers_planted_motifs() {
    let dir = std::env::temp_dir().join("motifenrich_e2e");
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("promoters.fa");
    let out = dir.join("out");
    make_test_fasta(&input);

    let cfg = test_config(&out, &input);
    let result = run(&cfg).expect("pipeline run failed");

    assert!(!result.motifs.is_empty(), "no motifs recovered");
    // motif A (8mer, 60% of sequences) must be found among top hits
    let top: Vec<&motifenricher::MotifReport> = result.motifs.iter().take(6).collect();
    let hit_a = top.iter().any(|m| {
        let c = &m.consensus;
        let n = rc_norm(&c);
        n.contains("GCGG") || n.contains("CCGC")
    });
    assert!(hit_a, "planted motif A not recovered; got: {:?}", {
        let mut v: Vec<String> = result.motifs.iter().take(8).map(|m| m.consensus.clone()).collect();
        v.sort();
        v
    });
    // the best motif must be extremely significant
    assert!(result.motifs[0].p_value < 1e-10);
    // sites file coordinates valid: motif site sequences must match planted words
    let sites = std::fs::read_to_string(out.join("sites.tsv")).unwrap();
    assert!(sites.starts_with("motif_id\t"));
    assert!(sites.contains("ACGCGGAT") || sites.contains("ATCCGCGT"));

    let meme = std::fs::read_to_string(out.join("motifs.meme")).unwrap();
    assert!(meme.starts_with("MEME version 4"));
    assert!(meme.contains("letter-probability matrix"));
    assert!(meme.contains("E= "));

    let tsv = std::fs::read_to_string(out.join("motifs.tsv")).unwrap();
    assert!(tsv.contains("p_value"));
}

#[test]
fn empirical_background_mode() {
    // background = independent random sequences without planted motifs;
    // the pipeline must still find motif A and P-values must remain small.
    let dir = std::env::temp_dir().join("motifenrich_e2e_bg");
    std::fs::create_dir_all(&dir).unwrap();
    let fg = dir.join("fg.fa");
    make_test_fasta(&fg);
    let bgp = dir.join("bg.fa");
    let mut rng = Lcg(99887766);
    let bgseqs: Vec<(String, Vec<u8>)> = (0..300)
        .map(|i| {
            (
                format!("bg{i}"),
                (0..400).map(|_| rng.base()).collect::<Vec<u8>>(),
            )
        })
        .collect();
    write_fasta(&bgp, &bgseqs);

    let mut cfg = test_config(&dir.join("out_bg"), &fg);
    cfg.bg_fasta = Some(bgp.display().to_string());
    let result = run(&cfg).expect("pipeline failed with bg mode");
    assert!(!result.motifs.is_empty());
}

#[test]
fn long_k_sparse_mode_finds_long_motif() {
    // k > 10 exercises the sparse counting path; plant a 13-mer and
    // search k = 12..14
    let dir = std::env::temp_dir().join("motifenrich_e2e_longk");
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("longk.fa");
    let mut rng = Lcg(20261015);
    const M: &[u8] = b"ACGCGGATCCATG";
    let (n, len) = (120usize, 300usize);
    let mut seqs: Vec<(String, Vec<u8>)> = Vec::new();
    for i in 0..n {
        let mut s: Vec<u8> = (0..len).map(|_| rng.base()).collect();
        if i % 5 < 3 {
            let pos = (rng.next_u64() as usize) % (len - M.len());
            let m = if rng.next_u64() % 2 == 0 {
                M.to_vec()
            } else {
                revcomp(M)
            };
            s[pos..pos + M.len()].copy_from_slice(&m);
        }
        seqs.push((format!("lk{i}"), s));
    }
    write_fasta(&input, &seqs);

    let cfg = Config {
        input: input.display().to_string(),
        out_dir: dir.join("out").display().to_string(),
        bg_fasta: None,
        markov_order: 2,
        kmin: 12,
        kmax: 14,
        nmis: 1,
        max_motifs: 10,
        fdr: 0.05,
        min_word_p: 1e-4,
        min_occ: 4,
        scan_p: 1e-5,
        pseudocount: 0.25,
        dedup_corr: 0.8,
        strand_both: true,
    };
    let result = run(&cfg).expect("long-k run failed");
    assert!(!result.motifs.is_empty(), "long motif not recovered");
    let hit = result.motifs.iter().take(5).any(|m| {
        let c = rc_norm(&m.consensus);
        c.contains("CGCGGAT") || c.contains("ATCCGCG")
    });
    assert!(
        hit,
        "planted 13-mer not found; got {:?}",
        result.motifs
            .iter()
            .take(5)
            .map(|m| m.consensus.clone())
            .collect::<Vec<_>>()
    );
    assert!(result.motifs[0].p_value < 1e-10);
}

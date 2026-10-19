//! FASTA parsing. Sequences are upper-cased; every encoded byte is either
//! A/C/G/T or a lower-case/iupac code that gets mapped to `N` upstream.

use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Clone)]
pub struct SeqRec {
    /// Full defline without the leading '>'.
    pub header: String,
    /// Upper-case sequence bytes (may contain non-ACGT letters, treated as N).
    pub seq: Vec<u8>,
    /// 2-bit encoded bases: 0..=3 for ACGT, 4 for anything else.
    pub enc: Vec<u8>,
}

impl SeqRec {
    pub fn new(header: String, seq: Vec<u8>) -> Self {
        let seq: Vec<u8> = seq.iter().map(|&b| b.to_ascii_uppercase()).collect();
        let enc = seq.iter().map(|&b| super::kmer::enc_base(b)).collect();
        SeqRec { header, seq, enc }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }
}

pub fn read_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<SeqRec>> {
    let path = path.as_ref();
    let file = File::open(path).with_context(|| format!("无法打开 FASTA 文件: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut recs = Vec::new();
    let mut header: Option<String> = None;
    let mut cur: Vec<u8> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(h) = header.take() {
                if !cur.is_empty() {
                    recs.push(SeqRec::new(h, std::mem::take(&mut cur)));
                }
            }
            header = Some(rest.trim().to_string());
        } else {
            if header.is_none() {
                anyhow::bail!("FASTA 格式错误: 在首个 '>' 之前出现序列行 ({})", path.display());
            }
            cur.extend_from_slice(line.as_bytes());
        }
    }
    if let Some(h) = header.take() {
        if !cur.is_empty() {
            recs.push(SeqRec::new(h, cur));
        }
    }
    Ok(recs)
}

/// Mononucleotide frequencies over all sequences (ACGT counts).
pub fn letter_freqs(seqs: &[SeqRec]) -> [f64; 4] {
    let mut counts = [0u64; 4];
    for r in seqs {
        for &e in &r.enc {
            if (e as usize) < 4 {
                counts[e as usize] += 1;
            }
        }
    }
    let total: u64 = counts.iter().sum();
    if total == 0 {
        return [0.25; 4];
    }
    let mut f = [0.0f64; 4];
    for i in 0..4 {
        f[i] = counts[i] as f64 / total as f64;
    }
    f
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn parse_roundtrip() {
        let dir = std::env::temp_dir().join("motifenrich_test");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("smoke.fa");
        let mut f = std::fs::File::create(&p).unwrap();
        writeln!(f, ">seq1 some desc").unwrap();
        writeln!(f, "acgtNN").unwrap();
        writeln!(f, "TTGA").unwrap();
        writeln!(f, ">seq2").unwrap();
        writeln!(f, "ACGTACGT").unwrap();
        drop(f);
        let recs = read_fasta(&p).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].header, "seq1 some desc");
        assert_eq!(recs[0].seq, b"ACGTNNTTGA");
        assert_eq!(recs[0].enc, vec![0, 1, 2, 3, 4, 4, 3, 3, 2, 0]);
        assert_eq!(recs[1].len(), 8);
        let freqs = letter_freqs(&recs);
        assert!((freqs.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }
}

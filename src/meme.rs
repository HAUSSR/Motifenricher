//! Output writers: minimal MEME format, motif summary TSV and site table.

use crate::fasta::SeqRec;
use crate::pwm::Site;
use crate::MotifReport;
use anyhow::Result;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Minimal MEME format (compatible with meme-suite tools such as
/// TOMTOM/FIMO/ceqlogo).
pub fn write_meme<P: AsRef<Path>>(
    path: P,
    motifs: &[MotifReport],
    bg: [f64; 4],
    both_strands: bool,
) -> Result<()> {
    let mut out = String::new();
    out.push_str("MEME version 4\n\nALPHABET= ACGT\n\n");
    if both_strands {
        out.push_str("strands: + -\n\n");
    } else {
        out.push_str("strands: +\n\n");
    }
    let _ = writeln!(
        out,
        "Background letter frequencies\nA {:.4} C {:.4} G {:.4} T {:.4}\n",
        bg[0], bg[1], bg[2], bg[3]
    );
    for m in motifs {
        let _ = writeln!(out, "MOTIF motif_{} {}\n", m.id, m.consensus);
        let _ = writeln!(
            out,
            "letter-probability matrix: alength= 4 w= {} nsites= {} E= {:.3e}",
            m.width, m.nsites, m.e_value
        );
        for col in &m.ppm {
            let _ = writeln!(
                out,
                "{:.6} {:.6} {:.6} {:.6}",
                col[0], col[1], col[2], col[3]
            );
        }
        out.push('\n');
    }
    Ok(File::create(path.as_ref())?.write_all(out.as_bytes())?)
}

pub fn write_motifs_tsv<P: AsRef<Path>>(path: P, motifs: &[MotifReport], total_seqs: usize) -> Result<()> {
    let mut out = String::new();
    out.push_str("motif_id\tconsensus\twidth\tnsites\tseqs_with_site\ttotal_seqs\tfrac_with_site\texpected_frac\tlog2_fold_enrichment\tp_value\tq_value_fdr\te_value\tthreshold_bits\tseed_word\tseed_p\n");
    for m in motifs {
        let _ = writeln!(
            out,
            "motif_{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.6}\t{:.3}\t{:.4e}\t{:.4e}\t{:.4e}\t{:.2}\t{}\t{:.4e}",
            m.id,
            m.consensus,
            m.width,
            m.nsites,
            m.n_seq_with,
            total_seqs,
            m.n_seq_with as f64 / total_seqs.max(1) as f64,
            m.exp_frac,
            m.log2_fe,
            m.p_value,
            m.q_value,
            m.e_value,
            m.threshold_bits,
            m.seed_word,
            m.seed_p,
        );
    }
    Ok(File::create(path.as_ref())?.write_all(out.as_bytes())?)
}

pub fn write_sites_tsv<P: AsRef<Path>>(path: P, motifs: &[MotifReport], seqs: &[SeqRec]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path.as_ref())?);
    writeln!(w, "motif_id\tseq_id\tstart\tend\tstrand\tscore_bits\tp_value\tsequence")?;
    for m in motifs {
        let empty: Vec<Site> = Vec::new();
        let sites: &Vec<Site> = m.sites.as_ref().unwrap_or(&empty);
        for s in sites {
            let rec = &seqs[s.seq];
            let (st, en) = s.start_end(rec.len(), m.width);
            let score_bits = s.units as f64 / crate::pwm::GRANULARITY as f64;
            writeln!(
                w,
                "motif_{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.4e}\t{}",
                m.id,
                rec.header,
                st,
                en,
                s.strand as char,
                score_bits,
                s.p,
                s.seq_string(&rec.enc, m.width),
            )?;
        }
    }
    Ok(w.flush()?)
}

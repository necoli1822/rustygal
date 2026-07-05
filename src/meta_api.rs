// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2007-2016 Doug Hyatt, Univ. of Tennessee / UT-Battelle (original Prodigal C)
// Copyright (C) 2026 Sunju Kim (Rust reimplementation)

//! Metagenomic gene-finding exposed as a plain library API (no `python` feature
//! required). `run_meta` reproduces the binary's `-p meta` path exactly; its
//! `gff`/`trans_faa`/`nuc` bytes are byte-identical to the `rustygal` binary's
//! `-f gff` / `-a` / `-d` output for the same contig.

use crate::metagenomic::MetagenomicBin;
use crate::sequence::Mask;
use crate::{bitmap, dprog, gene, metagenomic, node, sequence};

/// Per-contig metagenomic gene-finding output. Each field is the exact byte
/// stream the binary would write for that contig in `-p meta` mode.
pub struct MetaOutput {
    /// `-a` protein FASTA bytes.
    pub trans_faa: Vec<u8>,
    /// `-f gff` bytes (includes the per-sequence `# Sequence Data` / `# Model
    /// Data` comment header).
    pub gff: Vec<u8>,
    /// `-d` nucleotide FASTA bytes.
    pub nuc: Vec<u8>,
}

/// Build the 50 metagenomic training bins once. The result is read-only and may
/// be shared (`&[MetagenomicBin]`) across threads for parallel `run_meta` calls.
pub fn meta_bins() -> Vec<MetagenomicBin> {
    let mut meta = vec![MetagenomicBin::default(); metagenomic::NUM_META];
    metagenomic::initialize_metagenomic_bins(&mut meta);
    meta
}

/// Run metagenomic gene-finding on a single contig.
///
/// * `seq_num`     1-based contig index (the binary assigns these in input
///                 order). Used verbatim for the `ID={seq_num}_{n}` gene tag and
///                 the gff `seqnum=` field, so pass the same value the binary
///                 would to get byte-identical output.
/// * `header`      the FASTA header text *after* `>` (no leading `>`, no newline).
///                 Its first whitespace-delimited token becomes the contig name
///                 used in FASTA/gff records; the full string is the gff `seqhdr`.
/// * `dna`         the contig sequence bytes (A/C/G/T/N, upper or lower case;
///                 must contain only sequence letters — no newlines).
/// * `meta`        the bins from [`meta_bins`] (shared, read-only).
///
/// Defaults match the binary's meta mode: closed=0, no masking, per-bin
/// translation table. Each call allocates its own scratch and only reads `meta`,
/// so it is safe to call concurrently from multiple threads.
pub fn run_meta(seq_num: i32, header: &str, dna: &[u8], meta: &[MetagenomicBin]) -> MetaOutput {
    let closed = 0i32;
    let output = 3i32; // gff
    let is_meta = 1i32;
    let slen = dna.len() as i32;
    let slen_u = slen.max(0) as usize;

    // Bitmaps. build_bitmaps is byte-identical to next_seq_multi's encoding
    // (g/t/c/a -> same 2-bit codes, N -> code 01 + useq mask bit, gc = gc/len)
    // for clean DNA, so this matches what the binary feeds the scorer.
    let bytes_2bit = slen_u / 4 + 16;
    let bytes_1bit = slen_u / 8 + 16;
    let mut seq = vec![0u8; bytes_2bit];
    let mut rseq = vec![0u8; bytes_2bit];
    let mut useq = vec![0u8; bytes_1bit];
    let mut gc = 0.0f64;
    bitmap::build_bitmaps(dna, slen, &mut seq, &mut rseq, &mut useq, &mut gc);

    let cur_header_str = header;
    let mut short_header = String::new();
    sequence::calc_short_header(cur_header_str, &mut short_header, seq_num);

    // Meta default: no masking.
    let mlist: Vec<Mask> = Vec::new();
    let nmask = 0i32;

    // Same slen-proportional scratch caps as process_single_sequence.
    let node_cap = (slen_u / 8 + 65536).min(slen_u * 2 + 16).max(16);
    let gene_cap = (slen_u / 20 + 1024).min(gene::MAX_GENES).max(16);
    let mut nodes = vec![node::Node::default(); node_cap];
    let mut genes = vec![gene::Gene::default(); gene_cap];

    let mut nn: i32 = 0;
    let mut ng: i32 = 0;
    let mut ipath: i32;
    let mut max_phase: i32 = 0;
    let mut max_score: f64 = -100.0;

    // ---- META gene-finding (numeric logic verbatim from the binary) ----------
    let mut low = 0.88495 * gc - 0.0102337;
    if low > 0.65 {
        low = 0.65;
    }
    let mut high = 0.86596 * gc + 0.1131991;
    if high < 0.35 {
        high = 0.35;
    }

    for i in 0..metagenomic::NUM_META {
        if i == 0 || meta[i].tinf.trans_table != meta[i - 1].tinf.trans_table {
            for n in nodes[..].iter_mut() {
                *n = node::Node::default();
            }
            nn = node::add_nodes(
                &seq, &rseq, slen, &mut nodes, closed, &mlist, nmask, &meta[i].tinf,
            );
            nodes[..nn as usize].sort_unstable_by(|a, b| node::compare_nodes(a, b));
        }

        if meta[i].tinf.gc < low || meta[i].tinf.gc > high {
            continue;
        }

        node::reset_node_scores(&mut nodes, nn);
        node::score_nodes(&seq, &rseq, slen, &mut nodes, nn, &meta[i].tinf, closed, is_meta);
        node::record_overlapping_starts(&mut nodes, nn, &meta[i].tinf, 1);
        ipath = dprog::dprog(&mut nodes, nn, &meta[i].tinf, 1);

        if ipath >= 0 && nodes[ipath as usize].score > max_score {
            max_phase = i as i32;
            max_score = nodes[ipath as usize].score;
            dprog::eliminate_bad_genes(&mut nodes, ipath, &meta[i].tinf);
            ng = gene::add_genes(&mut genes, &nodes, ipath);
            gene::tweak_final_starts(&mut genes, ng, &mut nodes, nn, &meta[i].tinf);
            gene::record_gene_data(&mut genes, ng, &nodes, &meta[i].tinf, seq_num);
        }
    }

    // Node recovery with the winning bin (REQUIRED: the loop's
    // eliminate_bad_genes mutated `nodes`; rebuild clean nodes for output).
    for n in nodes[..].iter_mut() {
        *n = node::Node::default();
    }
    nn = node::add_nodes(
        &seq, &rseq, slen, &mut nodes, closed, &mlist, nmask,
        &meta[max_phase as usize].tinf,
    );
    nodes[..nn as usize].sort_unstable_by(|a, b| node::compare_nodes(a, b));
    node::score_nodes(
        &seq, &rseq, slen, &mut nodes, nn, &meta[max_phase as usize].tinf, closed, is_meta,
    );

    // ---- outputs via the same gene:: writers the binary uses -----------------
    let desc_str = std::str::from_utf8(&meta[max_phase as usize].desc)
        .unwrap()
        .trim_end_matches('\0');

    let mut gff = Vec::new();
    gene::print_genes(
        &mut gff, &genes, ng, &nodes, slen, output, seq_num, 1, desc_str,
        &meta[max_phase as usize].tinf, cur_header_str, &short_header, crate::VERSION,
    );

    let mut trans_faa = Vec::new();
    gene::write_translations(
        &mut trans_faa, &genes, ng, &nodes, &seq, &rseq, &useq, slen,
        &meta[max_phase as usize].tinf, seq_num, &short_header,
    );

    let mut nuc = Vec::new();
    gene::write_nucleotide_seqs(
        &mut nuc, &genes, ng, &nodes, &seq, &rseq, &useq, slen,
        &meta[max_phase as usize].tinf, seq_num, &short_header,
    );

    MetaOutput { trans_faa, gff, nuc }
}

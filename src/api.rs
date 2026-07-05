// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2007-2016 Doug Hyatt, Univ. of Tennessee / UT-Battelle (original Prodigal C)
// Copyright (C) 2026 Sunju Kim (Rust reimplementation)

use crate::*;

pub struct GeneFindingResult {

    pub genes: Vec<gene::Gene>,

    pub num_genes: i32,

    pub nodes: Vec<node::Node>,

    pub num_nodes: i32,

    pub training_info: Option<training::Training>,
}

pub struct GeneFinderConfig {
    pub meta: bool,
    pub closed: i32,
    pub do_mask: i32,
    pub min_gene: i32,
    pub min_edge_gene: i32,
    pub max_overlap: i32,
}

impl Default for GeneFinderConfig {
    fn default() -> Self {
        GeneFinderConfig {
            meta: false,
            closed: 0,
            do_mask: 0,
            min_gene: 90,
            min_edge_gene: 60,
            max_overlap: 60,
        }
    }
}

pub fn find_genes(
    sequence: &[u8],
    training_info: Option<&training::Training>,
    config: &GeneFinderConfig,
) -> Result<GeneFindingResult, String> {

    let mut seq = vec![0u8; sequence::MAX_SEQ / 4];
    let mut rseq = vec![0u8; sequence::MAX_SEQ / 4];
    let mut useq = vec![0u8; sequence::MAX_SEQ / 8];

    let slen = sequence.len() as i32;
    if slen > sequence::MAX_SEQ as i32 {
        return Err(format!("Sequence too long: {} > {}", slen, sequence::MAX_SEQ));
    }

    let mut gc = 0.0;
    bitmap::build_bitmaps(sequence, slen, &mut seq, &mut rseq, &mut useq, &mut gc);

    let max_slen = if slen > node::STT_NOD as i32 * 8 {
        (slen / 8) as usize
    } else {
        node::STT_NOD
    };

    let mut nodes = vec![node::Node::default(); max_slen];
    let mut genes = vec![gene::Gene::default(); gene::MAX_GENES];

    let mlist = vec![sequence::Mask { begin: 0, end: 0 }; sequence::MAX_MASKS];
    let nmask = 0;

    let tinf = training_info.ok_or("Training info required for single mode")?;

    let nn = node::add_nodes(&seq, &rseq, slen, &mut nodes, config.closed, &mlist, nmask, tinf);
    nodes[..nn as usize].sort_unstable_by(|a, b| node::compare_nodes(a, b));

    node::score_nodes(&seq, &rseq, slen, &mut nodes, nn, tinf, config.closed, if config.meta { 1 } else { 0 });

    node::record_overlapping_starts(&mut nodes, nn, tinf, 1);
    let ipath = dprog::dprog(&mut nodes, nn, tinf, 1);
    dprog::eliminate_bad_genes(&mut nodes, ipath, tinf);

    let ng = gene::add_genes(&mut genes, &nodes, ipath);

    gene::tweak_final_starts(&mut genes, ng, &mut nodes, nn, tinf);
    gene::record_gene_data(&mut genes, ng, &nodes, tinf, 1);

    Ok(GeneFindingResult {
        genes: genes[..ng as usize].to_vec(),
        num_genes: ng,
        nodes: nodes[..nn as usize].to_vec(),
        num_nodes: nn,
        training_info: Some(tinf.clone()),
    })
}

pub fn train_on_sequence(
    sequence: &[u8],
    translation_table: i32,
    force_nonsd: bool,
) -> Result<training::Training, String> {
    let slen = sequence.len() as i32;

    if slen < 20000 {
        return Err(format!(
            "Sequence too short for training: {} < 20000 bp. Use metagenomic mode instead.",
            slen
        ));
    }

    let mut seq = vec![0u8; sequence::MAX_SEQ / 4];
    let mut rseq = vec![0u8; sequence::MAX_SEQ / 4];
    let mut useq = vec![0u8; sequence::MAX_SEQ / 8];
    let mut gc = 0.0;

    bitmap::build_bitmaps(sequence, slen, &mut seq, &mut rseq, &mut useq, &mut gc);

    let mut tinf = training::Training::default();
    tinf.trans_table = translation_table;
    tinf.gc = gc;

    let max_slen = if slen > node::STT_NOD as i32 * 8 {
        (slen / 8) as usize
    } else {
        node::STT_NOD
    };
    let mut nodes = vec![node::Node::default(); max_slen];

    let mlist = vec![sequence::Mask { begin: 0, end: 0 }; sequence::MAX_MASKS];
    let nmask = 0;

    let nn = node::add_nodes(&seq, &rseq, slen, &mut nodes, 0, &mlist, nmask, &tinf);
    nodes[..nn as usize].sort_unstable_by(|a, b| node::compare_nodes(a, b));

    let gc_frame = sequence::calc_most_gc_frame(&seq, slen);
    node::record_gc_bias(&gc_frame, &mut nodes, nn, &mut tinf);

    node::record_overlapping_starts(&mut nodes, nn, &tinf, 0);
    let ipath = dprog::dprog(&mut nodes, nn, &tinf, 0);

    node::calc_dicodon_gene(&mut tinf, &seq, &rseq, slen, &nodes, ipath);
    node::raw_coding_score(&seq, &rseq, slen, &mut nodes, nn, &tinf);

    node::rbs_score(&seq, &rseq, slen, &mut nodes, nn, &tinf);
    node::train_starts_sd(&seq, &rseq, slen, &mut nodes, nn, &mut tinf);
    node::determine_sd_usage(&mut tinf);

    if force_nonsd {
        tinf.uses_sd = 0;
    }

    if tinf.uses_sd == 0 {
        node::train_starts_nonsd(&seq, &rseq, slen, &mut nodes, nn, &mut tinf);
    }

    let mut genes = vec![gene::Gene::default(); gene::MAX_GENES];
    let ng = gene::add_genes(&mut genes, &nodes, ipath);

    gene::tweak_final_starts(&mut genes, ng, &mut nodes, nn, &tinf);

    Ok(tinf)
}

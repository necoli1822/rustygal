// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2007-2016 Doug Hyatt, Univ. of Tennessee / UT-Battelle (original Prodigal C)
// Copyright (C) 2026 Sunju Kim (Rust reimplementation)

use crate::node::Node;
use crate::training::Training;

pub const MAX_GENES: usize = 30_000;

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Gene {

    pub begin: i32,

    pub end: i32,

    pub start_ndx: i32,

    pub stop_ndx: i32,

    pub gene_data: [u8; 500],

    pub score_data: [u8; 500],
}

impl Default for Gene {
    fn default() -> Self {
        Self {
            begin: 0,
            end: 0,
            start_ndx: 0,
            stop_ndx: 0,
            gene_data: [0; 500],
            score_data: [0; 500],
        }
    }
}

pub fn add_genes(genes: &mut [Gene], nodes: &[Node], dbeg: i32) -> i32 {
    use crate::sequence::STOP;

    if dbeg == -1 {
        return 0;
    }

    let mut path = dbeg;
    while nodes[path as usize].traceb != -1 {
        path = nodes[path as usize].traceb;
    }

    let mut ctr = 0i32;
    while path != -1 {

        if nodes[path as usize].elim == 1 {
            path = nodes[path as usize].tracef;
            continue;
        }

        if nodes[path as usize].strand == 1 && nodes[path as usize].type_ != STOP {
            genes[ctr as usize].begin = nodes[path as usize].ndx + 1;
            genes[ctr as usize].start_ndx = path;
        }

        if nodes[path as usize].strand == -1 && nodes[path as usize].type_ == STOP {
            genes[ctr as usize].begin = nodes[path as usize].ndx - 1;
            genes[ctr as usize].stop_ndx = path;
        }

        if nodes[path as usize].strand == 1 && nodes[path as usize].type_ == STOP {
            genes[ctr as usize].end = nodes[path as usize].ndx + 3;
            genes[ctr as usize].stop_ndx = path;
            ctr += 1;
        }

        if nodes[path as usize].strand == -1 && nodes[path as usize].type_ != STOP {
            genes[ctr as usize].end = nodes[path as usize].ndx + 1;
            genes[ctr as usize].start_ndx = path;
            ctr += 1;
        }

        path = nodes[path as usize].tracef;

        if ctr >= genes.len() as i32 {
            eprintln!("warning, max # of genes exceeded, truncating...");
            return ctr;
        }
    }

    ctr
}

pub fn record_gene_data(
    genes: &mut [Gene],
    ng: i32,
    nodes: &[Node],
    tinf: &Training,
    sctr: i32,
) {
    use crate::sequence::mer_text;

    let type_string = ["ATG", "GTG", "TTG", "Edge"];

    let sd_string = [
        "None",
        "GGA/GAG/AGG",
        "3Base/5BMM",
        "4Base/6BMM",
        "AGxAG",
        "AGxAG",
        "GGA/GAG/AGG",
        "GGxGG",
        "GGxGG",
        "AGxAG",
        "AGGAG(G)/GGAGG",
        "AGGA/GGAG/GAGG",
        "AGGA/GGAG/GAGG",
        "GGA/GAG/AGG",
        "GGxGG",
        "AGGA",
        "GGAG/GAGG",
        "AGxAGG/AGGxGG",
        "AGxAGG/AGGxGG",
        "AGxAGG/AGGxGG",
        "AGGAG/GGAGG",
        "AGGAG",
        "AGGAG",
        "GGAGG",
        "GGAGG",
        "AGGAGG",
        "AGGAGG",
        "AGGAGG",
    ];

    let sd_spacer = [
        "None",
        "3-4bp",
        "13-15bp",
        "13-15bp",
        "11-12bp",
        "3-4bp",
        "11-12bp",
        "11-12bp",
        "3-4bp",
        "5-10bp",
        "13-15bp",
        "3-4bp",
        "11-12bp",
        "5-10bp",
        "5-10bp",
        "5-10bp",
        "5-10bp",
        "11-12bp",
        "3-4bp",
        "5-10bp",
        "11-12bp",
        "3-4bp",
        "5-10bp",
        "3-4bp",
        "5-10bp",
        "11-12bp",
        "3-4bp",
        "5-10bp",
    ];

    for i in 0..ng as usize {

        genes[i].gene_data = [0; 500];
        genes[i].score_data = [0; 500];

        let ndx = genes[i].start_ndx as usize;
        let sndx = genes[i].stop_ndx as usize;

        let partial_left = if (nodes[ndx].edge == 1 && nodes[ndx].strand == 1)
            || (nodes[sndx].edge == 1 && nodes[ndx].strand == -1)
        {
            1
        } else {
            0
        };

        let partial_right = if (nodes[sndx].edge == 1 && nodes[ndx].strand == 1)
            || (nodes[ndx].edge == 1 && nodes[ndx].strand == -1)
        {
            1
        } else {
            0
        };

        let st_type = if nodes[ndx].edge == 1 {
            3
        } else {
            nodes[ndx].type_ as usize
        };

        let mut gene_data = format!(
            "ID={}_{};partial={}{};start_type={};",
            sctr,
            i + 1,
            partial_left,
            partial_right,
            type_string[st_type]
        );

        let rbs1 = tinf.rbs_wt[nodes[ndx].rbs[0] as usize] * tinf.st_wt;
        let rbs2 = tinf.rbs_wt[nodes[ndx].rbs[1] as usize] * tinf.st_wt;

        if tinf.uses_sd == 1 {

            if rbs1 > rbs2 {
                gene_data.push_str(&format!(
                    "rbs_motif={};rbs_spacer={}",
                    sd_string[nodes[ndx].rbs[0] as usize], sd_spacer[nodes[ndx].rbs[0] as usize]
                ));
            } else {
                gene_data.push_str(&format!(
                    "rbs_motif={};rbs_spacer={}",
                    sd_string[nodes[ndx].rbs[1] as usize], sd_spacer[nodes[ndx].rbs[1] as usize]
                ));
            }
        } else {

            let mut qt = vec![0u8; 10];
            mer_text(&mut qt, nodes[ndx].mot.len, nodes[ndx].mot.ndx);

            if tinf.no_mot > -0.5 && rbs1 > rbs2 && rbs1 > nodes[ndx].mot.score * tinf.st_wt {
                gene_data.push_str(&format!(
                    "rbs_motif={};rbs_spacer={}",
                    sd_string[nodes[ndx].rbs[0] as usize], sd_spacer[nodes[ndx].rbs[0] as usize]
                ));
            } else if tinf.no_mot > -0.5 && rbs2 >= rbs1 && rbs2 > nodes[ndx].mot.score * tinf.st_wt
            {
                gene_data.push_str(&format!(
                    "rbs_motif={};rbs_spacer={}",
                    sd_string[nodes[ndx].rbs[1] as usize], sd_spacer[nodes[ndx].rbs[1] as usize]
                ));
            } else if nodes[ndx].mot.len == 0 {
                gene_data.push_str("rbs_motif=None;rbs_spacer=None");
            } else {
                let motif_str =
                    std::str::from_utf8(&qt[..nodes[ndx].mot.len as usize]).unwrap_or("None");
                gene_data.push_str(&format!(
                    "rbs_motif={};rbs_spacer={}bp",
                    motif_str, nodes[ndx].mot.spacer
                ));
            }
        }

        gene_data.push_str(&format!(";gc_cont={:.3}", nodes[ndx].gc_cont));

        let gene_data_bytes = gene_data.as_bytes();
        let copy_len = gene_data_bytes.len().min(499);
        genes[i].gene_data[..copy_len].copy_from_slice(&gene_data_bytes[..copy_len]);
        genes[i].gene_data[copy_len] = 0;

        let confidence = calculate_confidence(nodes[ndx].cscore + nodes[ndx].sscore, tinf.st_wt);

        let score_data = format!(
            "conf={:.2};score={:.2};cscore={:.2};sscore={:.2};rscore={:.2};uscore={:.2};tscore={:.2};",
            confidence,
            nodes[ndx].cscore + nodes[ndx].sscore,
            nodes[ndx].cscore,
            nodes[ndx].sscore,
            nodes[ndx].rscore,
            nodes[ndx].uscore,
            nodes[ndx].tscore
        );

        let score_data_bytes = score_data.as_bytes();
        let copy_len = score_data_bytes.len().min(499);
        genes[i].score_data[..copy_len].copy_from_slice(&score_data_bytes[..copy_len]);
        genes[i].score_data[copy_len] = 0;
    }
}

pub fn tweak_final_starts(
    genes: &mut [Gene],
    ng: i32,
    nodes: &mut [Node],
    nn: i32,
    tinf: &Training,
) {
    use crate::node::{intergenic_mod, MAX_SAM_OVLP};
    use crate::sequence::STOP;

    for i in 0..ng as usize {
        let ndx = genes[i].start_ndx as usize;
        let sc = nodes[ndx].sscore + nodes[ndx].cscore;
        let mut igm = 0.0f64;

        if i > 0 && nodes[ndx].strand == 1 && nodes[genes[i - 1].start_ndx as usize].strand == 1 {
            igm = intergenic_mod(
                &nodes[genes[i - 1].stop_ndx as usize],
                &nodes[ndx],
                tinf,
            );
        }
        if i > 0 && nodes[ndx].strand == 1 && nodes[genes[i - 1].start_ndx as usize].strand == -1 {
            igm = intergenic_mod(
                &nodes[genes[i - 1].start_ndx as usize],
                &nodes[ndx],
                tinf,
            );
        }
        if i < ng as usize - 1
            && nodes[ndx].strand == -1
            && nodes[genes[i + 1].start_ndx as usize].strand == 1
        {
            igm = intergenic_mod(&nodes[ndx], &nodes[genes[i + 1].start_ndx as usize], tinf);
        }
        if i < ng as usize - 1
            && nodes[ndx].strand == -1
            && nodes[genes[i + 1].start_ndx as usize].strand == -1
        {
            igm = intergenic_mod(&nodes[ndx], &nodes[genes[i + 1].stop_ndx as usize], tinf);
        }

        let mut maxndx = [-1i32; 2];
        let mut maxsc = [0.0f64; 2];
        let mut maxigm = [0.0f64; 2];

        let search_start = if ndx >= 100 { ndx - 100 } else { 0 };
        let search_end = if ndx + 100 < nn as usize {
            ndx + 100
        } else {
            nn as usize
        };

        for j in search_start..search_end {
            if j == ndx {
                continue;
            }
            if nodes[j].type_ == STOP || nodes[j].stop_val != nodes[ndx].stop_val {
                continue;
            }

            let mut tigm = 0.0f64;

            if i > 0 && nodes[j].strand == 1 && nodes[genes[i - 1].start_ndx as usize].strand == 1
            {
                if nodes[genes[i - 1].stop_ndx as usize].ndx - nodes[j].ndx > MAX_SAM_OVLP {
                    continue;
                }
                tigm = intergenic_mod(&nodes[genes[i - 1].stop_ndx as usize], &nodes[j], tinf);
            }
            if i > 0 && nodes[j].strand == 1 && nodes[genes[i - 1].start_ndx as usize].strand == -1
            {
                if nodes[genes[i - 1].start_ndx as usize].ndx - nodes[j].ndx >= 0 {
                    continue;
                }
                tigm = intergenic_mod(&nodes[genes[i - 1].start_ndx as usize], &nodes[j], tinf);
            }
            if i < ng as usize - 1
                && nodes[j].strand == -1
                && nodes[genes[i + 1].start_ndx as usize].strand == 1
            {
                if nodes[j].ndx - nodes[genes[i + 1].start_ndx as usize].ndx >= 0 {
                    continue;
                }
                tigm = intergenic_mod(&nodes[j], &nodes[genes[i + 1].start_ndx as usize], tinf);
            }
            if i < ng as usize - 1
                && nodes[j].strand == -1
                && nodes[genes[i + 1].start_ndx as usize].strand == -1
            {
                if nodes[j].ndx - nodes[genes[i + 1].stop_ndx as usize].ndx > MAX_SAM_OVLP {
                    continue;
                }
                tigm = intergenic_mod(&nodes[j], &nodes[genes[i + 1].stop_ndx as usize], tinf);
            }

            if maxndx[0] == -1 {
                maxndx[0] = j as i32;
                maxsc[0] = nodes[j].cscore + nodes[j].sscore;
                maxigm[0] = tigm;
            } else if nodes[j].cscore + nodes[j].sscore + tigm > maxsc[0] {
                maxndx[1] = maxndx[0];
                maxsc[1] = maxsc[0];
                maxigm[1] = maxigm[0];
                maxndx[0] = j as i32;
                maxsc[0] = nodes[j].cscore + nodes[j].sscore;
                maxigm[0] = tigm;
            } else if maxndx[1] == -1 || nodes[j].cscore + nodes[j].sscore + tigm > maxsc[1] {
                maxndx[1] = j as i32;
                maxsc[1] = nodes[j].cscore + nodes[j].sscore;
                maxigm[1] = tigm;
            }
        }

        for j in 0..2 {
            let mndx = maxndx[j];
            if mndx == -1 {
                continue;
            }
            let mndx_usize = mndx as usize;

            if nodes[mndx_usize].tscore < nodes[ndx].tscore
                && maxsc[j] - nodes[mndx_usize].tscore >= sc - nodes[ndx].tscore + tinf.st_wt
                && nodes[mndx_usize].rscore > nodes[ndx].rscore
                && nodes[mndx_usize].uscore > nodes[ndx].uscore
                && nodes[mndx_usize].cscore > nodes[ndx].cscore
                && (nodes[mndx_usize].ndx - nodes[ndx].ndx).abs() > 15
            {
                maxsc[j] += nodes[ndx].tscore - nodes[mndx_usize].tscore;
            }

            else if (nodes[mndx_usize].ndx - nodes[ndx].ndx).abs() <= 15
                && nodes[mndx_usize].rscore + nodes[mndx_usize].tscore
                    > nodes[ndx].rscore + nodes[ndx].tscore
                && nodes[ndx].edge == 0
                && nodes[mndx_usize].edge == 0
            {
                if nodes[ndx].cscore > nodes[mndx_usize].cscore {
                    maxsc[j] += nodes[ndx].cscore - nodes[mndx_usize].cscore;
                }
                if nodes[ndx].uscore > nodes[mndx_usize].uscore {
                    maxsc[j] += nodes[ndx].uscore - nodes[mndx_usize].uscore;
                }
                if igm > maxigm[j] {
                    maxsc[j] += igm - maxigm[j];
                }
            } else {
                maxsc[j] = -1000.0;
            }
        }

        let mut mndx = -1i32;
        for j in 0..2 {
            if maxndx[j] == -1 {
                continue;
            }
            if mndx == -1 && maxsc[j] + maxigm[j] > sc + igm {
                mndx = j as i32;
            } else if mndx >= 0 && maxsc[j] + maxigm[j] > maxsc[mndx as usize] + maxigm[mndx as usize]
            {
                mndx = j as i32;
            }
        }

        if mndx != -1 {
            let best_ndx = maxndx[mndx as usize] as usize;
            if nodes[best_ndx].strand == 1 {
                genes[i].start_ndx = maxndx[mndx as usize];
                genes[i].begin = nodes[best_ndx].ndx + 1;
            } else if nodes[best_ndx].strand == -1 {
                genes[i].start_ndx = maxndx[mndx as usize];
                genes[i].end = nodes[best_ndx].ndx + 1;
            }
        }
    }
}

pub fn print_genes(
    output_ptr: &mut dyn std::io::Write,
    genes: &[Gene],
    ng: i32,
    nodes: &[Node],
    slen: i32,
    flag: i32,
    sctr: i32,
    is_meta: i32,
    mdesc: &str,
    tinf: &Training,
    header: &str,
    short_hdr: &str,
    version: &str,
) {

    let seq_data = format!("seqnum={};seqlen={};seqhdr=\"{}\"", sctr, slen, header);

    let mut run_data = if is_meta == 0 {
        format!("version=Rustygal.v{};run_type=Single;model=\"Ab initio\";", version)
    } else {
        format!(
            "version=Rustygal.v{};run_type=Metagenomic;model=\"{}\";",
            version, mdesc
        )
    };

    run_data.push_str(&format!(
        "gc_cont={:.2};transl_table={};uses_sd={}",
        tinf.gc * 100.0,
        tinf.trans_table,
        tinf.uses_sd
    ));

    if flag == 3 && sctr == 1 {
        let _ = writeln!(output_ptr, "##gff-version  3");
    }

    if flag == 0 {
        let _ = writeln!(output_ptr, "DEFINITION  {};{}", seq_data, run_data);
        let _ = writeln!(output_ptr, "FEATURES             Location/Qualifiers");
    } else if flag != 1 {
        let _ = writeln!(output_ptr, "# Sequence Data: {}", seq_data);
        let _ = writeln!(output_ptr, "# Model Data: {}", run_data);
    }

    for i in 0..ng as usize {
        let ndx = genes[i].start_ndx as usize;
        let sndx = genes[i].stop_ndx as usize;

        let gene_data_str = std::str::from_utf8(&genes[i].gene_data)
            .unwrap_or("")
            .trim_end_matches('\0');
        let score_data_str = std::str::from_utf8(&genes[i].score_data)
            .unwrap_or("")
            .trim_end_matches('\0');

        if nodes[ndx].strand == 1 {

            let left = if nodes[ndx].edge == 1 {
                format!("<{}", genes[i].begin)
            } else {
                format!("{}", genes[i].begin)
            };

            let right = if nodes[sndx].edge == 1 {
                format!(">{}", genes[i].end)
            } else {
                format!("{}", genes[i].end)
            };

            match flag {
                0 => {
                    let _ = writeln!(output_ptr, "     CDS             {}..{}", left, right);
                    let _ = writeln!(
                        output_ptr,
                        "                     /note=\"{};{}\"",
                        gene_data_str, score_data_str
                    );
                }
                1 => {
                    let _ = writeln!(
                        output_ptr,
                        "gene_prodigal={}|1|f|y|y|3|0|{}|{}|{}|{}|-1|-1|1.0",
                        i + 1,
                        genes[i].begin,
                        genes[i].end,
                        genes[i].begin,
                        genes[i].end
                    );
                }
                2 => {
                    let _ = writeln!(
                        output_ptr,
                        ">{}_{}_{}_ +",
                        i + 1,
                        genes[i].begin,
                        genes[i].end
                    );
                }
                3 => {
                    let _ = writeln!(
                        output_ptr,
                        "{}\tRustygal_v{}\tCDS\t{}\t{}\t{:.1}\t+\t0\t{};{}",
                        short_hdr,
                        version,
                        genes[i].begin,
                        genes[i].end,
                        nodes[ndx].cscore + nodes[ndx].sscore,
                        gene_data_str,
                        score_data_str
                    );
                }
                _ => {}
            }
        } else {

            let left = if nodes[sndx].edge == 1 {
                format!("<{}", genes[i].begin)
            } else {
                format!("{}", genes[i].begin)
            };

            let right = if nodes[ndx].edge == 1 {
                format!(">{}", genes[i].end)
            } else {
                format!("{}", genes[i].end)
            };

            match flag {
                0 => {
                    let _ = writeln!(
                        output_ptr,
                        "     CDS             complement({}..{})",
                        left, right
                    );
                    let _ = writeln!(
                        output_ptr,
                        "                     /note=\"{};{}\"",
                        gene_data_str, score_data_str
                    );
                }
                1 => {
                    let _ = writeln!(
                        output_ptr,
                        "gene_prodigal={}|1|r|y|y|3|0|{}|{}|{}|{}|-1|-1|1.0",
                        i + 1,
                        slen + 1 - genes[i].end,
                        slen + 1 - genes[i].begin,
                        slen + 1 - genes[i].end,
                        slen + 1 - genes[i].begin
                    );
                }
                2 => {
                    let _ = writeln!(
                        output_ptr,
                        ">{}_{}_{}_ -",
                        i + 1,
                        genes[i].begin,
                        genes[i].end
                    );
                }
                3 => {
                    let _ = writeln!(
                        output_ptr,
                        "{}\tRustygal_v{}\tCDS\t{}\t{}\t{:.1}\t-\t0\t{};{}",
                        short_hdr,
                        version,
                        genes[i].begin,
                        genes[i].end,
                        nodes[ndx].cscore + nodes[ndx].sscore,
                        gene_data_str,
                        score_data_str
                    );
                }
                _ => {}
            }
        }
    }

    if flag == 0 {
        let _ = writeln!(output_ptr, "//");
    }
}

pub fn write_translations(
    trans_ptr: &mut dyn std::io::Write,
    genes: &[Gene],
    ng: i32,
    nodes: &[Node],
    seq: &[u8],
    rseq: &[u8],
    useq: &[u8],
    slen: i32,
    tinf: &Training,
    _sctr: i32,
    short_hdr: &str,
) {
    use crate::sequence::{amino, is_n};

    for i in 0..ng as usize {
        let gene_data_str = std::str::from_utf8(&genes[i].gene_data)
            .unwrap_or("")
            .trim_end_matches('\0');

        if nodes[genes[i].start_ndx as usize].strand == 1 {

            let _ = writeln!(
                trans_ptr,
                ">{}_{} # {} # {} # 1 # {}",
                short_hdr,
                i + 1,
                genes[i].begin,
                genes[i].end,
                gene_data_str
            );

            let mut col = 0;
            for j in (genes[i].begin..genes[i].end).step_by(3) {
                if is_n(useq, j - 1) || is_n(useq, j) || is_n(useq, j + 1) {
                    let _ = write!(trans_ptr, "X");
                } else {
                    let is_start = j == genes[i].begin;
                    let is_complete = 1 - nodes[genes[i].start_ndx as usize].edge;
                    let _ = write!(
                        trans_ptr,
                        "{}",
                        amino(seq, j - 1, tinf, if is_start { is_complete } else { 0 }) as char
                    );
                }
                col += 1;
                if col % 60 == 0 {
                    let _ = writeln!(trans_ptr);
                }
            }
            if col % 60 != 0 {
                let _ = writeln!(trans_ptr);
            }
        } else {

            let _ = writeln!(
                trans_ptr,
                ">{}_{} # {} # {} # -1 # {}",
                short_hdr,
                i + 1,
                genes[i].begin,
                genes[i].end,
                gene_data_str
            );

            let mut col = 0;
            for j in (slen + 1 - genes[i].end..slen + 1 - genes[i].begin).step_by(3) {
                if is_n(useq, slen - j) || is_n(useq, slen - 1 - j) || is_n(useq, slen - 2 - j) {
                    let _ = write!(trans_ptr, "X");
                } else {
                    let is_start = j == slen + 1 - genes[i].end;
                    let is_complete = 1 - nodes[genes[i].start_ndx as usize].edge;
                    let _ = write!(
                        trans_ptr,
                        "{}",
                        amino(rseq, j - 1, tinf, if is_start { is_complete } else { 0 }) as char
                    );
                }
                col += 1;
                if col % 60 == 0 {
                    let _ = writeln!(trans_ptr);
                }
            }
            if col % 60 != 0 {
                let _ = writeln!(trans_ptr);
            }
        }
    }
}

pub fn write_nucleotide_seqs(
    nuc_ptr: &mut dyn std::io::Write,
    genes: &[Gene],
    ng: i32,
    nodes: &[Node],
    seq: &[u8],
    rseq: &[u8],
    useq: &[u8],
    slen: i32,
    _tinf: &Training,
    _sctr: i32,
    short_hdr: &str,
) {
    use crate::sequence::{is_a, is_c, is_g, is_n, is_t};

    for i in 0..ng as usize {
        let gene_data_str = std::str::from_utf8(&genes[i].gene_data)
            .unwrap_or("")
            .trim_end_matches('\0');

        if nodes[genes[i].start_ndx as usize].strand == 1 {

            let _ = writeln!(
                nuc_ptr,
                ">{}_{} # {} # {} # 1 # {}",
                short_hdr,
                i + 1,
                genes[i].begin,
                genes[i].end,
                gene_data_str
            );

            let mut col = 0;
            for j in genes[i].begin - 1..genes[i].end {
                let base = if is_a(seq, j) {
                    'A'
                } else if is_t(seq, j) {
                    'T'
                } else if is_g(seq, j) {
                    'G'
                } else if is_c(seq, j) && !is_n(useq, j) {
                    'C'
                } else {
                    'N'
                };
                let _ = write!(nuc_ptr, "{}", base);
                col += 1;
                if col % 70 == 0 {
                    let _ = writeln!(nuc_ptr);
                }
            }
            if col % 70 != 0 {
                let _ = writeln!(nuc_ptr);
            }
        } else {

            let _ = writeln!(
                nuc_ptr,
                ">{}_{} # {} # {} # -1 # {}",
                short_hdr,
                i + 1,
                genes[i].begin,
                genes[i].end,
                gene_data_str
            );

            let mut col = 0;
            for j in slen - genes[i].end..slen + 1 - genes[i].begin {
                let base = if is_a(rseq, j) {
                    'A'
                } else if is_t(rseq, j) {
                    'T'
                } else if is_g(rseq, j) {
                    'G'
                } else if is_c(rseq, j) && !is_n(useq, slen - 1 - j) {
                    'C'
                } else {
                    'N'
                };
                let _ = write!(nuc_ptr, "{}", base);
                col += 1;
                if col % 70 == 0 {
                    let _ = writeln!(nuc_ptr);
                }
            }
            if col % 70 != 0 {
                let _ = writeln!(nuc_ptr);
            }
        }
    }
}

pub fn calculate_confidence(score: f64, start_weight: f64) -> f64 {
    let conf = if score / start_weight < 41.0 {
        let exp_val = (score / start_weight).exp();
        (exp_val / (exp_val + 1.0)) * 100.0
    } else {
        99.99
    };

    if conf <= 50.0 {
        50.0
    } else {
        conf
    }
}

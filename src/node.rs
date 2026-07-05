// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2007-2016 Doug Hyatt, Univ. of Tennessee / UT-Battelle (original Prodigal C)
// Copyright (C) 2026 Sunju Kim (Rust reimplementation)

use crate::sequence::Mask;
use crate::training::Training;

pub const STT_NOD: usize = 200_000;

pub const MIN_GENE: i32 = 90;
pub const MIN_EDGE_GENE: i32 = 60;
pub const MAX_SAM_OVLP: i32 = 60;
pub const ST_WINDOW: i32 = 60;
pub const OPER_DIST: i32 = 60;
pub const EDGE_BONUS: f64 = 0.74;
pub const EDGE_UPS: f64 = -1.00;
pub const META_PEN: f64 = 7.5;

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Motif {

    pub ndx: i32,

    pub len: i32,

    pub spacer: i32,

    pub spacendx: i32,

    pub score: f64,
}

impl Default for Motif {
    fn default() -> Self {
        Self {
            ndx: 0,
            len: 0,
            spacer: 0,
            spacendx: 0,
            score: 0.0,
        }
    }
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Node {

    pub type_: i32,

    pub edge: i32,

    pub ndx: i32,

    pub strand: i32,

    pub stop_val: i32,

    pub star_ptr: [i32; 3],

    pub gc_bias: i32,

    pub gc_score: [f64; 3],

    pub cscore: f64,

    pub gc_cont: f64,

    pub rbs: [i32; 2],

    pub mot: Motif,

    pub uscore: f64,

    pub tscore: f64,

    pub rscore: f64,

    pub sscore: f64,

    pub traceb: i32,

    pub tracef: i32,

    pub ov_mark: i32,

    pub score: f64,

    pub elim: i32,
}

impl Default for Node {
    fn default() -> Self {
        Self {
            type_: 0,
            edge: 0,
            ndx: 0,
            strand: 0,
            stop_val: 0,
            star_ptr: [0; 3],
            gc_bias: 0,
            gc_score: [0.0; 3],
            cscore: 0.0,
            gc_cont: 0.0,
            rbs: [0; 2],
            mot: Motif::default(),
            uscore: 0.0,
            tscore: 0.0,
            rscore: 0.0,
            sscore: 0.0,
            traceb: 0,
            tracef: 0,
            ov_mark: 0,
            score: 0.0,
            elim: 0,
        }
    }
}

pub fn add_nodes(
    seq: &[u8],
    rseq: &[u8],
    slen: i32,
    nodes: &mut [Node],
    closed: i32,
    mlist: &[Mask],
    nmask: i32,
    tinf: &Training,
) -> i32 {
    use crate::sequence::{classify_start, is_stop, ATG, GTG, STOP, TTG};

    let mut nn = 0i32;
    let mut last = [0i32; 3];
    let mut saw_start = [0i32; 3];
    let mut min_dist = [0i32; 3];

    let slmod = slen % 3;
    for i in 0..3 {
        last[((i + slmod) % 3) as usize] = slen + i;
        saw_start[(i % 3) as usize] = 0;
        min_dist[(i % 3) as usize] = MIN_EDGE_GENE;
        if closed == 0 {
            while last[((i + slmod) % 3) as usize] + 2 > slen - 1 {
                last[((i + slmod) % 3) as usize] -= 3;
            }
        }
    }

    let mut i = slen - 3;
    while i >= 0 {

        if nn as usize + 4 >= nodes.len() {
            break;
        }
        let frame = (i % 3) as usize;

        if is_stop(seq, i, tinf) {
            if saw_start[frame] == 1 {
                if !is_stop(seq, last[frame], tinf) {
                    nodes[nn as usize].edge = 1;
                }
                nodes[nn as usize].ndx = last[frame];
                nodes[nn as usize].type_ = STOP;
                nodes[nn as usize].strand = 1;
                nodes[nn as usize].stop_val = i;
                nn += 1;
            }
            min_dist[frame] = MIN_GENE;
            last[frame] = i;
            saw_start[frame] = 0;
            i -= 1;
            continue;
        }

        if last[frame] >= slen {
            i -= 1;
            continue;
        }

        let stype_fwd = classify_start(seq, i, tinf);
        if stype_fwd == ATG
            && ((last[frame] - i + 3) >= min_dist[frame])
            && !cross_mask(i, last[frame], mlist, nmask)
        {
            nodes[nn as usize].ndx = i;
            nodes[nn as usize].type_ = ATG;
            saw_start[frame] = 1;
            nodes[nn as usize].stop_val = last[frame];
            nodes[nn as usize].strand = 1;
            nn += 1;
        } else if stype_fwd == GTG
            && ((last[frame] - i + 3) >= min_dist[frame])
            && !cross_mask(i, last[frame], mlist, nmask)
        {
            nodes[nn as usize].ndx = i;
            nodes[nn as usize].type_ = GTG;
            saw_start[frame] = 1;
            nodes[nn as usize].stop_val = last[frame];
            nodes[nn as usize].strand = 1;
            nn += 1;
        } else if stype_fwd == TTG
            && ((last[frame] - i + 3) >= min_dist[frame])
            && !cross_mask(i, last[frame], mlist, nmask)
        {
            nodes[nn as usize].ndx = i;
            nodes[nn as usize].type_ = TTG;
            saw_start[frame] = 1;
            nodes[nn as usize].stop_val = last[frame];
            nodes[nn as usize].strand = 1;
            nn += 1;
        } else if i <= 2
            && closed == 0
            && ((last[frame] - i) > MIN_EDGE_GENE)
            && !cross_mask(i, last[frame], mlist, nmask)
        {
            nodes[nn as usize].ndx = i;
            nodes[nn as usize].type_ = ATG;
            saw_start[frame] = 1;
            nodes[nn as usize].edge = 1;
            nodes[nn as usize].stop_val = last[frame];
            nodes[nn as usize].strand = 1;
            nn += 1;
        }

        i -= 1;
    }

    for i in 0..3 {
        let frame = (i % 3) as usize;
        if saw_start[frame] == 1 {
            if !is_stop(seq, last[frame], tinf) {
                nodes[nn as usize].edge = 1;
            }
            nodes[nn as usize].ndx = last[frame];
            nodes[nn as usize].type_ = STOP;
            nodes[nn as usize].strand = 1;
            nodes[nn as usize].stop_val = i - 6;
            nn += 1;
        }
    }

    for i in 0..3 {
        last[((i + slmod) % 3) as usize] = slen + i;
        saw_start[(i % 3) as usize] = 0;
        min_dist[(i % 3) as usize] = MIN_EDGE_GENE;
        if closed == 0 {
            while last[((i + slmod) % 3) as usize] + 2 > slen - 1 {
                last[((i + slmod) % 3) as usize] -= 3;
            }
        }
    }

    let mut i = slen - 3;
    while i >= 0 {

        if nn as usize + 4 >= nodes.len() {
            break;
        }
        let frame = (i % 3) as usize;

        if is_stop(rseq, i, tinf) {
            if saw_start[frame] == 1 {
                if !is_stop(rseq, last[frame], tinf) {
                    nodes[nn as usize].edge = 1;
                }
                nodes[nn as usize].ndx = slen - last[frame] - 1;
                nodes[nn as usize].type_ = STOP;
                nodes[nn as usize].strand = -1;
                nodes[nn as usize].stop_val = slen - i - 1;
                nn += 1;
            }
            min_dist[frame] = MIN_GENE;
            last[frame] = i;
            saw_start[frame] = 0;
            i -= 1;
            continue;
        }

        if last[frame] >= slen {
            i -= 1;
            continue;
        }

        let stype_rev = classify_start(rseq, i, tinf);
        if stype_rev == ATG
            && ((last[frame] - i + 3) >= min_dist[frame])
            && !cross_mask(slen - last[frame] - 1, slen - i - 1, mlist, nmask)
        {
            nodes[nn as usize].ndx = slen - i - 1;
            nodes[nn as usize].type_ = ATG;
            saw_start[frame] = 1;
            nodes[nn as usize].stop_val = slen - last[frame] - 1;
            nodes[nn as usize].strand = -1;
            nn += 1;
        } else if stype_rev == GTG
            && ((last[frame] - i + 3) >= min_dist[frame])
            && !cross_mask(slen - last[frame] - 1, slen - i - 1, mlist, nmask)
        {
            nodes[nn as usize].ndx = slen - i - 1;
            nodes[nn as usize].type_ = GTG;
            saw_start[frame] = 1;
            nodes[nn as usize].stop_val = slen - last[frame] - 1;
            nodes[nn as usize].strand = -1;
            nn += 1;
        } else if stype_rev == TTG
            && ((last[frame] - i + 3) >= min_dist[frame])
            && !cross_mask(slen - last[frame] - 1, slen - i - 1, mlist, nmask)
        {
            nodes[nn as usize].ndx = slen - i - 1;
            nodes[nn as usize].type_ = TTG;
            saw_start[frame] = 1;
            nodes[nn as usize].stop_val = slen - last[frame] - 1;
            nodes[nn as usize].strand = -1;
            nn += 1;
        } else if i <= 2
            && closed == 0
            && ((last[frame] - i) > MIN_EDGE_GENE)
            && !cross_mask(slen - last[frame] - 1, slen - i - 1, mlist, nmask)
        {
            nodes[nn as usize].ndx = slen - i - 1;
            nodes[nn as usize].type_ = ATG;
            saw_start[frame] = 1;
            nodes[nn as usize].edge = 1;
            nodes[nn as usize].stop_val = slen - last[frame] - 1;
            nodes[nn as usize].strand = -1;
            nn += 1;
        }

        i -= 1;
    }

    for i in 0..3 {
        let frame = (i % 3) as usize;
        if saw_start[frame] == 1 {
            if !is_stop(rseq, last[frame], tinf) {
                nodes[nn as usize].edge = 1;
            }
            nodes[nn as usize].ndx = slen - last[frame] - 1;
            nodes[nn as usize].type_ = STOP;
            nodes[nn as usize].strand = -1;
            nodes[nn as usize].stop_val = slen - i + 5;
            nn += 1;
        }
    }

    nn
}

fn add_nodes_forward(
    seq: &[u8],
    slen: i32,
    closed: i32,
    mlist: &[Mask],
    nmask: i32,
    tinf: &Training,
) -> (Vec<Node>, i32) {
    use crate::sequence::{classify_start, is_stop, ATG, GTG, STOP, TTG};

    let mut nodes = vec![Node::default(); STT_NOD];
    let mut nn = 0i32;
    let mut last = [0i32; 3];
    let mut saw_start = [0i32; 3];
    let mut min_dist = [0i32; 3];

    let slmod = slen % 3;
    for i in 0..3 {
        last[((i + slmod) % 3) as usize] = slen + i;
        saw_start[(i % 3) as usize] = 0;
        min_dist[(i % 3) as usize] = MIN_EDGE_GENE;
        if closed == 0 {
            while last[((i + slmod) % 3) as usize] + 2 > slen - 1 {
                last[((i + slmod) % 3) as usize] -= 3;
            }
        }
    }

    let mut i = slen - 3;
    while i >= 0 {

        if nn as usize + 4 >= nodes.len() {
            break;
        }
        let frame = (i % 3) as usize;

        if is_stop(seq, i, tinf) {
            if saw_start[frame] == 1 {
                if !is_stop(seq, last[frame], tinf) {
                    nodes[nn as usize].edge = 1;
                }
                nodes[nn as usize].ndx = last[frame];
                nodes[nn as usize].type_ = STOP;
                nodes[nn as usize].strand = 1;
                nodes[nn as usize].stop_val = i;
                nn += 1;
            }
            min_dist[frame] = MIN_GENE;
            last[frame] = i;
            saw_start[frame] = 0;
            i -= 1;
            continue;
        }

        if last[frame] >= slen {
            i -= 1;
            continue;
        }

        let stype_fwd = classify_start(seq, i, tinf);
        if stype_fwd == ATG
            && ((last[frame] - i + 3) >= min_dist[frame])
            && !cross_mask(i, last[frame], mlist, nmask)
        {
            nodes[nn as usize].ndx = i;
            nodes[nn as usize].type_ = ATG;
            saw_start[frame] = 1;
            nodes[nn as usize].stop_val = last[frame];
            nodes[nn as usize].strand = 1;
            nn += 1;
        } else if stype_fwd == GTG
            && ((last[frame] - i + 3) >= min_dist[frame])
            && !cross_mask(i, last[frame], mlist, nmask)
        {
            nodes[nn as usize].ndx = i;
            nodes[nn as usize].type_ = GTG;
            saw_start[frame] = 1;
            nodes[nn as usize].stop_val = last[frame];
            nodes[nn as usize].strand = 1;
            nn += 1;
        } else if stype_fwd == TTG
            && ((last[frame] - i + 3) >= min_dist[frame])
            && !cross_mask(i, last[frame], mlist, nmask)
        {
            nodes[nn as usize].ndx = i;
            nodes[nn as usize].type_ = TTG;
            saw_start[frame] = 1;
            nodes[nn as usize].stop_val = last[frame];
            nodes[nn as usize].strand = 1;
            nn += 1;
        } else if i <= 2
            && closed == 0
            && ((last[frame] - i) > MIN_EDGE_GENE)
            && !cross_mask(i, last[frame], mlist, nmask)
        {
            nodes[nn as usize].ndx = i;
            nodes[nn as usize].type_ = ATG;
            saw_start[frame] = 1;
            nodes[nn as usize].edge = 1;
            nodes[nn as usize].stop_val = last[frame];
            nodes[nn as usize].strand = 1;
            nn += 1;
        }

        i -= 1;
    }

    for i in 0..3 {
        let frame = (i % 3) as usize;
        if saw_start[frame] == 1 {
            if !is_stop(seq, last[frame], tinf) {
                nodes[nn as usize].edge = 1;
            }
            nodes[nn as usize].ndx = last[frame];
            nodes[nn as usize].type_ = STOP;
            nodes[nn as usize].strand = 1;
            nodes[nn as usize].stop_val = i - 6;
            nn += 1;
        }
    }

    (nodes, nn)
}

fn add_nodes_reverse(
    rseq: &[u8],
    slen: i32,
    closed: i32,
    mlist: &[Mask],
    nmask: i32,
    tinf: &Training,
) -> (Vec<Node>, i32) {
    use crate::sequence::{classify_start, is_stop, ATG, GTG, STOP, TTG};

    let mut nodes = vec![Node::default(); STT_NOD];
    let mut nn = 0i32;
    let mut last = [0i32; 3];
    let mut saw_start = [0i32; 3];
    let mut min_dist = [0i32; 3];

    let slmod = slen % 3;
    for i in 0..3 {
        last[((i + slmod) % 3) as usize] = slen + i;
        saw_start[(i % 3) as usize] = 0;
        min_dist[(i % 3) as usize] = MIN_EDGE_GENE;
        if closed == 0 {
            while last[((i + slmod) % 3) as usize] + 2 > slen - 1 {
                last[((i + slmod) % 3) as usize] -= 3;
            }
        }
    }

    let mut i = slen - 3;
    while i >= 0 {

        if nn as usize + 4 >= nodes.len() {
            break;
        }
        let frame = (i % 3) as usize;

        if is_stop(rseq, i, tinf) {
            if saw_start[frame] == 1 {
                if !is_stop(rseq, last[frame], tinf) {
                    nodes[nn as usize].edge = 1;
                }
                nodes[nn as usize].ndx = slen - last[frame] - 1;
                nodes[nn as usize].type_ = STOP;
                nodes[nn as usize].strand = -1;
                nodes[nn as usize].stop_val = slen - i - 1;
                nn += 1;
            }
            min_dist[frame] = MIN_GENE;
            last[frame] = i;
            saw_start[frame] = 0;
            i -= 1;
            continue;
        }

        if last[frame] >= slen {
            i -= 1;
            continue;
        }

        let stype_rev = classify_start(rseq, i, tinf);
        if stype_rev == ATG
            && ((last[frame] - i + 3) >= min_dist[frame])
            && !cross_mask(slen - last[frame] - 1, slen - i - 1, mlist, nmask)
        {
            nodes[nn as usize].ndx = slen - i - 1;
            nodes[nn as usize].type_ = ATG;
            saw_start[frame] = 1;
            nodes[nn as usize].stop_val = slen - last[frame] - 1;
            nodes[nn as usize].strand = -1;
            nn += 1;
        } else if stype_rev == GTG
            && ((last[frame] - i + 3) >= min_dist[frame])
            && !cross_mask(slen - last[frame] - 1, slen - i - 1, mlist, nmask)
        {
            nodes[nn as usize].ndx = slen - i - 1;
            nodes[nn as usize].type_ = GTG;
            saw_start[frame] = 1;
            nodes[nn as usize].stop_val = slen - last[frame] - 1;
            nodes[nn as usize].strand = -1;
            nn += 1;
        } else if stype_rev == TTG
            && ((last[frame] - i + 3) >= min_dist[frame])
            && !cross_mask(slen - last[frame] - 1, slen - i - 1, mlist, nmask)
        {
            nodes[nn as usize].ndx = slen - i - 1;
            nodes[nn as usize].type_ = TTG;
            saw_start[frame] = 1;
            nodes[nn as usize].stop_val = slen - last[frame] - 1;
            nodes[nn as usize].strand = -1;
            nn += 1;
        } else if i <= 2
            && closed == 0
            && ((last[frame] - i) > MIN_EDGE_GENE)
            && !cross_mask(slen - last[frame] - 1, slen - i - 1, mlist, nmask)
        {
            nodes[nn as usize].ndx = slen - i - 1;
            nodes[nn as usize].type_ = ATG;
            saw_start[frame] = 1;
            nodes[nn as usize].edge = 1;
            nodes[nn as usize].stop_val = slen - last[frame] - 1;
            nodes[nn as usize].strand = -1;
            nn += 1;
        }

        i -= 1;
    }

    for i in 0..3 {
        let frame = (i % 3) as usize;
        if saw_start[frame] == 1 {
            if !is_stop(rseq, last[frame], tinf) {
                nodes[nn as usize].edge = 1;
            }
            nodes[nn as usize].ndx = slen - last[frame] - 1;
            nodes[nn as usize].type_ = STOP;
            nodes[nn as usize].strand = -1;
            nodes[nn as usize].stop_val = slen - i + 5;
            nn += 1;
        }
    }

    (nodes, nn)
}

pub fn add_nodes_parallel(
    seq: &[u8],
    rseq: &[u8],
    slen: i32,
    nodes: &mut [Node],
    closed: i32,
    mlist: &[Mask],
    nmask: i32,
    tinf: &Training,
) -> i32 {

    let ((nodes_fwd, nn_fwd), (nodes_rev, nn_rev)) = rayon::join(
        || add_nodes_forward(seq, slen, closed, mlist, nmask, tinf),
        || add_nodes_reverse(rseq, slen, closed, mlist, nmask, tinf),
    );

    nodes[0..nn_fwd as usize].copy_from_slice(&nodes_fwd[0..nn_fwd as usize]);
    nodes[nn_fwd as usize..(nn_fwd + nn_rev) as usize]
        .copy_from_slice(&nodes_rev[0..nn_rev as usize]);

    nn_fwd + nn_rev
}

pub fn reset_node_scores(nodes: &mut [Node], nn: i32) {
    for i in 0..nn as usize {
        for j in 0..3 {
            nodes[i].star_ptr[j] = 0;
            nodes[i].gc_score[j] = 0.0;
        }
        for j in 0..2 {
            nodes[i].rbs[j] = 0;
        }
        nodes[i].score = 0.0;
        nodes[i].cscore = 0.0;
        nodes[i].sscore = 0.0;
        nodes[i].rscore = 0.0;
        nodes[i].tscore = 0.0;
        nodes[i].uscore = 0.0;
        nodes[i].traceb = -1;
        nodes[i].tracef = -1;
        nodes[i].ov_mark = -1;
        nodes[i].elim = 0;
        nodes[i].gc_bias = 0;
        nodes[i].mot = Motif::default();
    }
}

pub fn compare_nodes(a: &Node, b: &Node) -> std::cmp::Ordering {
    match a.ndx.cmp(&b.ndx) {
        std::cmp::Ordering::Equal => b.strand.cmp(&a.strand),
        other => other,
    }
}

pub fn stopcmp_nodes(a: &Node, b: &Node) -> std::cmp::Ordering {
    match a.stop_val.cmp(&b.stop_val) {
        std::cmp::Ordering::Equal => match b.strand.cmp(&a.strand) {

            std::cmp::Ordering::Equal => a.ndx.cmp(&b.ndx),
            other => other,
        },
        other => other,
    }
}

pub fn record_overlapping_starts(nodes: &mut [Node], nn: i32, tinf: &Training, flag: i32) {
    use crate::sequence::STOP;

    for i in 0..nn as usize {
        for j in 0..3 {
            nodes[i].star_ptr[j] = -1;
        }
        if nodes[i].type_ != STOP || nodes[i].edge == 1 {
            continue;
        }

        if nodes[i].strand == 1 {

            let mut max_sc = -100.0;
            let mut j = i as i32 + 3;
            while j >= 0 {
                let j_usize = j as usize;
                if j_usize >= nn as usize || nodes[j_usize].ndx > nodes[i].ndx + 2 {
                    j -= 1;
                    continue;
                }
                if nodes[j_usize].ndx + MAX_SAM_OVLP < nodes[i].ndx {
                    break;
                }
                if nodes[j_usize].strand == 1 && nodes[j_usize].type_ != STOP {
                    if nodes[j_usize].stop_val <= nodes[i].ndx {
                        j -= 1;
                        continue;
                    }
                    let fr = (nodes[j_usize].ndx % 3) as usize;
                    if flag == 0 && nodes[i].star_ptr[fr] == -1 {
                        nodes[i].star_ptr[fr] = j;
                    } else if flag == 1 {
                        let sc = nodes[j_usize].cscore
                            + nodes[j_usize].sscore
                            + intergenic_mod(&nodes[i], &nodes[j_usize], tinf);
                        if sc > max_sc {
                            nodes[i].star_ptr[fr] = j;
                            max_sc = sc;
                        }
                    }
                }
                j -= 1;
            }
        } else {

            let mut max_sc = -100.0;
            for j in (i as i32 - 3)..(nn as i32) {
                if j < 0 || nodes[j as usize].ndx < nodes[i].ndx - 2 {
                    continue;
                }
                if nodes[j as usize].ndx - MAX_SAM_OVLP > nodes[i].ndx {
                    break;
                }
                if nodes[j as usize].strand == -1 && nodes[j as usize].type_ != STOP {
                    if nodes[j as usize].stop_val >= nodes[i].ndx {
                        continue;
                    }
                    let fr = (nodes[j as usize].ndx % 3) as usize;
                    if flag == 0 && nodes[i].star_ptr[fr] == -1 {
                        nodes[i].star_ptr[fr] = j;
                    } else if flag == 1 {
                        let sc = nodes[j as usize].cscore
                            + nodes[j as usize].sscore
                            + intergenic_mod(&nodes[j as usize], &nodes[i], tinf);
                        if sc > max_sc {
                            nodes[i].star_ptr[fr] = j;
                            max_sc = sc;
                        }
                    }
                }
            }
        }
    }
}

pub fn record_gc_bias(gc: &[i32], nodes: &mut [Node], nn: i32, tinf: &mut Training) {
    use crate::sequence::{max_fr, STOP};

    if nn == 0 {
        return;
    }

    let mut ctr = [[0i32; 3]; 3];
    let mut last = [0i32; 3];

    for i in (0..nn).rev() {
        let i_usize = i as usize;
        let fr = (nodes[i_usize].ndx % 3) as usize;
        let frmod = 3 - fr as i32;

        if nodes[i_usize].strand == 1 && nodes[i_usize].type_ == STOP {
            for j in 0..3 {
                ctr[fr][j] = 0;
            }
            last[fr] = nodes[i_usize].ndx;
            ctr[fr][((gc[nodes[i_usize].ndx as usize] + frmod) % 3) as usize] = 1;
        } else if nodes[i_usize].strand == 1 {
            let mut j = last[fr] - 3;
            while j >= nodes[i_usize].ndx {
                ctr[fr][((gc[j as usize] + frmod) % 3) as usize] += 1;
                j -= 3;
            }
            let mfr = max_fr(ctr[fr][0], ctr[fr][1], ctr[fr][2]);
            nodes[i_usize].gc_bias = mfr;
            for j in 0..3 {
                nodes[i_usize].gc_score[j] = 3.0 * ctr[fr][j] as f64;
                nodes[i_usize].gc_score[j] /=
                    (nodes[i_usize].stop_val - nodes[i_usize].ndx + 3) as f64;
            }
            last[fr] = nodes[i_usize].ndx;
        }
    }

    for i in 0..nn as usize {
        let fr = (nodes[i].ndx % 3) as usize;
        let frmod = fr as i32;

        if nodes[i].strand == -1 && nodes[i].type_ == STOP {
            for j in 0..3 {
                ctr[fr][j] = 0;
            }
            last[fr] = nodes[i].ndx;
            ctr[fr][(((3 - gc[nodes[i].ndx as usize]) + frmod) % 3) as usize] = 1;
        } else if nodes[i].strand == -1 {
            let mut j = last[fr] + 3;
            while j <= nodes[i].ndx {
                ctr[fr][(((3 - gc[j as usize]) + frmod) % 3) as usize] += 1;
                j += 3;
            }
            let mfr = max_fr(ctr[fr][0], ctr[fr][1], ctr[fr][2]);
            nodes[i].gc_bias = mfr;
            for j in 0..3 {
                nodes[i].gc_score[j] = 3.0 * ctr[fr][j] as f64;
                nodes[i].gc_score[j] /= (nodes[i].ndx - nodes[i].stop_val + 3) as f64;
            }
            last[fr] = nodes[i].ndx;
        }
    }

    for i in 0..3 {
        tinf.bias[i] = 0.0;
    }
    for i in 0..nn as usize {
        if nodes[i].type_ != STOP {
            let len = (nodes[i].stop_val - nodes[i].ndx).abs() + 1;
            tinf.bias[nodes[i].gc_bias as usize] +=
                (nodes[i].gc_score[nodes[i].gc_bias as usize] * len as f64) / 1000.0;
        }
    }
    let tot = tinf.bias[0] + tinf.bias[1] + tinf.bias[2];
    for i in 0..3 {
        tinf.bias[i] *= 3.0 / tot;
    }
}

pub fn calc_dicodon_gene(
    tinf: &mut Training,
    seq: &[u8],
    rseq: &[u8],
    slen: i32,
    nodes: &[Node],
    dbeg: i32,
) {
    use crate::sequence::{calc_mer_bg, mer_ndx, STOP};

    let mut counts = [0i32; 4096];
    let mut prob = [0.0f64; 4096];
    let mut bg = [0.0f64; 4096];
    let mut glob = 0i32;

    let mut left = -1;
    let mut right = -1;
    calc_mer_bg(6, seq, rseq, slen, &mut bg);

    let mut path = dbeg;
    let mut in_gene = 0;

    while path != -1 {
        let path_usize = path as usize;

        if nodes[path_usize].strand == -1 && nodes[path_usize].type_ != STOP {
            in_gene = -1;
            left = slen - nodes[path_usize].ndx - 1;
        }
        if nodes[path_usize].strand == 1 && nodes[path_usize].type_ == STOP {
            in_gene = 1;
            right = nodes[path_usize].ndx + 2;
        }
        if in_gene == -1 && nodes[path_usize].strand == -1 && nodes[path_usize].type_ == STOP {
            right = slen - nodes[path_usize].ndx + 1;
            let _count_before = glob;
            let mut i = left;
            while i < right - 5 {
                counts[mer_ndx(6, rseq, i) as usize] += 1;
                glob += 1;
                i += 3;
            }
            in_gene = 0;
        }
        if in_gene == 1 && nodes[path_usize].strand == 1 && nodes[path_usize].type_ != STOP {
            left = nodes[path_usize].ndx;
            let _count_before = glob;
            let mut i = left;
            while i < right - 5 {
                counts[mer_ndx(6, seq, i) as usize] += 1;
                glob += 1;
                i += 3;
            }
            in_gene = 0;
        }
        path = nodes[path_usize].traceb;
    }

    for i in 0..4096 {
        prob[i] = counts[i] as f64 / glob as f64;
        if prob[i] == 0.0 && bg[i] != 0.0 {
            tinf.gene_dc[i] = -5.0;
        } else if bg[i] == 0.0 {
            tinf.gene_dc[i] = 0.0;
        } else {
            tinf.gene_dc[i] = (prob[i] / bg[i]).ln();
        }
        if tinf.gene_dc[i] > 5.0 {
            tinf.gene_dc[i] = 5.0;
        }
        if tinf.gene_dc[i] < -5.0 {
            tinf.gene_dc[i] = -5.0;
        }
    }
}

pub fn score_nodes(
    seq: &[u8],
    rseq: &[u8],
    slen: i32,
    nodes: &mut [Node],
    nn: i32,
    tinf: &Training,
    closed: i32,
    is_meta: i32,
) {
    use crate::sequence::{is_stop, STOP};

    calc_orf_gc(seq, rseq, slen, nodes, nn, tinf);
    raw_coding_score(seq, rseq, slen, nodes, nn, tinf);

    if tinf.uses_sd == 1 {
        rbs_score(seq, rseq, slen, nodes, nn, tinf);
    } else {
        for i in 0..nn as usize {
            if nodes[i].type_ == STOP || nodes[i].edge == 1 {
                continue;
            }
            find_best_upstream_motif(tinf, seq, rseq, slen, &mut nodes[i], 2);
        }
    }

    for i in 0..nn as usize {
        if nodes[i].type_ == STOP {
            continue;
        }

        let mut edge_gene = 0;
        if nodes[i].edge == 1 {
            edge_gene += 1;
        }
        if (nodes[i].strand == 1 && !is_stop(seq, nodes[i].stop_val, tinf))
            || (nodes[i].strand == -1 && !is_stop(rseq, slen - 1 - nodes[i].stop_val, tinf))
        {
            edge_gene += 1;
        }

        if nodes[i].edge == 1 {
            nodes[i].tscore = EDGE_BONUS * tinf.st_wt / edge_gene as f64;
            nodes[i].uscore = 0.0;
            nodes[i].rscore = 0.0;
        } else {

            nodes[i].tscore = tinf.type_wt[nodes[i].type_ as usize] * tinf.st_wt;

            let rbs1 = tinf.rbs_wt[nodes[i].rbs[0] as usize];
            let rbs2 = tinf.rbs_wt[nodes[i].rbs[1] as usize];
            let sd_score = dmax(rbs1, rbs2) * tinf.st_wt;
            if tinf.uses_sd == 1 {
                nodes[i].rscore = sd_score;
            } else {
                nodes[i].rscore = tinf.st_wt * nodes[i].mot.score;
                if nodes[i].rscore < sd_score && tinf.no_mot > -0.5 {
                    nodes[i].rscore = sd_score;
                }
            }

            if nodes[i].strand == 1 {
                score_upstream_composition(seq, slen, &mut nodes[i], tinf);
            } else {
                score_upstream_composition(rseq, slen, &mut nodes[i], tinf);
            }

            if closed == 0 && nodes[i].ndx <= 2 && nodes[i].strand == 1 {
                nodes[i].uscore += EDGE_UPS * tinf.st_wt;
            } else if closed == 0 && nodes[i].ndx >= slen - 3 && nodes[i].strand == -1 {
                nodes[i].uscore += EDGE_UPS * tinf.st_wt;
            } else if i < 500 && nodes[i].strand == 1 {
                for j in (0..i).rev() {
                    if nodes[j].edge == 1 && nodes[i].stop_val == nodes[j].stop_val {
                        nodes[i].uscore += EDGE_UPS * tinf.st_wt;
                        break;
                    }
                }
            } else if (i as i32) >= nn - 500 && nodes[i].strand == -1 {

                for j in (i + 1)..(nn as usize) {
                    if nodes[j].edge == 1 && nodes[i].stop_val == nodes[j].stop_val {
                        nodes[i].uscore += EDGE_UPS * tinf.st_wt;
                        break;
                    }
                }
            }
        }

        if ((nodes[i].ndx <= 2 && nodes[i].strand == 1)
            || (nodes[i].ndx >= slen - 3 && nodes[i].strand == -1))
            && nodes[i].edge == 0
            && closed == 0
        {
            edge_gene += 1;
            nodes[i].edge = 1;
            nodes[i].tscore = 0.0;
            nodes[i].uscore = EDGE_BONUS * tinf.st_wt / edge_gene as f64;
            nodes[i].rscore = 0.0;
        }

        if nodes[i].edge == 0 && edge_gene == 1 {
            nodes[i].uscore -= 0.5 * EDGE_BONUS * tinf.st_wt;
        }

        if edge_gene == 0 && (nodes[i].ndx - nodes[i].stop_val).abs() < 250 {

            let negf = 250.0 / ((nodes[i].ndx - nodes[i].stop_val).abs() as f32) as f64;
            let posf = ((nodes[i].ndx - nodes[i].stop_val).abs() as f32) as f64 / 250.0;
            if nodes[i].rscore < 0.0 {
                nodes[i].rscore *= negf;
            }
            if nodes[i].uscore < 0.0 {
                nodes[i].uscore *= negf;
            }
            if nodes[i].tscore < 0.0 {
                nodes[i].tscore *= negf;
            }
            if nodes[i].rscore > 0.0 {
                nodes[i].rscore *= posf;
            }
            if nodes[i].uscore > 0.0 {
                nodes[i].uscore *= posf;
            }
            if nodes[i].tscore > 0.0 {
                nodes[i].tscore *= posf;
            }
        }

        if is_meta == 1
            && slen < 3000
            && edge_gene == 0
            && (nodes[i].cscore < 5.0 || (nodes[i].ndx - nodes[i].stop_val).abs() < 120)
        {
            nodes[i].cscore -= META_PEN * dmax(0.0, (3000 - slen) as f64 / 2700.0);
        }

        nodes[i].sscore = nodes[i].tscore + nodes[i].rscore + nodes[i].uscore;

        if nodes[i].cscore < 0.0 {
            if edge_gene > 0 && nodes[i].edge == 0 {
                if is_meta == 0 || slen > 1500 {
                    nodes[i].sscore -= tinf.st_wt;
                } else {
                    nodes[i].sscore -= 10.31 - 0.004 * slen as f64;
                }
            } else if is_meta == 1 && slen < 3000 && nodes[i].edge == 1 {
                let min_meta_len = (slen as f64).sqrt() * 5.0;
                if (nodes[i].ndx - nodes[i].stop_val).abs() as f64 >= min_meta_len {
                    if nodes[i].cscore >= 0.0 {
                        nodes[i].cscore = -1.0;
                    }
                    nodes[i].sscore = 0.0;
                    nodes[i].uscore = 0.0;
                }
            } else {
                nodes[i].sscore -= 0.5;
            }
        } else if nodes[i].cscore < 5.0
            && is_meta == 1
            && (nodes[i].ndx - nodes[i].stop_val).abs() < 120
            && nodes[i].sscore < 0.0
        {
            nodes[i].sscore -= tinf.st_wt;
        }
    }

    for i in 0..nn as usize {
        if nodes[i].strand == -1
            && nodes[i].stop_val >= 361240
            && nodes[i].stop_val <= 361260
            && nodes[i].type_ != STOP
        {
            let _type_str = match nodes[i].type_ {
                0 => "ATG",
                1 => "GTG",
                2 => "TTG",
                _ => "UNK",
            };
        }
    }
}

pub fn raw_coding_score(
    seq: &[u8],
    rseq: &[u8],
    slen: i32,
    nodes: &mut [Node],
    nn: i32,
    tinf: &Training,
) {
    use crate::sequence::mer_ndx;

    let mut last = [0i32; 3];
    let mut score = [0.0f64; 3];

    let no_stop = if tinf.trans_table != 11 {

        let ns = ((1.0 - tinf.gc) * (1.0 - tinf.gc) * tinf.gc) / 8.0
            + ((1.0 - tinf.gc) * (1.0 - tinf.gc) * (1.0 - tinf.gc)) / 8.0;
        1.0 - ns
    } else {
        let ns = ((1.0 - tinf.gc) * (1.0 - tinf.gc) * tinf.gc) / 4.0
            + ((1.0 - tinf.gc) * (1.0 - tinf.gc) * (1.0 - tinf.gc)) / 8.0;
        1.0 - ns
    };

    for i in 0..3 {
        score[i] = 0.0;
    }
    for i in (0..nn).rev() {
        let fr = (nodes[i as usize].ndx % 3) as usize;
        if nodes[i as usize].strand == 1 && nodes[i as usize].type_ == crate::sequence::STOP {
            last[fr] = nodes[i as usize].ndx;
            score[fr] = 0.0;
        } else if nodes[i as usize].strand == 1 {
            let mut j = last[fr] - 3;
            while j >= nodes[i as usize].ndx {
                let idx = mer_ndx(6, seq, j) as usize;

                score[fr] += tinf.gene_dc[idx];
                j -= 3;
            }
            nodes[i as usize].cscore = score[fr];
            last[fr] = nodes[i as usize].ndx;
        }
    }

    for i in 0..3 {
        score[i] = 0.0;
    }
    for i in 0..nn {
        let fr = (nodes[i as usize].ndx % 3) as usize;
        if nodes[i as usize].strand == -1 && nodes[i as usize].type_ == crate::sequence::STOP {
            last[fr] = nodes[i as usize].ndx;
            score[fr] = 0.0;
        } else if nodes[i as usize].strand == -1 {

            let _is_debug_node = nodes[i as usize].strand == -1
                && nodes[i as usize].ndx == 361763
                && nodes[i as usize].stop_val == 361250;

            let mut j = last[fr] + 3;
            while j <= nodes[i as usize].ndx {
                let idx = mer_ndx(6, rseq, slen - j - 1) as usize;
                let gene_dc_val = tinf.gene_dc[idx];
                score[fr] += gene_dc_val;
                j += 3;
            }

            nodes[i as usize].cscore = score[fr];
            last[fr] = nodes[i as usize].ndx;

            if nodes[i as usize].strand == -1
                && nodes[i as usize].stop_val >= 361240
                && nodes[i as usize].stop_val <= 361260
            {
            }
        }
    }

    for i in 0..3 {
        score[i] = -10000.0;
    }
    for i in 0..nn {
        let fr = (nodes[i as usize].ndx % 3) as usize;
        if nodes[i as usize].strand == 1 && nodes[i as usize].type_ == crate::sequence::STOP {
            score[fr] = -10000.0;
        } else if nodes[i as usize].strand == 1 {
            if nodes[i as usize].cscore > score[fr] {
                score[fr] = nodes[i as usize].cscore;
            } else {
                nodes[i as usize].cscore -= score[fr] - nodes[i as usize].cscore;
            }
        }
    }

    for i in 0..3 {
        score[i] = -10000.0;
    }
    for i in (0..nn).rev() {
        let fr = (nodes[i as usize].ndx % 3) as usize;
        if nodes[i as usize].strand == -1 && nodes[i as usize].type_ == crate::sequence::STOP {
            score[fr] = -10000.0;
        } else if nodes[i as usize].strand == -1 {
            let _cscore_before_pass2 = nodes[i as usize].cscore;
            if nodes[i as usize].cscore > score[fr] {
                score[fr] = nodes[i as usize].cscore;
            } else {
                nodes[i as usize].cscore -= score[fr] - nodes[i as usize].cscore;
            }

            if nodes[i as usize].strand == -1
                && nodes[i as usize].stop_val >= 361240
                && nodes[i as usize].stop_val <= 361260
            {
            }
        }
    }

    for i in 0..nn {
        let fr = (nodes[i as usize].ndx % 3) as usize;
        if nodes[i as usize].strand == 1 && nodes[i as usize].type_ == crate::sequence::STOP {
            score[fr] = -10000.0;
        } else if nodes[i as usize].strand == 1 {

            let gsize = (((nodes[i as usize].stop_val - nodes[i as usize].ndx).abs() + 3) as f32) as f64 / 3.0;
            let mut lfac = if gsize > 1000.0 {
                let mut lf = ((1.0 - no_stop.powf(1000.0)) / no_stop.powf(1000.0)).ln();
                lf -= ((1.0 - no_stop.powf(80.0)) / no_stop.powf(80.0)).ln();
                lf * (gsize - 80.0) / 920.0
            } else {
                let mut lf = ((1.0 - no_stop.powf(gsize)) / no_stop.powf(gsize)).ln();
                lf -= ((1.0 - no_stop.powf(80.0)) / no_stop.powf(80.0)).ln();
                lf
            };

            if lfac > score[fr] {
                score[fr] = lfac;
            } else {
                lfac -= dmax(dmin(score[fr] - lfac, lfac), 0.0);
            }

            if lfac > 3.0 && nodes[i as usize].cscore < 0.5 * lfac {
                nodes[i as usize].cscore = 0.5 * lfac;
            }
            nodes[i as usize].cscore += lfac;
        }
    }

    for i in (0..nn).rev() {
        let fr = (nodes[i as usize].ndx % 3) as usize;
        if nodes[i as usize].strand == -1 && nodes[i as usize].type_ == crate::sequence::STOP {
            score[fr] = -10000.0;
        } else if nodes[i as usize].strand == -1 {

            let gsize = (((nodes[i as usize].stop_val - nodes[i as usize].ndx).abs() + 3) as f32) as f64 / 3.0;
            let _cscore_before = nodes[i as usize].cscore;
            let mut lfac = if gsize > 1000.0 {
                let mut lf = ((1.0 - no_stop.powf(1000.0)) / no_stop.powf(1000.0)).ln();
                lf -= ((1.0 - no_stop.powf(80.0)) / no_stop.powf(80.0)).ln();
                lf * (gsize - 80.0) / 920.0
            } else {
                let mut lf = ((1.0 - no_stop.powf(gsize)) / no_stop.powf(gsize)).ln();
                lf -= ((1.0 - no_stop.powf(80.0)) / no_stop.powf(80.0)).ln();
                lf
            };
            let _lfac_initial = lfac;

            if lfac > score[fr] {
                score[fr] = lfac;
            } else {
                lfac -= dmax(dmin(score[fr] - lfac, lfac), 0.0);
            }

            if lfac > 3.0 && nodes[i as usize].cscore < 0.5 * lfac {
                nodes[i as usize].cscore = 0.5 * lfac;
            }
            nodes[i as usize].cscore += lfac;

            if nodes[i as usize].strand == -1
                && nodes[i as usize].stop_val >= 361240
                && nodes[i as usize].stop_val <= 361260
                && nodes[i as usize].type_ != crate::sequence::STOP
            {
            }
        }
    }
}

pub fn calc_orf_gc(
    seq: &[u8],
    _rseq: &[u8],
    _slen: i32,
    nodes: &mut [Node],
    nn: i32,
    _tinf: &Training,
) {
    use crate::sequence::{is_gc, STOP};

    let mut last = [0i32; 3];
    let mut gc = [0.0f64; 3];

    for i in 0..3 {
        gc[i] = 0.0;
    }
    for i in (0..nn).rev() {
        let fr = (nodes[i as usize].ndx % 3) as usize;
        if nodes[i as usize].strand == 1 && nodes[i as usize].type_ == STOP {
            last[fr] = nodes[i as usize].ndx;
            gc[fr] = (is_gc(seq, nodes[i as usize].ndx) as u8) as f64
                + (is_gc(seq, nodes[i as usize].ndx + 1) as u8) as f64
                + (is_gc(seq, nodes[i as usize].ndx + 2) as u8) as f64;
        } else if nodes[i as usize].strand == 1 {
            let mut j = last[fr] - 3;
            while j >= nodes[i as usize].ndx {
                gc[fr] += (is_gc(seq, j) as u8) as f64
                    + (is_gc(seq, j + 1) as u8) as f64
                    + (is_gc(seq, j + 2) as u8) as f64;
                j -= 3;
            }
            let gsize = ((nodes[i as usize].stop_val - nodes[i as usize].ndx).abs() + 3) as f32;
            nodes[i as usize].gc_cont = gc[fr] / (gsize as f64);
            last[fr] = nodes[i as usize].ndx;
        }
    }

    for i in 0..3 {
        gc[i] = 0.0;
    }
    for i in 0..nn {
        let fr = (nodes[i as usize].ndx % 3) as usize;
        if nodes[i as usize].strand == -1 && nodes[i as usize].type_ == STOP {
            last[fr] = nodes[i as usize].ndx;
            gc[fr] = (is_gc(seq, nodes[i as usize].ndx) as u8) as f64
                + (is_gc(seq, nodes[i as usize].ndx - 1) as u8) as f64
                + (is_gc(seq, nodes[i as usize].ndx - 2) as u8) as f64;
        } else if nodes[i as usize].strand == -1 {
            let mut j = last[fr] + 3;
            while j <= nodes[i as usize].ndx {
                gc[fr] += (is_gc(seq, j) as u8) as f64
                    + (is_gc(seq, j + 1) as u8) as f64
                    + (is_gc(seq, j + 2) as u8) as f64;
                j += 3;
            }
            let gsize = ((nodes[i as usize].stop_val - nodes[i as usize].ndx).abs() + 3) as f32;
            nodes[i as usize].gc_cont = gc[fr] / (gsize as f64);
            last[fr] = nodes[i as usize].ndx;
        }
    }
}

pub fn rbs_score(
    seq: &[u8],
    rseq: &[u8],
    slen: i32,
    nodes: &mut [Node],
    nn: i32,
    tinf: &Training,
) {
    use crate::sequence::{shine_dalgarno_exact, shine_dalgarno_mm, STOP};
    use rayon::prelude::*;

    nodes[..nn as usize].par_iter_mut().for_each(|node| {
        if node.type_ == STOP || node.edge == 1 {
            return;
        }
        node.rbs[0] = 0;
        node.rbs[1] = 0;

        if node.strand == 1 {

            for j in (node.ndx - 20)..=(node.ndx - 6) {
                if j < 0 {
                    continue;
                }
                let cur_sc_exact = shine_dalgarno_exact(seq, j, node.ndx, &tinf.rbs_wt);
                let cur_sc_mm = shine_dalgarno_mm(seq, j, node.ndx, &tinf.rbs_wt);

                if cur_sc_exact > node.rbs[0] {
                    node.rbs[0] = cur_sc_exact;
                }
                if cur_sc_mm > node.rbs[1] {
                    node.rbs[1] = cur_sc_mm;
                }
            }
        } else if node.strand == -1 {

            for j in (slen - node.ndx - 21)..=(slen - node.ndx - 7) {
                if j > slen - 1 {
                    continue;
                }
                let cur_sc_exact =
                    shine_dalgarno_exact(rseq, j, slen - 1 - node.ndx, &tinf.rbs_wt);
                let cur_sc_mm = shine_dalgarno_mm(rseq, j, slen - 1 - node.ndx, &tinf.rbs_wt);

                if cur_sc_exact > node.rbs[0] {
                    node.rbs[0] = cur_sc_exact;
                }
                if cur_sc_mm > node.rbs[1] {
                    node.rbs[1] = cur_sc_mm;
                }
            }
        }
    });
}

pub fn score_upstream_composition(
    seq: &[u8],
    slen: i32,
    nod: &mut Node,
    tinf: &Training,
) {
    use crate::sequence::mer_ndx;

    let start = if nod.strand == 1 {
        nod.ndx
    } else {
        slen - 1 - nod.ndx
    };

    nod.uscore = 0.0;
    let mut count = 0;
    for i in 1..45 {
        if i > 2 && i < 15 {
            continue;
        }
        if start - i < 0 {
            continue;
        }
        nod.uscore += 0.4 * tinf.st_wt * tinf.ups_comp[count][mer_ndx(1, seq, start - i) as usize];
        count += 1;
    }
}

pub fn determine_sd_usage(tinf: &mut Training) {
    tinf.uses_sd = 1;
    if tinf.rbs_wt[0] >= 0.0 {
        tinf.uses_sd = 0;
    }
    if tinf.rbs_wt[16] < 1.0
        && tinf.rbs_wt[13] < 1.0
        && tinf.rbs_wt[15] < 1.0
        && (tinf.rbs_wt[0] >= -0.5
            || (tinf.rbs_wt[22] < 2.0 && tinf.rbs_wt[24] < 2.0 && tinf.rbs_wt[27] < 2.0))
    {
        tinf.uses_sd = 0;
    }
}

pub fn intergenic_mod(n1: &Node, n2: &Node, tinf: &Training) -> f64 {

    let s1 = n1.strand;
    let s2 = n2.strand;
    let x1 = n1.ndx;
    let x2 = n2.ndx;

    if s1 != s2 {
        return 0.0 - 0.15 * tinf.st_wt;
    }

    let mut rval = 0.0;

    if x1 + 2 == x2 || x1 - 1 == x2 {
        if s1 == 1 {
            if n2.rscore < 0.0 {
                rval -= n2.rscore;
            }
            if n2.uscore < 0.0 {
                rval -= n2.uscore;
            }
        } else {
            if n1.rscore < 0.0 {
                rval -= n1.rscore;
            }
            if n1.uscore < 0.0 {
                rval -= n1.uscore;
            }
        }
    }

    let ovlp = if s1 == 1 { x1 + 2 >= x2 } else { x1 >= x2 + 2 };

    let dist = (x1 - x2).abs();
    if dist > 3 * OPER_DIST {
        rval -= 0.15 * tinf.st_wt;
    } else if (dist <= OPER_DIST && !ovlp) || dist < (OPER_DIST / 4) {
        rval += (2.0 - dist as f64 / OPER_DIST as f64) * 0.15 * tinf.st_wt;
    }

    rval
}

pub fn train_starts_sd(
    seq: &[u8],
    rseq: &[u8],
    slen: i32,
    nodes: &mut [Node],
    nn: i32,
    tinf: &mut Training,
) {
    use crate::sequence::STOP;

    let wt = tinf.st_wt;
    let mut sthresh = 35.0;

    for j in 0..3 {
        tinf.type_wt[j] = 0.0;
    }
    for j in 0..28 {
        tinf.rbs_wt[j] = 0.0;
    }
    for i in 0..32 {
        for j in 0..4 {
            tinf.ups_comp[i][j] = 0.0;
        }
    }

    let mut tbg = [0.0f64; 3];
    for i in 0..nn as usize {
        if nodes[i].type_ == STOP {
            continue;
        }
        tbg[nodes[i].type_ as usize] += 1.0;
    }
    let mut sum = 0.0;
    for i in 0..3 {
        sum += tbg[i];
    }
    for i in 0..3 {
        tbg[i] /= sum;
    }

    for iter in 0..10 {

        let mut rbg = [0.0f64; 28];
        for j in 0..nn as usize {
            if nodes[j].type_ == STOP || nodes[j].edge == 1 {
                continue;
            }
            let max_rb = if tinf.rbs_wt[nodes[j].rbs[0] as usize]
                > tinf.rbs_wt[nodes[j].rbs[1] as usize] + 1.0
                || nodes[j].rbs[1] == 0
            {
                nodes[j].rbs[0]
            } else if tinf.rbs_wt[nodes[j].rbs[0] as usize]
                < tinf.rbs_wt[nodes[j].rbs[1] as usize] - 1.0
                || nodes[j].rbs[0] == 0
            {
                nodes[j].rbs[1]
            } else {
                dmax(nodes[j].rbs[0] as f64, nodes[j].rbs[1] as f64) as i32
            };
            rbg[max_rb as usize] += 1.0;
        }
        sum = 0.0;
        for j in 0..28 {
            sum += rbg[j];
        }
        for j in 0..28 {
            rbg[j] /= sum;
        }

        let mut rreal = [0.0f64; 28];
        let mut treal = [0.0f64; 3];

        let mut best = [0.0f64; 3];
        let mut bndx = [-1i32; 3];
        let mut rbs = [0i32; 3];
        let mut type_ = [0i32; 3];

        for j in 0..nn as usize {
            if nodes[j].type_ != STOP && nodes[j].edge == 1 {
                continue;
            }
            let fr = (nodes[j].ndx % 3) as usize;

            if nodes[j].type_ == STOP && nodes[j].strand == 1 {
                if best[fr] >= sthresh && bndx[fr] >= 0 && nodes[bndx[fr] as usize].ndx % 3 == fr as i32
                {
                    rreal[rbs[fr] as usize] += 1.0;
                    treal[type_[fr] as usize] += 1.0;
                    if iter == 9 {
                        count_upstream_composition(
                            seq,
                            slen,
                            1,
                            nodes[bndx[fr] as usize].ndx,
                            tinf,
                        );
                    }
                }
                best[fr] = 0.0;
                bndx[fr] = -1;
                rbs[fr] = 0;
                type_[fr] = 0;
            } else if nodes[j].strand == 1 {
                let max_rb = if tinf.rbs_wt[nodes[j].rbs[0] as usize]
                    > tinf.rbs_wt[nodes[j].rbs[1] as usize] + 1.0
                    || nodes[j].rbs[1] == 0
                {
                    nodes[j].rbs[0]
                } else if tinf.rbs_wt[nodes[j].rbs[0] as usize]
                    < tinf.rbs_wt[nodes[j].rbs[1] as usize] - 1.0
                    || nodes[j].rbs[0] == 0
                {
                    nodes[j].rbs[1]
                } else {
                    dmax(nodes[j].rbs[0] as f64, nodes[j].rbs[1] as f64) as i32
                };

                let score = nodes[j].cscore
                    + wt * tinf.rbs_wt[max_rb as usize]
                    + wt * tinf.type_wt[nodes[j].type_ as usize];

                if score >= best[fr] {
                    best[fr] = score;
                    bndx[fr] = j as i32;
                    type_[fr] = nodes[j].type_;
                    rbs[fr] = max_rb;
                }
            }
        }

        best = [0.0f64; 3];
        bndx = [-1i32; 3];
        rbs = [0i32; 3];
        type_ = [0i32; 3];

        for j in (0..nn as usize).rev() {
            if nodes[j].type_ != STOP && nodes[j].edge == 1 {
                continue;
            }
            let fr = (nodes[j].ndx % 3) as usize;

            if nodes[j].type_ == STOP && nodes[j].strand == -1 {
                if best[fr] >= sthresh && bndx[fr] >= 0 && nodes[bndx[fr] as usize].ndx % 3 == fr as i32
                {
                    rreal[rbs[fr] as usize] += 1.0;
                    treal[type_[fr] as usize] += 1.0;
                    if iter == 9 {
                        count_upstream_composition(
                            rseq,
                            slen,
                            -1,
                            nodes[bndx[fr] as usize].ndx,
                            tinf,
                        );
                    }
                }
                best[fr] = 0.0;
                bndx[fr] = -1;
                rbs[fr] = 0;
                type_[fr] = 0;
            } else if nodes[j].strand == -1 {
                let max_rb = if tinf.rbs_wt[nodes[j].rbs[0] as usize]
                    > tinf.rbs_wt[nodes[j].rbs[1] as usize] + 1.0
                    || nodes[j].rbs[1] == 0
                {
                    nodes[j].rbs[0]
                } else if tinf.rbs_wt[nodes[j].rbs[0] as usize]
                    < tinf.rbs_wt[nodes[j].rbs[1] as usize] - 1.0
                    || nodes[j].rbs[0] == 0
                {
                    nodes[j].rbs[1]
                } else {
                    dmax(nodes[j].rbs[0] as f64, nodes[j].rbs[1] as f64) as i32
                };

                let score = nodes[j].cscore
                    + wt * tinf.rbs_wt[max_rb as usize]
                    + wt * tinf.type_wt[nodes[j].type_ as usize];

                if score >= best[fr] {
                    best[fr] = score;
                    bndx[fr] = j as i32;
                    type_[fr] = nodes[j].type_;
                    rbs[fr] = max_rb;
                }
            }
        }

        sum = 0.0;
        for j in 0..28 {
            sum += rreal[j];
        }
        if sum == 0.0 {
            for j in 0..28 {
                tinf.rbs_wt[j] = 0.0;
            }
        } else {
            for j in 0..28 {
                rreal[j] /= sum;
                if rbg[j] != 0.0 {
                    tinf.rbs_wt[j] = (rreal[j] / rbg[j]).ln();
                } else {
                    tinf.rbs_wt[j] = -4.0;
                }
                if tinf.rbs_wt[j] > 4.0 {
                    tinf.rbs_wt[j] = 4.0;
                }
                if tinf.rbs_wt[j] < -4.0 {
                    tinf.rbs_wt[j] = -4.0;
                }
            }
        }

        sum = 0.0;
        for j in 0..3 {
            sum += treal[j];
        }
        if sum == 0.0 {
            for j in 0..3 {
                tinf.type_wt[j] = 0.0;
            }
        } else {
            for j in 0..3 {
                treal[j] /= sum;
                if tbg[j] != 0.0 {
                    tinf.type_wt[j] = (treal[j] / tbg[j]).ln();
                } else {
                    tinf.type_wt[j] = -4.0;
                }
                if tinf.type_wt[j] > 4.0 {
                    tinf.type_wt[j] = 4.0;
                }
                if tinf.type_wt[j] < -4.0 {
                    tinf.type_wt[j] = -4.0;
                }
            }
        }

        if sum <= (nn as f64) / 2000.0 {
            sthresh /= 2.0;
        }
    }

    for i in 0..32 {
        sum = 0.0;
        for j in 0..4 {
            sum += tinf.ups_comp[i][j];
        }
        if sum == 0.0 {
            for j in 0..4 {
                tinf.ups_comp[i][j] = 0.0;
            }
        } else {
            for j in 0..4 {
                tinf.ups_comp[i][j] /= sum;
                if tinf.gc > 0.1 && tinf.gc < 0.9 {
                    if j == 0 || j == 3 {
                        tinf.ups_comp[i][j] =
                            (tinf.ups_comp[i][j] * 2.0 / (1.0 - tinf.gc)).ln();
                    } else {
                        tinf.ups_comp[i][j] = (tinf.ups_comp[i][j] * 2.0 / tinf.gc).ln();
                    }
                } else if tinf.gc <= 0.1 {
                    if j == 0 || j == 3 {
                        tinf.ups_comp[i][j] = (tinf.ups_comp[i][j] * 2.0 / 0.90).ln();
                    } else {
                        tinf.ups_comp[i][j] = (tinf.ups_comp[i][j] * 2.0 / 0.10).ln();
                    }
                } else {
                    if j == 0 || j == 3 {
                        tinf.ups_comp[i][j] = (tinf.ups_comp[i][j] * 2.0 / 0.10).ln();
                    } else {
                        tinf.ups_comp[i][j] = (tinf.ups_comp[i][j] * 2.0 / 0.90).ln();
                    }
                }
                if tinf.ups_comp[i][j] > 4.0 {
                    tinf.ups_comp[i][j] = 4.0;
                }
                if tinf.ups_comp[i][j] < -4.0 {
                    tinf.ups_comp[i][j] = -4.0;
                }
            }
        }
    }
}

pub fn train_starts_nonsd(
    seq: &[u8],
    rseq: &[u8],
    slen: i32,
    nodes: &mut [Node],
    nn: i32,
    tinf: &mut Training,
) {
    use crate::sequence::STOP;

    let wt = tinf.st_wt;
    let mut sthresh = 35.0;

    for i in 0..32 {
        for j in 0..4 {
            tinf.ups_comp[i][j] = 0.0;
        }
    }

    for i in 0..3 {
        tinf.type_wt[i] = 0.0;
    }
    let mut tbg = [0.0f64; 3];
    for i in 0..nn as usize {
        if nodes[i].type_ == STOP {
            continue;
        }
        tbg[nodes[i].type_ as usize] += 1.0;
    }
    let mut sum = 0.0;
    for i in 0..3 {
        sum += tbg[i];
    }
    for i in 0..3 {
        tbg[i] /= sum;
    }

    let mut mgood = [[[0i32; 4096]; 4]; 4];

    for iter in 0..20 {

        let stage = if iter < 4 {
            0
        } else if iter < 12 {
            1
        } else {
            2
        };

        let mut mbg = [[[0.0f64; 4096]; 4]; 4];
        let mut zbg = 0.0f64;

        for j in 0..nn as usize {
            if nodes[j].type_ == STOP || nodes[j].edge == 1 {
                continue;
            }
            find_best_upstream_motif(tinf, seq, rseq, slen, &mut nodes[j], stage);
            update_motif_counts(&mut mbg, &mut zbg, seq, rseq, slen, &nodes[j], stage);
        }

        sum = 0.0;
        for j in 0..4 {
            for k in 0..4 {
                for l in 0..4096 {
                    sum += mbg[j][k][l];
                }
            }
        }
        sum += zbg;
        for j in 0..4 {
            for k in 0..4 {
                for l in 0..4096 {
                    mbg[j][k][l] /= sum;
                }
            }
        }
        zbg /= sum;

        let mut mreal = [[[0.0f64; 4096]; 4]; 4];
        let mut zreal = 0.0f64;
        let mut treal = [0.0f64; 3];
        let mut ngenes = 0.0f64;

        let mut best = [0.0f64; 3];
        let mut bndx = [-1i32; 3];

        for j in 0..nn as usize {
            if nodes[j].type_ != STOP && nodes[j].edge == 1 {
                continue;
            }
            let fr = (nodes[j].ndx % 3) as usize;

            if nodes[j].type_ == STOP && nodes[j].strand == 1 {
                if best[fr] >= sthresh && bndx[fr] >= 0 {
                    ngenes += 1.0;
                    treal[nodes[bndx[fr] as usize].type_ as usize] += 1.0;
                    update_motif_counts(
                        &mut mreal,
                        &mut zreal,
                        seq,
                        rseq,
                        slen,
                        &nodes[bndx[fr] as usize],
                        stage,
                    );
                    if iter == 19 {
                        count_upstream_composition(
                            seq,
                            slen,
                            1,
                            nodes[bndx[fr] as usize].ndx,
                            tinf,
                        );
                    }
                }
                best[fr] = 0.0;
                bndx[fr] = -1;
            } else if nodes[j].strand == 1 {
                let score = nodes[j].cscore + wt * nodes[j].mot.score + wt * tinf.type_wt[nodes[j].type_ as usize];

                if score >= best[fr] {
                    best[fr] = nodes[j].cscore + wt * nodes[j].mot.score;
                    best[fr] += wt * tinf.type_wt[nodes[j].type_ as usize];
                    bndx[fr] = j as i32;
                }
            }
        }

        best = [0.0f64; 3];
        bndx = [-1i32; 3];

        for j in (0..nn as usize).rev() {
            if nodes[j].type_ != STOP && nodes[j].edge == 1 {
                continue;
            }
            let fr = (nodes[j].ndx % 3) as usize;

            if nodes[j].type_ == STOP && nodes[j].strand == -1 {
                if best[fr] >= sthresh && bndx[fr] >= 0 {
                    ngenes += 1.0;
                    treal[nodes[bndx[fr] as usize].type_ as usize] += 1.0;
                    update_motif_counts(
                        &mut mreal,
                        &mut zreal,
                        seq,
                        rseq,
                        slen,
                        &nodes[bndx[fr] as usize],
                        stage,
                    );
                    if iter == 19 {
                        count_upstream_composition(
                            rseq,
                            slen,
                            -1,
                            nodes[bndx[fr] as usize].ndx,
                            tinf,
                        );
                    }
                }
                best[fr] = 0.0;
                bndx[fr] = -1;
            } else if nodes[j].strand == -1 {
                let score = nodes[j].cscore + wt * nodes[j].mot.score + wt * tinf.type_wt[nodes[j].type_ as usize];

                if score >= best[fr] {
                    best[fr] = nodes[j].cscore + wt * nodes[j].mot.score;
                    best[fr] += wt * tinf.type_wt[nodes[j].type_ as usize];
                    bndx[fr] = j as i32;
                }
            }
        }

        if stage < 2 {
            build_coverage_map(&mut mreal, &mut mgood, ngenes, stage);
        }

        sum = 0.0;
        for j in 0..4 {
            for k in 0..4 {
                for l in 0..4096 {
                    sum += mreal[j][k][l];
                }
            }
        }
        sum += zreal;

        if sum == 0.0 {
            for j in 0..4 {
                for k in 0..4 {
                    for l in 0..4096 {
                        tinf.mot_wt[j][k][l] = 0.0;
                    }
                }
            }
            tinf.no_mot = 0.0;
        } else {
            for j in 0..4 {
                for k in 0..4 {
                    for l in 0..4096 {
                        if mgood[j][k][l] == 0 {
                            zreal += mreal[j][k][l];
                            zbg += mreal[j][k][l];
                            mreal[j][k][l] = 0.0;
                            mbg[j][k][l] = 0.0;
                        }
                        mreal[j][k][l] /= sum;
                        if mbg[j][k][l] != 0.0 {
                            tinf.mot_wt[j][k][l] = (mreal[j][k][l] / mbg[j][k][l]).ln();
                        } else {
                            tinf.mot_wt[j][k][l] = -4.0;
                        }
                        if tinf.mot_wt[j][k][l] > 4.0 {
                            tinf.mot_wt[j][k][l] = 4.0;
                        }
                        if tinf.mot_wt[j][k][l] < -4.0 {
                            tinf.mot_wt[j][k][l] = -4.0;
                        }
                    }
                }
            }
        }

        zreal /= sum;
        if zbg != 0.0 {
            tinf.no_mot = (zreal / zbg).ln();
        } else {
            tinf.no_mot = -4.0;
        }
        if tinf.no_mot > 4.0 {
            tinf.no_mot = 4.0;
        }
        if tinf.no_mot < -4.0 {
            tinf.no_mot = -4.0;
        }

        sum = 0.0;
        for j in 0..3 {
            sum += treal[j];
        }
        if sum == 0.0 {
            for j in 0..3 {
                tinf.type_wt[j] = 0.0;
            }
        } else {
            for j in 0..3 {
                treal[j] /= sum;
                if tbg[j] != 0.0 {
                    tinf.type_wt[j] = (treal[j] / tbg[j]).ln();
                } else {
                    tinf.type_wt[j] = -4.0;
                }
                if tinf.type_wt[j] > 4.0 {
                    tinf.type_wt[j] = 4.0;
                }
                if tinf.type_wt[j] < -4.0 {
                    tinf.type_wt[j] = -4.0;
                }
            }
        }

        if sum <= (nn as f64) / 2000.0 {
            sthresh /= 2.0;
        }
    }

    for i in 0..32 {
        sum = 0.0;
        for j in 0..4 {
            sum += tinf.ups_comp[i][j];
        }
        if sum == 0.0 {
            for j in 0..4 {
                tinf.ups_comp[i][j] = 0.0;
            }
        } else {
            for j in 0..4 {
                tinf.ups_comp[i][j] /= sum;
                if tinf.gc > 0.1 && tinf.gc < 0.9 {
                    if j == 0 || j == 3 {
                        tinf.ups_comp[i][j] =
                            (tinf.ups_comp[i][j] * 2.0 / (1.0 - tinf.gc)).ln();
                    } else {
                        tinf.ups_comp[i][j] = (tinf.ups_comp[i][j] * 2.0 / tinf.gc).ln();
                    }
                } else if tinf.gc <= 0.1 {
                    if j == 0 || j == 3 {
                        tinf.ups_comp[i][j] = (tinf.ups_comp[i][j] * 2.0 / 0.90).ln();
                    } else {
                        tinf.ups_comp[i][j] = (tinf.ups_comp[i][j] * 2.0 / 0.10).ln();
                    }
                } else {
                    if j == 0 || j == 3 {
                        tinf.ups_comp[i][j] = (tinf.ups_comp[i][j] * 2.0 / 0.10).ln();
                    } else {
                        tinf.ups_comp[i][j] = (tinf.ups_comp[i][j] * 2.0 / 0.90).ln();
                    }
                }
                if tinf.ups_comp[i][j] > 4.0 {
                    tinf.ups_comp[i][j] = 4.0;
                }
                if tinf.ups_comp[i][j] < -4.0 {
                    tinf.ups_comp[i][j] = -4.0;
                }
            }
        }
    }
}

pub fn count_upstream_composition(
    seq: &[u8],
    slen: i32,
    strand: i32,
    pos: i32,
    tinf: &mut Training,
) {
    use crate::sequence::mer_ndx;

    let start = if strand == 1 { pos } else { slen - 1 - pos };

    let mut count = 0;
    for i in 1..45 {
        if i > 2 && i < 15 {
            continue;
        }
        if start - i >= 0 {
            tinf.ups_comp[count][mer_ndx(1, seq, start - i) as usize] += 1.0;
        }
        count += 1;
    }
}

pub fn build_coverage_map(
    cvg: &mut [[[f64; 4096]; 4]; 4],
    cnts: &mut [[[i32; 4096]; 4]; 4],
    gc: f64,
    _stage: i32,
) {
    let thresh = 0.2;

    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4096 {
                cnts[i][j][k] = 0;
            }
        }
    }

    for i in 0..4 {
        for j in 0..64 {
            if cvg[0][i][j] / gc >= thresh {
                for k in 0..4 {
                    cnts[0][k][j] = 1;
                }
            }
        }
    }

    for i in 0..4 {
        for j in 0..256 {
            let decomp0 = (j & 252) >> 2;
            let decomp1 = j & 63;
            if cnts[0][i][decomp0] == 0 || cnts[0][i][decomp1] == 0 {
                continue;
            }
            cnts[1][i][j] = 1;
        }
    }

    for i in 0..4 {
        for j in 0..1024 {
            let decomp0 = (j & 1008) >> 4;
            let decomp1 = (j & 252) >> 2;
            let decomp2 = j & 63;
            if cnts[0][i][decomp0] == 0 || cnts[0][i][decomp1] == 0 || cnts[0][i][decomp2] == 0 {
                continue;
            }
            cnts[2][i][j] = 1;
            let mut tmp = j;
            for k in (0..=16).step_by(16) {
                tmp = tmp ^ k;
                for l in (0..=32).step_by(32) {
                    tmp = tmp ^ l;
                    if cnts[2][i][tmp] == 0 {
                        cnts[2][i][tmp] = 2;
                    }
                }
            }
        }
    }

    for i in 0..4 {
        for j in 0..4096 {
            let decomp0 = (j & 4092) >> 2;
            let decomp1 = j & 1023;
            if cnts[2][i][decomp0] == 0 || cnts[2][i][decomp1] == 0 {
                continue;
            }
            if cnts[2][i][decomp0] == 1 && cnts[2][i][decomp1] == 1 {
                cnts[3][i][j] = 1;
            } else {
                cnts[3][i][j] = 2;
            }
        }
    }
}

pub fn find_best_upstream_motif(
    tinf: &Training,
    seq: &[u8],
    rseq: &[u8],
    slen: i32,
    nod: &mut Node,
    stage: i32,
) {
    use crate::sequence::{mer_ndx, STOP};

    if nod.type_ == STOP || nod.edge == 1 {
        return;
    }

    let (wseq, start) = if nod.strand == 1 {
        (seq, nod.ndx)
    } else {
        (rseq, slen - 1 - nod.ndx)
    };

    let mut max_sc = -100.0;
    let mut max_spacer = 0i32;
    let mut max_spacendx = 0i32;
    let mut max_len = 0i32;
    let mut max_ndx = 0i32;

    for i in (0..=3).rev() {
        for j in (start - 18 - i)..=(start - 6 - i) {
            if j < 0 {
                continue;
            }
            let spacer = start - j - i - 3;
            let spacendx: i32 = if j <= start - 16 - i {
                3
            } else if j <= start - 14 - i {
                2
            } else if j >= start - 7 - i {
                1
            } else {
                0
            };
            let index = mer_ndx(i + 3, wseq, j);
            let score = tinf.mot_wt[i as usize][spacendx as usize][index as usize];
            if score > max_sc {
                max_sc = score;
                max_spacendx = spacendx;
                max_spacer = spacer;
                max_ndx = index;
                max_len = i + 3;
            }
        }
    }

    if stage == 2 && (max_sc == -4.0 || max_sc < tinf.no_mot + 0.69) {
        nod.mot.ndx = 0;
        nod.mot.len = 0;
        nod.mot.spacendx = 0;
        nod.mot.spacer = 0;
        nod.mot.score = tinf.no_mot;
    } else {
        nod.mot.ndx = max_ndx;
        nod.mot.len = max_len;
        nod.mot.spacendx = max_spacendx;
        nod.mot.spacer = max_spacer;
        nod.mot.score = max_sc;
    }
}

pub fn update_motif_counts(
    cvg: &mut [[[f64; 4096]; 4]; 4],
    no_mot: &mut f64,
    seq: &[u8],
    rseq: &[u8],
    slen: i32,
    nod: &Node,
    stage: i32,
) {
    use crate::sequence::{mer_ndx, STOP};

    if nod.type_ == STOP || nod.edge == 1 {
        return;
    }
    if nod.mot.len == 0 {
        *no_mot += 1.0;
        return;
    }

    let (wseq, start) = if nod.strand == 1 {
        (seq, nod.ndx)
    } else {
        (rseq, slen - 1 - nod.ndx)
    };

    if stage == 0 {
        for i in (0..=3).rev() {
            for j in (start - 18 - i)..=(start - 6 - i) {
                if j < 0 {
                    continue;
                }
                let _spacendx = if j <= start - 16 - i {
                    3
                } else if j <= start - 14 - i {
                    2
                } else if j >= start - 7 - i {
                    1
                } else {
                    0
                };
                for k in 0..4 {
                    cvg[i as usize][k][mer_ndx(i + 3, wseq, j) as usize] += 1.0;
                }
            }
        }
    }

    else if stage == 1 {
        cvg[(nod.mot.len - 3) as usize][nod.mot.spacendx as usize][nod.mot.ndx as usize] += 1.0;
        for i in 0..(nod.mot.len - 3) {
            for j in (start - nod.mot.spacer - nod.mot.len)..=(start - nod.mot.spacer - (i + 3)) {
                if j < 0 {
                    continue;
                }
                let spacendx = if j <= start - 16 - i {
                    3
                } else if j <= start - 14 - i {
                    2
                } else if j >= start - 7 - i {
                    1
                } else {
                    0
                };
                cvg[i as usize][spacendx][mer_ndx(i + 3, wseq, j) as usize] += 1.0;
            }
        }
    }

    else if stage == 2 {
        cvg[(nod.mot.len - 3) as usize][nod.mot.spacendx as usize][nod.mot.ndx as usize] += 1.0;
    }
}

pub fn write_start_file(
    start_ptr: &mut dyn std::io::Write,
    nodes: &[Node],
    nn: i32,
    tinf: &Training,
    num_seq: i32,
    slen: i32,
    is_meta: i32,
    header: &str,
    short_header: &str,
    vers: &str,
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

    let seq_data = format!("seqnum={};seqlen={};seqhdr=\"{}\"", num_seq, slen, header);

    let mut run_data = if is_meta == 0 {
        format!(
            "version=Rustygal.v{};run_type=Single;model=\"Ab initio\";",
            vers
        )
    } else {
        format!(
            "version=Rustygal.v{};run_type=Metagenomic;model=\"{}\";",
            vers, short_header
        )
    };

    run_data.push_str(&format!(
        "gc_cont={:.2};transl_table={};uses_sd={}",
        tinf.gc * 100.0,
        tinf.trans_table,
        tinf.uses_sd
    ));

    let mut sorted_nodes: Vec<Node> = nodes[..nn as usize].to_vec();

    sorted_nodes.sort_by(stopcmp_nodes);

    let _ = writeln!(start_ptr, "# Sequence Data: {}", seq_data);
    let _ = writeln!(start_ptr, "# Run Data: {}\n", run_data);

    let _ = write!(
        start_ptr,
        "Beg\tEnd\tStd\tTotal\tCodPot\tStrtSc\tCodon\tRBSMot\t"
    );
    let _ = writeln!(
        start_ptr,
        "Spacer\tRBSScr\tUpsScr\tTypeScr\tGCCont"
    );

    let mut prev_stop = -1i32;
    let mut prev_strand = 0i32;

    for node in &sorted_nodes {
        use crate::sequence::STOP;

        if node.type_ == STOP {
            continue;
        }

        let st_type = if node.edge == 1 { 3 } else { node.type_ as usize };

        if node.stop_val != prev_stop || node.strand != prev_strand {
            prev_stop = node.stop_val;
            prev_strand = node.strand;
            let _ = writeln!(start_ptr);
        }

        if node.strand == 1 {
            let _ = write!(
                start_ptr,
                "{}\t{}\t+\t{:.2}\t{:.2}\t{:.2}\t{}\t",
                node.ndx + 1,
                node.stop_val + 3,
                node.cscore + node.sscore,
                node.cscore,
                node.sscore,
                type_string[st_type]
            );
        } else {
            let _ = write!(
                start_ptr,
                "{}\t{}\t-\t{:.2}\t{:.2}\t{:.2}\t{}\t",
                node.stop_val - 1,
                node.ndx + 1,
                node.cscore + node.sscore,
                node.cscore,
                node.sscore,
                type_string[st_type]
            );
        }

        let rbs1 = tinf.rbs_wt[node.rbs[0] as usize] * tinf.st_wt;
        let rbs2 = tinf.rbs_wt[node.rbs[1] as usize] * tinf.st_wt;

        if tinf.uses_sd == 1 {

            if rbs1 > rbs2 {
                let _ = write!(
                    start_ptr,
                    "{}\t{}\t{:.2}\t",
                    sd_string[node.rbs[0] as usize],
                    sd_spacer[node.rbs[0] as usize],
                    node.rscore
                );
            } else {
                let _ = write!(
                    start_ptr,
                    "{}\t{}\t{:.2}\t",
                    sd_string[node.rbs[1] as usize],
                    sd_spacer[node.rbs[1] as usize],
                    node.rscore
                );
            }
        } else {

            let mut qt = vec![0u8; 10];
            mer_text(&mut qt, node.mot.len, node.mot.ndx);

            if tinf.no_mot > -0.5
                && rbs1 > rbs2
                && rbs1 > node.mot.score * tinf.st_wt
            {
                let _ = write!(
                    start_ptr,
                    "{}\t{}\t{:.2}\t",
                    sd_string[node.rbs[0] as usize],
                    sd_spacer[node.rbs[0] as usize],
                    node.rscore
                );
            } else if tinf.no_mot > -0.5
                && rbs2 >= rbs1
                && rbs2 > node.mot.score * tinf.st_wt
            {
                let _ = write!(
                    start_ptr,
                    "{}\t{}\t{:.2}\t",
                    sd_string[node.rbs[1] as usize],
                    sd_spacer[node.rbs[1] as usize],
                    node.rscore
                );
            } else {
                if node.mot.len == 0 {
                    let _ = write!(start_ptr, "None\tNone\t{:.2}\t", node.rscore);
                } else {
                    let motif_str = std::str::from_utf8(&qt[..node.mot.len as usize])
                        .unwrap_or("None");
                    let _ = write!(
                        start_ptr,
                        "{}\t{}bp\t{:.2}\t",
                        motif_str, node.mot.spacer, node.rscore
                    );
                }
            }
        }

        let _ = writeln!(
            start_ptr,
            "{:.2}\t{:.2}\t{:.3}",
            node.uscore, node.tscore, node.gc_cont
        );
    }

    let _ = writeln!(start_ptr);

}

pub fn cross_mask(left: i32, right: i32, mlist: &[Mask], nmask: i32) -> bool {
    for i in 0..nmask as usize {
        if i >= mlist.len() {
            break;
        }

        if right < mlist[i].begin || left > mlist[i].end {
            continue;
        }

        return true;
    }
    false
}

pub fn dmax(a: f64, b: f64) -> f64 {
    if a > b {
        a
    } else {
        b
    }
}

pub fn dmin(a: f64, b: f64) -> f64 {
    if a < b {
        a
    } else {
        b
    }
}

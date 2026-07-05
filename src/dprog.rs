// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2007-2016 Doug Hyatt, Univ. of Tennessee / UT-Battelle (original Prodigal C)
// Copyright (C) 2026 Sunju Kim (Rust reimplementation)

use crate::node::Node;
use crate::training::Training;

pub const MAX_OPP_OVLP: i32 = 200;
pub const MAX_NODE_DIST: i32 = 500;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum ConnectionKind {
    ForwardStart = 0,
    ForwardStop = 1,
    BackwardStart = 2,
    BackwardStop = 3,
}

impl ConnectionKind {

    #[inline(always)]
    pub fn from_node(strand: i32, type_: i32, stop_const: i32) -> Self {
        let is_backward = (strand != 1) as u8;
        let is_stop = (type_ == stop_const) as u8;
        let kind = 2 * is_backward + is_stop;

        match kind {
            0 => ConnectionKind::ForwardStart,
            1 => ConnectionKind::ForwardStop,
            2 => ConnectionKind::BackwardStart,
            _ => ConnectionKind::BackwardStop,
        }
    }
}

#[inline]
pub fn dprog_auto(nodes: &mut [Node], nn: i32, tinf: &Training, flag: i32) -> i32 {
    dprog(nodes, nn, tinf, flag)
}

#[inline(always)]
fn skip_one(n1s: i8, n1_stop: bool, n1f: u8, n2s: i8, n2_stop: bool, n2f: u8) -> u8 {
    let n1_fwd = n1s == 1;
    let n2_fwd = n2s == 1;
    let s = (!n1_stop && !n2_stop && n1s == n2s)
        || (n1_fwd && !n1_stop && !n2_fwd)
        || (!n1_fwd && n1_stop && n2_fwd)
        || (!n1_fwd && !n1_stop && n2_fwd && n2_stop)
        || (n1s == n2s && n1_fwd && !n1_stop && n2_stop && n1f != n2f)
        || (n1s == n2s && !n1_fwd && n1_stop && !n2_stop && n1f != n2f);
    s as u8
}

fn compute_skippable(
    strand: &[i8],
    stop: &[u8],
    frame: &[u8],
    min: i32,
    i: i32,
    skip: &mut [u8],
) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx2") {

            unsafe { compute_skippable_avx2(strand, stop, frame, min, i, skip) };
            return;
        }
    }
    compute_skippable_scalar(strand, stop, frame, min, i, skip);
}

#[inline]
fn compute_skippable_scalar(
    strand: &[i8],
    stop: &[u8],
    frame: &[u8],
    min: i32,
    i: i32,
    skip: &mut [u8],
) {
    let iu = i as usize;
    let n2s = strand[iu];
    let n2_stop = stop[iu] != 0;
    let n2f = frame[iu];
    for j in min as usize..iu {
        skip[j] = skip_one(strand[j], stop[j] != 0, frame[j], n2s, n2_stop, n2f);
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn compute_skippable_avx2(
    strand: &[i8],
    stop: &[u8],
    frame: &[u8],
    min: i32,
    i: i32,
    skip: &mut [u8],
) {
    use std::arch::x86_64::*;

    let iu = i as usize;
    let n2s = strand[iu];
    let n2_stopb = stop[iu] != 0;
    let n2f = frame[iu];

    let all1 = _mm256_set1_epi8(-1);
    let ones = _mm256_set1_epi8(1);
    let fwd_splat = _mm256_set1_epi8(1);
    let stop_splat = _mm256_set1_epi8(1);
    let n2s_splat = _mm256_set1_epi8(n2s);
    let n2f_splat = _mm256_set1_epi8(n2f as i8);

    let n2_fwd = if n2s == 1 { all1 } else { _mm256_setzero_si256() };
    let n2_stop = if n2_stopb { all1 } else { _mm256_setzero_si256() };

    let strand_p = strand.as_ptr();
    let stop_p = stop.as_ptr();
    let frame_p = frame.as_ptr();
    let skip_p = skip.as_mut_ptr();

    let mut j = min as usize;
    while j + 32 <= iu {
        let s_j = _mm256_loadu_si256(strand_p.add(j) as *const __m256i);
        let st_j = _mm256_loadu_si256(stop_p.add(j) as *const __m256i);
        let f_j = _mm256_loadu_si256(frame_p.add(j) as *const __m256i);

        let n1_fwd = _mm256_cmpeq_epi8(s_j, fwd_splat);
        let n1_stop = _mm256_cmpeq_epi8(st_j, stop_splat);
        let samestrand = _mm256_cmpeq_epi8(s_j, n2s_splat);
        let framematch = _mm256_cmpeq_epi8(f_j, n2f_splat);
        let framediff = _mm256_xor_si256(framematch, all1);

        let t1 = _mm256_andnot_si256(n1_stop, _mm256_andnot_si256(n2_stop, samestrand));

        let t2 = _mm256_andnot_si256(n2_fwd, _mm256_andnot_si256(n1_stop, n1_fwd));

        let t3 = _mm256_and_si256(_mm256_andnot_si256(n1_fwd, n1_stop), n2_fwd);

        let t4 = _mm256_andnot_si256(
            n1_fwd,
            _mm256_andnot_si256(n1_stop, _mm256_and_si256(n2_fwd, n2_stop)),
        );

        let t5 = _mm256_and_si256(
            samestrand,
            _mm256_andnot_si256(
                n1_stop,
                _mm256_and_si256(n1_fwd, _mm256_and_si256(n2_stop, framediff)),
            ),
        );

        let t6 = _mm256_and_si256(
            samestrand,
            _mm256_andnot_si256(
                n1_fwd,
                _mm256_and_si256(n1_stop, _mm256_andnot_si256(n2_stop, framediff)),
            ),
        );

        let mut s = _mm256_or_si256(t1, t2);
        s = _mm256_or_si256(s, t3);
        s = _mm256_or_si256(s, t4);
        s = _mm256_or_si256(s, t5);
        s = _mm256_or_si256(s, t6);

        s = _mm256_and_si256(s, ones);
        _mm256_storeu_si256(skip_p.add(j) as *mut __m256i, s);
        j += 32;
    }

    while j < iu {
        *skip_p.add(j) = skip_one(*strand_p.add(j), *stop_p.add(j) != 0, *frame_p.add(j), n2s, n2_stopb, n2f);
        j += 1;
    }
}

pub fn dprog(nodes: &mut [Node], nn: i32, tinf: &Training, flag: i32) -> i32 {
    use crate::sequence::STOP;

    if nn == 0 {
        return -1;
    }

    for i in 0..nn as usize {
        nodes[i].score = 0.0;
        nodes[i].traceb = -1;
        nodes[i].tracef = -1;
    }

    let nnu = nn as usize;
    let mut soa_strand = vec![0i8; nnu];
    let mut soa_stop = vec![0u8; nnu];
    let mut soa_frame = vec![0u8; nnu];
    for k in 0..nnu {
        soa_strand[k] = nodes[k].strand as i8;
        soa_stop[k] = (nodes[k].type_ == STOP) as u8;
        soa_frame[k] = (nodes[k].ndx % 3) as u8;
    }
    let mut skip = vec![0u8; nnu];

    for i in 0..nn {

        let mut min = if i < MAX_NODE_DIST {
            0
        } else {
            i - MAX_NODE_DIST
        };

        if nodes[i as usize].strand == -1
            && nodes[i as usize].type_ != STOP
            && nodes[min as usize].ndx >= nodes[i as usize].stop_val
        {
            let _min_before = min;
            while min >= 0 && nodes[i as usize].ndx != nodes[i as usize].stop_val {
                min -= 1;
            }
        }
        if nodes[i as usize].strand == 1
            && nodes[i as usize].type_ == STOP
            && nodes[min as usize].ndx >= nodes[i as usize].stop_val
        {
            let _min_before = min;
            while min >= 0 && nodes[i as usize].ndx != nodes[i as usize].stop_val {
                min -= 1;
            }
        }
        if min < MAX_NODE_DIST {
            min = 0;
        } else {
            min = min - MAX_NODE_DIST;
        }

        let kind = ConnectionKind::from_node(
            nodes[i as usize].strand,
            nodes[i as usize].type_,
            STOP,
        );

        compute_skippable(&soa_strand, &soa_stop, &soa_frame, min, i, &mut skip);

        for j in min..i {

            if skip[j as usize] != 0 {
                continue;
            }

            match kind {
                ConnectionKind::ForwardStart => {
                    score_connection_forward_start(nodes, j, i, tinf, flag);
                }
                ConnectionKind::ForwardStop => {
                    score_connection_forward_stop(nodes, j, i, tinf, flag);
                }
                ConnectionKind::BackwardStart => {
                    score_connection_backward_start(nodes, j, i, tinf, flag);
                }
                ConnectionKind::BackwardStop => {
                    score_connection_backward_stop(nodes, j, i, tinf, flag);
                }
            }
        }
    }

    let mut max_sc = -1.0f64;
    let mut max_ndx = -1i32;

    for i in (0..nn).rev() {
        if nodes[i as usize].strand == 1 && nodes[i as usize].type_ != STOP {
            continue;
        }
        if nodes[i as usize].strand == -1 && nodes[i as usize].type_ == STOP {
            continue;
        }
        if nodes[i as usize].score > max_sc {
            max_sc = nodes[i as usize].score;
            max_ndx = i;
        }
    }

    let mut path = max_ndx;
    while path != -1 && nodes[path as usize].traceb != -1 {
        let nxt = nodes[path as usize].traceb;
        if nodes[path as usize].strand == -1
            && nodes[path as usize].type_ == STOP
            && nodes[nxt as usize].strand == 1
            && nodes[nxt as usize].type_ == STOP
            && nodes[path as usize].ov_mark != -1
            && nodes[path as usize].ndx > nodes[nxt as usize].ndx
        {
            let tmp = nodes[path as usize].star_ptr[nodes[path as usize].ov_mark as usize];
            let mut i = tmp;
            while nodes[i as usize].ndx != nodes[tmp as usize].stop_val {
                i -= 1;
            }
            nodes[path as usize].traceb = tmp;
            nodes[tmp as usize].traceb = i;
            nodes[i as usize].ov_mark = -1;
            nodes[i as usize].traceb = nxt;
        }
        path = nodes[path as usize].traceb;
    }

    path = max_ndx;
    while path != -1 && nodes[path as usize].traceb != -1 {
        let nxt = nodes[path as usize].traceb;

        if nodes[path as usize].strand == -1
            && nodes[path as usize].type_ != STOP
            && nodes[nxt as usize].strand == 1
            && nodes[nxt as usize].type_ == STOP
        {
            let mut i = path;
            while nodes[i as usize].ndx != nodes[path as usize].stop_val {
                i -= 1;
            }
            nodes[path as usize].traceb = i;
            nodes[i as usize].traceb = nxt;
        }
        if nodes[path as usize].strand == 1
            && nodes[path as usize].type_ == STOP
            && nodes[nxt as usize].strand == 1
            && nodes[nxt as usize].type_ == STOP
        {
            nodes[path as usize].traceb =
                nodes[nxt as usize].star_ptr[(nodes[path as usize].ndx % 3) as usize];
            let new_traceb = nodes[path as usize].traceb;
            nodes[new_traceb as usize].traceb = nxt;
        }
        if nodes[path as usize].strand == -1
            && nodes[path as usize].type_ == STOP
            && nodes[nxt as usize].strand == -1
            && nodes[nxt as usize].type_ == STOP
        {
            nodes[path as usize].traceb =
                nodes[path as usize].star_ptr[(nodes[nxt as usize].ndx % 3) as usize];
            let new_traceb = nodes[path as usize].traceb;
            nodes[new_traceb as usize].traceb = nxt;
        }
        path = nodes[path as usize].traceb;
    }

    path = max_ndx;
    while path != -1 && nodes[path as usize].traceb != -1 {
        let traceb = nodes[path as usize].traceb;
        nodes[traceb as usize].tracef = path;
        path = traceb;
    }

    if max_ndx == -1 || nodes[max_ndx as usize].traceb == -1 {
        return -1;
    } else {
        return max_ndx;
    }
}

#[inline]
fn score_connection_forward_start(
    nodes: &mut [Node],
    p1: i32,
    p2: i32,
    tinf: &Training,
    flag: i32,
) {
    use crate::node::intergenic_mod;
    use crate::sequence::STOP;

    let left;
    let right = nodes[p2 as usize].ndx;
    let mut score = 0.0f64;

    let p1_strand = nodes[p1 as usize].strand;
    let p1_type = nodes[p1 as usize].type_;

    if p1_strand == 1 && p1_type != STOP {
        return;
    }

    if p1_strand == -1 && p1_type == STOP {
        return;
    }

    if nodes[p1 as usize].traceb == -1 && p1_strand == 1 && p1_type == STOP {
        return;
    }
    if nodes[p1 as usize].traceb == -1 && p1_strand == -1 && p1_type != STOP {
        return;
    }

    if p1_strand == 1 && p1_type == STOP {
        left = nodes[p1 as usize].ndx + 2;
        if left >= right {
            return;
        }
        if flag == 1 {
            score = intergenic_mod(&nodes[p1 as usize], &nodes[p2 as usize], tinf);
        }

    }

    else if p1_strand == -1 && p1_type != STOP {
        left = nodes[p1 as usize].ndx;
        if left >= right {
            return;
        }
        if flag == 1 {
            score = intergenic_mod(&nodes[p1 as usize], &nodes[p2 as usize], tinf);
        }
    } else {

        return;
    }

    if nodes[p1 as usize].score + score >= nodes[p2 as usize].score {
        nodes[p2 as usize].score = nodes[p1 as usize].score + score;
        nodes[p2 as usize].traceb = p1;
        nodes[p2 as usize].ov_mark = -1;
    }
}

#[inline]
fn score_connection_forward_stop(
    nodes: &mut [Node],
    p1: i32,
    p2: i32,
    tinf: &Training,
    flag: i32,
) {
    use crate::node::intergenic_mod;
    use crate::sequence::STOP;

    let mut left = nodes[p1 as usize].ndx;
    let mut right = nodes[p2 as usize].ndx;
    let ovlp = 0i32;
    let mut score = 0.0f64;
    let mut scr_mod = 0.0f64;

    let p1_strand = nodes[p1 as usize].strand;
    let p1_type = nodes[p1 as usize].type_;

    if p1_strand == -1 && p1_type == STOP {
        return;
    }

    if p1_strand == -1 && p1_type != STOP {
        return;
    }

    if nodes[p1 as usize].traceb == -1 && p1_strand == 1 && p1_type == STOP {
        return;
    }

    if p1_strand == 1 && p1_type != STOP {
        if nodes[p2 as usize].stop_val >= nodes[p1 as usize].ndx {
            return;
        }
        if nodes[p1 as usize].ndx % 3 != nodes[p2 as usize].ndx % 3 {
            return;
        }
        right += 2;
        if flag == 0 {
            scr_mod = tinf.bias[0] * nodes[p1 as usize].gc_score[0]
                + tinf.bias[1] * nodes[p1 as usize].gc_score[1]
                + tinf.bias[2] * nodes[p1 as usize].gc_score[2];
        } else if flag == 1 {
            score = nodes[p1 as usize].cscore + nodes[p1 as usize].sscore;
        }
    }

    else if p1_strand == 1 && p1_type == STOP {
        if nodes[p2 as usize].stop_val >= nodes[p1 as usize].ndx {
            return;
        }
        if nodes[p1 as usize].star_ptr[(nodes[p2 as usize].ndx % 3) as usize] == -1 {
            return;
        }
        let n3_idx = nodes[p1 as usize].star_ptr[(nodes[p2 as usize].ndx % 3) as usize] as usize;
        left = nodes[n3_idx].ndx;
        right += 2;
        if flag == 0 {
            scr_mod = tinf.bias[0] * nodes[n3_idx].gc_score[0]
                + tinf.bias[1] * nodes[n3_idx].gc_score[1]
                + tinf.bias[2] * nodes[n3_idx].gc_score[2];
        } else if flag == 1 {
            score = nodes[n3_idx].cscore
                + nodes[n3_idx].sscore
                + intergenic_mod(&nodes[p1 as usize], &nodes[n3_idx], tinf);
        }
    } else {
        return;
    }

    if flag == 0 {
        score = ((right - left + 1 - (ovlp * 2)) as f64) * scr_mod;
    }

    if nodes[p1 as usize].score + score >= nodes[p2 as usize].score {
        nodes[p2 as usize].score = nodes[p1 as usize].score + score;
        nodes[p2 as usize].traceb = p1;
        nodes[p2 as usize].ov_mark = -1;
    }
}

#[inline]
fn score_connection_backward_start(
    nodes: &mut [Node],
    p1: i32,
    p2: i32,
    tinf: &Training,
    flag: i32,
) {
    use crate::sequence::STOP;

    let mut left = nodes[p1 as usize].ndx;
    let right = nodes[p2 as usize].ndx;
    let ovlp: i32;
    let mut score = 0.0f64;
    let mut scr_mod = 0.0f64;

    let p1_strand = nodes[p1 as usize].strand;
    let p1_type = nodes[p1 as usize].type_;

    if p1_strand == -1 && p1_type != STOP {
        return;
    }

    if p1_strand == 1 && p1_type != STOP {
        return;
    }

    if nodes[p1 as usize].traceb == -1 && p1_strand == 1 && p1_type == STOP {
        return;
    }
    if nodes[p1 as usize].traceb == -1 && p1_strand == -1 && p1_type != STOP {
        return;
    }

    if p1_strand == 1 && p1_type == STOP {
        if nodes[p2 as usize].stop_val - 2 >= nodes[p1 as usize].ndx + 2 {
            return;
        }
        ovlp = (nodes[p1 as usize].ndx + 2) - (nodes[p2 as usize].stop_val - 2) + 1;
        if ovlp >= MAX_OPP_OVLP {
            return;
        }
        if (nodes[p1 as usize].ndx + 2 - nodes[p2 as usize].stop_val - 2 + 1)
            >= (nodes[p2 as usize].ndx - nodes[p1 as usize].ndx + 3 + 1)
        {
            return;
        }
        let bnd = if nodes[p1 as usize].traceb == -1 {
            0
        } else {
            nodes[nodes[p1 as usize].traceb as usize].ndx
        };
        if (nodes[p1 as usize].ndx + 2 - nodes[p2 as usize].stop_val - 2 + 1)
            >= (nodes[p2 as usize].stop_val - 3 - bnd + 1)
        {
            return;
        }
        left = nodes[p2 as usize].stop_val - 2;
        if flag == 0 {
            scr_mod = tinf.bias[0] * nodes[p2 as usize].gc_score[0]
                + tinf.bias[1] * nodes[p2 as usize].gc_score[1]
                + tinf.bias[2] * nodes[p2 as usize].gc_score[2];
        } else if flag == 1 {
            score = nodes[p2 as usize].cscore + nodes[p2 as usize].sscore - 0.15 * tinf.st_wt;
        }

        if flag == 0 {
            score = ((right - left + 1 - (ovlp * 2)) as f64) * scr_mod;
        }
    }

    else if p1_strand == -1 && p1_type == STOP {

        if nodes[p1 as usize].stop_val <= nodes[p2 as usize].ndx {
            return;
        }
        if nodes[p1 as usize].ndx % 3 != nodes[p2 as usize].ndx % 3 {
            return;
        }
        left -= 2;
        if flag == 0 {
            scr_mod = tinf.bias[0] * nodes[p2 as usize].gc_score[0]
                + tinf.bias[1] * nodes[p2 as usize].gc_score[1]
                + tinf.bias[2] * nodes[p2 as usize].gc_score[2];
        } else if flag == 1 {
            score = nodes[p2 as usize].cscore + nodes[p2 as usize].sscore;
        }

        if flag == 0 {
            score = ((right - left + 1) as f64) * scr_mod;
        }
    } else {
        return;
    }

    if nodes[p1 as usize].score + score >= nodes[p2 as usize].score {
        nodes[p2 as usize].score = nodes[p1 as usize].score + score;
        nodes[p2 as usize].traceb = p1;
        nodes[p2 as usize].ov_mark = -1;
    }
}

#[inline]
fn score_connection_backward_stop(
    nodes: &mut [Node],
    p1: i32,
    p2: i32,
    tinf: &Training,
    flag: i32,
) {
    use crate::node::intergenic_mod;
    use crate::sequence::STOP;

    let mut left = nodes[p1 as usize].ndx;
    let mut right = nodes[p2 as usize].ndx;
    let mut ovlp = 0i32;
    let mut maxfr = -1i32;
    let mut score = 0.0f64;
    let mut scr_mod = 0.0f64;

    let p1_strand = nodes[p1 as usize].strand;
    let p1_type = nodes[p1 as usize].type_;

    if p1_strand == 1 && p1_type != STOP {
        return;
    }

    if nodes[p1 as usize].traceb == -1 && p1_strand == 1 && p1_type == STOP {
        return;
    }
    if nodes[p1 as usize].traceb == -1 && p1_strand == -1 && p1_type != STOP {
        return;
    }

    if p1_strand == 1 && p1_type == STOP {
        left += 2;
        right -= 2;
        if left >= right {
            return;
        }

        maxfr = -1;
        let mut maxval = 0.0f64;
        for i in 0..3 {
            if nodes[p2 as usize].star_ptr[i] == -1 {
                continue;
            }
            let n3_idx = nodes[p2 as usize].star_ptr[i] as usize;
            ovlp = left - nodes[n3_idx].stop_val + 3;
            if ovlp <= 0 || ovlp >= MAX_OPP_OVLP {
                continue;
            }
            if ovlp >= nodes[n3_idx].ndx - left {
                continue;
            }
            if nodes[p1 as usize].traceb == -1 {
                continue;
            }
            if ovlp >= nodes[n3_idx].stop_val - nodes[nodes[p1 as usize].traceb as usize].ndx - 2 {
                continue;
            }
            if (flag == 1
                && nodes[n3_idx].cscore + nodes[n3_idx].sscore
                    + intergenic_mod(&nodes[n3_idx], &nodes[p2 as usize], tinf)
                    > maxval)
                || (flag == 0
                    && tinf.bias[0] * nodes[n3_idx].gc_score[0]
                        + tinf.bias[1] * nodes[n3_idx].gc_score[1]
                        + tinf.bias[2] * nodes[n3_idx].gc_score[2]
                        > maxval)
            {
                maxfr = i as i32;
                maxval = nodes[n3_idx].cscore
                    + nodes[n3_idx].sscore
                    + intergenic_mod(&nodes[n3_idx], &nodes[p2 as usize], tinf);
            }
        }
        if maxfr != -1 {
            let n3_idx = nodes[p2 as usize].star_ptr[maxfr as usize] as usize;
            if flag == 0 {
                scr_mod = tinf.bias[0] * nodes[n3_idx].gc_score[0]
                    + tinf.bias[1] * nodes[n3_idx].gc_score[1]
                    + tinf.bias[2] * nodes[n3_idx].gc_score[2];
            } else if flag == 1 {
                score = nodes[n3_idx].cscore
                    + nodes[n3_idx].sscore
                    + intergenic_mod(&nodes[n3_idx], &nodes[p2 as usize], tinf);
            }
        } else if flag == 1 {
            score = intergenic_mod(&nodes[p1 as usize], &nodes[p2 as usize], tinf);
        }
    }

    else if p1_strand == -1 && p1_type != STOP {
        right -= 2;
        if left >= right {
            return;
        }
        if flag == 1 {
            score = intergenic_mod(&nodes[p1 as usize], &nodes[p2 as usize], tinf);
        }
    }

    else if p1_strand == -1 && p1_type == STOP {
        if nodes[p1 as usize].stop_val <= nodes[p2 as usize].ndx {
            return;
        }
        if nodes[p2 as usize].star_ptr[(nodes[p1 as usize].ndx % 3) as usize] == -1 {
            return;
        }
        let n3_idx = nodes[p2 as usize].star_ptr[(nodes[p1 as usize].ndx % 3) as usize] as usize;
        left -= 2;
        right = nodes[n3_idx].ndx;
        if flag == 0 {
            scr_mod = tinf.bias[0] * nodes[n3_idx].gc_score[0]
                + tinf.bias[1] * nodes[n3_idx].gc_score[1]
                + tinf.bias[2] * nodes[n3_idx].gc_score[2];
        } else if flag == 1 {
            score = nodes[n3_idx].cscore
                + nodes[n3_idx].sscore
                + intergenic_mod(&nodes[n3_idx], &nodes[p2 as usize], tinf);
        }
    } else {
        return;
    }

    if flag == 0 && scr_mod != 0.0 {
        score = ((right - left + 1 - (ovlp * 2)) as f64) * scr_mod;
    }

    if nodes[p1 as usize].score + score >= nodes[p2 as usize].score {
        nodes[p2 as usize].score = nodes[p1 as usize].score + score;
        nodes[p2 as usize].traceb = p1;
        nodes[p2 as usize].ov_mark = maxfr;
    }
}

pub fn score_connection(
    nodes: &mut [Node],
    p1: i32,
    p2: i32,
    tinf: &Training,
    flag: i32,
) {
    use crate::node::intergenic_mod;
    use crate::sequence::STOP;

    let mut left = nodes[p1 as usize].ndx;
    let mut right = nodes[p2 as usize].ndx;
    let mut ovlp = 0i32;
    let mut maxfr = -1i32;
    let mut score = 0.0f64;
    let mut scr_mod = 0.0f64;

    if nodes[p1 as usize].type_ != STOP
        && nodes[p2 as usize].type_ != STOP
        && nodes[p1 as usize].strand == nodes[p2 as usize].strand
    {
        return;
    }

    else if nodes[p1 as usize].strand == 1
        && nodes[p1 as usize].type_ != STOP
        && nodes[p2 as usize].strand == -1
    {
        return;
    }

    else if nodes[p1 as usize].strand == -1
        && nodes[p1 as usize].type_ == STOP
        && nodes[p2 as usize].strand == 1
    {
        return;
    }

    else if nodes[p1 as usize].strand == -1
        && nodes[p1 as usize].type_ != STOP
        && nodes[p2 as usize].strand == 1
        && nodes[p2 as usize].type_ == STOP
    {
        return;
    }

    if nodes[p1 as usize].traceb == -1
        && nodes[p1 as usize].strand == 1
        && nodes[p1 as usize].type_ == STOP
    {
        return;
    }
    if nodes[p1 as usize].traceb == -1
        && nodes[p1 as usize].strand == -1
        && nodes[p1 as usize].type_ != STOP
    {
        return;
    }

    if nodes[p1 as usize].strand == nodes[p2 as usize].strand
        && nodes[p1 as usize].strand == 1
        && nodes[p1 as usize].type_ != STOP
        && nodes[p2 as usize].type_ == STOP
    {
        if nodes[p2 as usize].stop_val >= nodes[p1 as usize].ndx {
            return;
        }
        if nodes[p1 as usize].ndx % 3 != nodes[p2 as usize].ndx % 3 {
            return;
        }
        right += 2;
        if flag == 0 {
            scr_mod = tinf.bias[0] * nodes[p1 as usize].gc_score[0]
                + tinf.bias[1] * nodes[p1 as usize].gc_score[1]
                + tinf.bias[2] * nodes[p1 as usize].gc_score[2];
        } else if flag == 1 {
            score = nodes[p1 as usize].cscore + nodes[p1 as usize].sscore;
        }
    }

    else if nodes[p1 as usize].strand == nodes[p2 as usize].strand
        && nodes[p1 as usize].strand == -1
        && nodes[p1 as usize].type_ == STOP
        && nodes[p2 as usize].type_ != STOP
    {
        if nodes[p1 as usize].stop_val <= nodes[p2 as usize].ndx {
            return;
        }
        if nodes[p1 as usize].ndx % 3 != nodes[p2 as usize].ndx % 3 {
            return;
        }
        left -= 2;
        if flag == 0 {
            scr_mod = tinf.bias[0] * nodes[p2 as usize].gc_score[0]
                + tinf.bias[1] * nodes[p2 as usize].gc_score[1]
                + tinf.bias[2] * nodes[p2 as usize].gc_score[2];
        } else if flag == 1 {
            score = nodes[p2 as usize].cscore + nodes[p2 as usize].sscore;
        }
    }

    else if nodes[p1 as usize].strand == 1
        && nodes[p1 as usize].type_ == STOP
        && nodes[p2 as usize].strand == 1
        && nodes[p2 as usize].type_ != STOP
    {
        left += 2;
        if left >= right {
            return;
        }
        if flag == 1 {
            score = intergenic_mod(&nodes[p1 as usize], &nodes[p2 as usize], tinf);
        }
    }

    else if nodes[p1 as usize].strand == 1
        && nodes[p1 as usize].type_ == STOP
        && nodes[p2 as usize].strand == -1
        && nodes[p2 as usize].type_ == STOP
    {
        left += 2;
        right -= 2;
        if left >= right {
            return;
        }

        maxfr = -1;
        let mut maxval = 0.0f64;
        for i in 0..3 {
            if nodes[p2 as usize].star_ptr[i] == -1 {
                continue;
            }
            let n3_idx = nodes[p2 as usize].star_ptr[i] as usize;
            ovlp = left - nodes[n3_idx].stop_val + 3;
            if ovlp <= 0 || ovlp >= MAX_OPP_OVLP {
                continue;
            }
            if ovlp >= nodes[n3_idx].ndx - left {
                continue;
            }
            if nodes[p1 as usize].traceb == -1 {
                continue;
            }
            if ovlp >= nodes[n3_idx].stop_val - nodes[nodes[p1 as usize].traceb as usize].ndx - 2 {
                continue;
            }
            if (flag == 1
                && nodes[n3_idx].cscore + nodes[n3_idx].sscore
                    + intergenic_mod(&nodes[n3_idx], &nodes[p2 as usize], tinf)
                    > maxval)
                || (flag == 0
                    && tinf.bias[0] * nodes[n3_idx].gc_score[0]
                        + tinf.bias[1] * nodes[n3_idx].gc_score[1]
                        + tinf.bias[2] * nodes[n3_idx].gc_score[2]
                        > maxval)
            {
                maxfr = i as i32;
                maxval = nodes[n3_idx].cscore
                    + nodes[n3_idx].sscore
                    + intergenic_mod(&nodes[n3_idx], &nodes[p2 as usize], tinf);
            }
        }
        if maxfr != -1 {
            let n3_idx = nodes[p2 as usize].star_ptr[maxfr as usize] as usize;
            if flag == 0 {
                scr_mod = tinf.bias[0] * nodes[n3_idx].gc_score[0]
                    + tinf.bias[1] * nodes[n3_idx].gc_score[1]
                    + tinf.bias[2] * nodes[n3_idx].gc_score[2];
            } else if flag == 1 {
                score = nodes[n3_idx].cscore
                    + nodes[n3_idx].sscore
                    + intergenic_mod(&nodes[n3_idx], &nodes[p2 as usize], tinf);
            }
        } else if flag == 1 {
            score = intergenic_mod(&nodes[p1 as usize], &nodes[p2 as usize], tinf);
        }
    }

    else if nodes[p1 as usize].strand == -1
        && nodes[p1 as usize].type_ != STOP
        && nodes[p2 as usize].strand == -1
        && nodes[p2 as usize].type_ == STOP
    {
        right -= 2;
        if left >= right {
            return;
        }
        if flag == 1 {
            score = intergenic_mod(&nodes[p1 as usize], &nodes[p2 as usize], tinf);
        }
    }

    else if nodes[p1 as usize].strand == -1
        && nodes[p1 as usize].type_ != STOP
        && nodes[p2 as usize].strand == 1
        && nodes[p2 as usize].type_ != STOP
    {
        if left >= right {
            return;
        }
        if flag == 1 {
            score = intergenic_mod(&nodes[p1 as usize], &nodes[p2 as usize], tinf);
        }
    }

    else if nodes[p1 as usize].strand == 1
        && nodes[p2 as usize].strand == 1
        && nodes[p1 as usize].type_ == STOP
        && nodes[p2 as usize].type_ == STOP
    {
        if nodes[p2 as usize].stop_val >= nodes[p1 as usize].ndx {
            return;
        }
        if nodes[p1 as usize].star_ptr[(nodes[p2 as usize].ndx % 3) as usize] == -1 {
            return;
        }
        let n3_idx = nodes[p1 as usize].star_ptr[(nodes[p2 as usize].ndx % 3) as usize] as usize;
        left = nodes[n3_idx].ndx;
        right += 2;
        if flag == 0 {
            scr_mod = tinf.bias[0] * nodes[n3_idx].gc_score[0]
                + tinf.bias[1] * nodes[n3_idx].gc_score[1]
                + tinf.bias[2] * nodes[n3_idx].gc_score[2];
        } else if flag == 1 {
            score = nodes[n3_idx].cscore
                + nodes[n3_idx].sscore
                + intergenic_mod(&nodes[p1 as usize], &nodes[n3_idx], tinf);
        }
    }

    else if nodes[p1 as usize].strand == -1
        && nodes[p1 as usize].type_ == STOP
        && nodes[p2 as usize].strand == -1
        && nodes[p2 as usize].type_ == STOP
    {
        if nodes[p1 as usize].stop_val <= nodes[p2 as usize].ndx {
            return;
        }
        if nodes[p2 as usize].star_ptr[(nodes[p1 as usize].ndx % 3) as usize] == -1 {
            return;
        }
        let n3_idx = nodes[p2 as usize].star_ptr[(nodes[p1 as usize].ndx % 3) as usize] as usize;
        left -= 2;
        right = nodes[n3_idx].ndx;
        if flag == 0 {
            scr_mod = tinf.bias[0] * nodes[n3_idx].gc_score[0]
                + tinf.bias[1] * nodes[n3_idx].gc_score[1]
                + tinf.bias[2] * nodes[n3_idx].gc_score[2];
        } else if flag == 1 {
            score = nodes[n3_idx].cscore
                + nodes[n3_idx].sscore
                + intergenic_mod(&nodes[n3_idx], &nodes[p2 as usize], tinf);
        }
    }

    else if nodes[p1 as usize].strand == 1
        && nodes[p1 as usize].type_ == STOP
        && nodes[p2 as usize].strand == -1
        && nodes[p2 as usize].type_ != STOP
    {
        if nodes[p2 as usize].stop_val - 2 >= nodes[p1 as usize].ndx + 2 {
            return;
        }
        ovlp = (nodes[p1 as usize].ndx + 2) - (nodes[p2 as usize].stop_val - 2) + 1;
        if ovlp >= MAX_OPP_OVLP {
            return;
        }
        if (nodes[p1 as usize].ndx + 2 - nodes[p2 as usize].stop_val - 2 + 1)
            >= (nodes[p2 as usize].ndx - nodes[p1 as usize].ndx + 3 + 1)
        {
            return;
        }
        let bnd = if nodes[p1 as usize].traceb == -1 {
            0
        } else {
            nodes[nodes[p1 as usize].traceb as usize].ndx
        };
        if (nodes[p1 as usize].ndx + 2 - nodes[p2 as usize].stop_val - 2 + 1)
            >= (nodes[p2 as usize].stop_val - 3 - bnd + 1)
        {
            return;
        }
        left = nodes[p2 as usize].stop_val - 2;
        if flag == 0 {
            scr_mod = tinf.bias[0] * nodes[p2 as usize].gc_score[0]
                + tinf.bias[1] * nodes[p2 as usize].gc_score[1]
                + tinf.bias[2] * nodes[p2 as usize].gc_score[2];
        } else if flag == 1 {
            score = nodes[p2 as usize].cscore + nodes[p2 as usize].sscore - 0.15 * tinf.st_wt;
        }
    }

    if flag == 0 {
        score = ((right - left + 1 - (ovlp * 2)) as f64) * scr_mod;
    }

    if nodes[p1 as usize].score + score >= nodes[p2 as usize].score {
        let _old_score = nodes[p2 as usize].score;
        nodes[p2 as usize].score = nodes[p1 as usize].score + score;
        nodes[p2 as usize].traceb = p1;
        nodes[p2 as usize].ov_mark = maxfr;
    }
}

pub fn eliminate_bad_genes(nodes: &mut [Node], dbeg: i32, tinf: &Training) {
    use crate::node::intergenic_mod;
    use crate::sequence::STOP;

    if dbeg == -1 {
        return;
    }

    let mut path = dbeg;
    while nodes[path as usize].traceb != -1 {
        path = nodes[path as usize].traceb;
    }

    while nodes[path as usize].tracef != -1 {
        if nodes[path as usize].strand == 1 && nodes[path as usize].type_ == STOP {
            let tracef = nodes[path as usize].tracef;
            let modifier = intergenic_mod(&nodes[path as usize], &nodes[tracef as usize], tinf);
            nodes[tracef as usize].sscore += modifier;
        }
        if nodes[path as usize].strand == -1 && nodes[path as usize].type_ != STOP {
            let tracef = nodes[path as usize].tracef;
            let modifier = intergenic_mod(&nodes[path as usize], &nodes[tracef as usize], tinf);
            nodes[path as usize].sscore += modifier;
        }
        path = nodes[path as usize].tracef;
    }

    path = dbeg;
    while nodes[path as usize].traceb != -1 {
        path = nodes[path as usize].traceb;
    }

    while nodes[path as usize].tracef != -1 {
        if nodes[path as usize].strand == 1
            && nodes[path as usize].type_ != STOP
            && nodes[path as usize].cscore + nodes[path as usize].sscore < 0.0
        {
            nodes[path as usize].elim = 1;
            let tracef = nodes[path as usize].tracef;
            nodes[tracef as usize].elim = 1;
        }
        if nodes[path as usize].strand == -1 && nodes[path as usize].type_ == STOP {
            let tracef = nodes[path as usize].tracef;
            if nodes[tracef as usize].cscore + nodes[tracef as usize].sscore < 0.0 {
                nodes[path as usize].elim = 1;
                nodes[tracef as usize].elim = 1;
            }
        }
        path = nodes[path as usize].tracef;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_connection_kind_forward_start() {
        use crate::sequence::STOP;
        let kind = ConnectionKind::from_node(1, 0, STOP);
        assert_eq!(kind, ConnectionKind::ForwardStart);
    }

    #[test]
    fn test_connection_kind_forward_stop() {
        use crate::sequence::STOP;
        let kind = ConnectionKind::from_node(1, STOP, STOP);
        assert_eq!(kind, ConnectionKind::ForwardStop);
    }

    #[test]
    fn test_connection_kind_backward_start() {
        use crate::sequence::STOP;
        let kind = ConnectionKind::from_node(-1, 0, STOP);
        assert_eq!(kind, ConnectionKind::BackwardStart);
    }

    #[test]
    fn test_connection_kind_backward_stop() {
        use crate::sequence::STOP;
        let kind = ConnectionKind::from_node(-1, STOP, STOP);
        assert_eq!(kind, ConnectionKind::BackwardStop);
    }

    fn create_test_nodes(n: usize) -> Vec<Node> {
        use crate::sequence::{ATG, GTG, STOP, TTG};

        let mut nodes = Vec::with_capacity(n);
        for i in 0..n {
            let mut node = Node::default();
            node.ndx = (i * 30) as i32;

            match i % 8 {
                0 => {

                    node.strand = 1;
                    node.type_ = STOP;
                    node.stop_val = (i.saturating_sub(4) * 30) as i32;
                }
                1 => {

                    node.strand = 1;
                    node.type_ = ATG;
                    node.stop_val = ((i + 3) * 30) as i32;
                }
                2 => {

                    node.strand = 1;
                    node.type_ = GTG;
                    node.stop_val = ((i + 2) * 30) as i32;
                }
                3 => {

                    node.strand = 1;
                    node.type_ = STOP;
                    node.stop_val = (i.saturating_sub(2) * 30) as i32;
                }
                4 => {

                    node.strand = -1;
                    node.type_ = STOP;
                    node.stop_val = ((i + 4) * 30) as i32;
                }
                5 => {

                    node.strand = -1;
                    node.type_ = ATG;
                    node.stop_val = (i.saturating_sub(1) * 30) as i32;
                }
                6 => {

                    node.strand = -1;
                    node.type_ = TTG;
                    node.stop_val = (i.saturating_sub(2) * 30) as i32;
                }
                7 => {

                    node.strand = -1;
                    node.type_ = STOP;
                    node.stop_val = ((i + 2) * 30) as i32;
                }
                _ => unreachable!(),
            }

            node.star_ptr = [-1, -1, -1];
            if node.type_ == STOP {

                let frame = (node.ndx % 3) as usize;
                if i > 0 {
                    node.star_ptr[frame] = (i - 1) as i32;
                }
            }

            node.cscore = ((i % 10) as f64) * 0.1;
            node.sscore = ((i % 5) as f64) * 0.05;
            node.gc_score = [0.3, 0.3, 0.4];

            nodes.push(node);
        }
        nodes
    }

    fn create_test_training() -> Training {
        let mut tinf = Training::default();
        tinf.bias = [0.0, 0.0, 0.0];
        tinf.st_wt = 4.35;
        tinf
    }

    #[test]
    fn test_dprog_basic() {

        let mut nodes = create_test_nodes(100);
        let tinf = create_test_training();

        let result = dprog(&mut nodes, 100, &tinf, 0);

        assert!(result >= -1);
    }

    #[test]
    fn test_dprog_empty() {
        let mut nodes: Vec<Node> = Vec::new();
        let tinf = create_test_training();

        assert_eq!(dprog(&mut nodes, 0, &tinf, 0), -1);
    }

    #[test]
    fn test_dprog_single_node() {
        let mut nodes = create_test_nodes(1);
        let tinf = create_test_training();

        let result = dprog(&mut nodes, 1, &tinf, 0);

        assert_eq!(result, -1);
    }

    #[test]
    fn test_dprog_flag1() {

        let mut nodes = create_test_nodes(80);
        let tinf = create_test_training();

        let result = dprog(&mut nodes, 80, &tinf, 1);

        assert!(result >= -1);
    }

    #[test]
    fn test_dprog_auto_matches_dprog() {

        let mut nodes1 = create_test_nodes(100);
        let mut nodes2 = nodes1.clone();
        let tinf = create_test_training();

        let result1 = dprog(&mut nodes1, 100, &tinf, 0);
        let result2 = dprog_auto(&mut nodes2, 100, &tinf, 0);

        assert_eq!(result1, result2);

        for i in 0..100 {
            assert_eq!(nodes1[i].score, nodes2[i].score);
            assert_eq!(nodes1[i].traceb, nodes2[i].traceb);
        }
    }
}

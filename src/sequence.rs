// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2007-2016 Doug Hyatt, Univ. of Tennessee / UT-Battelle (original Prodigal C)
// Copyright (C) 2026 Sunju Kim (Rust reimplementation)

use crate::fptr::FilePtr;
use crate::training::Training;
use crate::translation;

// MAX_SEQ (300 Mb) is the single-genome TRAINING window: read_seq_training still trains on at
// most the first 300 Mb of input (largest human chromosome, chr1, is 248 Mb), preserving the
// original training-window semantics exactly. It is also used by the api.rs library helpers as
// their allocation size / overflow guard. Prodigal's original cap was 32 Mb.
//
// The standalone binary's gene-finding read path (next_seq_multi) no longer pre-allocates a
// fixed MAX_SEQ buffer: it collects each record's bytes into a right-sized Vec and builds the
// bitmaps per record (matching the meta_api / library path). It only enforces MAX_SEQ_GUARD as
// an overflow ceiling. NOTE: the metagenomic library path (meta_api::run_meta, used by rust-ise)
// allocates per-sequence dynamically and never reads MAX_SEQ, so this does not affect isscan.
pub const MAX_SEQ: usize = 300_000_000;

// Overflow ceiling for the standalone binary's per-record read path. Bitmap positions are i32
// and 2*len must stay < 2^31 (~2.147e9), so a 1 Gb ceiling keeps 2*len < 2^31 with headroom.
// This replaces the old 300 Mb functional cap: sequences up to this limit are now processed,
// and anything larger is rejected cleanly instead of causing an i32 index overflow.
pub const MAX_SEQ_GUARD: usize = 1_000_000_000;
pub const MAX_LINE: usize = 10_000;
pub const WINDOW: usize = 120;
pub const MASK_SIZE: usize = 50;
pub const MAX_MASKS: usize = 5_000;

pub const ATG: i32 = 0;
pub const GTG: i32 = 1;
pub const TTG: i32 = 2;
pub const STOP: i32 = 3;

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Nucleotide {
    A = 0b000,
    G = 0b001,
    C = 0b010,
    T = 0b011,
    N = 0b110,
}

impl Nucleotide {
    pub const A_VAL: u8 = 0b000;
    pub const G_VAL: u8 = 0b001;
    pub const C_VAL: u8 = 0b010;
    pub const T_VAL: u8 = 0b011;
    pub const N_VAL: u8 = 0b110;
}

const ASCII_TO_NUCLEOTIDE: [u8; 256] = {
    let mut table = [Nucleotide::N_VAL; 256];
    table[b'A' as usize] = Nucleotide::A_VAL;
    table[b'a' as usize] = Nucleotide::A_VAL;
    table[b'G' as usize] = Nucleotide::G_VAL;
    table[b'g' as usize] = Nucleotide::G_VAL;
    table[b'C' as usize] = Nucleotide::C_VAL;
    table[b'c' as usize] = Nucleotide::C_VAL;
    table[b'T' as usize] = Nucleotide::T_VAL;
    table[b't' as usize] = Nucleotide::T_VAL;
    table
};

#[inline(always)]
pub fn ascii_to_nucleotide(ch: u8) -> u8 {
    ASCII_TO_NUCLEOTIDE[ch as usize]
}

#[inline(always)]
pub fn complement(nuc: u8) -> u8 {
    nuc ^ 0b011
}

const TAA_IS_STOP: [u8; 34] = [
    0,
    1,
    1,
    1,
    1,
    1,
    0,
    0,
    0,
    1,
    1,
    1,
    1,
    1,
    1,
    0,
    1,
    0,
    0,
    0,
    0,
    1,
    0,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
];

const TAG_IS_STOP: [u8; 34] = [
    0,
    1,
    1,
    1,
    1,
    1,
    0,
    0,
    0,
    1,
    1,
    1,
    1,
    1,
    1,
    0,
    1,
    0,
    0,
    0,
    0,
    1,
    0,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
];

const TGA_IS_STOP: [u8; 34] = [
    0,
    1,
    0,
    0,
    0,
    0,
    1,
    0,
    0,
    0,
    1,
    1,
    1,
    0,
    0,
    1,
    1,
    0,
    0,
    0,
    0,
    0,
    1,
    1,
    0,
    0,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
];

#[allow(dead_code)]
const AGA_IS_STOP: [u8; 34] = [
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

#[allow(dead_code)]
const AGG_IS_STOP: [u8; 34] = [
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

#[allow(dead_code)]
const ATG_IS_START: [u8; 34] = [1; 34];

const GTG_IS_START: [u8; 34] = [
    1,
    0,
    1,
    0,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    0,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    0,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
];

const TTG_IS_START: [u8; 34] = [
    0,
    0,
    0,
    0,
    1,
    1,
    1,
    1,
    1,
    0,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    0,
    0,
    0,
    0,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
];

#[inline(always)]
pub fn is_stop_fast(digits: &[u8], slen: usize, pos: usize, tt: usize, strand: i32) -> bool {
    let tt = tt.min(33);

    let (x0, x1, x2) = if strand == 1 {

        (digits[pos], digits[pos + 1], digits[pos + 2])
    } else {

        (
            complement(digits[slen - 1 - pos]),
            complement(digits[slen - 2 - pos]),
            complement(digits[slen - 3 - pos]),
        )
    };

    if x0 == Nucleotide::T_VAL && x1 == Nucleotide::A_VAL && x2 == Nucleotide::A_VAL {
        return TAA_IS_STOP[tt] == 1;
    }

    if x0 == Nucleotide::T_VAL && x1 == Nucleotide::A_VAL && x2 == Nucleotide::G_VAL {
        return TAG_IS_STOP[tt] == 1;
    }

    if x0 == Nucleotide::T_VAL && x1 == Nucleotide::G_VAL && x2 == Nucleotide::A_VAL {
        return TGA_IS_STOP[tt] == 1;
    }

    if tt == 2 {

        if x0 == Nucleotide::A_VAL && x1 == Nucleotide::G_VAL && x2 == Nucleotide::A_VAL {
            return true;
        }

        if x0 == Nucleotide::A_VAL && x1 == Nucleotide::G_VAL && x2 == Nucleotide::G_VAL {
            return true;
        }
    }

    false
}

#[inline(always)]
pub fn is_start_fast(digits: &[u8], slen: usize, pos: usize, tt: usize, strand: i32) -> bool {
    let tt = tt.min(33);

    let (x0, x1, x2) = if strand == 1 {
        (digits[pos], digits[pos + 1], digits[pos + 2])
    } else {
        (
            complement(digits[slen - 1 - pos]),
            complement(digits[slen - 2 - pos]),
            complement(digits[slen - 3 - pos]),
        )
    };

    if x0 == Nucleotide::A_VAL && x1 == Nucleotide::T_VAL && x2 == Nucleotide::G_VAL {
        return true;
    }

    if x0 == Nucleotide::G_VAL && x1 == Nucleotide::T_VAL && x2 == Nucleotide::G_VAL {
        return GTG_IS_START[tt] == 1;
    }

    if x0 == Nucleotide::T_VAL && x1 == Nucleotide::T_VAL && x2 == Nucleotide::G_VAL {
        return TTG_IS_START[tt] == 1;
    }

    false
}

#[inline(always)]
pub fn start_type_fast(digits: &[u8], slen: usize, pos: usize, strand: i32) -> i32 {
    let (x0, x1, x2) = if strand == 1 {
        (digits[pos], digits[pos + 1], digits[pos + 2])
    } else {
        (
            complement(digits[slen - 1 - pos]),
            complement(digits[slen - 2 - pos]),
            complement(digits[slen - 3 - pos]),
        )
    };

    if x2 != Nucleotide::G_VAL {
        return -1;
    }

    if x1 != Nucleotide::T_VAL {
        return -1;
    }

    if x0 == Nucleotide::A_VAL {
        return ATG;
    }

    if x0 == Nucleotide::G_VAL {
        return GTG;
    }

    if x0 == Nucleotide::T_VAL {
        return TTG;
    }

    -1
}

pub fn bitmap_to_digits(seq: &[u8], useq: &[u8], slen: i32, digits: &mut [u8]) {
    for i in 0..slen as usize {
        if crate::bitmap::test(useq, i as i32) == 1 {
            digits[i] = Nucleotide::N_VAL;
        } else {
            let bit0 = crate::bitmap::test(seq, (i * 2) as i32);
            let bit1 = crate::bitmap::test(seq, (i * 2 + 1) as i32);

            digits[i] = match (bit0, bit1) {
                (0, 0) => Nucleotide::A_VAL,
                (0, 1) => Nucleotide::C_VAL,
                (1, 0) => Nucleotide::G_VAL,
                (1, 1) => Nucleotide::T_VAL,
                _ => Nucleotide::N_VAL,
            };
        }
    }
}

#[inline(always)]
pub fn is_gc_fast(digits: &[u8], pos: usize) -> bool {
    let nuc = digits[pos];
    nuc == Nucleotide::G_VAL || nuc == Nucleotide::C_VAL
}

pub const ACCEPT: &str = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.:^*$@!+_?-|";

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Mask {
    pub begin: i32,
    pub end: i32,
}

/// Reads the single-genome training sequence, collecting the processed nucleotide bytes into
/// `dna` (a growable buffer sized to the sequence's true length by the caller, which then builds
/// the bitmaps via `bitmap::build_bitmaps`). Multiple records are concatenated with the same
/// 12-base "TTAATTAATTAA" spacer the original bit-level code inserted. Training still stops after
/// the first `MAX_SEQ` (300 Mb) bases, preserving the original training-window semantics exactly.
/// Returns the number of bases collected (`dna.len()`).
pub fn read_seq_training(
    fptr: &mut FilePtr,
    do_mask: i32,
    mlist: &mut [Mask],
    nmask: &mut i32,
    dna: &mut Vec<u8>,
) -> i32 {
    use std::io::BufRead;

    dna.clear();

    let mut line = String::with_capacity(MAX_LINE + 1);
    let mut hdr = 0i32;
    let mut fhdr = 0i32;
    let mut len = 0i32;
    let mut wrn = 0i32;
    let mut mask_beg = -1i32;

    loop {
        line.clear();
        match fptr.read_line(&mut line) {
            Ok(0) => break,
            Ok(_) => {},
            Err(_) => break,
        }

        if hdr == 0 && !line.ends_with('\n') && wrn == 0 {
            wrn = 1;
            eprintln!("\n\nWarning: saw non-sequence line longer than {} chars, sequence might not be read correctly.\n", MAX_LINE);
        }

        let line_bytes = line.as_bytes();
        let is_fasta_header = line_bytes.first() == Some(&b'>');
        let is_sq_header = line_bytes.len() >= 2 && line_bytes[0] == b'S' && line_bytes[1] == b'Q';
        let is_origin_header = line.len() > 6 && line.starts_with("ORIGIN");

        if is_fasta_header || is_sq_header || is_origin_header {
            hdr = 1;
            if fhdr > 0 {

                for i in 0..12 {

                    // i%4 in {0,1} -> 'T' (both 2-bit bits set), else -> 'A' (no bits): the exact
                    // 12-base spacer the original bit-level code inserted between concatenated
                    // records. build_bitmaps re-derives the same bits from these bytes.
                    dna.push(if i % 4 == 0 || i % 4 == 1 { b'T' } else { b'A' });
                    len += 1;
                }
            }
            fhdr += 1;
        } else if hdr == 1 && line_bytes.len() >= 2 && line_bytes[0] == b'/' && line_bytes[1] == b'/' {

            hdr = 0;
        } else if hdr == 1 {

            if line.contains("Expand") && line.contains("gap") {
                if let Some(gap_pos) = line.find("gap") {
                    let gap_str = &line[gap_pos + 4..];
                    if let Ok(gapsize) = gap_str.trim_start().split_whitespace().next()
                        .unwrap_or("0").parse::<usize>() {
                        if gapsize >= 1 && gapsize <= MAX_LINE {

                            line.clear();
                            for _ in 0..gapsize {
                                line.push('n');
                            }
                        } else {
                            eprintln!("Error: gap size in gbk file can't exceed line size.");
                            std::process::exit(51);
                        }
                    }
                }
            }

            for ch in line.bytes() {

                if ch < b'A' || ch > b'z' {
                    continue;
                }

                if do_mask == 1 && mask_beg != -1 && ch != b'N' && ch != b'n' {
                    if len - mask_beg >= MASK_SIZE as i32 {
                        if *nmask == MAX_MASKS as i32 {
                            eprintln!("Error: saw too many regions of 'N''s in the sequence.");
                            std::process::exit(52);
                        }
                        mlist[*nmask as usize].begin = mask_beg;
                        mlist[*nmask as usize].end = len - 1;
                        *nmask += 1;
                    }
                    mask_beg = -1;
                }

                if do_mask == 1 && mask_beg == -1 && (ch == b'N' || ch == b'n') {
                    mask_beg = len;
                }

                // Collect the raw byte; build_bitmaps applies the identical 2-bit/N encoding
                // (g/G,t/T,c/C,a/A -> same codes; anything else -> N code + useq mask bit).
                dna.push(ch);
                len += 1;
            }
        }

        if len + MAX_LINE as i32 >= MAX_SEQ as i32 {
            eprintln!("\n\nWarning: Sequence is long (max {} for training).", MAX_SEQ);
            eprintln!("Training on the first {} bases.\n", MAX_SEQ);
            break;
        }
    }

    if fhdr > 1 {
        for i in 0..12 {
            dna.push(if i % 4 == 0 || i % 4 == 1 { b'T' } else { b'A' });
            len += 1;
        }
    }

    len
}

/// Reads the next record from a multi-FASTA / Genbank / EMBL stream, collecting the processed
/// nucleotide bytes into `dna` (a growable buffer). The caller sizes the 2-bit/1-bit bitmaps to
/// `dna.len()` and builds them via `bitmap::build_bitmaps`, so each record gets right-sized
/// buffers instead of a shared MAX_SEQ buffer. Returns the number of bases read (`dna.len()`),
/// or -1 at end of input. Sequences larger than `MAX_SEQ_GUARD` are rejected cleanly.
pub fn next_seq_multi(
    fptr: &mut FilePtr,
    dna: &mut Vec<u8>,
    sctr: &mut i32,
    do_mask: i32,
    mlist: &mut [Mask],
    nmask: &mut i32,
    cur_hdr: &mut [u8],
    new_hdr: &mut [u8],
) -> i32 {
    use std::io::BufRead;

    dna.clear();

    let mut line = String::with_capacity(MAX_LINE + 1);
    let mut reading_seq = 0i32;
    let mut genbank_end = 0i32;
    let mut len = 0i32;
    let mut wrn = 0i32;
    let mut mask_beg = -1i32;

    let default_new = format!("Prodigal_Seq_{}", *sctr + 2);
    let default_bytes = default_new.as_bytes();
    let copy_len = default_bytes.len().min(new_hdr.len() - 1);
    new_hdr[..copy_len].copy_from_slice(&default_bytes[..copy_len]);
    new_hdr[copy_len] = 0;

    if *sctr > 0 {
        reading_seq = 1;
    }

    loop {
        line.clear();
        match fptr.read_line(&mut line) {
            Ok(0) => break,
            Ok(_) => {},
            Err(_) => break,
        }

        let line_bytes = line.as_bytes();

        if reading_seq == 0 && !line.ends_with('\n') && wrn == 0 {
            wrn = 1;
            eprintln!("\n\nWarning: saw non-sequence line longer than {} chars, sequence might not be read correctly.\n", MAX_LINE);
        }

        if line.len() > 10 && line.starts_with("DEFINITION") {
            let def_content = line.get(12..).unwrap_or("").trim_end_matches('\n');
            if genbank_end == 0 {

                let def_bytes = def_content.as_bytes();
                let copy_len = def_bytes.len().min(cur_hdr.len() - 1);
                cur_hdr[..copy_len].copy_from_slice(&def_bytes[..copy_len]);
                cur_hdr[copy_len] = 0;
            } else {

                let def_bytes = def_content.as_bytes();
                let copy_len = def_bytes.len().min(new_hdr.len() - 1);
                new_hdr[..copy_len].copy_from_slice(&def_bytes[..copy_len]);
                new_hdr[copy_len] = 0;
            }
        }

        let is_fasta_header = line_bytes.first() == Some(&b'>');
        let is_sq_header = line_bytes.len() >= 2 && line_bytes[0] == b'S' && line_bytes[1] == b'Q';
        let is_origin_header = line.len() > 6 && line.starts_with("ORIGIN");

        if is_fasta_header || is_sq_header || is_origin_header {

            if reading_seq == 1 || genbank_end == 1 || *sctr > 0 {
                if is_fasta_header {

                    let hdr_content = line.get(1..).unwrap_or("").trim_end_matches('\n');
                    let hdr_bytes = hdr_content.as_bytes();
                    let copy_len = hdr_bytes.len().min(new_hdr.len() - 1);
                    new_hdr[..copy_len].copy_from_slice(&hdr_bytes[..copy_len]);
                    new_hdr[copy_len] = 0;
                }
                break;
            }

            if is_fasta_header {
                let hdr_content = line.get(1..).unwrap_or("").trim_end_matches('\n');
                let hdr_bytes = hdr_content.as_bytes();
                let copy_len = hdr_bytes.len().min(cur_hdr.len() - 1);
                cur_hdr[..copy_len].copy_from_slice(&hdr_bytes[..copy_len]);
                cur_hdr[copy_len] = 0;
            }
            reading_seq = 1;
        } else if reading_seq == 1 && line_bytes.len() >= 2 && line_bytes[0] == b'/' && line_bytes[1] == b'/' {

            reading_seq = 0;
            genbank_end = 1;
        } else if reading_seq == 1 {

            if line.contains("Expand") && line.contains("gap") {
                if let Some(gap_pos) = line.find("gap") {
                    let gap_str = &line[gap_pos + 4..];
                    if let Ok(gapsize) = gap_str.trim_start().split_whitespace().next()
                        .unwrap_or("0").parse::<usize>() {
                        if gapsize >= 1 && gapsize <= MAX_LINE {
                            line.clear();
                            for _ in 0..gapsize {
                                line.push('n');
                            }
                        } else {
                            eprintln!("Error: gap size in gbk file can't exceed line size.");
                            std::process::exit(54);
                        }
                    }
                }
            }

            for ch in line.bytes() {
                if ch < b'A' || ch > b'z' {
                    continue;
                }

                if do_mask == 1 && mask_beg != -1 && ch != b'N' && ch != b'n' {
                    if len - mask_beg >= MASK_SIZE as i32 {
                        if *nmask == MAX_MASKS as i32 {
                            eprintln!("Error: saw too many regions of 'N''s in the sequence.");
                            std::process::exit(55);
                        }
                        mlist[*nmask as usize].begin = mask_beg;
                        mlist[*nmask as usize].end = len - 1;
                        *nmask += 1;
                    }
                    mask_beg = -1;
                }

                if do_mask == 1 && mask_beg == -1 && (ch == b'N' || ch == b'n') {
                    mask_beg = len;
                }

                // Collect the raw byte; build_bitmaps applies the identical 2-bit/N encoding.
                dna.push(ch);
                len += 1;
            }
        }

        if dna.len() + MAX_LINE >= MAX_SEQ_GUARD {
            eprintln!("\n\nError: sequence exceeds the {} bp limit supported by this build.", MAX_SEQ_GUARD);
            std::process::exit(56);
        }
    }

    if len == 0 {
        return -1;
    }

    *sctr += 1;
    len
}

pub fn rcom_seq(seq: &[u8], rseq: &mut [u8], useq: &[u8], slen: i32) {
    let bit_len = slen * 2;

    for i in 0..(bit_len as usize / 8 + 1) {
        rseq[i] = 0;
    }

    for i in 0..bit_len {
        if crate::bitmap::test(seq, i) == 0 {
            crate::bitmap::set(rseq, bit_len - i - 1 + if i % 2 == 0 { -1 } else { 1 });
        }
    }

    for i in 0..slen {
        if crate::bitmap::test(useq, i) == 1 {
            crate::bitmap::toggle(rseq, bit_len - 1 - i * 2);
            crate::bitmap::toggle(rseq, bit_len - 2 - i * 2);
        }
    }
}

pub fn calc_short_header(header: &str, short_header: &mut String, sctr: i32) {

    short_header.clear();

    let word = match header.find(|c: char| c == ' ' || c == '\t' || c == '\r' || c == '\n') {
        Some(pos) => &header[..pos],
        None => header,
    };

    if word.is_empty() {
        short_header.push_str(&format!("Prodigal_Seq_{}", sctr));
    } else {
        short_header.push_str(word);
    }
}

#[inline(always)]
fn base2(seq: &[u8], n: i32) -> u8 {
    let p = (n as usize) * 2;
    (seq[p >> 3] >> (p & 7)) & 0b11
}

#[inline]
pub fn is_a(seq: &[u8], n: i32) -> bool {
    base2(seq, n) == 0
}

#[inline]
pub fn is_c(seq: &[u8], n: i32) -> bool {
    base2(seq, n) == 2
}

#[inline]
pub fn is_g(seq: &[u8], n: i32) -> bool {
    base2(seq, n) == 1
}

#[inline]
pub fn is_t(seq: &[u8], n: i32) -> bool {
    base2(seq, n) == 3
}

#[inline]
pub fn is_n(useq: &[u8], n: i32) -> bool {

    crate::bitmap::test(useq, n) == 1
}

#[inline]
pub fn is_gc(seq: &[u8], n: i32) -> bool {

    let b = base2(seq, n);
    b == 1 || b == 2
}

#[inline]
pub fn is_atg(seq: &[u8], n: i32) -> bool {
    is_a(seq, n) && is_t(seq, n + 1) && is_g(seq, n + 2)
}

#[inline]
pub fn is_gtg(seq: &[u8], n: i32) -> bool {
    is_g(seq, n) && is_t(seq, n + 1) && is_g(seq, n + 2)
}

#[inline]
pub fn is_ttg(seq: &[u8], n: i32) -> bool {
    is_t(seq, n) && is_t(seq, n + 1) && is_g(seq, n + 2)
}

pub fn is_start(seq: &[u8], n: i32, tinf: &Training) -> bool {

    if is_atg(seq, n) {
        return true;
    }

    if tinf.trans_table == 6
        || tinf.trans_table == 10
        || tinf.trans_table == 14
        || tinf.trans_table == 15
        || tinf.trans_table == 16
        || tinf.trans_table == 22
    {
        return false;
    }

    if is_gtg(seq, n) {

        if tinf.trans_table == 1
            || tinf.trans_table == 3
            || tinf.trans_table == 12
            || tinf.trans_table == 22
        {
            return false;
        }
        return true;
    }

    if is_ttg(seq, n) {

        if tinf.trans_table < 4
            || tinf.trans_table == 9
            || (tinf.trans_table >= 21 && tinf.trans_table < 25)
        {
            return false;
        }
        return true;
    }

    false
}

#[inline]
pub fn classify_start(seq: &[u8], n: i32, tinf: &Training) -> i32 {
    if !is_start(seq, n, tinf) {
        return -1;
    }
    if is_atg(seq, n) {
        ATG
    } else if is_gtg(seq, n) {
        GTG
    } else if is_ttg(seq, n) {
        TTG
    } else {
        -1
    }
}

pub fn is_stop(seq: &[u8], n: i32, tinf: &Training) -> bool {
    if tinf.trans_table == 4 {

        (is_t(seq, n) && is_a(seq, n + 1) && is_a(seq, n + 2))
            || (is_t(seq, n) && is_a(seq, n + 1) && is_g(seq, n + 2))
    } else if tinf.trans_table == 11 {

        (is_t(seq, n) && is_a(seq, n + 1) && is_a(seq, n + 2))
            || (is_t(seq, n) && is_a(seq, n + 1) && is_g(seq, n + 2))
            || (is_t(seq, n) && is_g(seq, n + 1) && is_a(seq, n + 2))
    } else {

        (is_t(seq, n) && is_a(seq, n + 1) && is_a(seq, n + 2))
            || (is_t(seq, n) && is_a(seq, n + 1) && is_g(seq, n + 2))
    }
}

pub fn gc_content(seq: &[u8], a: i32, b: i32) -> f64 {
    let mut gc = 0.0;
    let mut sum = 0.0;
    for i in a..=b {
        if is_g(seq, i) || is_c(seq, i) {
            gc += 1.0;
        }
        sum += 1.0;
    }
    gc / sum
}

pub fn amino(seq: &[u8], n: i32, tinf: &Training, is_init: i32) -> u8 {

    let x0 = get_nucleotide_3bit(seq, n);
    let x1 = get_nucleotide_3bit(seq, n + 1);
    let x2 = get_nucleotide_3bit(seq, n + 2);

    let aa = translation::translate_codon(tinf.trans_table as usize, x0, x1, x2);
    if aa == b'*' {
        return b'*';
    }

    if is_init == 1 && is_start(seq, n, tinf) {
        return b'M';
    }

    aa
}

#[inline(always)]
pub fn amino_fast(
    digits: &[u8],
    slen: usize,
    pos: usize,
    tt: usize,
    strand: i32,
    is_init: bool,
) -> u8 {
    translation::translate_codon_init(digits, slen, pos, tt, strand, is_init)
}

#[inline(always)]
fn get_nucleotide_3bit(seq: &[u8], n: i32) -> u8 {
    let bit0 = crate::bitmap::test(seq, n * 2);
    let bit1 = crate::bitmap::test(seq, n * 2 + 1);

    match (bit0, bit1) {
        (0, 0) => Nucleotide::A_VAL,
        (0, 1) => Nucleotide::C_VAL,
        (1, 0) => Nucleotide::G_VAL,
        (1, 1) => Nucleotide::T_VAL,
        _ => Nucleotide::N_VAL,
    }
}

pub fn amino_num(letter: u8) -> i32 {
    match letter {
        b'A' => 0, b'C' => 1, b'D' => 2, b'E' => 3,
        b'F' => 4, b'G' => 5, b'H' => 6, b'I' => 7,
        b'K' => 8, b'L' => 9, b'M' => 10, b'N' => 11,
        b'P' => 12, b'Q' => 13, b'R' => 14, b'S' => 15,
        b'T' => 16, b'V' => 17, b'W' => 18, b'Y' => 19,
        _ => 20,
    }
}

pub fn amino_letter(num: i32) -> u8 {
    let letters = b"ACDEFGHIKLMNPQRSTVWY";
    if num >= 0 && num < 20 {
        letters[num as usize]
    } else {
        b'X'
    }
}

pub fn rframe(a: i32, b: i32) -> i32 {
    (b - a + 3) % 3
}

pub fn max_fr(a: i32, b: i32, c: i32) -> i32 {

    if a > b {
        if a > c {
            0
        } else {
            2
        }
    } else {
        if b > c {
            1
        } else {
            2
        }
    }
}

pub fn calc_most_gc_frame(seq: &[u8], slen: i32) -> Vec<i32> {
    let mut gp = vec![-1i32; slen as usize];
    let mut fwd = vec![0i32; slen as usize];
    let mut bwd = vec![0i32; slen as usize];
    let mut tot = vec![0i32; slen as usize];

    for i in 0..3 {
        for j in (i..slen).step_by(1) {
            if j < 3 {
                fwd[j as usize] = is_gc(seq, j) as i32;
            } else {
                fwd[j as usize] = fwd[(j - 3) as usize] + is_gc(seq, j) as i32;
            }
            if j < 3 {
                bwd[(slen - j - 1) as usize] = is_gc(seq, slen - j - 1) as i32;
            } else {
                bwd[(slen - j - 1) as usize] = bwd[(slen - j + 2) as usize] + is_gc(seq, slen - j - 1) as i32;
            }
        }
    }

    for i in 0..slen {
        tot[i as usize] = fwd[i as usize] + bwd[i as usize] - is_gc(seq, i) as i32;
        if i - WINDOW as i32 / 2 >= 0 {
            tot[i as usize] -= fwd[(i - WINDOW as i32 / 2) as usize];
        }
        if i + WINDOW as i32 / 2 < slen {
            tot[i as usize] -= bwd[(i + WINDOW as i32 / 2) as usize];
        }
    }

    let mut i = 0;
    while i < slen - 2 {
        let win = max_fr(tot[i as usize], tot[(i + 1) as usize], tot[(i + 2) as usize]);
        for j in 0..3 {
            gp[(i + j) as usize] = win;
        }
        i += 3;
    }

    gp
}

#[inline(always)]
pub fn mer_ndx(length: i32, seq: &[u8], pos: i32) -> i32 {
    let mut ndx = 0;
    for i in 0..(2 * length) {
        ndx |= (crate::bitmap::test(seq, pos * 2 + i) as i32) << i;
    }
    ndx
}

pub fn mer_text(qt: &mut [u8], length: i32, ndx: i32) {
    let letters = [b'A', b'G', b'C', b'T'];
    if length == 0 {
        qt[0] = b'N';
        qt[1] = b'o';
        qt[2] = b'n';
        qt[3] = b'e';
        qt[4] = 0;
    } else {
        for i in 0..length as usize {
            let val = ((ndx & (1 << (2 * i))) + (ndx & (1 << (2 * i + 1)))) >> (i * 2);
            qt[i] = letters[val as usize];
        }
        qt[length as usize] = 0;
    }
}

pub fn calc_mer_bg(mer_len: i32, seq: &[u8], rseq: &[u8], slen: i32, bg: &mut [f64]) {
    let mut size = 1;
    for _ in 0..mer_len {
        size *= 4;
    }

    let mut counts = vec![0i32; size];
    let mut glob = 0;

    for i in 0..=(slen - mer_len) {
        counts[mer_ndx(mer_len, seq, i) as usize] += 1;
        counts[mer_ndx(mer_len, rseq, i) as usize] += 1;
        glob += 2;
    }

    for i in 0..size {
        bg[i] = counts[i] as f64 / glob as f64;
    }
}

pub fn shine_dalgarno_exact(seq: &[u8], pos: i32, start: i32, rwt: &[f64]) -> i32 {
    let limit = imin(6, start - 4 - pos);
    let mut match_scores = [-10.0; 6];

    for i in 0..limit {
        if pos + i >= 0 {
            if i % 3 == 0 && is_a(seq, pos + i) {
                match_scores[i as usize] = 2.0;
            } else if i % 3 != 0 && is_g(seq, pos + i) {
                match_scores[i as usize] = 3.0;
            }
        }
    }

    let mut max_val = 0;
    for i in (3..=limit).rev() {
        for j in 0..=(limit - i) {
            let mut cur_ctr = -2.0;
            let mut mism = 0;

            for k in j..(j + i) {
                cur_ctr += match_scores[k as usize];
                if match_scores[k as usize] < 0.0 {
                    mism += 1;
                }
            }

            if mism > 0 {
                continue;
            }

            let rdis = start - (pos + j + i);
            let dis_flag = if rdis < 5 && i < 5 {
                2
            } else if rdis < 5 && i >= 5 {
                1
            } else if rdis > 10 && rdis <= 12 && i < 5 {
                1
            } else if rdis > 10 && rdis <= 12 && i >= 5 {
                2
            } else if rdis >= 13 {
                3
            } else {
                0
            };

            if rdis > 15 || cur_ctr < 6.0 {
                continue;
            }

            let cur_val = match (cur_ctr as i32, dis_flag) {
                (6, 2) => 1,
                (6, 3) => 2,
                (8, 3) | (9, 3) => 3,
                (6, 1) => 6,
                (11, 3) | (12, 3) | (14, 3) => 10,
                (8, 2) | (9, 2) => 11,
                (8, 1) | (9, 1) => 12,
                (6, 0) => 13,
                (8, 0) => 15,
                (9, 0) => 16,
                (11, 2) | (12, 2) => 20,
                (11, 1) => 21,
                (11, 0) => 22,
                (12, 1) => 23,
                (12, 0) => 24,
                (14, 2) => 25,
                (14, 1) => 26,
                (14, 0) => 27,
                _ => 0,
            };

            if rwt[cur_val as usize] < rwt[max_val as usize] {
                continue;
            }
            if rwt[cur_val as usize] == rwt[max_val as usize] && cur_val < max_val {
                continue;
            }
            max_val = cur_val;
        }
    }

    max_val
}

pub fn shine_dalgarno_mm(seq: &[u8], pos: i32, start: i32, rwt: &[f64]) -> i32 {
    let limit = imin(6, start - 4 - pos);
    let mut match_scores = [-10.0; 6];

    for i in 0..limit {
        if pos + i >= 0 {
            if i % 3 == 0 {

                match_scores[i as usize] = if is_a(seq, pos + i) { 2.0 } else { -3.0 };
            } else {

                match_scores[i as usize] = if is_g(seq, pos + i) { 3.0 } else { -2.0 };
            }
        }
    }

    let mut max_val = 0i32;

    for i in (5..=limit).rev() {
        for j in 0..=(limit - i) {
            let mut cur_ctr = -2.0;
            let mut mism = 0;

            for k in j..(j + i) {
                cur_ctr += match_scores[k as usize];
                if match_scores[k as usize] < 0.0 {
                    mism += 1;
                }

                if match_scores[k as usize] < 0.0 && (k <= j + 1 || k >= j + i - 2) {
                    cur_ctr -= 10.0;
                }
            }

            if mism != 1 {
                continue;
            }

            let rdis = start - (pos + j + i);
            let dis_flag = if rdis < 5 {
                1
            } else if rdis > 10 && rdis <= 12 {
                2
            } else if rdis >= 13 {
                3
            } else {
                0
            };

            if rdis > 15 || cur_ctr < 6.0 {
                continue;
            }

            let cur_val = if cur_ctr < 6.0 {
                0
            } else if cur_ctr == 6.0 && dis_flag == 3 {
                2
            } else if cur_ctr == 7.0 && dis_flag == 3 {
                2
            } else if cur_ctr == 9.0 && dis_flag == 3 {
                3
            } else if cur_ctr == 6.0 && dis_flag == 2 {
                4
            } else if cur_ctr == 6.0 && dis_flag == 1 {
                5
            } else if cur_ctr == 6.0 && dis_flag == 0 {
                9
            } else if cur_ctr == 7.0 && dis_flag == 2 {
                7
            } else if cur_ctr == 7.0 && dis_flag == 1 {
                8
            } else if cur_ctr == 7.0 && dis_flag == 0 {
                14
            } else if cur_ctr == 9.0 && dis_flag == 2 {
                17
            } else if cur_ctr == 9.0 && dis_flag == 1 {
                18
            } else if cur_ctr == 9.0 && dis_flag == 0 {
                19
            } else {
                0
            };

            if cur_val >= rwt.len() as i32 || max_val >= rwt.len() as i32 {
                continue;
            }

            if rwt[cur_val as usize] < rwt[max_val as usize] {
                continue;
            }
            if rwt[cur_val as usize] == rwt[max_val as usize] && cur_val < max_val {
                continue;
            }
            max_val = cur_val;
        }
    }

    max_val
}

pub fn imin(a: i32, b: i32) -> i32 {
    a.min(b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bitmap;

    #[test]
    fn test_nucleotide_functions() {
        let mut seq = vec![0u8; 100];

        let mut bctr = 0;

        bctr += 2;

        bitmap::set(&mut seq, bctr);
        bitmap::set(&mut seq, bctr + 1);
        bctr += 2;

        bitmap::set(&mut seq, bctr);
        let _ = bctr;

        assert!(is_a(&seq, 0));
        assert!(is_t(&seq, 1));
        assert!(is_g(&seq, 2));
        assert!(is_atg(&seq, 0));
    }

    #[test]
    fn test_gc_content() {
        let mut seq = vec![0u8; 100];

        bitmap::set(&mut seq, 3);
        bitmap::set(&mut seq, 4);
        bitmap::set(&mut seq, 6);
        bitmap::set(&mut seq, 7);

        let gc = gc_content(&seq, 0, 3);
        assert!((gc - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_mer_ndx() {
        let mut seq = vec![0u8; 100];

        let mut bctr = 0;

        bctr += 2;

        bitmap::set(&mut seq, bctr + 1);
        bctr += 2;

        bitmap::set(&mut seq, bctr);
        bctr += 2;

        bitmap::set(&mut seq, bctr);
        bitmap::set(&mut seq, bctr + 1);
        bctr += 2;

        bctr += 2;

        bitmap::set(&mut seq, bctr + 1);
        let _ = bctr;

        let idx = mer_ndx(6, &seq, 0);

        assert!(idx >= 0 && idx < 4096);
    }

    #[test]
    fn test_is_gc() {
        let mut seq = vec![0u8; 100];

        bitmap::set(&mut seq, 0);
        assert!(is_gc(&seq, 0));

        let mut seq2 = vec![0u8; 100];
        bitmap::set(&mut seq2, 1);
        assert!(is_gc(&seq2, 0));

        let seq3 = vec![0u8; 100];
        assert!(!is_gc(&seq3, 0));
    }

    #[test]
    fn test_ascii_to_nucleotide() {

        assert_eq!(ascii_to_nucleotide(b'A'), Nucleotide::A_VAL);
        assert_eq!(ascii_to_nucleotide(b'a'), Nucleotide::A_VAL);
        assert_eq!(ascii_to_nucleotide(b'G'), Nucleotide::G_VAL);
        assert_eq!(ascii_to_nucleotide(b'g'), Nucleotide::G_VAL);
        assert_eq!(ascii_to_nucleotide(b'C'), Nucleotide::C_VAL);
        assert_eq!(ascii_to_nucleotide(b'c'), Nucleotide::C_VAL);
        assert_eq!(ascii_to_nucleotide(b'T'), Nucleotide::T_VAL);
        assert_eq!(ascii_to_nucleotide(b't'), Nucleotide::T_VAL);

        assert_eq!(ascii_to_nucleotide(b'N'), Nucleotide::N_VAL);
        assert_eq!(ascii_to_nucleotide(b'n'), Nucleotide::N_VAL);
        assert_eq!(ascii_to_nucleotide(b'X'), Nucleotide::N_VAL);
        assert_eq!(ascii_to_nucleotide(b' '), Nucleotide::N_VAL);
    }

    #[test]
    fn test_complement_xor() {

        assert_eq!(complement(Nucleotide::A_VAL), Nucleotide::T_VAL);

        assert_eq!(complement(Nucleotide::T_VAL), Nucleotide::A_VAL);

        assert_eq!(complement(Nucleotide::G_VAL), Nucleotide::C_VAL);

        assert_eq!(complement(Nucleotide::C_VAL), Nucleotide::G_VAL);

        assert_eq!(complement(complement(Nucleotide::A_VAL)), Nucleotide::A_VAL);
        assert_eq!(complement(complement(Nucleotide::G_VAL)), Nucleotide::G_VAL);
    }

    #[test]
    fn test_is_stop_fast_standard_table() {

        let taa: Vec<u8> = vec![Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::A_VAL];
        assert!(is_stop_fast(&taa, 3, 0, 11, 1));

        let tag: Vec<u8> = vec![Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::G_VAL];
        assert!(is_stop_fast(&tag, 3, 0, 11, 1));

        let tga: Vec<u8> = vec![Nucleotide::T_VAL, Nucleotide::G_VAL, Nucleotide::A_VAL];
        assert!(is_stop_fast(&tga, 3, 0, 11, 1));

        let atg: Vec<u8> = vec![Nucleotide::A_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL];
        assert!(!is_stop_fast(&atg, 3, 0, 11, 1));
    }

    #[test]
    fn test_is_stop_fast_mycoplasma_table() {

        let taa: Vec<u8> = vec![Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::A_VAL];
        assert!(is_stop_fast(&taa, 3, 0, 4, 1));

        let tag: Vec<u8> = vec![Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::G_VAL];
        assert!(is_stop_fast(&tag, 3, 0, 4, 1));

        let tga: Vec<u8> = vec![Nucleotide::T_VAL, Nucleotide::G_VAL, Nucleotide::A_VAL];
        assert!(!is_stop_fast(&tga, 3, 0, 4, 1));
    }

    #[test]
    fn test_is_start_fast() {

        let atg: Vec<u8> = vec![Nucleotide::A_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL];
        assert!(is_start_fast(&atg, 3, 0, 11, 1));
        assert!(is_start_fast(&atg, 3, 0, 4, 1));
        assert!(is_start_fast(&atg, 3, 0, 1, 1));

        let gtg: Vec<u8> = vec![Nucleotide::G_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL];
        assert!(is_start_fast(&gtg, 3, 0, 11, 1));
        assert!(!is_start_fast(&gtg, 3, 0, 1, 1));

        let ttg: Vec<u8> = vec![Nucleotide::T_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL];
        assert!(is_start_fast(&ttg, 3, 0, 11, 1));
        assert!(!is_start_fast(&ttg, 3, 0, 1, 1));
    }

    #[test]
    fn test_start_type_fast() {
        let atg: Vec<u8> = vec![Nucleotide::A_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL];
        assert_eq!(start_type_fast(&atg, 3, 0, 1), ATG);

        let gtg: Vec<u8> = vec![Nucleotide::G_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL];
        assert_eq!(start_type_fast(&gtg, 3, 0, 1), GTG);

        let ttg: Vec<u8> = vec![Nucleotide::T_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL];
        assert_eq!(start_type_fast(&ttg, 3, 0, 1), TTG);

        let taa: Vec<u8> = vec![Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::A_VAL];
        assert_eq!(start_type_fast(&taa, 3, 0, 1), -1);
    }

    #[test]
    fn test_reverse_strand_stop_detection() {

        let seq: Vec<u8> = vec![Nucleotide::T_VAL, Nucleotide::T_VAL, Nucleotide::A_VAL];

        assert!(is_stop_fast(&seq, 3, 0, 11, -1));
    }

    #[test]
    fn test_bitmap_to_digits_conversion() {

        let mut seq = vec![0u8; 100];
        let useq = vec![0u8; 100];

        bitmap::set(&mut seq, 2);
        bitmap::set(&mut seq, 3);

        bitmap::set(&mut seq, 4);

        bitmap::set(&mut seq, 6);
        bitmap::set(&mut seq, 7);

        let mut digits = vec![0u8; 6];
        bitmap_to_digits(&seq, &useq, 6, &mut digits);

        assert_eq!(digits[0], Nucleotide::A_VAL);
        assert_eq!(digits[1], Nucleotide::T_VAL);
        assert_eq!(digits[2], Nucleotide::G_VAL);
        assert_eq!(digits[3], Nucleotide::T_VAL);
        assert_eq!(digits[4], Nucleotide::A_VAL);
        assert_eq!(digits[5], Nucleotide::A_VAL);

        assert!(is_start_fast(&digits, 6, 0, 11, 1));
        assert!(is_stop_fast(&digits, 6, 3, 11, 1));
    }

    #[test]
    fn test_is_gc_fast() {
        let digits: Vec<u8> = vec![
            Nucleotide::A_VAL,
            Nucleotide::G_VAL,
            Nucleotide::C_VAL,
            Nucleotide::T_VAL,
        ];

        assert!(!is_gc_fast(&digits, 0));
        assert!(is_gc_fast(&digits, 1));
        assert!(is_gc_fast(&digits, 2));
        assert!(!is_gc_fast(&digits, 3));
    }

    #[test]
    fn test_equivalence_with_original_is_stop() {

        let test_codons: [(u8, u8, u8); 4] = [
            (Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::A_VAL),
            (Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::G_VAL),
            (Nucleotide::T_VAL, Nucleotide::G_VAL, Nucleotide::A_VAL),
            (Nucleotide::A_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL),
        ];

        for (x0, x1, x2) in test_codons.iter() {
            let mut seq = vec![0u8; 100];
            let mut bctr = 0;

            for &nuc in [*x0, *x1, *x2].iter() {
                let (bit0, bit1) = match nuc {
                    0b000 => (0, 0),
                    0b001 => (1, 0),
                    0b010 => (0, 1),
                    0b011 => (1, 1),
                    _ => (0, 1),
                };
                if bit0 == 1 { bitmap::set(&mut seq, bctr); }
                if bit1 == 1 { bitmap::set(&mut seq, bctr + 1); }
                bctr += 2;
            }

            let digits: Vec<u8> = vec![*x0, *x1, *x2];

            for tt in [4, 11].iter() {
                let tinf = Training { trans_table: *tt, ..Default::default() };

                let orig_result = is_stop(&seq, 0, &tinf);
                let fast_result = is_stop_fast(&digits, 3, 0, *tt as usize, 1);

                assert_eq!(
                    orig_result, fast_result,
                    "Mismatch for codon ({:?},{:?},{:?}) table {}: orig={}, fast={}",
                    x0, x1, x2, tt, orig_result, fast_result
                );
            }
        }
    }

    #[test]
    fn test_equivalence_with_original_is_start() {

        let test_codons: [(u8, u8, u8); 4] = [
            (Nucleotide::A_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL),
            (Nucleotide::G_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL),
            (Nucleotide::T_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL),
            (Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::A_VAL),
        ];

        for (x0, x1, x2) in test_codons.iter() {
            let mut seq = vec![0u8; 100];
            let mut bctr = 0;

            for &nuc in [*x0, *x1, *x2].iter() {
                let (bit0, bit1) = match nuc {
                    0b000 => (0, 0),
                    0b001 => (1, 0),
                    0b010 => (0, 1),
                    0b011 => (1, 1),
                    _ => (0, 1),
                };
                if bit0 == 1 { bitmap::set(&mut seq, bctr); }
                if bit1 == 1 { bitmap::set(&mut seq, bctr + 1); }
                bctr += 2;
            }

            let digits: Vec<u8> = vec![*x0, *x1, *x2];

            for tt in [1, 4, 11].iter() {
                let tinf = Training { trans_table: *tt, ..Default::default() };

                let orig_result = is_start(&seq, 0, &tinf);
                let fast_result = is_start_fast(&digits, 3, 0, *tt as usize, 1);

                assert_eq!(
                    orig_result, fast_result,
                    "Mismatch for codon ({:?},{:?},{:?}) table {}: orig={}, fast={}",
                    x0, x1, x2, tt, orig_result, fast_result
                );
            }
        }
    }

    fn create_bitmap_seq(nucs: &[u8]) -> Vec<u8> {
        let mut seq = vec![0u8; 100];
        let mut bctr = 0i32;

        for &nuc in nucs {
            let (bit0, bit1) = match nuc {
                0b000 => (0, 0),
                0b001 => (1, 0),
                0b010 => (0, 1),
                0b011 => (1, 1),
                _ => (0, 1),
            };
            if bit0 == 1 { bitmap::set(&mut seq, bctr); }
            if bit1 == 1 { bitmap::set(&mut seq, bctr + 1); }
            bctr += 2;
        }

        seq
    }

    #[test]
    fn test_amino_atg() {

        let seq = create_bitmap_seq(&[
            Nucleotide::A_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL
        ]);

        let mut tinf = Box::new(Training::default());
        tinf.trans_table = 11;

        assert_eq!(amino(&seq, 0, &tinf, 0), b'M');
        assert_eq!(amino(&seq, 0, &tinf, 1), b'M');
    }

    #[test]
    fn test_amino_stop_codons() {
        let mut tinf = Box::new(Training::default());
        tinf.trans_table = 11;

        let taa = create_bitmap_seq(&[
            Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::A_VAL
        ]);
        assert_eq!(amino(&taa, 0, &tinf, 0), b'*');

        let tag = create_bitmap_seq(&[
            Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::G_VAL
        ]);
        assert_eq!(amino(&tag, 0, &tinf, 0), b'*');

        let tga = create_bitmap_seq(&[
            Nucleotide::T_VAL, Nucleotide::G_VAL, Nucleotide::A_VAL
        ]);
        assert_eq!(amino(&tga, 0, &tinf, 0), b'*');
    }

    #[test]
    fn test_amino_alternative_starts() {
        let mut tinf = Box::new(Training::default());
        tinf.trans_table = 11;

        let gtg = create_bitmap_seq(&[
            Nucleotide::G_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL
        ]);
        assert_eq!(amino(&gtg, 0, &tinf, 0), b'V');
        assert_eq!(amino(&gtg, 0, &tinf, 1), b'M');

        let ttg = create_bitmap_seq(&[
            Nucleotide::T_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL
        ]);
        assert_eq!(amino(&ttg, 0, &tinf, 0), b'L');
        assert_eq!(amino(&ttg, 0, &tinf, 1), b'M');
    }

    #[test]
    fn test_amino_common_codons() {
        let mut tinf = Box::new(Training::default());
        tinf.trans_table = 11;

        let gct = create_bitmap_seq(&[
            Nucleotide::G_VAL, Nucleotide::C_VAL, Nucleotide::T_VAL
        ]);
        assert_eq!(amino(&gct, 0, &tinf, 0), b'A');

        let ttt = create_bitmap_seq(&[
            Nucleotide::T_VAL, Nucleotide::T_VAL, Nucleotide::T_VAL
        ]);
        assert_eq!(amino(&ttt, 0, &tinf, 0), b'F');

        let aaa = create_bitmap_seq(&[
            Nucleotide::A_VAL, Nucleotide::A_VAL, Nucleotide::A_VAL
        ]);
        assert_eq!(amino(&aaa, 0, &tinf, 0), b'K');

        let ggg = create_bitmap_seq(&[
            Nucleotide::G_VAL, Nucleotide::G_VAL, Nucleotide::G_VAL
        ]);
        assert_eq!(amino(&ggg, 0, &tinf, 0), b'G');
    }

    #[test]
    fn test_amino_table_differences() {

        let tga = create_bitmap_seq(&[
            Nucleotide::T_VAL, Nucleotide::G_VAL, Nucleotide::A_VAL
        ]);

        let mut tinf_11 = Box::new(Training::default());
        tinf_11.trans_table = 11;

        let mut tinf_4 = Box::new(Training::default());
        tinf_4.trans_table = 4;

        assert_eq!(amino(&tga, 0, &tinf_11, 0), b'*');
        assert_eq!(amino(&tga, 0, &tinf_4, 0), b'W');
    }

    #[test]
    fn test_amino_fast() {

        let digits = vec![
            Nucleotide::A_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL,
            Nucleotide::G_VAL, Nucleotide::C_VAL, Nucleotide::T_VAL,
            Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::A_VAL,
        ];

        assert_eq!(amino_fast(&digits, 9, 0, 11, 1, false), b'M');
        assert_eq!(amino_fast(&digits, 9, 3, 11, 1, false), b'A');
        assert_eq!(amino_fast(&digits, 9, 6, 11, 1, false), b'*');
    }

    #[test]
    fn test_get_nucleotide_3bit() {

        let seq = create_bitmap_seq(&[
            Nucleotide::A_VAL, Nucleotide::G_VAL, Nucleotide::C_VAL, Nucleotide::T_VAL
        ]);

        assert_eq!(get_nucleotide_3bit(&seq, 0), Nucleotide::A_VAL);
        assert_eq!(get_nucleotide_3bit(&seq, 1), Nucleotide::G_VAL);
        assert_eq!(get_nucleotide_3bit(&seq, 2), Nucleotide::C_VAL);
        assert_eq!(get_nucleotide_3bit(&seq, 3), Nucleotide::T_VAL);
    }
}

// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2007-2016 Doug Hyatt, Univ. of Tennessee / UT-Battelle (original Prodigal C)
// Copyright (C) 2026 Sunju Kim (Rust reimplementation)

#[inline(always)]
pub fn test(bitmap: &[u8], pos: i32) -> u8 {
    let ndx = pos as usize;
    if (bitmap[ndx >> 3] & (1 << (ndx & 0x07))) != 0 {
        1
    } else {
        0
    }
}

#[inline(always)]
pub fn clear(bitmap: &mut [u8], pos: i32) {
    let ndx = pos as usize;
    bitmap[ndx >> 3] &= !(1 << (ndx & 0x07));
}

#[inline(always)]
pub fn set(bitmap: &mut [u8], pos: i32) {
    let ndx = pos as usize;
    bitmap[ndx >> 3] |= 1 << (ndx & 0x07);
}

#[inline(always)]
pub fn toggle(bitmap: &mut [u8], pos: i32) {
    let ndx = pos as usize;
    bitmap[ndx >> 3] ^= 1 << (ndx & 0x07);
}

pub fn build_bitmaps(
    sequence: &[u8],
    slen: i32,
    seq: &mut [u8],
    rseq: &mut [u8],
    useq: &mut [u8],
    gc: &mut f64,
) {
    let mut bctr = 0i32;
    let mut gc_count = 0i32;

    let bitmap_bytes = (slen as usize * 2 + 7) / 8;
    seq[..bitmap_bytes].fill(0);
    useq[..(slen as usize + 7) / 8].fill(0);

    for (i, &ch) in sequence.iter().enumerate().take(slen as usize) {

        match ch {
            b'g' | b'G' => {

                set(seq, bctr);
                gc_count += 1;
            }
            b't' | b'T' => {

                set(seq, bctr);
                set(seq, bctr + 1);
            }
            b'c' | b'C' => {

                set(seq, bctr + 1);
                gc_count += 1;
            }
            b'a' | b'A' => {

            }
            _ => {

                set(seq, bctr + 1);
                set(useq, i as i32);
            }
        }
        bctr += 2;
    }

    *gc = gc_count as f64 / slen as f64;

    crate::sequence::rcom_seq(seq, rseq, useq, slen);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bitmap_operations() {

        let mut bitmap = vec![0u8; 4];

        assert_eq!(test(&bitmap, 0), 0);
        assert_eq!(test(&bitmap, 7), 0);
        assert_eq!(test(&bitmap, 15), 0);
        assert_eq!(test(&bitmap, 31), 0);

        set(&mut bitmap, 5);
        assert_eq!(test(&bitmap, 5), 1);
        assert_eq!(bitmap[0], 0b00100000);

        set(&mut bitmap, 13);
        assert_eq!(test(&bitmap, 13), 1);
        assert_eq!(bitmap[1], 0b00100000);

        clear(&mut bitmap, 5);
        assert_eq!(test(&bitmap, 5), 0);
        assert_eq!(bitmap[0], 0);

        toggle(&mut bitmap, 7);
        assert_eq!(test(&bitmap, 7), 1);
        assert_eq!(bitmap[0], 0b10000000);

        toggle(&mut bitmap, 7);
        assert_eq!(test(&bitmap, 7), 0);
        assert_eq!(bitmap[0], 0);
    }

    #[test]
    fn test_bitmap_byte_boundaries() {
        let mut bitmap = vec![0u8; 4];

        set(&mut bitmap, 0);
        set(&mut bitmap, 7);
        set(&mut bitmap, 8);
        set(&mut bitmap, 15);

        assert_eq!(test(&bitmap, 0), 1);
        assert_eq!(test(&bitmap, 7), 1);
        assert_eq!(test(&bitmap, 8), 1);
        assert_eq!(test(&bitmap, 15), 1);

        assert_eq!(bitmap[0], 0b10000001);
        assert_eq!(bitmap[1], 0b10000001);
    }

    #[test]
    fn test_bitmap_multiple_operations() {
        let mut bitmap = vec![0u8; 2];

        set(&mut bitmap, 3);
        set(&mut bitmap, 5);
        set(&mut bitmap, 11);

        assert_eq!(test(&bitmap, 3), 1);
        assert_eq!(test(&bitmap, 5), 1);
        assert_eq!(test(&bitmap, 11), 1);
        assert_eq!(test(&bitmap, 4), 0);

        clear(&mut bitmap, 5);
        assert_eq!(test(&bitmap, 5), 0);
        assert_eq!(test(&bitmap, 3), 1);

        toggle(&mut bitmap, 3);
        toggle(&mut bitmap, 4);
        assert_eq!(test(&bitmap, 3), 0);
        assert_eq!(test(&bitmap, 4), 1);
    }
}

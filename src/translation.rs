// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2007-2016 Doug Hyatt, Univ. of Tennessee / UT-Battelle (original Prodigal C)
// Copyright (C) 2026 Sunju Kim (Rust reimplementation)

use crate::sequence::Nucleotide;

#[inline(always)]
pub fn codon_index(x0: u8, x1: u8, x2: u8) -> usize {
    ((x0 as usize) << 4) + ((x1 as usize) << 2) + (x2 as usize)
}

#[allow(dead_code)]
#[inline(always)]
fn to_ncbi_index(x0: u8, x1: u8, x2: u8) -> usize {

    const MAP: [usize; 4] = [2, 3, 1, 0];
    let n0 = MAP[(x0 & 0b11) as usize];
    let n1 = MAP[(x1 & 0b11) as usize];
    let n2 = MAP[(x2 & 0b11) as usize];
    (n0 * 16) + (n1 * 4) + n2
}

const fn build_table_from_ncbi(ncbi_str: &[u8; 64]) -> [u8; 64] {
    let mut table = [0u8; 64];

    const OUR_TO_NCBI: [usize; 4] = [2, 3, 1, 0];

    let mut our_idx = 0;
    while our_idx < 64 {

        let our_x0 = (our_idx >> 4) & 0x3;
        let our_x1 = (our_idx >> 2) & 0x3;
        let our_x2 = our_idx & 0x3;

        let ncbi_x0 = OUR_TO_NCBI[our_x0];
        let ncbi_x1 = OUR_TO_NCBI[our_x1];
        let ncbi_x2 = OUR_TO_NCBI[our_x2];

        let ncbi_idx = (ncbi_x0 * 16) + (ncbi_x1 * 4) + ncbi_x2;

        table[our_idx] = ncbi_str[ncbi_idx];

        our_idx += 1;
    }

    table
}

const NCBI_TABLE_1: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

const NCBI_TABLE_2: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";

const NCBI_TABLE_3: &[u8; 64] = b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

const NCBI_TABLE_4: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

const NCBI_TABLE_5: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG";

const NCBI_TABLE_6: &[u8; 64] = b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

const NCBI_TABLE_9: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";

const NCBI_TABLE_10: &[u8; 64] = b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

const NCBI_TABLE_11: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

const NCBI_TABLE_12: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

const NCBI_TABLE_13: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG";

const NCBI_TABLE_14: &[u8; 64] = b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";

const NCBI_TABLE_15: &[u8; 64] = b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

const NCBI_TABLE_16: &[u8; 64] = b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

const NCBI_TABLE_21: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG";

const NCBI_TABLE_22: &[u8; 64] = b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

const NCBI_TABLE_23: &[u8; 64] = b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

const NCBI_TABLE_24: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG";

const NCBI_TABLE_25: &[u8; 64] = b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

const TABLE_1: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_1);
const TABLE_2: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_2);
const TABLE_3: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_3);
const TABLE_4: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_4);
const TABLE_5: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_5);
const TABLE_6: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_6);
const TABLE_9: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_9);
const TABLE_10: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_10);
const TABLE_11: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_11);
const TABLE_12: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_12);
const TABLE_13: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_13);
const TABLE_14: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_14);
const TABLE_15: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_15);
const TABLE_16: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_16);
const TABLE_21: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_21);
const TABLE_22: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_22);
const TABLE_23: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_23);
const TABLE_24: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_24);
const TABLE_25: [u8; 64] = build_table_from_ncbi(NCBI_TABLE_25);

const TABLE_DEFAULT: [u8; 64] = TABLE_1;

pub const TRANSLATION_TABLES: [[u8; 64]; 34] = [
    TABLE_DEFAULT,
    TABLE_1,
    TABLE_2,
    TABLE_3,
    TABLE_4,
    TABLE_5,
    TABLE_6,
    TABLE_DEFAULT,
    TABLE_DEFAULT,
    TABLE_9,
    TABLE_10,
    TABLE_11,
    TABLE_12,
    TABLE_13,
    TABLE_14,
    TABLE_15,
    TABLE_16,
    TABLE_DEFAULT,
    TABLE_DEFAULT,
    TABLE_DEFAULT,
    TABLE_DEFAULT,
    TABLE_21,
    TABLE_22,
    TABLE_23,
    TABLE_24,
    TABLE_25,
    TABLE_DEFAULT,
    TABLE_DEFAULT,
    TABLE_DEFAULT,
    TABLE_DEFAULT,
    TABLE_DEFAULT,
    TABLE_DEFAULT,
    TABLE_DEFAULT,
    TABLE_DEFAULT,
];

#[inline(always)]
pub fn translate_codon(tt: usize, x0: u8, x1: u8, x2: u8) -> u8 {
    let tt = if tt < 34 { tt } else { 11 };
    let idx = codon_index(x0, x1, x2);
    TRANSLATION_TABLES[tt][idx]
}

#[inline(always)]
pub fn translate_codon_at(
    digits: &[u8],
    slen: usize,
    pos: usize,
    tt: usize,
    strand: i32,
) -> u8 {
    let (x0, x1, x2) = if strand == 1 {

        (digits[pos], digits[pos + 1], digits[pos + 2])
    } else {

        (
            crate::sequence::complement(digits[slen - 1 - pos]),
            crate::sequence::complement(digits[slen - 2 - pos]),
            crate::sequence::complement(digits[slen - 3 - pos]),
        )
    };

    translate_codon(tt, x0, x1, x2)
}

#[inline(always)]
pub fn translate_codon_init(
    digits: &[u8],
    slen: usize,
    pos: usize,
    tt: usize,
    strand: i32,
    is_init: bool,
) -> u8 {
    let (x0, x1, x2) = if strand == 1 {
        (digits[pos], digits[pos + 1], digits[pos + 2])
    } else {
        (
            crate::sequence::complement(digits[slen - 1 - pos]),
            crate::sequence::complement(digits[slen - 2 - pos]),
            crate::sequence::complement(digits[slen - 3 - pos]),
        )
    };

    let aa = translate_codon(tt, x0, x1, x2);
    if aa == b'*' {
        return b'*';
    }

    if is_init {

        if x2 == Nucleotide::G_VAL && x1 == Nucleotide::T_VAL {
            if x0 == Nucleotide::A_VAL || x0 == Nucleotide::G_VAL || x0 == Nucleotide::T_VAL {

                return b'M';
            }
        }
    }

    aa
}

#[inline(always)]
pub fn is_stop_codon(tt: usize, x0: u8, x1: u8, x2: u8) -> bool {
    translate_codon(tt, x0, x1, x2) == b'*'
}

pub fn translate_codon_str(tt: usize, codon: &str) -> char {
    if codon.len() != 3 {
        return 'X';
    }

    let bytes = codon.as_bytes();
    let x0 = crate::sequence::ascii_to_nucleotide(bytes[0]);
    let x1 = crate::sequence::ascii_to_nucleotide(bytes[1]);
    let x2 = crate::sequence::ascii_to_nucleotide(bytes[2]);

    translate_codon(tt, x0, x1, x2) as char
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::Nucleotide;

    #[test]
    fn test_codon_index() {

        assert_eq!(codon_index(0, 0, 0), 0);

        assert_eq!(codon_index(0, 0, 1), 1);

        assert_eq!(codon_index(3, 3, 3), 63);

        assert_eq!(codon_index(
            Nucleotide::A_VAL,
            Nucleotide::T_VAL,
            Nucleotide::G_VAL
        ), 13);
    }

    #[test]
    fn test_translate_atg() {

        let x0 = Nucleotide::A_VAL;
        let x1 = Nucleotide::T_VAL;
        let x2 = Nucleotide::G_VAL;

        assert_eq!(translate_codon(1, x0, x1, x2), b'M');
        assert_eq!(translate_codon(4, x0, x1, x2), b'M');
        assert_eq!(translate_codon(11, x0, x1, x2), b'M');
    }

    #[test]
    fn test_translate_stop_codons_table_11() {

        assert_eq!(translate_codon(11,
            Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::A_VAL
        ), b'*');

        assert_eq!(translate_codon(11,
            Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::G_VAL
        ), b'*');

        assert_eq!(translate_codon(11,
            Nucleotide::T_VAL, Nucleotide::G_VAL, Nucleotide::A_VAL
        ), b'*');
    }

    #[test]
    fn test_translate_stop_codons_table_4() {

        assert_eq!(translate_codon(4,
            Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::A_VAL
        ), b'*');

        assert_eq!(translate_codon(4,
            Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::G_VAL
        ), b'*');

        assert_eq!(translate_codon(4,
            Nucleotide::T_VAL, Nucleotide::G_VAL, Nucleotide::A_VAL
        ), b'W');
    }

    #[test]
    fn test_translate_common_codons() {

        assert_eq!(translate_codon(11,
            Nucleotide::G_VAL, Nucleotide::C_VAL, Nucleotide::T_VAL
        ), b'A');
        assert_eq!(translate_codon(11,
            Nucleotide::G_VAL, Nucleotide::C_VAL, Nucleotide::C_VAL
        ), b'A');
        assert_eq!(translate_codon(11,
            Nucleotide::G_VAL, Nucleotide::C_VAL, Nucleotide::A_VAL
        ), b'A');
        assert_eq!(translate_codon(11,
            Nucleotide::G_VAL, Nucleotide::C_VAL, Nucleotide::G_VAL
        ), b'A');

        assert_eq!(translate_codon(11,
            Nucleotide::T_VAL, Nucleotide::T_VAL, Nucleotide::T_VAL
        ), b'F');
        assert_eq!(translate_codon(11,
            Nucleotide::T_VAL, Nucleotide::T_VAL, Nucleotide::C_VAL
        ), b'F');

        assert_eq!(translate_codon(11,
            Nucleotide::A_VAL, Nucleotide::A_VAL, Nucleotide::A_VAL
        ), b'K');
        assert_eq!(translate_codon(11,
            Nucleotide::A_VAL, Nucleotide::A_VAL, Nucleotide::G_VAL
        ), b'K');

        assert_eq!(translate_codon(11,
            Nucleotide::G_VAL, Nucleotide::G_VAL, Nucleotide::T_VAL
        ), b'G');
        assert_eq!(translate_codon(11,
            Nucleotide::G_VAL, Nucleotide::G_VAL, Nucleotide::C_VAL
        ), b'G');
        assert_eq!(translate_codon(11,
            Nucleotide::G_VAL, Nucleotide::G_VAL, Nucleotide::A_VAL
        ), b'G');
        assert_eq!(translate_codon(11,
            Nucleotide::G_VAL, Nucleotide::G_VAL, Nucleotide::G_VAL
        ), b'G');
    }

    #[test]
    fn test_translate_codon_str() {

        assert_eq!(translate_codon_str(11, "ATG"), 'M');
        assert_eq!(translate_codon_str(11, "TAA"), '*');
        assert_eq!(translate_codon_str(11, "TAG"), '*');
        assert_eq!(translate_codon_str(11, "TGA"), '*');
        assert_eq!(translate_codon_str(11, "GCT"), 'A');
        assert_eq!(translate_codon_str(11, "TTT"), 'F');
        assert_eq!(translate_codon_str(11, "AAA"), 'K');
        assert_eq!(translate_codon_str(11, "GGG"), 'G');

        assert_eq!(translate_codon_str(11, "atg"), 'M');
        assert_eq!(translate_codon_str(11, "taa"), '*');
    }

    #[test]
    fn test_translate_codon_init() {

        let atg = vec![
            Nucleotide::A_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL
        ];
        let gtg = vec![
            Nucleotide::G_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL
        ];
        let ttg = vec![
            Nucleotide::T_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL
        ];

        assert_eq!(translate_codon_init(&atg, 3, 0, 11, 1, true), b'M');
        assert_eq!(translate_codon_init(&atg, 3, 0, 11, 1, false), b'M');

        assert_eq!(translate_codon_init(&gtg, 3, 0, 11, 1, true), b'M');
        assert_eq!(translate_codon_init(&gtg, 3, 0, 11, 1, false), b'V');

        assert_eq!(translate_codon_init(&ttg, 3, 0, 11, 1, true), b'M');
        assert_eq!(translate_codon_init(&ttg, 3, 0, 11, 1, false), b'L');
    }

    #[test]
    fn test_is_stop_codon() {

        assert!(is_stop_codon(11,
            Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::A_VAL
        ));
        assert!(is_stop_codon(11,
            Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::G_VAL
        ));
        assert!(is_stop_codon(11,
            Nucleotide::T_VAL, Nucleotide::G_VAL, Nucleotide::A_VAL
        ));

        assert!(!is_stop_codon(11,
            Nucleotide::A_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL
        ));

        assert!(!is_stop_codon(4,
            Nucleotide::T_VAL, Nucleotide::G_VAL, Nucleotide::A_VAL
        ));
    }

    #[test]
    fn test_all_64_codons_table_11() {

        let expected = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

        let ncbi_nucs: [u8; 4] = [3, 2, 0, 1];

        let mut idx = 0;
        for &n0 in &ncbi_nucs {
            for &n1 in &ncbi_nucs {
                for &n2 in &ncbi_nucs {
                    let aa = translate_codon(11, n0, n1, n2);
                    assert_eq!(aa, expected[idx],
                        "Codon {} (idx={}) expected '{}' got '{}'",
                        idx, idx, expected[idx] as char, aa as char);
                    idx += 1;
                }
            }
        }
    }

    #[test]
    fn test_all_64_codons_table_1() {

        let expected = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

        let ncbi_nucs: [u8; 4] = [3, 2, 0, 1];

        let mut idx = 0;
        for &n0 in &ncbi_nucs {
            for &n1 in &ncbi_nucs {
                for &n2 in &ncbi_nucs {
                    let aa = translate_codon(1, n0, n1, n2);
                    assert_eq!(aa, expected[idx],
                        "Codon {} (idx={}) expected '{}' got '{}'",
                        idx, idx, expected[idx] as char, aa as char);
                    idx += 1;
                }
            }
        }
    }

    #[test]
    fn test_all_64_codons_table_4() {

        let expected = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

        let ncbi_nucs: [u8; 4] = [3, 2, 0, 1];

        let mut idx = 0;
        for &n0 in &ncbi_nucs {
            for &n1 in &ncbi_nucs {
                for &n2 in &ncbi_nucs {
                    let aa = translate_codon(4, n0, n1, n2);
                    assert_eq!(aa, expected[idx],
                        "Table 4 codon {} (idx={}) expected '{}' got '{}'",
                        idx, idx, expected[idx] as char, aa as char);
                    idx += 1;
                }
            }
        }
    }

    #[test]
    fn test_translate_codon_at_forward() {

        let seq = vec![
            Nucleotide::A_VAL, Nucleotide::T_VAL, Nucleotide::G_VAL,
            Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::A_VAL,
        ];

        assert_eq!(translate_codon_at(&seq, 6, 0, 11, 1), b'M');
        assert_eq!(translate_codon_at(&seq, 6, 3, 11, 1), b'*');
    }

    #[test]
    fn test_translate_codon_at_reverse() {

        let seq = vec![
            Nucleotide::C_VAL, Nucleotide::A_VAL, Nucleotide::T_VAL,
        ];

        assert_eq!(translate_codon_at(&seq, 3, 0, 11, -1), b'M');
    }

    #[test]
    fn test_table_differences() {

        let aa_11 = translate_codon(11,
            Nucleotide::T_VAL, Nucleotide::G_VAL, Nucleotide::A_VAL
        );
        let aa_4 = translate_codon(4,
            Nucleotide::T_VAL, Nucleotide::G_VAL, Nucleotide::A_VAL
        );
        assert_eq!(aa_11, b'*');
        assert_eq!(aa_4, b'W');

        let aa_11 = translate_codon(11,
            Nucleotide::A_VAL, Nucleotide::G_VAL, Nucleotide::A_VAL
        );
        let aa_2 = translate_codon(2,
            Nucleotide::A_VAL, Nucleotide::G_VAL, Nucleotide::A_VAL
        );
        assert_eq!(aa_11, b'R');
        assert_eq!(aa_2, b'*');

        let aa_11 = translate_codon(11,
            Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::A_VAL
        );
        let aa_6 = translate_codon(6,
            Nucleotide::T_VAL, Nucleotide::A_VAL, Nucleotide::A_VAL
        );
        assert_eq!(aa_11, b'*');
        assert_eq!(aa_6, b'Q');
    }

    #[test]
    fn test_memory_size() {

        let table_size = std::mem::size_of::<[[u8; 64]; 34]>();
        assert_eq!(table_size, 34 * 64);
    }
}

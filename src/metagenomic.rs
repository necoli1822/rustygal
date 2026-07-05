// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2007-2016 Doug Hyatt, Univ. of Tennessee / UT-Battelle (original Prodigal C)
// Copyright (C) 2026 Sunju Kim (Rust reimplementation)

use crate::training::Training;

pub const NUM_BIN: usize = 6;
pub const NUM_META: usize = 50;
pub const SAMPLE_LEN: i32 = 120;
pub const MAX_SAMPLE: i32 = 200;

#[repr(C)]
#[derive(Debug, Clone)]
pub struct MetagenomicBin {

    pub index: i32,

    pub clusnum: i32,

    pub desc: [u8; 500],

    pub weight: f64,

    pub gc: f64,

    pub tinf: Box<Training>,
}

impl Default for MetagenomicBin {
    fn default() -> Self {
        Self {
            index: 0,
            clusnum: 0,
            desc: [0; 500],
            weight: 0.0,
            gc: 0.0,
            tinf: Box::new(Training::default()),
        }
    }
}

pub fn initialize_metagenomic_bins(meta: &mut [MetagenomicBin]) {
    use crate::training::*;

    initialize_metagenome_0(&mut meta[0].tinf);
    initialize_metagenome_1(&mut meta[1].tinf);
    initialize_metagenome_2(&mut meta[2].tinf);
    initialize_metagenome_3(&mut meta[3].tinf);
    initialize_metagenome_4(&mut meta[4].tinf);
    initialize_metagenome_5(&mut meta[5].tinf);
    initialize_metagenome_6(&mut meta[6].tinf);
    initialize_metagenome_7(&mut meta[7].tinf);
    initialize_metagenome_8(&mut meta[8].tinf);
    initialize_metagenome_9(&mut meta[9].tinf);
    initialize_metagenome_10(&mut meta[10].tinf);
    initialize_metagenome_11(&mut meta[11].tinf);
    initialize_metagenome_12(&mut meta[12].tinf);
    initialize_metagenome_13(&mut meta[13].tinf);
    initialize_metagenome_14(&mut meta[14].tinf);
    initialize_metagenome_15(&mut meta[15].tinf);
    initialize_metagenome_16(&mut meta[16].tinf);
    initialize_metagenome_17(&mut meta[17].tinf);
    initialize_metagenome_18(&mut meta[18].tinf);
    initialize_metagenome_19(&mut meta[19].tinf);
    initialize_metagenome_20(&mut meta[20].tinf);
    initialize_metagenome_21(&mut meta[21].tinf);
    initialize_metagenome_22(&mut meta[22].tinf);
    initialize_metagenome_23(&mut meta[23].tinf);
    initialize_metagenome_24(&mut meta[24].tinf);
    initialize_metagenome_25(&mut meta[25].tinf);
    initialize_metagenome_26(&mut meta[26].tinf);
    initialize_metagenome_27(&mut meta[27].tinf);
    initialize_metagenome_28(&mut meta[28].tinf);
    initialize_metagenome_29(&mut meta[29].tinf);
    initialize_metagenome_30(&mut meta[30].tinf);
    initialize_metagenome_31(&mut meta[31].tinf);
    initialize_metagenome_32(&mut meta[32].tinf);
    initialize_metagenome_33(&mut meta[33].tinf);
    initialize_metagenome_34(&mut meta[34].tinf);
    initialize_metagenome_35(&mut meta[35].tinf);
    initialize_metagenome_36(&mut meta[36].tinf);
    initialize_metagenome_37(&mut meta[37].tinf);
    initialize_metagenome_38(&mut meta[38].tinf);
    initialize_metagenome_39(&mut meta[39].tinf);
    initialize_metagenome_40(&mut meta[40].tinf);
    initialize_metagenome_41(&mut meta[41].tinf);
    initialize_metagenome_42(&mut meta[42].tinf);
    initialize_metagenome_43(&mut meta[43].tinf);
    initialize_metagenome_44(&mut meta[44].tinf);
    initialize_metagenome_45(&mut meta[45].tinf);
    initialize_metagenome_46(&mut meta[46].tinf);
    initialize_metagenome_47(&mut meta[47].tinf);
    initialize_metagenome_48(&mut meta[48].tinf);
    initialize_metagenome_49(&mut meta[49].tinf);

    write_desc(&mut meta[0].desc, 0, "Mycoplasma_bovis_PG45", "B", 29.31,
               meta[0].tinf.trans_table, meta[0].tinf.uses_sd);
    write_desc(&mut meta[1].desc, 1, "Mycoplasma_pneumoniae_M129", "B", 40.01,
               meta[1].tinf.trans_table, meta[1].tinf.uses_sd);
    write_desc(&mut meta[2].desc, 2, "Mycoplasma_suis_Illinois", "B", 31.08,
               meta[2].tinf.trans_table, meta[2].tinf.uses_sd);
    write_desc(&mut meta[3].desc, 3, "Aeropyrum_pernix_K1", "A", 56.31,
               meta[3].tinf.trans_table, meta[3].tinf.uses_sd);
    write_desc(&mut meta[4].desc, 4, "Akkermansia_muciniphila_ATCC_BAA_835", "B", 55.76,
               meta[4].tinf.trans_table, meta[4].tinf.uses_sd);
    write_desc(&mut meta[5].desc, 5, "Anaplasma_marginale_Maries", "B", 49.76,
               meta[5].tinf.trans_table, meta[5].tinf.uses_sd);
    write_desc(&mut meta[6].desc, 6, "Anaplasma_phagocytophilum_HZ", "B", 41.64,
               meta[6].tinf.trans_table, meta[6].tinf.uses_sd);
    write_desc(&mut meta[7].desc, 7, "Archaeoglobus_fulgidus_DSM_4304", "A", 48.58,
               meta[7].tinf.trans_table, meta[7].tinf.uses_sd);
    write_desc(&mut meta[8].desc, 8, "Bacteroides_fragilis_NCTC_9343", "B", 43.19,
               meta[8].tinf.trans_table, meta[8].tinf.uses_sd);
    write_desc(&mut meta[9].desc, 9, "Brucella_canis_ATCC_23365", "B", 57.21,
               meta[9].tinf.trans_table, meta[9].tinf.uses_sd);
    write_desc(&mut meta[10].desc, 10, "Burkholderia_rhizoxinica_HKI_454", "B", 59.70,
               meta[10].tinf.trans_table, meta[10].tinf.uses_sd);
    write_desc(&mut meta[11].desc, 11, "Candidatus_Amoebophilus_asiaticus_5a2", "B", 35.05,
               meta[11].tinf.trans_table, meta[11].tinf.uses_sd);
    write_desc(&mut meta[12].desc, 12, "Candidatus_Korarchaeum_cryptofilum_OPF8", "A", 49.00,
               meta[12].tinf.trans_table, meta[12].tinf.uses_sd);
    write_desc(&mut meta[13].desc, 13, "Catenulispora_acidiphila_DSM_44928", "B", 69.77,
               meta[13].tinf.trans_table, meta[13].tinf.uses_sd);
    write_desc(&mut meta[14].desc, 14, "Cenarchaeum_symbiosum_B", "A", 57.19,
               meta[14].tinf.trans_table, meta[14].tinf.uses_sd);
    write_desc(&mut meta[15].desc, 15, "Chlorobium_phaeobacteroides_BS1", "B", 48.93,
               meta[15].tinf.trans_table, meta[15].tinf.uses_sd);
    write_desc(&mut meta[16].desc, 16, "Chlorobium_tepidum_TLS", "B", 56.53,
               meta[16].tinf.trans_table, meta[16].tinf.uses_sd);
    write_desc(&mut meta[17].desc, 17, "Desulfotomaculum_acetoxidans_DSM_771", "B", 41.55,
               meta[17].tinf.trans_table, meta[17].tinf.uses_sd);
    write_desc(&mut meta[18].desc, 18, "Desulfurococcus_kamchatkensis_1221n", "B", 45.34,
               meta[18].tinf.trans_table, meta[18].tinf.uses_sd);
    write_desc(&mut meta[19].desc, 19, "Erythrobacter_litoralis_HTCC2594", "B", 63.07,
               meta[19].tinf.trans_table, meta[19].tinf.uses_sd);
    write_desc(&mut meta[20].desc, 20, "Escherichia_coli_UMN026", "B", 50.72,
               meta[20].tinf.trans_table, meta[20].tinf.uses_sd);
    write_desc(&mut meta[21].desc, 21, "Haloquadratum_walsbyi_DSM_16790", "A", 47.86,
               meta[21].tinf.trans_table, meta[21].tinf.uses_sd);
    write_desc(&mut meta[22].desc, 22, "Halorubrum_lacusprofundi_ATCC_49239", "A", 57.14,
               meta[22].tinf.trans_table, meta[22].tinf.uses_sd);
    write_desc(&mut meta[23].desc, 23, "Hyperthermus_butylicus_DSM_5456", "A", 53.74,
               meta[23].tinf.trans_table, meta[23].tinf.uses_sd);
    write_desc(&mut meta[24].desc, 24, "Ignisphaera_aggregans_DSM_17230", "A", 35.69,
               meta[24].tinf.trans_table, meta[24].tinf.uses_sd);
    write_desc(&mut meta[25].desc, 25, "Marinobacter_aquaeolei_VT8", "B", 57.27,
               meta[25].tinf.trans_table, meta[25].tinf.uses_sd);
    write_desc(&mut meta[26].desc, 26, "Methanopyrus_kandleri_AV19", "A", 61.16,
               meta[26].tinf.trans_table, meta[26].tinf.uses_sd);
    write_desc(&mut meta[27].desc, 27, "Methanosphaerula_palustris_E1_9c", "A", 55.35,
               meta[27].tinf.trans_table, meta[27].tinf.uses_sd);
    write_desc(&mut meta[28].desc, 28, "Methanothermobacter_thermautotrophicus_Delta_H", "B", 49.54,
               meta[28].tinf.trans_table, meta[28].tinf.uses_sd);
    write_desc(&mut meta[29].desc, 29, "Methylacidiphilum_infernorum_V4", "B", 45.48,
               meta[29].tinf.trans_table, meta[29].tinf.uses_sd);
    write_desc(&mut meta[30].desc, 30, "Mycobacterium_leprae_TN", "B", 57.80,
               meta[30].tinf.trans_table, meta[30].tinf.uses_sd);
    write_desc(&mut meta[31].desc, 31, "Natrialba_magadii_ATCC_43099", "A", 61.42,
               meta[31].tinf.trans_table, meta[31].tinf.uses_sd);
    write_desc(&mut meta[32].desc, 32, "Orientia_tsutsugamushi_Boryong", "B", 30.53,
               meta[32].tinf.trans_table, meta[32].tinf.uses_sd);
    write_desc(&mut meta[33].desc, 33, "Pelotomaculum_thermopropionicum_SI", "B", 52.96,
               meta[33].tinf.trans_table, meta[33].tinf.uses_sd);
    write_desc(&mut meta[34].desc, 34, "Prochlorococcus_marinus_MIT_9313", "B", 50.74,
               meta[34].tinf.trans_table, meta[34].tinf.uses_sd);
    write_desc(&mut meta[35].desc, 35, "Pyrobaculum_aerophilum_IM2", "A", 51.36,
               meta[35].tinf.trans_table, meta[35].tinf.uses_sd);
    write_desc(&mut meta[36].desc, 36, "Ralstonia_solanacearum_PSI07", "B", 66.13,
               meta[36].tinf.trans_table, meta[36].tinf.uses_sd);
    write_desc(&mut meta[37].desc, 37, "Rhizobium_NGR234", "B", 58.49,
               meta[37].tinf.trans_table, meta[37].tinf.uses_sd);
    write_desc(&mut meta[38].desc, 38, "Rhodococcus_jostii_RHA1", "B", 65.05,
               meta[38].tinf.trans_table, meta[38].tinf.uses_sd);
    write_desc(&mut meta[39].desc, 39, "Rickettsia_conorii_Malish_7", "B", 32.44,
               meta[39].tinf.trans_table, meta[39].tinf.uses_sd);
    write_desc(&mut meta[40].desc, 40, "Rothia_dentocariosa_ATCC_17931", "B", 53.69,
               meta[40].tinf.trans_table, meta[40].tinf.uses_sd);
    write_desc(&mut meta[41].desc, 41, "Shigella_dysenteriae_Sd197", "B", 51.25,
               meta[41].tinf.trans_table, meta[41].tinf.uses_sd);
    write_desc(&mut meta[42].desc, 42, "Synechococcus_CC9605", "B", 59.22,
               meta[42].tinf.trans_table, meta[42].tinf.uses_sd);
    write_desc(&mut meta[43].desc, 43, "Synechococcus_JA_2_3B_a_2_13_", "B", 58.45,
               meta[43].tinf.trans_table, meta[43].tinf.uses_sd);
    write_desc(&mut meta[44].desc, 44, "Thermoplasma_volcanium_GSS1", "A", 39.92,
               meta[44].tinf.trans_table, meta[44].tinf.uses_sd);
    write_desc(&mut meta[45].desc, 45, "Treponema_pallidum_Nichols", "B", 52.77,
               meta[45].tinf.trans_table, meta[45].tinf.uses_sd);
    write_desc(&mut meta[46].desc, 46, "Tropheryma_whipplei_TW08_27", "B", 46.31,
               meta[46].tinf.trans_table, meta[46].tinf.uses_sd);
    write_desc(&mut meta[47].desc, 47, "Xenorhabdus_nematophila_ATCC_19061", "B", 44.15,
               meta[47].tinf.trans_table, meta[47].tinf.uses_sd);
    write_desc(&mut meta[48].desc, 48, "Xylella_fastidiosa_Temecula1", "B", 51.78,
               meta[48].tinf.trans_table, meta[48].tinf.uses_sd);
    write_desc(&mut meta[49].desc, 49, "_Nostoc_azollae__0708", "B", 38.45,
               meta[49].tinf.trans_table, meta[49].tinf.uses_sd);
}

fn write_desc(desc: &mut [u8; 500], index: i32, name: &str, domain: &str,
              gc: f64, trans_table: i32, uses_sd: i32) {
    let formatted = format!("{}|{}|{}|{:.1}|{}|{}", index, name, domain, gc, trans_table, uses_sd);
    let bytes = formatted.as_bytes();
    let len = bytes.len().min(499);
    desc[..len].copy_from_slice(&bytes[..len]);
    desc[len] = 0;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_initialize_metagenomic_bins() {
        let mut meta = vec![MetagenomicBin::default(); NUM_META];

        initialize_metagenomic_bins(&mut meta);

        assert_eq!(meta[0].tinf.gc, 0.293);
        assert_eq!(meta[0].tinf.trans_table, 4);
        assert_eq!(meta[0].tinf.uses_sd, 1);

        let desc_str = std::str::from_utf8(&meta[0].desc)
            .unwrap()
            .trim_end_matches('\0');
        assert!(desc_str.starts_with("0|Mycoplasma_bovis_PG45|B|29.3|"));

        assert_eq!(meta[49].tinf.trans_table, 11);
        let desc_str = std::str::from_utf8(&meta[49].desc)
            .unwrap()
            .trim_end_matches('\0');
        assert!(desc_str.starts_with("49|_Nostoc_azollae__0708|B|38.5|"));

        for i in 0..NUM_META {
            assert!(meta[i].tinf.gc > 0.0);
            assert!(meta[i].desc[0] != 0);
        }
    }

    #[test]
    fn test_write_desc() {
        let mut desc = [0u8; 500];
        write_desc(&mut desc, 0, "Test_organism", "B", 50.5, 11, 1);

        let desc_str = std::str::from_utf8(&desc)
            .unwrap()
            .trim_end_matches('\0');
        assert_eq!(desc_str, "0|Test_organism|B|50.5|11|1");
    }

    #[test]
    fn test_metagenomic_bin_default() {
        let bin = MetagenomicBin::default();
        assert_eq!(bin.index, 0);
        assert_eq!(bin.clusnum, 0);
        assert_eq!(bin.weight, 0.0);
        assert_eq!(bin.gc, 0.0);
    }
}

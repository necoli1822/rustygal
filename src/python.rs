// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2007-2016 Doug Hyatt, Univ. of Tennessee / UT-Battelle (original Prodigal C)
// Copyright (C) 2026 Sunju Kim (Rust reimplementation)

use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3::exceptions::{PyValueError, PyRuntimeError};
use std::io::Write;

use crate::*;

#[pymodule]
pub fn rustygal(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyGeneFinder>()?;
    m.add_class::<PyGenes>()?;
    m.add_class::<PyGene>()?;
    m.add_class::<PyTrainingInfo>()?;
    m.add_class::<PySequence>()?;

    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("PRODIGAL_VERSION", "2.6.3")?;
    m.add("MIN_SINGLE_GENOME", 20000)?;
    m.add("IDEAL_SINGLE_GENOME", 100000)?;

    Ok(())
}

#[pyclass(name = "GeneFinder")]
pub struct PyGeneFinder {
    meta: bool,
    closed: i32,
    do_mask: i32,
    min_mask: i32,
    min_gene: i32,
    min_edge_gene: i32,
    max_overlap: i32,
    training_info: Option<training::Training>,
}

#[pymethods]
impl PyGeneFinder {
    #[new]
    #[pyo3(signature = (
        training_info=None,
        *,
        meta=false,
        closed=false,
        mask=false,
        min_mask=50,
        min_gene=90,
        min_edge_gene=60,
        max_overlap=60
    ))]
    fn new(
        training_info: Option<PyTrainingInfo>,
        meta: bool,
        closed: bool,
        mask: bool,
        min_mask: i32,
        min_gene: i32,
        min_edge_gene: i32,
        max_overlap: i32,
    ) -> PyResult<Self> {
        Ok(PyGeneFinder {
            meta,
            closed: if closed { 1 } else { 0 },
            do_mask: if mask { 1 } else { 0 },
            min_mask,
            min_gene,
            min_edge_gene,
            max_overlap,
            training_info: training_info.map(|ti| ti.inner),
        })
    }

    fn find_genes(&mut self, _py: Python<'_>, sequence: &Bound<'_, PyAny>) -> PyResult<PyGenes> {

        let seq_bytes: Vec<u8> = if let Ok(s) = sequence.extract::<String>() {
            s.into_bytes()
        } else if let Ok(b) = sequence.extract::<Vec<u8>>() {
            b
        } else {
            return Err(PyValueError::new_err("sequence must be str or bytes"));
        };

        let config = crate::api::GeneFinderConfig {
            meta: self.meta,
            closed: self.closed,
            do_mask: self.do_mask,
            min_gene: self.min_gene,
            min_edge_gene: self.min_edge_gene,
            max_overlap: self.max_overlap,
        };

        let result = crate::api::find_genes(
            &seq_bytes,
            self.training_info.as_ref(),
            &config,
        ).map_err(|e| PyRuntimeError::new_err(e))?;

        let py_genes: Vec<PyGene> = result.genes.iter().map(|g| {
            PyGene::from_gene(g, &result.nodes, &seq_bytes)
        }).collect();

        Ok(PyGenes {
            genes: py_genes,
            score: 0.0,
            sequence: seq_bytes,
            training_info: result.training_info,
        })
    }

    #[pyo3(signature = (sequence, *_sequences, force_nonsd=false, _start_weight=4.35, translation_table=11))]
    fn train(
        &mut self,
        _py: Python<'_>,
        sequence: &Bound<'_, PyAny>,
        _sequences: &Bound<'_, pyo3::types::PyTuple>,
        force_nonsd: bool,
        _start_weight: f64,
        translation_table: i32,
    ) -> PyResult<PyTrainingInfo> {

        let seq_bytes: Vec<u8> = if let Ok(s) = sequence.extract::<String>() {
            s.into_bytes()
        } else if let Ok(b) = sequence.extract::<Vec<u8>>() {
            b
        } else {
            return Err(PyValueError::new_err("sequence must be str or bytes"));
        };

        let tinf = crate::api::train_on_sequence(
            &seq_bytes,
            translation_table,
            force_nonsd,
        ).map_err(|e| PyRuntimeError::new_err(e))?;

        self.training_info = Some(tinf.clone());

        Ok(PyTrainingInfo { inner: tinf })
    }

    #[getter]
    fn meta(&self) -> bool {
        self.meta
    }

    #[getter]
    fn closed(&self) -> bool {
        self.closed != 0
    }

    #[getter]
    fn mask(&self) -> bool {
        self.do_mask != 0
    }

    #[getter]
    fn min_gene(&self) -> i32 {
        self.min_gene
    }

    #[getter]
    fn min_edge_gene(&self) -> i32 {
        self.min_edge_gene
    }

    #[getter]
    fn max_overlap(&self) -> i32 {
        self.max_overlap
    }

    #[getter]
    fn training_info(&self) -> Option<PyTrainingInfo> {
        self.training_info.as_ref().map(|tinf| PyTrainingInfo { inner: tinf.clone() })
    }
}

#[pyclass(name = "Genes")]
pub struct PyGenes {
    genes: Vec<PyGene>,
    score: f64,
    sequence: Vec<u8>,
    training_info: Option<training::Training>,
}

#[pymethods]
impl PyGenes {
    fn __len__(&self) -> usize {
        self.genes.len()
    }

    fn __getitem__(&self, index: isize) -> PyResult<PyGene> {
        let idx = if index < 0 {
            (self.genes.len() as isize + index) as usize
        } else {
            index as usize
        };

        self.genes.get(idx)
            .cloned()
            .ok_or_else(|| PyValueError::new_err("index out of range"))
    }

    #[getter]
    fn score(&self) -> f64 {
        self.score
    }

    #[pyo3(signature = (file, sequence_id, header=true, _include_translation_table=false))]
    fn write_gff(
        &self,
        py: Python<'_>,
        file: &Bound<'_, PyAny>,
        sequence_id: &str,
        header: bool,
        _include_translation_table: bool,
    ) -> PyResult<usize> {
        let mut output = Vec::new();

        if header {
            writeln!(&mut output, "##gff-version  3").unwrap();
        }

        for (i, gene) in self.genes.iter().enumerate() {
            let strand_char = if gene.strand == 1 { '+' } else { '-' };

            writeln!(
                &mut output,
                "{}\tProdigal\tCDS\t{}\t{}\t{:.2}\t{}\t0\tID={}",
                sequence_id,
                gene.begin,
                gene.end,
                gene.score,
                strand_char,
                i + 1
            ).unwrap();
        }

        let output_str = String::from_utf8(output).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let result = file.call_method1("write", (output_str,))?;
        let bytes_written: usize = result.extract()?;

        Ok(bytes_written)
    }

    #[pyo3(signature = (file, sequence_id))]
    fn write_genbank(
        &self,
        py: Python<'_>,
        file: &Bound<'_, PyAny>,
        sequence_id: &str,
    ) -> PyResult<usize> {
        let mut output = Vec::new();

        for gene in &self.genes {
            let (left, right) = if gene.strand == 1 {
                (gene.begin, gene.end)
            } else {
                (gene.end, gene.begin)
            };

            if gene.strand == 1 {
                writeln!(&mut output, "     CDS             {}..{}", left, right).unwrap();
            } else {
                writeln!(&mut output, "     CDS             complement({}..{})", left, right).unwrap();
            }
        }

        let output_str = String::from_utf8(output).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let result = file.call_method1("write", (output_str,))?;
        let bytes_written: usize = result.extract()?;

        Ok(bytes_written)
    }

    #[pyo3(signature = (file, sequence_id, width=70))]
    fn write_genes(
        &self,
        py: Python<'_>,
        file: &Bound<'_, PyAny>,
        sequence_id: &str,
        width: usize,
    ) -> PyResult<usize> {
        let mut output = Vec::new();

        for (i, gene) in self.genes.iter().enumerate() {
            let strand_num = if gene.strand == 1 { 1 } else { -1 };

            writeln!(
                &mut output,
                ">{}_{} # {} # {} # {} # ID={};partial={}{};",
                sequence_id,
                i + 1,
                gene.begin,
                gene.end,
                strand_num,
                i + 1,
                if gene.partial_begin { "1" } else { "0" },
                if gene.partial_end { "1" } else { "0" }
            ).unwrap();

            for chunk in gene.gene_sequence.chunks(width) {
                if let Ok(s) = std::str::from_utf8(chunk) {
                    writeln!(&mut output, "{}", s).unwrap();
                }
            }
        }

        let output_str = String::from_utf8(output).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let result = file.call_method1("write", (output_str,))?;
        let bytes_written: usize = result.extract()?;

        Ok(bytes_written)
    }

    #[pyo3(signature = (file, sequence_id, width=70))]
    fn write_translations(
        &self,
        py: Python<'_>,
        file: &Bound<'_, PyAny>,
        sequence_id: &str,
        width: usize,
    ) -> PyResult<usize> {
        let mut output = Vec::new();

        let trans_table = self.training_info.as_ref().map(|t| t.trans_table).unwrap_or(11);

        for (i, gene) in self.genes.iter().enumerate() {
            let strand_num = if gene.strand == 1 { 1 } else { -1 };

            writeln!(
                &mut output,
                ">{}_{} # {} # {} # {} # ID={};partial={}{};",
                sequence_id,
                i + 1,
                gene.begin,
                gene.end,
                strand_num,
                i + 1,
                if gene.partial_begin { "1" } else { "0" },
                if gene.partial_end { "1" } else { "0" }
            ).unwrap();

            let protein = Self::simple_translate(&gene.gene_sequence, trans_table);

            for chunk in protein.as_bytes().chunks(width) {
                if let Ok(s) = std::str::from_utf8(chunk) {
                    writeln!(&mut output, "{}", s).unwrap();
                }
            }
        }

        let output_str = String::from_utf8(output).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let result = file.call_method1("write", (output_str,))?;
        let bytes_written: usize = result.extract()?;

        Ok(bytes_written)
    }
}

impl PyGenes {

    fn simple_translate(seq: &[u8], _trans_table: i32) -> String {
        let mut protein = String::new();

        for codon_bytes in seq.chunks(3) {
            if codon_bytes.len() < 3 {
                break;
            }

            let codon = std::str::from_utf8(codon_bytes).unwrap_or("NNN").to_uppercase();
            let aa = match codon.as_str() {
                "TTT" | "TTC" => 'F',
                "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => 'L',
                "ATT" | "ATC" | "ATA" => 'I',
                "ATG" => 'M',
                "GTT" | "GTC" | "GTA" | "GTG" => 'V',
                "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => 'S',
                "CCT" | "CCC" | "CCA" | "CCG" => 'P',
                "ACT" | "ACC" | "ACA" | "ACG" => 'T',
                "GCT" | "GCC" | "GCA" | "GCG" => 'A',
                "TAT" | "TAC" => 'Y',
                "TAA" | "TAG" | "TGA" => '*',
                "CAT" | "CAC" => 'H',
                "CAA" | "CAG" => 'Q',
                "AAT" | "AAC" => 'N',
                "AAA" | "AAG" => 'K',
                "GAT" | "GAC" => 'D',
                "GAA" | "GAG" => 'E',
                "TGT" | "TGC" => 'C',
                "TGG" => 'W',
                "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => 'R',
                "GGT" | "GGC" | "GGA" | "GGG" => 'G',
                _ => 'X',
            };
            protein.push(aa);
        }

        protein
    }
}

#[pyclass(name = "Gene")]
#[derive(Clone)]
pub struct PyGene {
    begin: i32,
    end: i32,
    strand: i32,
    start_type: String,
    rbs_motif: Option<String>,
    rbs_spacer: Option<String>,
    gc_cont: f64,
    cscore: f64,
    rscore: f64,
    sscore: f64,
    tscore: f64,
    uscore: f64,
    score: f64,
    partial_begin: bool,
    partial_end: bool,
    gene_sequence: Vec<u8>,
}

impl PyGene {

    fn parse_data_string(data: &str) -> std::collections::HashMap<String, String> {
        let mut map = std::collections::HashMap::new();
        for pair in data.split(';') {
            if let Some(eq_pos) = pair.find('=') {
                let key = pair[..eq_pos].trim().to_string();
                let value = pair[eq_pos + 1..].trim().to_string();
                map.insert(key, value);
            }
        }
        map
    }

    fn from_gene(gene: &gene::Gene, nodes: &[node::Node], sequence: &[u8]) -> Self {
        let gene_data_str = std::str::from_utf8(&gene.gene_data)
            .unwrap_or("")
            .trim_end_matches('\0');

        let score_data_str = std::str::from_utf8(&gene.score_data)
            .unwrap_or("")
            .trim_end_matches('\0');

        let gene_map = Self::parse_data_string(gene_data_str);
        let score_map = Self::parse_data_string(score_data_str);

        let partial_str = gene_map.get("partial").map(|s| s.as_str()).unwrap_or("00");
        let partial_begin = partial_str.chars().next() == Some('1');
        let partial_end = partial_str.chars().nth(1) == Some('1');

        let start_type = gene_map
            .get("start_type")
            .map(|s| s.clone())
            .unwrap_or_else(|| "ATG".to_string());

        let rbs_motif = gene_map.get("rbs_motif").and_then(|s| {
            if s == "None" { None } else { Some(s.clone()) }
        });

        let rbs_spacer = gene_map.get("rbs_spacer").and_then(|s| {
            if s == "None" { None } else { Some(s.clone()) }
        });

        let gc_cont = gene_map
            .get("gc_cont")
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.5);

        let cscore = score_map.get("cscore").and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        let rscore = score_map.get("rscore").and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        let sscore = score_map.get("sscore").and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        let tscore = score_map.get("tscore").and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        let uscore = score_map.get("uscore").and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        let score = score_map.get("score").and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);

        let strand = if gene.start_ndx < nodes.len() as i32 {
            nodes[gene.start_ndx as usize].strand
        } else {
            if gene.begin < gene.end { 1 } else { -1 }
        };

        let gene_sequence = if strand == 1 {
            sequence[(gene.begin - 1) as usize..gene.end as usize].to_vec()
        } else {

            let subseq = &sequence[(gene.begin - 1) as usize..gene.end as usize];
            subseq.iter().rev().map(|&b| match b {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'G' | b'g' => b'C',
                b'C' | b'c' => b'G',
                _ => b'N',
            }).collect()
        };

        PyGene {
            begin: gene.begin,
            end: gene.end,
            strand,
            start_type,
            rbs_motif,
            rbs_spacer,
            gc_cont,
            cscore,
            rscore,
            sscore,
            tscore,
            uscore,
            score,
            partial_begin,
            partial_end,
            gene_sequence,
        }
    }
}

#[pymethods]
impl PyGene {
    #[getter]
    fn begin(&self) -> i32 {
        self.begin
    }

    #[getter]
    fn end(&self) -> i32 {
        self.end
    }

    #[getter]
    fn strand(&self) -> i32 {
        self.strand
    }

    #[getter]
    fn start_type(&self) -> &str {
        &self.start_type
    }

    #[getter]
    fn rbs_motif(&self) -> Option<&str> {
        self.rbs_motif.as_deref()
    }

    #[getter]
    fn rbs_spacer(&self) -> Option<&str> {
        self.rbs_spacer.as_deref()
    }

    #[getter]
    fn gc_cont(&self) -> f64 {
        self.gc_cont
    }

    #[getter]
    fn cscore(&self) -> f64 {
        self.cscore
    }

    #[getter]
    fn rscore(&self) -> f64 {
        self.rscore
    }

    #[getter]
    fn sscore(&self) -> f64 {
        self.sscore
    }

    #[getter]
    fn tscore(&self) -> f64 {
        self.tscore
    }

    #[getter]
    fn uscore(&self) -> f64 {
        self.uscore
    }

    #[getter]
    fn score(&self) -> f64 {
        self.score
    }

    #[getter]
    fn partial_begin(&self) -> bool {
        self.partial_begin
    }

    #[getter]
    fn partial_end(&self) -> bool {
        self.partial_end
    }

    fn sequence(&self) -> PyResult<String> {
        String::from_utf8(self.gene_sequence.clone())
            .map_err(|e| PyValueError::new_err(format!("Invalid UTF-8: {}", e)))
    }

    #[pyo3(signature = (translation_table=11))]
    fn translate(&self, translation_table: i32) -> PyResult<String> {
        Ok(PyGenes::simple_translate(&self.gene_sequence, translation_table))
    }

    fn confidence(&self) -> f64 {
        if self.score > 100.0 {
            99.99
        } else if self.score < 0.0 {
            0.0
        } else {
            self.score
        }
    }
}

#[pyclass(name = "TrainingInfo")]
#[derive(Clone)]
pub struct PyTrainingInfo {
    inner: training::Training,
}

#[pymethods]
impl PyTrainingInfo {
    #[getter]
    fn gc(&self) -> f64 {
        self.inner.gc
    }

    #[getter]
    fn translation_table(&self) -> i32 {
        self.inner.trans_table
    }
}

#[pyclass(name = "Sequence")]
pub struct PySequence {
    gc: f64,
    length: usize,
}

#[pymethods]
impl PySequence {
    #[getter]
    fn gc(&self) -> f64 {
        self.gc
    }

    fn __len__(&self) -> usize {
        self.length
    }
}

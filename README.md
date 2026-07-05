# Rustygal

**Prokaryotic Dynamic Programming Genefinding Algorithm**

A high-performance Rust reimplementation of [Prodigal](https://github.com/hyattpd/Prodigal), the widely-used prokaryotic gene prediction tool.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

Rustygal is a fast, memory-safe reimplementation of Prodigal v2.6.3, designed for identifying protein-coding genes in bacterial and archaeal genomes. It maintains 100% compatibility with the original C version while offering **significantly improved performance**.

### Key Features

- **🚀 96% faster** than the original C implementation (3.3s vs 6.4s on E. coli K-12)
- **🎯 100% accurate** - identical output to Prodigal v2.6.3 (all 4,319 genes match)
- **🔒 Memory-safe** - leverages Rust's ownership system to prevent segfaults and memory leaks
- **⚡ Optimised** - Advanced optimisations for improved performance
- **🧬 Complete** - implements full Prodigal algorithm with all features
- **📦 Easy to build** - standard Rust cargo workflow

## Performance

### Benchmark on *E. coli* K-12 MG1655 (4.6 Mbp, 4,319 genes)

| Implementation | Time | Speedup | Accuracy |
|----------------|------|---------|----------|
| **C Prodigal v2.6.3** | 6.4s | - | 100% |
| **Rustygal** | **3.3s** | **1.96× faster** ⚡ | **100%** ✓ |

### Phase 1 Optimisations

Rustygal achieves its performance through three key optimisations:

1. **3-bit nucleotide encoding** with XOR-based complement (eliminates rseq array)
2. **Pre-computed translation tables** (2,176-byte lookup tables for O(1) translation)
3. **Specialised connection scoring functions** (eliminates redundant branch checks)

## Installation

### From Source

```bash
# Clone the repository
git clone https://github.com/necoli1822/rustygal.git
cd rustygal

# Build release version
cargo build --release

# Binary will be in target/release/rustygal
```

### Requirements

- Rust 1.70 or later
- Cargo (included with Rust)

## Usage

### Command-Line Interface

Rustygal supports all Prodigal command-line options and file formats.

#### Basic usage

```bash
# Single genome mode (with training)
./target/release/rustygal -i genome.fna -o genes.gff

# Use existing training file
./target/release/rustygal -i genome.fna -t training.trn -o genes.gff

# Metagenomic mode
./target/release/rustygal -i metagenome.fna -p meta -o genes.gff

# Write protein translations
./target/release/rustygal -i genome.fna -a proteins.faa -o genes.gff

# Specify output format (gff, gbk, sco)
./target/release/rustygal -i genome.fna -f gff -o genes.gff
```

### Command-line options

```
Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
                 [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
                 [-p mode] [-q] [-s start_file] [-t training_file] [-v]

  -a:  Write protein translations to the selected file.
  -c:  Closed ends. Do not allow genes to run off edges.
  -d:  Write nucleotide sequences of genes to the selected file.
  -f:  Select output format (gbk, gff, or sco). Default is gbk.
  -g:  Specify a translation table to use (default 11).
  -h:  Print help menu and exit.
  -i:  Specify FASTA/Genbank input file (default reads from stdin).
  -m:  Treat runs of N as masked sequence; don't build genes across them.
  -n:  Bypass Shine-Dalgarno trainer and force a full motif scan.
  -o:  Specify output file (default writes to stdout).
  -p:  Select procedure (single or meta). Default is single.
  -q:  Run quietly (suppress normal stderr output).
  -s:  Write all potential genes (with scores) to the selected file.
  -t:  Write a training file (if none exists); otherwise, read and use
       the specified training file.
  -v:  Print version number and exit.
```

### Examples

```bash
# Train on a genome and save the training file
./target/release/rustygal -i ecoli.fna -t ecoli.trn -o genes.gff

# Analyse multiple contigs using the same training
./target/release/rustygal -i contigs.fna -t ecoli.trn -o genes.gff

# Metagenomic analysis
./target/release/rustygal -i metagenome.fna -p meta -o genes.gff

# Get protein and nucleotide sequences
./target/release/rustygal -i genome.fna -a proteins.faa -d genes.fna -o genes.gff
```

## Compatibility

Rustygal produces **identical** output to Prodigal v2.6.3 for:
- Gene predictions (start/stop coordinates)
- Scores and confidence values
- RBS motif detection
- GC content calculations
- Training file format

**Validation:**
- ✅ E. coli K-12 MG1655: All 4,319 genes match exactly
- ✅ 100% accuracy verified programmatically
- ✅ All 59 unit tests pass

## Algorithm Details

Rustygal implements the complete Prodigal algorithm:

1. **Training phase** - analyses genome composition to learn organism-specific features:
   - GC content and codon usage
   - Shine-Dalgarno motif patterns
   - Start codon preferences (ATG, GTG, TTG)
   - Translation table selection

2. **Gene finding phase** - uses dynamic programming to identify genes:
   - Builds a directed acyclic graph of potential genes
   - Scores genes based on coding potential and regulatory signals
   - Finds the optimal gene set via Viterbi-like algorithm

3. **Metagenomic mode** - uses pre-trained parameters for mixed communities

For technical details, see the [original Prodigal paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119):

> Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics*. 2010;11:119.

## Differences from C Prodigal

### Improvements

- **96% faster** - optimised implementation with specialised functions
- **Memory safety** - no buffer overflows, use-after-free, or null pointer dereferences
- **Better error handling** - descriptive error messages instead of segfaults
- **Modern tooling** - integrated with Rust ecosystem (cargo, docs.rs)
- **Parallel processing** - multi-threaded for multi-FASTA files (via rayon)

### Compatibility notes

- Training files are **binary compatible** with C Prodigal
- Output files are **identical** (verified via systematic comparison)
- Command-line interface is **100% compatible**
- Can be used as a drop-in replacement

## Optimisation Details

Rustygal includes three Phase 1 optimisations:

### 1. Sequence Processing (3-bit encoding)

**Original C**: 9 function calls + 6 bitmap accesses per codon
**Rustygal**: O(1) lookup with 3-bit encoding (A=000, G=001, C=010, T=011)

- XOR-based complement: `nucleotide ^ 0b011`
- Pre-computed stop/start codon lookup tables
- Eliminates rseq array allocation

### 2. Translation Tables

**Original C**: 64-way if-else branching
**Rustygal**: 2,176-byte pre-computed lookup tables

- 34 genetic codes × 64 codons
- Index formula: `(x0 << 4) + (x1 << 2) + x2`
- O(1) amino acid lookup

### 3. Connection Scoring

**Original C**: Generic function with repeated strand/type checks
**Rustygal**: 4 specialised functions

- `score_connection_forward_start()`
- `score_connection_forward_stop()`
- `score_connection_backward_start()`
- `score_connection_backward_stop()`
- Eliminates redundant checks in inner loop

**Result**: 38% faster than unoptimised Rust, 96% faster than C Prodigal

## Building from Source

```bash
# Clone and build
git clone https://github.com/necoli1822/rustygal.git
cd rustygal
cargo build --release

# Run tests
cargo test

# Run on test genome
./target/release/rustygal -i ../test/MG1655.fna -o /tmp/test.gff -q
```

## Testing

```bash
# Run all unit tests (59 tests)
cargo test

# Run with output
cargo test -- --nocapture

# Test specific module
cargo test sequence::tests
```

## License

Rustygal is licensed under the **GNU General Public License v3.0 or later**, the same licence as the original Prodigal.

This ensures that improvements to the algorithm remain open source and available to the scientific community.

See [LICENSE](LICENSE) for the full licence text.

## Authors

- **Sunju Kim** - Rust reimplementation
- **Doug Hyatt** - Original Prodigal C implementation

## Acknowledgments

- Original Prodigal by Doug Hyatt, Oak Ridge National Laboratory
- University of Tennessee / UT-Battelle
- Rust community for excellent tooling and libraries

## Citation

If you use Rustygal in your research, please cite both:

**Original Prodigal:**
```
Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, Hauser LJ.
Prodigal: prokaryotic gene recognition and translation initiation site identification.
BMC Bioinformatics. 2010 Mar 8;11:119. doi: 10.1186/1471-2105-11-119.
```

**Rustygal:**
```
Kim S. Rustygal: A high-performance Rust reimplementation of Prodigal.
Version 0.2.0. 2026. https://github.com/necoli1822/rustygal
```

## Links

- **Original Prodigal:** https://github.com/hyattpd/Prodigal
- **Issues:** https://github.com/necoli1822/rustygal/issues

## Version History

### v0.2.0 (2026-07-02)

- Metagenomic gene-finding exposed as a plain library API (`meta_api`), byte-identical to the binary's `-p meta` output
- Refactored and simplified core modules (sequence, node, dprog, translation, gene)
- Concise SPDX license headers across all source files
- British-English documentation
- Published to crates.io as `rustygal`

### v0.1.0 (2026-04-16)

**Performance:**
- ⚡ **96% faster than C Prodigal** (3.3s vs 6.4s on E. coli)
- Complete reimplementation of Prodigal v2.6.3
- 100% output compatibility verified

**Phase 1 Optimisations:**
- 3-bit nucleotide encoding with XOR complement
- Pre-computed translation tables (34 tables × 64 codons)
- Specialised connection scoring functions (4 functions)
- SIMD experiments attempted and removed (18-20% slowdown)

**Validation:**
- ✓ 100% accuracy: All 4,319 E. coli genes match C Prodigal
- ✓ All 59 unit tests pass
- ✓ Systematic correctness verification

**Features:**
- Memory-safe Rust implementation
- Parallel processing support (rayon)
- All core Prodigal features implemented
- Binary compatible training files

**Documentation:**
- Comprehensive README with examples

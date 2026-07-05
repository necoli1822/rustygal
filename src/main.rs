// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2007-2016 Doug Hyatt, Univ. of Tennessee / UT-Battelle (original Prodigal C)
// Copyright (C) 2026 Sunju Kim (Rust reimplementation)

use rustygal::*;
use rustygal::fptr::FilePtr;
use std::env;
use std::process;
use std::fs::File;
use std::io::{self, Write, BufRead, BufReader};
use rayon::prelude::*;

#[derive(Clone)]
struct SequenceData {

    seq: Vec<u8>,

    rseq: Vec<u8>,

    useq: Vec<u8>,

    slen: i32,

    gc: f64,

    header: Vec<u8>,

    short_header: String,

    mlist: Vec<sequence::Mask>,

    nmask: i32,

    seq_num: i32,
}

struct SequenceResult {

    output_main: Vec<u8>,

    output_trans: Option<Vec<u8>>,

    output_nuc: Option<Vec<u8>>,

    output_start: Option<Vec<u8>>,
}

const VERSION: &str = "0.2.0";
const DATE: &str = "June, 2026";
const MIN_SINGLE_GENOME: i32 = 20_000;
const IDEAL_SINGLE_GENOME: i32 = 100_000;

fn version() {
    eprintln!("\nRustygal V{}: {}\n", VERSION, DATE);
    eprintln!("Rust port of Prodigal by Sunju Kim <n.e.coli.1822@gmail.com>");
    eprintln!("Original Prodigal by Doug Hyatt\n");
    process::exit(0);
}

fn usage(msg: &str) {
    eprintln!("\n{}\n", msg);
    eprintln!("Usage:  rustygal [-a trans_file] [-c] [-d nuc_file] [-f output_type]");
    eprintln!("                 [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]");
    eprintln!("                 [-p mode] [-q] [-s start_file] [-t training_file] [-v] [-w threads]");
    eprintln!("\nDo 'rustygal -h' for more information.\n");
    process::exit(15);
}

fn help() {
    eprintln!("\nUsage:  rustygal [-a trans_file] [-c] [-d nuc_file] [-f output_type]");
    eprintln!("                 [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]");
    eprintln!("                 [-p mode] [-q] [-s start_file] [-t training_file] [-v] [-w threads]\n");
    eprintln!("         -a:  Write protein translations to the selected file.");
    eprintln!("         -c:  Closed ends.  Do not allow genes to run off edges.");
    eprintln!("         -d:  Write nucleotide sequences of genes to the selected file.");
    eprintln!("         -f:  Select output format (gbk, gff, or sco).  Default is gbk.");
    eprintln!("         -g:  Specify a translation table to use (default 11).");
    eprintln!("         -h:  Print help menu and exit.");
    eprintln!("         -i:  Specify FASTA/Genbank input file (default reads from stdin).");
    eprintln!("         -m:  Treat runs of N as masked sequence; don't build genes across them.");
    eprintln!("         -n:  Bypass Shine-Dalgarno trainer and force a full motif scan.");
    eprintln!("         -o:  Specify output file (default writes to stdout).");
    eprintln!("         -p:  Select procedure (single or meta).  Default is single.");
    eprintln!("         -q:  Run quietly (suppress normal stderr output).");
    eprintln!("         -s:  Write all potential genes (with scores) to the selected file.");
    eprintln!("         -t:  Write a training file (if none exists); otherwise, read and use");
    eprintln!("              the specified training file.");
    eprintln!("         -v:  Print version number and exit.");
    eprintln!("         -w:  Number of worker threads (default: auto-detect CPU cores).\n");
    process::exit(0);
}

fn copy_standard_input_to_file(path: &str, quiet: bool) -> io::Result<()> {
    if !quiet {
        eprintln!("Piped input detected, copying stdin to a tmp file...");
    }

    let mut file = File::create(path)?;
    let stdin = io::stdin();
    let reader = BufReader::new(stdin);

    for line in reader.lines() {
        let line = line?;
        writeln!(file, "{}", line)?;
    }

    if !quiet {
        eprintln!("done!");
        eprintln!("-------------------------------------");
    }

    Ok(())
}

fn read_all_sequences(
    input_ptr: &mut FilePtr,
    do_mask: i32,
    quiet: bool
) -> io::Result<Vec<SequenceData>> {
    let mut sequences = Vec::new();

    let mut mlist = vec![sequence::Mask { begin: 0, end: 0 }; sequence::MAX_MASKS];
    let mut cur_header = vec![0u8; 4000];
    let mut new_header = vec![0u8; 4000];
    let mut dna: Vec<u8> = Vec::new();

    let mut num_seq = 0;

    loop {
        let mut nmask = 0;

        let slen = sequence::next_seq_multi(
            input_ptr,
            &mut dna,
            &mut num_seq,
            do_mask,
            &mut mlist,
            &mut nmask,
            &mut cur_header,
            &mut new_header
        );

        if slen == -1 {
            break;
        }

        if slen == 0 {
            eprintln!("\nSequence read failed (file must be Fasta, Genbank, or EMBL format).\n");
            process::exit(14);
        }

        // Right-size the bitmaps to this record's true length and build them through the same
        // build_bitmaps path the library uses, so the encoding (and gc) is byte-identical.
        let slen_usize = slen as usize;
        let bytes_2bit = slen_usize / 4 + 16;
        let bytes_1bit = slen_usize / 8 + 16;
        let mut seq = vec![0u8; bytes_2bit];
        let mut rseq = vec![0u8; bytes_2bit];
        let mut useq = vec![0u8; bytes_1bit];
        let mut gc = 0.0;
        bitmap::build_bitmaps(&dna, slen, &mut seq, &mut rseq, &mut useq, &mut gc);

        let header_end = cur_header.iter().position(|&b| b == 0).unwrap_or(cur_header.len());
        let header_trimmed = cur_header[..header_end].to_vec();
        let cur_header_str = std::str::from_utf8(&header_trimmed).unwrap_or("");

        let mut short_header = String::new();
        sequence::calc_short_header(cur_header_str, &mut short_header, num_seq);

        if !quiet {
            eprintln!("Read sequence #{} ({} bp)...", num_seq, slen);
        }

        sequences.push(SequenceData {
            seq,
            rseq,
            useq,
            slen,
            gc,
            header: header_trimmed,
            short_header,
            mlist: mlist[..nmask as usize].to_vec(),
            nmask,
            seq_num: num_seq,
        });

        cur_header.copy_from_slice(&new_header);
    }

    Ok(sequences)
}

fn process_single_sequence(
    seq_data: &SequenceData,
    tinf: &training::Training,
    meta: &[metagenomic::MetagenomicBin],
    is_meta: i32,
    closed: i32,
    output: i32,
    quiet: bool,
    write_trans: bool,
    write_nuc: bool,
    write_start: bool,
) -> SequenceResult {

    let slen_u = seq_data.slen.max(0) as usize;
    let node_cap = (slen_u / 8 + 65536).min(slen_u * 2 + 16).max(16);
    let gene_cap = (slen_u / 20 + 1024).min(gene::MAX_GENES).max(16);
    let mut nodes = vec![node::Node::default(); node_cap];
    let mut genes = vec![gene::Gene::default(); gene_cap];

    let cur_header_str = std::str::from_utf8(&seq_data.header)
        .unwrap_or("")
        .trim_end_matches('\0');

    let mut nn: i32 = 0;
    let mut ng: i32 = 0;
    let mut ipath: i32;
    let mut max_phase: i32 = 0;
    let mut max_score: f64 = -100.0;

    if !quiet {
        eprint!("Processing sequence #{} ({} bp)...", seq_data.seq_num, seq_data.slen);
    }

    if is_meta == 0 {

        nn = node::add_nodes(
            &seq_data.seq,
            &seq_data.rseq,
            seq_data.slen,
            &mut nodes,
            closed,
            &seq_data.mlist,
            seq_data.nmask,
            tinf
        );
        nodes[..nn as usize].sort_unstable_by(|a, b| node::compare_nodes(a, b));

        node::score_nodes(
            &seq_data.seq,
            &seq_data.rseq,
            seq_data.slen,
            &mut nodes,
            nn,
            tinf,
            closed,
            is_meta
        );

        node::record_overlapping_starts(&mut nodes, nn, tinf, 1);
        ipath = dprog::dprog(&mut nodes, nn, tinf, 1);
        dprog::eliminate_bad_genes(&mut nodes, ipath, tinf);

        ng = gene::add_genes(&mut genes, &nodes, ipath);
        gene::tweak_final_starts(&mut genes, ng, &mut nodes, nn, tinf);
        gene::record_gene_data(&mut genes, ng, &nodes, tinf, seq_data.seq_num);

    } else {

        let mut low = 0.88495 * seq_data.gc - 0.0102337;
        if low > 0.65 {
            low = 0.65;
        }
        let mut high = 0.86596 * seq_data.gc + 0.1131991;
        if high < 0.35 {
            high = 0.35;
        }

        for i in 0..metagenomic::NUM_META {

            if i == 0 || meta[i].tinf.trans_table != meta[i - 1].tinf.trans_table {
                for n in nodes[..].iter_mut() {
                    *n = node::Node::default();
                }
                nn = node::add_nodes(
                    &seq_data.seq,
                    &seq_data.rseq,
                    seq_data.slen,
                    &mut nodes,
                    closed,
                    &seq_data.mlist,
                    seq_data.nmask,
                    &meta[i].tinf
                );
                nodes[..nn as usize].sort_unstable_by(|a, b| node::compare_nodes(a, b));
            }

            if meta[i].tinf.gc < low || meta[i].tinf.gc > high {
                continue;
            }

            node::reset_node_scores(&mut nodes, nn);
            node::score_nodes(
                &seq_data.seq,
                &seq_data.rseq,
                seq_data.slen,
                &mut nodes,
                nn,
                &meta[i].tinf,
                closed,
                is_meta
            );
            node::record_overlapping_starts(&mut nodes, nn, &meta[i].tinf, 1);
            ipath = dprog::dprog(&mut nodes, nn, &meta[i].tinf, 1);

            if ipath >= 0 && nodes[ipath as usize].score > max_score {
                max_phase = i as i32;
                max_score = nodes[ipath as usize].score;
                dprog::eliminate_bad_genes(&mut nodes, ipath, &meta[i].tinf);
                ng = gene::add_genes(&mut genes, &nodes, ipath);
                gene::tweak_final_starts(&mut genes, ng, &mut nodes, nn, &meta[i].tinf);
                gene::record_gene_data(&mut genes, ng, &nodes, &meta[i].tinf, seq_data.seq_num);
            }
        }

        for n in nodes[..].iter_mut() {
            *n = node::Node::default();
        }
        nn = node::add_nodes(
            &seq_data.seq,
            &seq_data.rseq,
            seq_data.slen,
            &mut nodes,
            closed,
            &seq_data.mlist,
            seq_data.nmask,
            &meta[max_phase as usize].tinf
        );
        nodes[..nn as usize].sort_unstable_by(|a, b| node::compare_nodes(a, b));
        node::score_nodes(
            &seq_data.seq,
            &seq_data.rseq,
            seq_data.slen,
            &mut nodes,
            nn,
            &meta[max_phase as usize].tinf,
            closed,
            is_meta
        );
    }

    if !quiet {
        eprintln!("done!");
    }

    let mut output_main = Vec::new();
    let mut output_trans = None;
    let mut output_nuc = None;
    let mut output_start = None;

    if is_meta == 0 {
        gene::print_genes(
            &mut output_main,
            &genes,
            ng,
            &nodes,
            seq_data.slen,
            output,
            seq_data.seq_num,
            0,
            "",
            tinf,
            cur_header_str,
            &seq_data.short_header,
            VERSION
        );
    } else {
        let desc_str = std::str::from_utf8(&meta[max_phase as usize].desc)
            .unwrap()
            .trim_end_matches('\0');
        gene::print_genes(
            &mut output_main,
            &genes,
            ng,
            &nodes,
            seq_data.slen,
            output,
            seq_data.seq_num,
            1,
            desc_str,
            if is_meta == 0 { tinf } else { &meta[max_phase as usize].tinf },
            cur_header_str,
            &seq_data.short_header,
            VERSION
        );
    }

    if write_trans {
        let mut buf = Vec::new();
        gene::write_translations(
            &mut buf,
            &genes,
            ng,
            &nodes,
            &seq_data.seq,
            &seq_data.rseq,
            &seq_data.useq,
            seq_data.slen,
            if is_meta == 0 { tinf } else { &meta[max_phase as usize].tinf },
            seq_data.seq_num,
            &seq_data.short_header
        );
        output_trans = Some(buf);
    }

    if write_nuc {
        let mut buf = Vec::new();
        gene::write_nucleotide_seqs(
            &mut buf,
            &genes,
            ng,
            &nodes,
            &seq_data.seq,
            &seq_data.rseq,
            &seq_data.useq,
            seq_data.slen,
            if is_meta == 0 { tinf } else { &meta[max_phase as usize].tinf },
            seq_data.seq_num,
            &seq_data.short_header
        );
        output_nuc = Some(buf);
    }

    if write_start {
        let mut buf = Vec::new();
        if is_meta == 0 {
            node::write_start_file(
                &mut buf,
                &nodes,
                nn,
                tinf,
                seq_data.seq_num,
                seq_data.slen,
                0,
                "",
                VERSION,
                cur_header_str
            );
        } else {
            let desc_str = std::str::from_utf8(&meta[max_phase as usize].desc)
                .unwrap()
                .trim_end_matches('\0');
            node::write_start_file(
                &mut buf,
                &nodes,
                nn,
                &meta[max_phase as usize].tinf,
                seq_data.seq_num,
                seq_data.slen,
                1,
                desc_str,
                VERSION,
                cur_header_str
            );
        }
        output_start = Some(buf);
    }

    SequenceResult {
        output_main,
        output_trans,
        output_nuc,
        output_start,
    }
}

fn write_all_results(
    results: &[SequenceResult],
    output_ptr: &mut dyn Write,
    trans_ptr: &mut Option<File>,
    nuc_ptr: &mut Option<File>,
    start_ptr: &mut Option<File>,
) -> io::Result<()> {
    for result in results {

        output_ptr.write_all(&result.output_main)?;
        output_ptr.flush()?;

        if let (Some(ref data), Some(ref mut file)) = (&result.output_trans, trans_ptr.as_mut()) {
            file.write_all(data)?;
        }

        if let (Some(ref data), Some(ref mut file)) = (&result.output_nuc, nuc_ptr.as_mut()) {
            file.write_all(data)?;
        }

        if let (Some(ref data), Some(ref mut file)) = (&result.output_start, start_ptr.as_mut()) {
            file.write_all(data)?;
        }
    }

    Ok(())
}

fn main() {

    let mut nodes = vec![node::Node::default(); node::STT_NOD];
    let _genes = vec![gene::Gene::default(); gene::MAX_GENES];
    let mut tinf = training::Training::default();
    let mut meta = vec![metagenomic::MetagenomicBin::default(); metagenomic::NUM_META];
    let mut mlist = vec![sequence::Mask { begin: 0, end: 0 }; sequence::MAX_MASKS];

    let nn: i32;
    let slen: i32;
    let ipath: i32;
    let _ng: i32 = 0;
    let mut nmask: i32 = 0;
    let mut user_tt: i32 = 0;
    let mut is_meta: i32 = 0;
    let num_seq: i32;
    let mut quiet: bool = false;
    let _max_phase: i32 = 0;
    let _max_score: f64 = -100.0;
    let mut do_training: bool = false;
    let mut output: i32 = 0;
    let mut closed: i32 = 0;
    let mut do_mask: i32 = 0;
    let mut force_nonsd: i32 = 0;
    let mut piped: bool = false;
    let max_slen: i32 = 0;
    let mut num_threads: Option<usize> = None;

    let mut train_file: Option<String> = None;
    let mut start_file: Option<String> = None;
    let mut trans_file: Option<String> = None;
    let mut nuc_file: Option<String> = None;
    let mut input_file: Option<String> = None;
    let mut output_file: Option<String> = None;

    let pid = process::id();
    let input_copy = format!("tmp.prodigal.stdin.{}", pid);

    tinf.st_wt = 4.35;
    tinf.trans_table = 11;

    let args: Vec<String> = env::args().collect();

    if args.len() == 1 {
        usage("No arguments provided");
    }

    let mut i = 1;
    while i < args.len() {
        let arg = &args[i];

        if i == args.len() - 1 && matches!(arg.as_str(),
            "-t" | "-T" | "-a" | "-A" | "-g" | "-G" | "-f" | "-F" |
            "-s" | "-S" | "-i" | "-I" | "-o" | "-O" | "-p" | "-P" | "-d" | "-D" | "-w" | "-W")
        {
            usage("-a/-d/-f/-g/-i/-o/-p/-s/-t/-w options require parameters.");
        }

        match arg.as_str() {
            "-c" | "-C" => closed = 1,
            "-q" | "-Q" => quiet = true,
            "-m" | "-M" => do_mask = 1,
            "-n" | "-N" => force_nonsd = 1,
            "-h" | "-H" => help(),
            "-v" | "-V" => version(),
            "-a" | "-A" => {
                trans_file = Some(args[i + 1].clone());
                i += 1;
            }
            "-d" | "-D" => {
                nuc_file = Some(args[i + 1].clone());
                i += 1;
            }
            "-i" | "-I" => {
                input_file = Some(args[i + 1].clone());
                i += 1;
            }
            "-o" | "-O" => {
                output_file = Some(args[i + 1].clone());
                i += 1;
            }
            "-s" | "-S" => {
                start_file = Some(args[i + 1].clone());
                i += 1;
            }
            "-t" | "-T" => {
                train_file = Some(args[i + 1].clone());
                i += 1;
            }
            "-g" | "-G" => {
                tinf.trans_table = args[i + 1].parse::<i32>().unwrap_or_else(|_| {
                    usage("Invalid translation table specified.");
                    0
                });
                if tinf.trans_table < 1 || tinf.trans_table > 25 ||
                   tinf.trans_table == 7 || tinf.trans_table == 8 ||
                   (tinf.trans_table >= 17 && tinf.trans_table <= 20)
                {
                    usage("Invalid translation table specified.");
                }
                user_tt = tinf.trans_table;
                i += 1;
            }
            "-p" | "-P" => {
                let mode = &args[i + 1];
                let first_char = mode.chars().next().unwrap_or('x');
                if first_char == '0' || first_char == 's' || first_char == 'S' {
                    is_meta = 0;
                } else if first_char == '1' || first_char == 'm' || first_char == 'M' {
                    is_meta = 1;
                } else {
                    usage("Invalid meta/single genome type specified.");
                }
                i += 1;
            }
            "-f" | "-F" => {
                let fmt = &args[i + 1];
                if fmt.starts_with('0') || fmt.eq_ignore_ascii_case("gbk") {
                    output = 0;
                } else if fmt.starts_with('1') || fmt.eq_ignore_ascii_case("gca") {
                    output = 1;
                } else if fmt.starts_with('2') || fmt.eq_ignore_ascii_case("sco") {
                    output = 2;
                } else if fmt.starts_with('3') || fmt.eq_ignore_ascii_case("gff") {
                    output = 3;
                } else {
                    usage("Invalid output format specified.");
                }
                i += 1;
            }
            "-w" | "-W" => {
                num_threads = Some(args[i + 1].parse::<usize>().unwrap_or_else(|_| {
                    usage("Invalid number of threads specified.");
                    0
                }));
                i += 1;
            }
            _ => usage("Unknown option."),
        }

        i += 1;
    }

    if let Some(n) = num_threads {

        if is_meta == 0 && n > 4 && !quiet {
            eprintln!(
                "\nNote: single mode ('-p single') gains little beyond ~4 threads \
                 (gene-finding is limited by the serial training phase); -w {} will \
                 not run much faster than -w 4. For metagenomic / many-short-contig \
                 input, use '-p meta', which parallelizes well.\n",
                n
            );
        }
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .unwrap_or_else(|e| {
                eprintln!("\nError: failed to initialize thread pool with {} threads: {}\n", n, e);
                process::exit(5);
            });
    }

    if !quiet {
        eprintln!("-------------------------------------");
        eprintln!("RUSTYGAL v{} (Rust port of Prodigal)", VERSION);
        eprintln!("Univ of Tenn / Oak Ridge National Lab");
        eprintln!("Doug Hyatt, Loren Hauser, et al.     ");
        eprintln!("-------------------------------------");
    }

    if let Some(ref tf) = train_file {
        if is_meta == 1 {
            eprintln!("\nError: cannot specify metagenomic sequence with a training file.\n");
            process::exit(2);
        }

        let rv = training::read_training_file(tf, &mut tinf);
        if rv == 1 {
            do_training = true;
        } else {
            if force_nonsd == 1 {
                eprintln!("\nError: cannot force non-SD finder with a training file already created!\n");
                process::exit(3);
            }
            if !quiet {
                eprint!("Reading in training data from file {}...", tf);
            }
            if user_tt > 0 && user_tt != tinf.trans_table {
                eprintln!("\n\nWarning: user-specified translation table does not match");
                eprintln!("the one in the specified training file!\n");
            }
            if rv == -1 {
                eprintln!("\n\nError: training file did not read correctly!\n");
                process::exit(4);
            }
            if !quiet {
                eprintln!("done!");
                eprintln!("-------------------------------------");
            }
        }
    }

    if is_meta == 0 && train_file.is_none() && input_file.is_none() {
        use std::os::unix::io::AsRawFd;
        use std::os::unix::fs::MetadataExt;

        let stdin_fd = io::stdin().as_raw_fd();
        let metadata = std::fs::metadata(format!("/proc/self/fd/{}", stdin_fd));

        match metadata {
            Ok(meta) => {
                let mode = meta.mode();

                if (mode & 0o170000) == 0o020000 {
                    help();
                }

                else if (mode & 0o170000) == 0o010000 {
                    piped = true;
                    if copy_standard_input_to_file(&input_copy, quiet).is_err() {
                        eprintln!("\nError: can't copy stdin to file.\n");
                        process::exit(5);
                    }
                    input_file = Some(input_copy.clone());
                }

            }
            Err(_) => {
                eprintln!("\nError: can't fstat standard input.\n");
                process::exit(5);
            }
        }
    }

    let mut input_ptr: Option<FilePtr> = None;
    if let Some(ref inf) = input_file {
        match FilePtr::open(inf) {
            Ok(f) => input_ptr = Some(f),
            Err(_) => {
                eprintln!("\nError: can't open input file {}.\n", inf);
                process::exit(5);
            }
        }
    }
    if input_ptr.is_none() {
        match FilePtr::open("/dev/stdin") {
            Ok(f) => input_ptr = Some(f),
            Err(_) => {
                eprintln!("\nError: can't open stdin.\n");
                process::exit(5);
            }
        }
    }

    let mut output_ptr: File = if let Some(ref outf) = output_file {
        match File::create(outf) {
            Ok(f) => f,
            Err(_) => {
                eprintln!("\nError: can't open output file {}.\n", outf);
                process::exit(6);
            }
        }
    } else {

        match File::create("/dev/stdout") {
            Ok(f) => f,
            Err(_) => {
                eprintln!("\nError: can't open stdout.\n");
                process::exit(6);
            }
        }
    };

    let mut start_ptr: Option<File> = start_file.as_ref().and_then(|sf| {
        File::create(sf).ok().or_else(|| {
            eprintln!("\nError: can't open start file {}.\n", sf);
            process::exit(7);
        })
    });

    let mut trans_ptr: Option<File> = trans_file.as_ref().and_then(|tf| {
        File::create(tf).ok().or_else(|| {
            eprintln!("\nError: can't open translation file {}.\n", tf);
            process::exit(8);
        })
    });

    let mut nuc_ptr: Option<File> = nuc_file.as_ref().and_then(|nf| {
        File::create(nf).ok().or_else(|| {
            eprintln!("\nError: can't open gene nucleotide file {}.\n", nf);
            process::exit(16);
        })
    });

    if is_meta == 0 && (do_training || train_file.is_none()) {
        if !quiet {
            eprintln!("Request:  Single Genome, Phase:  Training");
            eprint!("Reading in the sequence(s) to train...");
        }

        let mut dna: Vec<u8> = Vec::new();
        slen = sequence::read_seq_training(
            input_ptr.as_mut().expect("input reader initialized above (opened or exited at startup)"),
            do_mask,
            &mut mlist,
            &mut nmask,
            &mut dna
        );

        if slen == 0 {
            eprintln!("\n\nSequence read failed (file must be Fasta, Genbank, or EMBL format).\n");
            process::exit(9);
        }

        if slen < MIN_SINGLE_GENOME {
            eprintln!("\n\nError:  Sequence must be {} characters (only {} read).",
                     MIN_SINGLE_GENOME, slen);
            eprintln!("(Consider running with the -p meta option or finding more contigs from the same genome.)\n");
            process::exit(10);
        }

        if slen < IDEAL_SINGLE_GENOME {
            eprintln!("\n\nWarning:  ideally Prodigal should be given at least {} bases for training.",
                     IDEAL_SINGLE_GENOME);
            eprintln!("You may get better results with the -p meta option.\n");
        }

        // Right-size the bitmaps to the training sequence's true length and build them through
        // the same build_bitmaps path (which also computes gc and the reverse-complement strand),
        // instead of a fixed MAX_SEQ buffer.
        let slen_usize = slen as usize;
        let mut seq = vec![0u8; slen_usize / 4 + 16];
        let mut rseq = vec![0u8; slen_usize / 4 + 16];
        let mut useq = vec![0u8; slen_usize / 8 + 16];
        let mut gc: f64 = 0.0;
        bitmap::build_bitmaps(&dna, slen, &mut seq, &mut rseq, &mut useq, &mut gc);
        tinf.gc = gc;

        if !quiet {
            eprintln!("{} bp seq created, {:.2} pct GC", slen, tinf.gc * 100.0);
        }

        if !quiet {
            eprint!("Locating all potential starts and stops...");
        }

        if slen > max_slen && slen > node::STT_NOD as i32 * 8 {
            let needed = (slen / 8) as usize;
            if needed > nodes.len() {
                let additional = needed - nodes.len();
                nodes.reserve_exact(additional);
                nodes.resize(needed, node::Node::default());
            }
        }

        nn = node::add_nodes(&seq, &rseq, slen, &mut nodes, closed, &mlist, nmask, &tinf);
        nodes[..nn as usize].sort_unstable_by(|a, b| node::compare_nodes(a, b));

        if !quiet {
            eprintln!("{} nodes", nn);
        }

        if !quiet {
            eprint!("Looking for GC bias in different frames...");
        }

        let gc_frame = sequence::calc_most_gc_frame(&seq, slen);
        node::record_gc_bias(&gc_frame, &mut nodes, nn, &mut tinf);

        if !quiet {
            eprintln!("frame bias scores: {:.2} {:.2} {:.2}",
                     tinf.bias[0], tinf.bias[1], tinf.bias[2]);
        }

        if !quiet {
            eprint!("Building initial set of genes to train from...");
        }

        node::record_overlapping_starts(&mut nodes, nn, &tinf, 0);
        ipath = dprog::dprog(&mut nodes, nn, &tinf, 0);

        if !quiet {
            eprintln!("done!");
        }

        if std::env::var("DEBUG_TRAINING").is_ok() {
            use crate::sequence::STOP;
            let mut path = ipath;
            let mut count = 0;
            eprintln!("TRAINING_GENES_START");
            while path != -1 {

                if nodes[path as usize].strand == 1 && nodes[path as usize].type_ != STOP {
                    eprintln!("{} {} {} {:.10}",
                              nodes[path as usize].ndx,
                              nodes[path as usize].stop_val,
                              nodes[path as usize].strand,
                              nodes[path as usize].cscore + nodes[path as usize].sscore);
                    count += 1;
                }

                if nodes[path as usize].strand == -1 && nodes[path as usize].type_ == STOP {
                    let start_node = nodes[path as usize].tracef;
                    eprintln!("{} {} {} {:.10}",
                              nodes[start_node as usize].ndx,
                              nodes[path as usize].ndx,
                              nodes[path as usize].strand,
                              nodes[start_node as usize].cscore + nodes[start_node as usize].sscore);
                    count += 1;
                }

                path = nodes[path as usize].traceb;
            }
            eprintln!("TRAINING_GENES_END total={}", count);
        }

        if !quiet {
            eprint!("Creating coding model and scoring nodes...");
        }

        node::calc_dicodon_gene(&mut tinf, &seq, &rseq, slen, &nodes, ipath);
        node::raw_coding_score(&seq, &rseq, slen, &mut nodes, nn, &tinf);

        if !quiet {
            eprintln!("done!");
        }

        if !quiet {
            eprint!("Examining upstream regions and training starts...");
        }

        node::rbs_score(&seq, &rseq, slen, &mut nodes, nn, &tinf);
        node::train_starts_sd(&seq, &rseq, slen, &mut nodes, nn, &mut tinf);
        node::determine_sd_usage(&mut tinf);

        if force_nonsd == 1 {
            tinf.uses_sd = 0;
        }

        if tinf.uses_sd == 0 {
            node::train_starts_nonsd(&seq, &rseq, slen, &mut nodes, nn, &mut tinf);
        }

        if !quiet {
            eprintln!("done!");
        }

        if do_training {
            if !quiet {
                eprint!("Writing data to training file {}...", train_file.as_ref().unwrap());
            }

            let rv = training::write_training_file(train_file.as_ref().unwrap(), &tinf);
            if rv != 0 {
                eprintln!("\nError: could not write training file!\n");
                process::exit(12);
            } else {
                if !quiet {
                    eprintln!("done!");
                }
                process::exit(0);
            }
        }

        if !quiet {
            eprintln!("-------------------------------------");
        }

        if let Some(ref inf) = input_file {
            match FilePtr::open(inf) {
                Ok(f) => input_ptr = Some(f),
                Err(_) => {
                    eprintln!("\nError: can't reopen input file {}.\n", inf);
                    process::exit(5);
                }
            }
        } else {
            match FilePtr::open("/dev/stdin") {
                Ok(f) => input_ptr = Some(f),
                Err(_) => {
                    eprintln!("\nError: can't reopen stdin.\n");
                    process::exit(5);
                }
            }
        }

        // The per-training bitmap buffers (seq/rseq/useq) are dropped at the end of this block;
        // gene finding below re-reads the input and allocates its own per-record buffers.
        for n in nodes.iter_mut() {
            *n = node::Node::default();
        }
    }

    if std::env::var("DUMP_GENE_DC").is_ok() {
        eprintln!("GENE_DC_DUMP_START");
        eprintln!("st_wt={:.10}", tinf.st_wt);
        eprintln!("type_wt[ATG]={:.10}", tinf.type_wt[0]);
        eprintln!("type_wt[GTG]={:.10}", tinf.type_wt[1]);
        eprintln!("type_wt[TTG]={:.10}", tinf.type_wt[2]);
        for i in 0..4096 {
            eprintln!("{} {:.10}", i, tinf.gene_dc[i]);
        }
        eprintln!("GENE_DC_DUMP_END");
    }

    if is_meta == 1 {
        if !quiet {
            eprintln!("Request:  Metagenomic, Phase:  Training");
            eprint!("Initializing training files...");
        }

        metagenomic::initialize_metagenomic_bins(&mut meta);

        if !quiet {
            eprintln!("done!");
            eprintln!("-------------------------------------");
        }
    }

    if !quiet {
        if is_meta == 1 {
            eprintln!("Request:  Metagenomic, Phase:  Gene Finding");
        } else {
            eprintln!("Request:  Single Genome, Phase:  Gene Finding");
        }
    }

    if !quiet {
        eprintln!("Reading sequences...");
    }

    let sequences = match read_all_sequences(input_ptr.as_mut().expect("input reader initialized above (opened or exited at startup)"), do_mask, quiet) {
        Ok(seqs) => seqs,
        Err(e) => {
            eprintln!("\nError reading sequences: {}\n", e);
            if piped {
                let _ = std::fs::remove_file(&input_copy);
            }
            process::exit(14);
        }
    };

    let num_sequences = sequences.len();
    if num_sequences == 0 {
        eprintln!("\nError: no sequences found in input file.\n");
        if piped {
            let _ = std::fs::remove_file(&input_copy);
        }
        process::exit(13);
    }

    if !quiet {
        eprintln!("Read {} sequence(s)", num_sequences);
        eprintln!("-------------------------------------");
    }

    if !quiet {
        eprintln!("Processing sequences in parallel...");
    }

    let results: Vec<SequenceResult> = sequences.par_iter()
        .map(|seq_data| {
            process_single_sequence(
                seq_data,
                &tinf,
                &meta,
                is_meta,
                closed,
                output,
                quiet,
                trans_file.is_some(),
                nuc_file.is_some(),
                start_file.is_some(),
            )
        })
        .collect();

    if !quiet {
        eprintln!("All sequences processed");
        eprintln!("-------------------------------------");
    }

    if !quiet {
        eprintln!("Writing output...");
    }

    match write_all_results(
        &results,
        &mut output_ptr,
        &mut trans_ptr,
        &mut nuc_ptr,
        &mut start_ptr,
    ) {
        Ok(_) => {},
        Err(e) => {
            eprintln!("\nError writing output: {}\n", e);
            if piped {
                let _ = std::fs::remove_file(&input_copy);
            }
            process::exit(16);
        }
    }

    if !quiet {
        eprintln!("done!");
    }

    num_seq = num_sequences as i32;

    if num_seq == 0 {
        eprintln!("\nError:  sequence read failed.\n");
        if piped {
            let _ = std::fs::remove_file(&input_copy);
        }
        process::exit(13);
    }

    if piped {
        let _ = std::fs::remove_file(&input_copy);
    }

    process::exit(0);
}

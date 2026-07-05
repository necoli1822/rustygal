// Verifies the library META API (`rustygal::meta_api::run_meta`) produces output
// byte-identical to the binary's `-p meta -a/-f gff`.
//
// Usage: meta_verify <input.fasta> <out.faa> <out.gff>
// Parses each contig, calls run_meta in input order (1-based seq_num), and writes
// the concatenated protein FASTA / gff. Also re-runs in parallel (rayon) and
// asserts the bytes match the serial run (thread-safety + determinism).

use rayon::prelude::*;
use rustygal::meta_api::{meta_bins, run_meta, MetaOutput};
use std::fs;

fn parse_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
    let data = fs::read(path).expect("read input");
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    let mut header: Option<String> = None;
    let mut dna: Vec<u8> = Vec::new();
    for line in data.split(|&b| b == b'\n') {
        if line.first() == Some(&b'>') {
            if let Some(h) = header.take() {
                out.push((h, std::mem::take(&mut dna)));
            }
            // header text after '>', trailing CR stripped
            let mut h = &line[1..];
            if h.last() == Some(&b'\r') {
                h = &h[..h.len() - 1];
            }
            header = Some(String::from_utf8_lossy(h).into_owned());
        } else {
            // keep only sequence letters, matching next_seq_multi's
            // `if ch < 'A' || ch > 'z' { continue }` filter
            for &c in line {
                if (b'A'..=b'z').contains(&c) {
                    dna.push(c);
                }
            }
        }
    }
    if let Some(h) = header.take() {
        out.push((h, dna));
    }
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: {} <input.fasta> <out.faa> <out.gff>", args[0]);
        std::process::exit(2);
    }
    let contigs = parse_fasta(&args[1]);
    let meta = meta_bins();

    // serial, input order
    let mut faa = Vec::new();
    let mut gff = Vec::new();
    let mut nuc = Vec::new();
    for (i, (h, d)) in contigs.iter().enumerate() {
        let o = run_meta(i as i32 + 1, h, d, &meta);
        faa.extend_from_slice(&o.trans_faa);
        gff.extend_from_slice(&o.gff);
        nuc.extend_from_slice(&o.nuc);
    }

    // parallel (proves Send+Sync + deterministic), then reassemble in order
    let par: Vec<MetaOutput> = contigs
        .par_iter()
        .enumerate()
        .map(|(i, (h, d))| run_meta(i as i32 + 1, h, d, &meta))
        .collect();
    let mut faa_par = Vec::new();
    let mut gff_par = Vec::new();
    for o in &par {
        faa_par.extend_from_slice(&o.trans_faa);
        gff_par.extend_from_slice(&o.gff);
    }
    assert!(faa == faa_par, "parallel faa != serial faa");
    assert!(gff == gff_par, "parallel gff != serial gff");
    eprintln!("parallel == serial OK ({} contigs)", contigs.len());

    fs::write(&args[2], &faa).expect("write faa");
    fs::write(&args[3], &gff).expect("write gff");
    fs::write(format!("{}.nuc", &args[2]), &nuc).expect("write nuc");
}

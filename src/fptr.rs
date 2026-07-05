// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2007-2016 Doug Hyatt, Univ. of Tennessee / UT-Battelle (original Prodigal C)
// Copyright (C) 2026 Sunju Kim (Rust reimplementation)

use std::fs::File;
use std::io::{BufRead, BufReader, Read};

#[cfg(feature = "gzip")]
use flate2::read::GzDecoder;

pub enum FilePtr {
    Regular(BufReader<File>),
    #[cfg(feature = "gzip")]
    Gzip(BufReader<GzDecoder<File>>),
}

impl FilePtr {

    pub fn open(path: &str) -> std::io::Result<Self> {
        let file = File::open(path)?;

        #[cfg(feature = "gzip")]
        {

            let mut reader = BufReader::new(file);
            let mut magic = [0u8; 2];

            {
                let buf = reader.fill_buf()?;
                if buf.len() >= 2 {
                    magic.copy_from_slice(&buf[0..2]);
                }
            }

            if magic[0] == 0x1f && magic[1] == 0x8b {

                let file = File::open(path)?;
                let gz = GzDecoder::new(file);
                Ok(FilePtr::Gzip(BufReader::new(gz)))
            } else {
                Ok(FilePtr::Regular(reader))
            }
        }

        #[cfg(not(feature = "gzip"))]
        {
            Ok(FilePtr::Regular(BufReader::new(file)))
        }
    }

    pub fn gets(&mut self, buf: &mut String) -> std::io::Result<usize> {
        self.read_line(buf)
    }
}

impl Read for FilePtr {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            FilePtr::Regular(f) => f.read(buf),
            #[cfg(feature = "gzip")]
            FilePtr::Gzip(f) => f.read(buf),
        }
    }
}

impl BufRead for FilePtr {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        match self {
            FilePtr::Regular(f) => f.fill_buf(),
            #[cfg(feature = "gzip")]
            FilePtr::Gzip(f) => f.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            FilePtr::Regular(f) => f.consume(amt),
            #[cfg(feature = "gzip")]
            FilePtr::Gzip(f) => f.consume(amt),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_open_regular_file() {

        let test_file = "/tmp/prodigal_test_regular.txt";
        let mut file = File::create(test_file).unwrap();
        writeln!(file, "Line 1").unwrap();
        writeln!(file, "Line 2").unwrap();
        writeln!(file, "Line 3").unwrap();
        drop(file);

        let mut fp = FilePtr::open(test_file).unwrap();

        let mut line1 = String::new();
        let mut line2 = String::new();
        let mut line3 = String::new();

        assert!(fp.gets(&mut line1).unwrap() > 0);
        assert_eq!(line1, "Line 1\n");

        assert!(fp.gets(&mut line2).unwrap() > 0);
        assert_eq!(line2, "Line 2\n");

        assert!(fp.gets(&mut line3).unwrap() > 0);
        assert_eq!(line3, "Line 3\n");

        let mut eof = String::new();
        assert_eq!(fp.gets(&mut eof).unwrap(), 0);

        std::fs::remove_file(test_file).unwrap();
    }

    #[test]
    #[cfg(feature = "gzip")]
    fn test_open_gzip_file() {
        use flate2::write::GzEncoder;
        use flate2::Compression;

        let test_file = "/tmp/prodigal_test_gzip.txt.gz";
        let file = File::create(test_file).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        writeln!(encoder, "Compressed Line 1").unwrap();
        writeln!(encoder, "Compressed Line 2").unwrap();
        encoder.finish().unwrap();

        let mut fp = FilePtr::open(test_file).unwrap();

        let mut line1 = String::new();
        let mut line2 = String::new();

        assert!(fp.gets(&mut line1).unwrap() > 0);
        assert_eq!(line1, "Compressed Line 1\n");

        assert!(fp.gets(&mut line2).unwrap() > 0);
        assert_eq!(line2, "Compressed Line 2\n");

        std::fs::remove_file(test_file).unwrap();
    }

    #[test]
    fn test_bufread_interface() {

        let test_file = "/tmp/prodigal_test_bufread.txt";
        let mut file = File::create(test_file).unwrap();
        writeln!(file, "Test line for BufRead").unwrap();
        drop(file);

        let mut fp = FilePtr::open(test_file).unwrap();

        let mut line = String::new();
        let bytes_read = fp.read_line(&mut line).unwrap();
        assert!(bytes_read > 0);
        assert_eq!(line, "Test line for BufRead\n");

        std::fs::remove_file(test_file).unwrap();
    }
}

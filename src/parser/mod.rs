//! Handles all the FASTA/FASTQ parsing
use std::fs::File;
use std::io::{stdin, Cursor, Read};
use std::path::Path;

#[cfg(feature = "bzip2")]
use bzip2::read::BzDecoder;
#[cfg(feature = "flate2")]
use flate2::read::MultiGzDecoder;
#[cfg(feature = "xz2")]
use liblzma::read::XzDecoder;
#[cfg(feature = "zstd")]
use zstd::stream::read::Decoder as ZstdDecoder;

use crate::errors::ParseError;
pub use crate::parser::fasta::Reader as FastaReader;
pub use crate::parser::fastq::Reader as FastqReader;

mod record;
mod utils;

mod fasta;
mod fastq;

pub use crate::parser::utils::FastxReader;

// Magic bytes for each compression format
#[cfg(feature = "flate2")]
const GZ_MAGIC: [u8; 2] = [0x1F, 0x8B];
#[cfg(feature = "bzip2")]
const BZ_MAGIC: [u8; 2] = [0x42, 0x5A];
#[cfg(feature = "xz2")]
const XZ_MAGIC: [u8; 2] = [0xFD, 0x37];
#[cfg(feature = "zstd")]
const ZST_MAGIC: [u8; 2] = [0x28, 0xB5];

fn get_fastx_reader<'a, R: 'a + io::Read + Send>(
    reader: R,
    first_byte: u8,
) -> Result<Box<dyn FastxReader + 'a>, ParseError> {
    match first_byte {
        b'>' => Ok(Box::new(FastaReader::new(reader))),
        b'@' => Ok(Box::new(FastqReader::new(reader))),
        _ => Err(ParseError::new_unknown_format(first_byte)),
    }
}

/// The main entry point of needletail if you're reading from something that implements [`std::io::Read`].
/// This automatically detects whether the file is:
/// 1. compressed: [`gzip`][gzip], [`bz`][bz], [`xz`][xz], and [`zstd`][zstd] are supported and will use the appropriate decoder
/// 2. FASTA or FASTQ: the right parser will be automatically instantiated
///
/// Option 1 is only available if the `compression` feature is enabled.
///
/// # Errors
///
/// If the object you're reading from has less than 2 bytes then a [`ParserError`](needletail::errors::ParserError) of the kind
/// [`ParseErrorKind::EmptyFile`](needletail::errors::ParseErrorKind::EmptyFile) is returned.
///
/// If the first byte in the object is unknown, then a `ParserError` of the kind
/// [`ParseErrorKind::UnknownFormat`](needletail::errors::ParseErrorKind::UnknownFormat) is returned.
///
/// # Examples
///
/// ```
/// use needletail::parse_fastx_reader;
///
/// let reader = ">read1\nACGT\nread2\nGGGG".as_bytes();
/// let mut fastx_reader = parse_fastx_reader(reader).expect("invalid reader");
/// let mut idx = 0;
/// let read_ids = [b"read1", b"read2"];
///
/// while let Some(r) = fastx_reader.next() {
///     let record = r.expect("invalid record");
///     assert_eq!(record.id(), read_ids[idx]);
///     idx += 1;
/// }
/// ```
///
/// [gzip]: https://www.gnu.org/software/gzip/
/// [bz]: https://sourceware.org/bzip2/
/// [xz]: https://tukaani.org/xz/format.html
/// [zstd]: https://facebook.github.io/zstd/
///
pub fn parse_fastx_reader<'a, R: 'a + io::Read + Send>(
    mut reader: R,
) -> Result<Box<dyn FastxReader + 'a>, ParseError> {
    let mut first_two_bytes = [0; 2];
    reader
        .read_exact(&mut first_two_bytes)
        .map_err(|_| ParseError::new_empty_file())?;
    let first_two_cursor = Cursor::new(first_two_bytes);
    let new_reader = first_two_cursor.chain(reader);

    match first_two_bytes {
        #[cfg(feature = "flate2")]
        GZ_MAGIC => {
            let mut gz_reader = MultiGzDecoder::new(new_reader);
            let mut first = [0; 1];
            gz_reader
                .read_exact(&mut first)
                .map_err(|e| match e.kind() {
                    io::ErrorKind::UnexpectedEof => ParseError::new_empty_file(),
                    _ => e.into(),
                })?;
            let r = Cursor::new(first).chain(gz_reader);
            get_fastx_reader(r, first[0])
        }
        #[cfg(feature = "bzip2")]
        BZ_MAGIC => {
            let mut bz_reader = BzDecoder::new(new_reader);
            let mut first = [0; 1];
            bz_reader
                .read_exact(&mut first)
                .map_err(|e| match e.kind() {
                    io::ErrorKind::UnexpectedEof => ParseError::new_empty_file(),
                    _ => e.into(),
                })?;
            let r = Cursor::new(first).chain(bz_reader);
            get_fastx_reader(r, first[0])
        }
        #[cfg(feature = "xz2")]
        XZ_MAGIC => {
            let mut xz_reader = XzDecoder::new(new_reader);
            let mut first = [0; 1];
            xz_reader
                .read_exact(&mut first)
                .map_err(|e| match e.kind() {
                    io::ErrorKind::UnexpectedEof => ParseError::new_empty_file(),
                    _ => e.into(),
                })?;
            let r = Cursor::new(first).chain(xz_reader);
            get_fastx_reader(r, first[0])
        }
        #[cfg(feature = "zstd")]
        ZST_MAGIC => {
            let mut zst_reader = ZstdDecoder::new(new_reader)?;
            let mut first = [0; 1];
            zst_reader
                .read_exact(&mut first)
                .map_err(|e| match e.kind() {
                    io::ErrorKind::UnexpectedEof => ParseError::new_empty_file(),
                    _ => e.into(),
                })?;
            let r = Cursor::new(first).chain(zst_reader);
            get_fastx_reader(r, first[0])
        }
        _ => get_fastx_reader(new_reader, first_two_bytes[0]),
    }
}

/// The main entry point of needletail if you're reading from stdin.
/// Shortcut to calling `parse_fastx_reader` with `stdin()`
pub fn parse_fastx_stdin() -> Result<Box<dyn FastxReader>, ParseError> {
    let stdin = stdin();
    parse_fastx_reader(stdin)
}

/// The main entry point of needletail if you're reading from a file.
/// Shortcut to calling `parse_fastx_reader` with a file
pub fn parse_fastx_file<P: AsRef<Path>>(path: P) -> Result<Box<dyn FastxReader>, ParseError> {
    parse_fastx_reader(File::open(&path)?)
}

pub use record::{mask_header_tabs, mask_header_utf8, write_fasta, write_fastq, SequenceRecord};
use std::io;
pub use utils::{Format, LineEnding};

#[cfg(test)]
mod test {
    use crate::errors::ParseErrorKind;
    use crate::parse_fastx_reader;
    use bzip2::read::BzEncoder;
    use bzip2::Compression as BzCompressionn;
    use flate2::write::GzEncoder;
    use flate2::Compression as GzCompression;
    use liblzma::write::XzEncoder;
    use zstd::stream::write::Encoder as ZstdEncoder;

    #[test]
    fn test_empty_file_raises_parser_error_of_same_kind() {
        let reader = "".as_bytes();
        let actual = parse_fastx_reader(reader);
        assert!(actual.is_err());

        let actual_err = actual.err().unwrap().kind;
        let expected_err = ParseErrorKind::EmptyFile;
        assert_eq!(actual_err, expected_err);
    }

    #[test]
    fn test_only_one_byte_in_file_raises_empty_file_error() {
        let reader = "@".as_bytes();
        let actual = parse_fastx_reader(reader);
        assert!(actual.is_err());

        let actual_err = actual.err().unwrap().kind;
        let expected_err = ParseErrorKind::EmptyFile;
        assert_eq!(actual_err, expected_err);
    }

    #[test]
    fn test_empty_gz_raises_empty_file_error() {
        let encoder = GzEncoder::new(Vec::new(), GzCompression::default());
        let compressed_bytes = encoder.finish().unwrap();
        let actual = parse_fastx_reader(compressed_bytes.as_slice());
        assert!(actual.is_err());

        let actual_err = actual.err().unwrap().kind;
        let expected_err = ParseErrorKind::EmptyFile;
        assert_eq!(actual_err, expected_err);
    }

    #[test]
    fn test_empty_bz_raises_empty_file_error() {
        let encoder = BzEncoder::new("".as_bytes(), BzCompressionn::default());
        let actual = parse_fastx_reader(encoder);
        assert!(actual.is_err());

        let actual_err = actual.err().unwrap().kind;
        let expected_err = ParseErrorKind::EmptyFile;
        assert_eq!(actual_err, expected_err);
    }

    #[test]
    fn test_empty_xz_raises_empty_file_error() {
        let encoder = XzEncoder::new(Vec::new(), 9);
        let compressed_bytes = encoder.finish().unwrap();
        let actual = parse_fastx_reader(compressed_bytes.as_slice());
        assert!(actual.is_err());

        let actual_err = actual.err().unwrap().kind;
        let expected_err = ParseErrorKind::EmptyFile;
        assert_eq!(actual_err, expected_err);
    }

    #[test]
    fn test_empty_zstd_raises_empty_file_error() {
        let encoder = ZstdEncoder::new(Vec::new(), zstd::DEFAULT_COMPRESSION_LEVEL).unwrap();
        let compressed_bytes = encoder.finish().unwrap();
        let actual = parse_fastx_reader(compressed_bytes.as_slice());
        assert!(actual.is_err());

        let actual_err = actual.err().unwrap().kind;
        let expected_err = ParseErrorKind::EmptyFile;
        assert_eq!(actual_err, expected_err);
    }
}

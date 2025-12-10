# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.7.1] - 2025-12-10

### Added
- Build Python 3.14 wheels (#109)

### Fixed
- Add missing `id`, `seq`, and `qual` attributes of the `Record` class to the Python stub file (#106)
- Fix multiline comment in `python.rs` module (#107)

## [0.7.0] - 2025-04-18

### Added

- Add the `phred_quality_score` method to the Python `Record` class (#98)
- Add Phred decoding functionality (#102)

## [0.4.1] - 2020-03-12

## Added
- `ParseErrorKind::EmptyFile` variant to handle cases where there are less than two bytes in a file [[#51][51]]

## [0.4.0] - 2020-07-08

## Changed
- Added `parse_fastx_file` which replaces `parse_sequence_reader` and offers an iterator like usage and is faster
than 0.3. Also adds `parse_fastx_reader` and `parse_fastx_stdin`.
- `SequenceRecord` now offers more information about the file such as line ending, which allows writing a file identical
to the input one.


## [0.3.0] - 2019-09-12
### Added
- Improved error reporting (i.e., a parse failure now gives the record it failed on).
- Significant code cleanup and additional linting (`cargo clippy`).
- Significant additional test coverage, including via fuzzing.
- Significant improvements to library documentation.

### Changed
- The `.kmers` method has been simplified and a new `.canonical_kmers` method has been introduced with much of the original's functionality.
- Added `parse_sequence_reader`, which replaces `fastx_stream` and `fastx_bytes`.
- `fastx_cli` updated and renamed to `parse_sequence_path`.
- `SeqRecord` is now `SequenceRecord` and many of its methods are now in the `Sequence` trait (e.g., working on byte slices).
- Automatic decompression now takes `Read` instead of `Read + Seek` so we can handle e.g. gzip files piped in through `stdin`.
- See [this link](https://github.com/onecodex/needletail/pull/26#issuecomment-530982670) for additional details on updating code to `v0.3.0`.

### Removed
- Single-file zip handling (zip requires `Seek`) 😞

## [0.3.1] - 2019-09-18
### Fixed
- Needletail no longer runs out of memory when parsing large, compressed files.

[51]: https://github.com/onecodex/needletail/issues/51

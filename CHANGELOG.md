# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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



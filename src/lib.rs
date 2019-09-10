#![crate_name = "needletail"]
pub mod bitkmer;
pub mod formats;
pub mod kmer;
pub mod sequence;
pub mod sequence_record;
mod util;

pub use formats::parse_sequences;
pub use sequence::Sequence;
pub use sequence_record::SequenceRecord;
pub use util::{ParseError, ParseErrorType};

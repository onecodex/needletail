#![no_main]
#[macro_use] extern crate libfuzzer_sys;
extern crate needletail;

use std::io::Cursor;

fuzz_target!(|data: &[u8]| {
    let cursor = Cursor::new([b"@", data].concat());
    let _ = needletail::parse_sequence_reader(cursor, |_ftype| {}, |_seq| {});
});

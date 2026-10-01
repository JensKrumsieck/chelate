#![no_main]

use chelate::FileType;
use std::io::BufReader;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let _ = chelate::parse(BufReader::new(data), FileType::MOL);
});

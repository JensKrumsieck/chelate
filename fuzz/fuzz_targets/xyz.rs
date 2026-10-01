#![no_main]

use chelate::FileType;
use libfuzzer_sys::fuzz_target;
use std::io::BufReader;

fuzz_target!(|data: &[u8]| {
    let _ = chelate::parse(BufReader::new(data), FileType::XYZ);
});

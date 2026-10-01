#![no_main]

use chelate::{FileType, atom::ToMolecule};
use libfuzzer_sys::fuzz_target;
use std::io::BufReader;

fuzz_target!(|data: &[u8]| {
    let Some((&selector, input)) = data.split_first() else {
        return;
    };
    let file_type = match selector % 5 {
        0 => FileType::CIF,
        1 => FileType::MOL,
        2 => FileType::MOL2,
        3 => FileType::PDB,
        _ => FileType::XYZ,
    };
    if let Ok(parsed) = chelate::parse(BufReader::new(input), file_type) {
        let _ = parsed.to_molecule();
    }
});

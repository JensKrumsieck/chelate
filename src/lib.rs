use error::FileError;
use std::{
    ffi::OsStr,
    fs::File,
    io::{self, BufReader, Read},
    path::Path,
    vec,
};
use types::{Atom, Bond};

pub mod cif;
pub mod error;
pub mod mol;
pub mod mol2;
pub mod pdb;
pub mod types;
pub mod xyz;

/// Parses a file based on the FileType and returns a vector of `Atom` and a vector of `Bond` objects.
/// # Examples
/// ```
/// use chelate;
/// let (atoms, bonds) = chelate::from_file("data/147288.cif").unwrap();
///
/// assert_eq!(atoms.len(), 206);
/// assert_eq!(bonds.len(), 230);
/// ```
pub fn from_file(filename: impl AsRef<Path>) -> Result<(Vec<Atom>, Vec<Bond>), FileError> {
    let file = File::open(&filename)?;
    let reader = BufReader::new(file);

    match filename.as_ref().extension().and_then(OsStr::to_str) {
        Some("cif") => parse(reader, FileType::CIF),
        Some("mol") => parse(reader, FileType::MOL),
        Some("mol2") => parse(reader, FileType::MOL2),
        Some("pdb") => parse(reader, FileType::PDB),
        Some("xyz") => parse(reader, FileType::XYZ),
        _ => Err(io::Error::new(io::ErrorKind::InvalidInput, "Unsupported file extension").into()),
    }
}

/// Enum to declare one of the supported chemical filetypes
#[non_exhaustive]
pub enum FileType {
    CIF,
    MOL,
    MOL2,
    PDB,
    XYZ,
}

/// Parses a file based on the FileType and returns a vector of `Atom` and a vector of `Bond` objects.
/// # Examples
/// ```
/// use chelate;
/// use chelate::FileType;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let file = File::open("data/147288.cif").unwrap();
/// let reader = BufReader::new(file);
/// let (atoms, bonds) = chelate::parse(reader, FileType::CIF).unwrap();
///
/// assert_eq!(atoms.len(), 206);
/// assert_eq!(bonds.len(), 230);
/// ```
pub fn parse<P: Read>(
    reader: BufReader<P>,
    type_: FileType,
) -> Result<(Vec<Atom>, Vec<Bond>), FileError> {
    match type_ {
        FileType::CIF => cif::parse(reader),
        FileType::MOL => mol::parse(reader),
        FileType::MOL2 => mol2::parse(reader),
        FileType::PDB => Ok((pdb::parse(reader)?, vec![])),
        FileType::XYZ => Ok((xyz::parse(reader)?, vec![])),
    }
}

/// Returns the trimmed content of the zero-based column range `start..end` of a fixed-width line.
/// Returns the available part if the line ends within the range and an empty string if it ends before.
fn column(line: &str, start: usize, end: usize) -> &str {
    line.get(start..end.min(line.len()))
        .unwrap_or_default()
        .trim()
}

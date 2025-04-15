use atom::{Atom, Bond, Molecule, ToMolecule};
use std::{
    ffi::OsStr,
    fs::File,
    io::{self, BufReader, Read},
    path::Path,
    vec,
};

pub mod atom;
pub mod cif;
pub mod mol;
pub mod mol2;
pub mod pdb;
pub mod xyz;


/// Parses a file based on the FileType and returns a Molecule type.
/// # Examples
/// ```
/// use chelate;
/// let mol = chelate::molecule_from_file("data/147288.cif").unwrap();
///
/// assert_eq!(mol.node_count(), 206);
/// assert_eq!(mol.edge_count(), 230);
/// ```
#[cfg(feature = "petgraph")]
pub fn molecule_from_file(filename: impl AsRef<Path>) -> io::Result<Molecule> {
    Ok(from_file(filename)?.to_molecule())
}

/// Parses a file based on the FileType and returns a vector of `Atom` and a vector of `Bond` objects.
/// # Examples
/// ```
/// use chelate;
/// let (atoms, bonds) = chelate::from_file("data/147288.cif").unwrap();
///
/// assert_eq!(atoms.len(), 206);
/// assert_eq!(bonds.len(), 230);
/// ```
pub fn from_file(filename: impl AsRef<Path>) -> io::Result<(Vec<Atom>, Vec<Bond>)> {
    let file = File::open(&filename)?;
    let reader = BufReader::new(file);

    match filename.as_ref().extension().and_then(OsStr::to_str) {
        Some("cif") => parse(reader, FileType::CIF),
        Some("mol") => parse(reader, FileType::MOL),
        Some("mol2") => parse(reader, FileType::MOL2),
        Some("pdb") => parse(reader, FileType::PDB),
        Some("xyz") => parse(reader, FileType::XYZ),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Unsupported file extension",
        )),
    }
}

/// Enum to declare one of the supported chemical filetypes
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
pub fn parse<P: Read>(reader: BufReader<P>, type_: FileType) -> io::Result<(Vec<Atom>, Vec<Bond>)> {
    match type_ {
        FileType::CIF => cif::parse(reader),
        FileType::MOL => mol::parse(reader),
        FileType::MOL2 => mol2::parse(reader),
        FileType::PDB => Ok((pdb::parse(reader)?, vec![])),
        FileType::XYZ => Ok((xyz::parse(reader)?, vec![])),
    }
}

fn normalize_symbol(symbol: &str) -> String {
    let normalized_symbol = if let Some(first_char) = symbol.chars().next() {
        first_char.to_uppercase().collect::<String>() + &symbol[1..].to_lowercase()
    } else {
        String::new()
    };
    normalized_symbol
}

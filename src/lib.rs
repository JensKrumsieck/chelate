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

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case("data/4n4n.cif", 15450, 14968)]
    #[case("data/4r21.cif", 6752, 6902)]
    #[case("data/147288.cif", 206, 230)]
    #[case("data/1484829.cif", 466, 528)]
    #[case("data/cif_noTrim.cif", 79, 89)]
    #[case("data/cif.cif", 79, 89)]
    #[case("data/CuHETMP.cif", 85, 92)]
    #[case("data/ligand.cif", 44, 46)]
    #[case("data/mmcif.cif", 1291, 1256)]
    #[case("data/benzene_3d.mol", 12, 12)]
    #[case("data/benzene_arom.mol", 12, 12)]
    #[case("data/benzene.mol", 6, 6)]
    #[case("data/tep.mol", 46, 50)]
    #[case("data/corrole.mol", 37, 41)]
    #[case("data/0001.mol2", 15450, 14898)]
    #[case("data/benzene.mol2", 12, 12)]
    #[case("data/myo.mol2", 1437, 1312)]
    #[case("data/ptcor.mol2", 129, 127)]
    #[case("data/tep.mol2", 46, 50)]
    #[case("data/VATTOC.mol2", 130, 146)]
    #[case("data/oriluy.pdb", 130, 151)]
    #[case("data/2spl.pdb", 1437, 1314)]
    #[case("data/1hv4.pdb", 9288, 9562)]
    #[case("data/0001.pdb", 15450, 14968)]
    #[case("data/cif.xyz", 102, 155)]
    #[case("data/mescho.xyz", 23, 23)]
    #[case("data/porphyrin.xyz", 37, 0)]
    fn test_molecule(
        #[case] filename: &str,
        #[case] atoms_count: usize,
        #[case] bonds_count: usize,
    ) {
        let mol = molecule_from_file(filename).unwrap();
        assert_eq!(mol.node_count(), atoms_count);
        assert_eq!(mol.edge_count(), bonds_count);
    }
}

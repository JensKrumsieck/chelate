use std::{
    ffi::OsStr,
    fs::File,
    io::{self, BufReader, Read},
    path::Path,
    vec,
};

pub mod cif;
pub mod mol;
pub mod mol2;
pub mod pdb;
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

pub(crate) fn normalize_symbol(symbol: &str) -> String {
    let normalized_symbol = if let Some(first_char) = symbol.chars().next() {
        first_char.to_uppercase().collect::<String>() + &symbol[1..].to_lowercase()
    } else {
        String::new()
    };
    normalized_symbol
}

#[derive(Debug)]
pub struct Atom {
    pub id: usize,
    pub atomic_number: u8,
    pub is_disordered: bool,
    pub coord: [f32; 3],
}

impl Atom {
    pub fn new(id: usize, atomic_number: u8, x: f32, y: f32, z: f32) -> Self {
        Atom {
            id,
            is_disordered: false,
            atomic_number,
            coord: [x, y, z],
        }
    }
}

static ATOMIC_SYMBOLS: [&str; 118] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
    "Fl", "Mc", "Lv", "Ts", "Og",
];

#[derive(Debug)]
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub order: u8,
    pub is_aromatic: bool,
}

impl Bond {
    pub fn new(atom1: Atom, atom2: Atom, order: u8, is_aromatic: bool) -> Self {
        Bond {
            atom1: atom1.id,
            atom2: atom2.id,
            order,
            is_aromatic,
        }
    }
}

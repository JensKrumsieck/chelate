use crate::types::Atom;

#[derive(Debug, PartialEq, Clone, Copy, PartialOrd)]
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub order: BondOrder,
}

impl Bond {
    pub fn new(atom1: &Atom, atom2: &Atom, order: BondOrder) -> Self {
        Bond {
            atom1: atom1.id,
            atom2: atom2.id,
            order,
        }
    }

    pub fn from_atoms(atoms: &[Atom]) -> Vec<Self> {
        #[cfg(feature = "rayon")]
        if atoms.len() < 600 {
            bond_from_atoms(atoms)
        } else {
            bond_from_atoms_parallel(atoms)
        }

        #[cfg(not(feature = "rayon"))]
        bond_from_atoms(atoms)
    }
}

fn bond_from_atoms(atoms: &[Atom]) -> Vec<Bond> {
    let mut bonds = Vec::with_capacity(atoms.len() * 3);

    for i in 0..atoms.len() {
        let atom_i = &atoms[i];
        for atom_j in &atoms[i + 1..] {
            if atom_i.bond_to_by_covalent_radii(atom_j, 25.0) {
                bonds.push(Bond::new(atom_i, atom_j, BondOrder::Single));
            }
        }
    }

    bonds
}

#[cfg(feature = "rayon")]
fn bond_from_atoms_parallel(atoms: &[Atom]) -> Vec<Bond> {
    use rayon::iter::{IntoParallelIterator, ParallelIterator};

    let n = atoms.len();
    (0..n)
        .into_par_iter()
        .flat_map_iter(|i| {
            let atom_i = &atoms[i];
            (i + 1..n).filter_map(move |j| {
                let atom_j = &atoms[j];
                if atom_i.bond_to_by_covalent_radii(atom_j, 25.0) {
                    Some(Bond::new(atom_i, atom_j, BondOrder::Single))
                } else {
                    None
                }
            })
        })
        .collect()
}

#[derive(Debug, PartialEq, Clone, Copy, PartialOrd)]
#[repr(u8)]
pub enum BondOrder {
    // 1 - 8 SDL Types
    Single = 1,
    Double = 2,
    Triple = 3,
    Aromatic = 4,
    SingleOrDouble = 5,
    SingleOrAromatic = 6,
    DoubleOrAromatic = 7,
    Any = 8,
    Coordination = 9,
    Hydrogen = 10,

    // additional sybyl types
    Amide = 11,
    Dummy = 12,
    Unknown = 13,
    NotConnected = 15,
}

impl BondOrder {
    pub fn from_sdf(value: u8) -> Option<Self> {
        match value {
            1 => Some(BondOrder::Single),
            2 => Some(BondOrder::Double),
            3 => Some(BondOrder::Triple),
            4 => Some(BondOrder::Aromatic),
            5 => Some(BondOrder::SingleOrDouble),
            6 => Some(BondOrder::SingleOrAromatic),
            7 => Some(BondOrder::DoubleOrAromatic),
            8 => Some(BondOrder::Any),
            9 => Some(BondOrder::Coordination),
            10 => Some(BondOrder::Hydrogen),
            _ => None,
        }
    }

    pub fn to_sdf(&self) -> u8 {
        let num = *self as u8;
        if num <= 10 { num } else { 8 }
    }

    pub fn from_sybyl(value: &str) -> Option<Self> {
        match value {
            "1" => Some(BondOrder::Single),
            "2" => Some(BondOrder::Double),
            "3" => Some(BondOrder::Triple),
            "am" => Some(BondOrder::Amide),
            "ar" => Some(BondOrder::Aromatic),
            "du" => Some(BondOrder::Dummy),
            "un" => Some(BondOrder::Unknown),
            "nc" => Some(BondOrder::NotConnected),
            _ => None,
        }
    }

    pub fn to_sybyl(&self) -> &'static str {
        match &self {
            Self::Single => "1",
            Self::Double => "2",
            Self::Triple => "3",
            Self::Aromatic => "ar",
            Self::SingleOrDouble => "2",
            Self::SingleOrAromatic | Self::DoubleOrAromatic => "ar",
            Self::Amide => "am",
            Self::Dummy => "du",
            Self::NotConnected => "nc",
            _ => "un",
        }
    }
}

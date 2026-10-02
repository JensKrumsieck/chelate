#[cfg(feature = "rayon")]
use crate::types::Atom;

#[derive(Debug)]
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub order: u8,
    pub is_aromatic: bool,
}

impl Bond {
    pub fn new(atom1: &Atom, atom2: &Atom, order: u8, is_aromatic: bool) -> Self {
        Bond {
            atom1: atom1.id,
            atom2: atom2.id,
            order,
            is_aromatic,
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
                bonds.push(Bond::new(atom_i, atom_j, 1, false));
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
                    Some(Bond::new(atom_i, atom_j, 1, false))
                } else {
                    None
                }
            })
        })
        .collect()
}

use nalgebra::{Point3, point};
#[cfg(feature = "petgraph")]
use petgraph::{Graph, Undirected};
use std::ops::{Deref, DerefMut};

#[cfg(feature = "petgraph")]
pub type Molecule = Graph<Atom, Bond, Undirected>;

#[cfg(feature = "petgraph")]
pub trait ToMolecule<T> {
    fn to_molecule(self) -> Molecule;
}

#[cfg(feature = "petgraph")]
impl ToMolecule<Self> for (Vec<Atom>, Vec<Bond>) {
    /// Takes a Tuple of atoms and bonds and converts it into a Molecule (UnGraph)
    /// Be aware: If no bonds can be found, this method will try to generate them!
    fn to_molecule(self) -> Molecule {
        let (atoms, mut bonds) = self;

        //generate bonds if none
        if bonds.is_empty() {
            bonds = Bond::from_atoms(&atoms);
        }

        let mut mol = Molecule::with_capacity(atoms.len(), bonds.len());
        for atom in atoms {
            mol.add_node(atom);
        }
        mol.extend_with_edges(bonds);

        mol
    }
}

#[cfg(feature = "petgraph")]
impl ToMolecule<Self> for Vec<Atom> {
    /// Converts a tuple of atoms and bonds into a `Molecule`.
    ///
    /// If no bonds are provided, they will be generated automatically.
    fn to_molecule(self) -> Molecule {
        let bonds = Bond::from_atoms(&self);
        (self, bonds).to_molecule()
    }
}

#[derive(Debug, Default, PartialEq)]
pub struct Atom {
    pub id: usize,
    pub atomic_number: u8,
    pub data: Box<AtomData>,
    pub coord: Point3<f32>,
}

impl Atom {
    pub fn new(id: usize, atomic_number: u8, x: f32, y: f32, z: f32) -> Self {
        Atom {
            id,
            atomic_number,
            coord: point![x, y, z],
            data: Default::default(),
        }
    }

    /// Checks whether atom is bond to `rhs` by their covalent radii while allowing a delta
    pub fn bond_to_by_covalent_radii(&self, rhs: &Atom, delta: f32) -> bool {
        let dist = nalgebra::distance_squared(&self.coord, &rhs.coord);

        //fast check -> take highest covalent radius ~((260+260+100)/100)²
        const FAST_FAIL: f32 = 38.0;
        if dist > FAST_FAIL {
            return false;
        }

        let check = (COVALENT_RADII_PM[self.atomic_number as usize - 1] as f32
            + COVALENT_RADII_PM[rhs.atomic_number as usize - 1] as f32
            + delta)
            / 100.0;
        dist < check * check
    }
}

impl Deref for Atom {
    type Target = AtomData;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl DerefMut for Atom {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}

#[derive(Debug, PartialEq)]
pub struct AtomData {
    pub name: String,
    pub resname: String,
    pub resid: i32,
    pub chain: char,
    pub disorder_group: usize,
    pub occupancy: f32,
}

impl Default for AtomData {
    fn default() -> Self {
        Self {
            name: Default::default(),
            resname: "UNK".to_string(),
            resid: Default::default(),
            chain: Default::default(),
            disorder_group: Default::default(),
            occupancy: 1.0,
        }
    }
}

pub(crate) static ATOMIC_SYMBOLS: [&str; 118] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
    "Fl", "Mc", "Lv", "Ts", "Og",
];

pub(crate) static COVALENT_RADII_PM: [u32; 118] = [
    31, 28, 128, 96, 84, 77, 71, 66, 64, 58, 166, 141, 121, 111, 107, 105, 102, 106, 203, 176, 170,
    160, 153, 139, 139, 132, 126, 124, 132, 122, 122, 122, 119, 120, 120, 116, 220, 195, 190, 175,
    164, 154, 147, 146, 142, 139, 145, 144, 142, 139, 139, 138, 139, 140, 244, 215, 207, 204, 203,
    201, 199, 198, 198, 196, 194, 192, 192, 189, 190, 187, 187, 175, 170, 162, 151, 144, 141, 136,
    136, 132, 145, 146, 148, 140, 150, 150, 260, 221, 215, 206, 200, 196, 190, 187, 180, 169, 168,
    168, 165, 167, 173, 176, 161, 157, 149, 143, 141, 134, 129, 128, 121, 122, 172, 171, 156, 162,
    156, 157,
];

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

#[cfg(feature = "petgraph")]
impl petgraph::IntoWeightedEdge<Bond> for Bond {
    type NodeId = u32;

    fn into_weighted_edge(self) -> (Self::NodeId, Self::NodeId, Bond) {
        //subtract 1 as atom ids start at 1 and index expects 0
        (self.atom1 as u32 - 1, self.atom2 as u32 - 1, self)
    }
}

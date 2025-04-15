use nalgebra::Vector3;
#[cfg(feature = "petgraph")]
use petgraph::{Graph, Undirected};
use std::ops::{Deref, DerefMut};

#[cfg(feature = "petgraph")]
pub type Molecule = Graph<Atom, Bond, Undirected>;

#[derive(Debug, Default, PartialEq)]
pub struct Atom {
    pub id: usize,
    pub atomic_number: u8,
    pub data: Box<AtomData>,
    pub coord: Vector3<f32>,
}

impl Atom {
    pub fn new(id: usize, atomic_number: u8, x: f32, y: f32, z: f32) -> Self {
        Atom {
            id,
            atomic_number,
            coord: Vector3::new(x, y, z),
            data: Default::default(),
        }
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

#[cfg(feature = "petgraph")]
impl petgraph::IntoWeightedEdge<Bond> for Bond {
    type NodeId = u32;

    fn into_weighted_edge(self) -> (Self::NodeId, Self::NodeId, Bond) {
        (self.atom1 as u32, self.atom2 as u32, self)
    }
}
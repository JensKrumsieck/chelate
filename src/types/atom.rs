use nalgebra::{Point3, point};
use smol_str::SmolStr;

#[derive(Debug, Default, PartialEq)]
pub struct Atom {
    pub id: usize,
    pub atomic_number: u8,
    pub data: AtomData,
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

#[derive(Debug, PartialEq)]
pub struct AtomData {
    pub name: SmolStr,
    pub resname: SmolStr,
    pub resid: i32,
    pub chain: SmolStr,
    pub disorder_group: usize,
    pub occupancy: f32,
}

impl Default for AtomData {
    fn default() -> Self {
        Self {
            name: Default::default(),
            resname: "UNK".into(),
            resid: Default::default(),
            chain: Default::default(),
            disorder_group: Default::default(),
            occupancy: 1.0,
        }
    }
}

pub static ATOMIC_SYMBOLS: [&str; 118] = [
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

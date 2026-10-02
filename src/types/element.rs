#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct Element(u8);

impl Element {
    pub fn new(atomic_number: u8) -> Self {
        Self(atomic_number)
    }

    pub fn atomic_number(&self) -> u8 {
        self.0
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

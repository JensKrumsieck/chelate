use nalgebra::{Point3, point};
use smol_str::SmolStr;

use crate::types::{Element, element::COVALENT_RADII_PM};

#[derive(Debug, Default, PartialEq)]
pub struct Atom {
    pub id: usize,
    pub symbol: Element,
    pub data: AtomData,
    pub coord: Point3<f32>,
}

impl Atom {
    pub fn new(id: usize, atomic_number: u8, x: f32, y: f32, z: f32) -> Self {
        Atom {
            id,
            symbol: Element::new(atomic_number),
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

        let check = (COVALENT_RADII_PM[self.symbol.atomic_number() as usize - 1] as f32
            + COVALENT_RADII_PM[rhs.symbol.atomic_number() as usize - 1] as f32
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

use crate::types::Cell;
use smol_str::SmolStr;

/// Crystallographic information of a structure, to be attached to a molecule.
/// Contains no atoms or bonds.
#[derive(Debug, PartialEq, Clone)]
pub struct Crystal {
    pub cell: Cell,
    pub space_group: SpaceGroup,
}

/// Space group of a crystal, any of the fields can be missing in a file
#[derive(Debug, Default, PartialEq, Clone)]
pub struct SpaceGroup {
    /// Hermann-Mauguin symbol, e.g. `P 21/c`
    pub symbol: Option<SmolStr>,
    /// International Tables number, 1 to 230
    pub number: Option<u8>,
    /// Symmetry operations as in the file, e.g. `-x, y+1/2, -z+1/2`
    pub symops: Vec<SmolStr>,
}

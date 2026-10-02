mod atom;
mod bond;
mod cell;
mod crystal;
mod element;

use smol_str::SmolStr;
use thiserror::Error;

pub use atom::Atom;
pub use bond::{Bond, BondOrder};
pub use cell::Cell;
pub use crystal::{Crystal, SpaceGroup};
pub use element::{ATOMIC_SYMBOLS, Element};

#[derive(Error, Debug)]
#[error("Could not parse Element: {0}")]
pub struct InvalidElementError(SmolStr);

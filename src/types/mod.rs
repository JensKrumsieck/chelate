mod atom;
mod bond;
mod element;

use smol_str::SmolStr;
use thiserror::Error;

pub use atom::Atom;
pub use bond::Bond;
pub use element::{ATOMIC_SYMBOLS, Element};

#[derive(Error, Debug)]
#[error("Could not parse Element: {0}")]
pub struct InvalidElementError(SmolStr);

//! Functions for parsing XYZ files (chemical/x-xyz)
//! XYZ format is the simplest file format just containing an atom symbol and XYZ cartesian coordinates
//! Documentation can be found here <https://en.wikipedia.org/wiki/XYZ_file_format>
use super::atomic_number;
use crate::atom::Atom;
use crate::error::{FileError, ParseError};
use std::io::{BufRead, BufReader, Read};

/// Parses a single line of an XYZ file and returns an `Atom` object, `None` for unknown elements.
/// The line should contain the atomic symbol followed by the x, y, and z coordinates.
/// Example line: `C 1.0 2.0 3.0`
fn parse_atom_line(line: &str, atom_count: &mut usize) -> Result<Option<Atom>, ParseError> {
    let mut iter = line.split_whitespace();
    let symbol = iter
        .next()
        .ok_or_else(|| ParseError::new("Missing atom"))?;
    let mut coord = || {
        iter.next()
            .and_then(|s| s.parse::<f32>().ok())
            .ok_or_else(|| ParseError::new("Invalid coordinates"))
    };
    let (x, y, z) = (coord()?, coord()?, coord()?);

    let Some(atomic_number) = atomic_number(symbol) else {
        return Ok(None);
    };
    *atom_count += 1;

    let mut atom = Atom::new(*atom_count, atomic_number, x, y, z);
    atom.name = symbol.into();

    Ok(Some(atom))
}

/// Parses an XYZ file and returns a vector of `Atom` objects.
/// # Examples
/// ```
/// use chelate::xyz;
/// use std::fs::File;
/// use nalgebra::Point3;
/// use std::io::BufReader;
///
/// let file = File::open("data/mescho.xyz").unwrap();
/// let reader = BufReader::new(file);
/// let atoms = xyz::parse(reader).unwrap();
///
/// assert_eq!(atoms.len(), 23);
/// assert_eq!(atoms[0].atomic_number, 6);
/// assert_eq!(atoms[0].coord, Point3::new(0.85246046633891, -1.08114766176821, 0.02536743820348));
/// assert_eq!(atoms[0].resname, "UNK");
/// assert_eq!(atoms[0].resid, 0);
/// assert_eq!(atoms[0].chain, char::default());
/// assert_eq!(atoms[0].occupancy, 1.0);
/// assert_eq!(atoms[0].name, "C");
/// ```
pub fn parse<P: Read>(reader: BufReader<P>) -> Result<Vec<Atom>, FileError> {
    let mut lines = reader.lines();
    // the number of atoms (after a byte order mark in some files) is followed by a comment line
    let count_line = lines.next().transpose()?.unwrap_or_default();
    let atom_len = count_line
        .trim_start_matches('\u{feff}')
        .split_whitespace()
        .next()
        .and_then(|s| s.parse::<usize>().ok())
        .ok_or_else(|| ParseError::new("Invalid atom count").with_line(1))?;
    lines.next().transpose()?;

    let mut atom_count = 0;
    let mut atoms = Vec::new();
    let mut line_count = 0;
    // only the first frame of files with multiple frames is read
    for (i, line) in lines.take(atom_len).enumerate() {
        line_count += 1;
        let atom = parse_atom_line(&line?, &mut atom_count).map_err(|e| e.with_line(i + 3))?;
        atoms.extend(atom);
    }
    if line_count < atom_len {
        let message = format!("Expected {atom_len} atoms, found {line_count}");
        return Err(ParseError::new(message).into());
    }
    Ok(atoms)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::fs::File;

    #[rstest]
    #[case("data/cif.xyz", 102)]
    #[case("data/mescho.xyz", 23)]
    #[case("data/porphyrin.xyz", 37)]
    fn test_xyz_files(#[case] filename: &str, #[case] len: usize) {
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let atoms = parse(reader).unwrap();

        assert_eq!(atoms.len(), len);
    }
}

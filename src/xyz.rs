use crate::{ATOMIC_SYMBOLS, Atom};
use std::io::{self, BufRead, BufReader, Read};

/// Parses a single line of an XYZ file and returns an `Atom` object.
/// The line should contain the atomic symbol followed by the x, y, and z coordinates.
/// Example line: `C 1.0 2.0 3.0`
fn parse_atom_line(line: &str, atom_count: &mut usize) -> Option<Atom> {
    let mut iter = line.split_whitespace();
    let symbol = iter.next()?;
    let atomic_number = ATOMIC_SYMBOLS.iter().position(|&s| s == symbol)? + 1;

    let x = iter.next()?.parse().ok()?;
    let y = iter.next()?.parse().ok()?;
    let z = iter.next()?.parse().ok()?;
    *atom_count += 1;
    Some(Atom::new(*atom_count, atomic_number as u8, x, y, z))
}

/// Parses an XYZ file and returns a vector of `Atom` objects.
/// # Examples
/// ```
/// use chelate::xyz;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let file = File::open("data/mescho.xyz").unwrap();
/// let reader = BufReader::new(file);
/// let atoms = xyz::parse(reader).unwrap();
///
/// assert_eq!(atoms.len(), 23);
/// ```
pub fn parse<P: Read>(reader: BufReader<P>) -> io::Result<Vec<Atom>> {
    let mut atom_count = 0;

    let mut atoms = Vec::new();
    for line in reader.lines().skip(2) {
        let line = line?;
        if let Some(atom) = parse_atom_line(&line, &mut atom_count) {
            atoms.push(atom);
        }
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

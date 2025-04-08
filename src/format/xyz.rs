use crate::{ATOMIC_SYMBOLS, Atom};
use std::io::{self, BufRead, BufReader, Read};

/// Parses a single line of an XYZ file and returns an `Atom` object.
/// The line should contain the atomic symbol followed by the x, y, and z coordinates.
/// Example line: "C 1.0 2.0 3.0"
/// # Examples
/// ```
/// use chelate::format::xyz::parse_line;
///
/// let line = "C 1.0 2.0 3.0";
/// let atom = parse_line(line).unwrap();
///
/// assert_eq!(atom.atomic_number, 6); // Carbon
/// assert_eq!(atom.coord[0], 1.0);
/// assert_eq!(atom.coord[1], 2.0);
/// assert_eq!(atom.coord[2], 3.0);
/// ```
///
pub fn parse_line(line: &str) -> Option<Atom> {
    let mut iter = line.split_whitespace();
    let symbol = iter.next()?;
    let atomic_number = ATOMIC_SYMBOLS.iter().position(|&s| s == symbol)? + 1;

    let x = iter.next()?.parse().ok()?;
    let y = iter.next()?.parse().ok()?;
    let z = iter.next()?.parse().ok()?;

    Some(Atom::new(atomic_number as u8, x, y, z))
}

/// Parses an XYZ file and returns a vector of `Atom` objects.
/// # Examples
/// ```
/// use chelate::format::xyz;
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
    let mut atoms = Vec::new();
    for line in reader.lines().skip(2) {
        let line = line?;
        if let Some(atom) = parse_line(&line) {
            atoms.push(atom);
        }
    }
    Ok(atoms)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use rstest::rstest;

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

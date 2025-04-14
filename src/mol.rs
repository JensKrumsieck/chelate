use super::normalize_symbol;
use crate::{ATOMIC_SYMBOLS, Atom, Bond};
use std::io::{self, BufRead, BufReader, Read};

/// Parses a single line of an MOL file and returns an `Atom` object.
/// The line should contain the x, y, and z coordinates followed by the atomic symbol.
/// Example line: "    1.3194   -1.2220   -0.8506 N   0  0  0  0  0  0  0  0  0  0  0  0"
fn parse_atom_line(line: &str, atom_count: &mut usize) -> Option<Atom> {
    let mut iter = line.split_whitespace();

    let x = iter.next()?.parse().ok()?;
    let y = iter.next()?.parse().ok()?;
    let z = iter.next()?.parse().ok()?;

    let symbol = iter.next()?;
    let atomic_number = ATOMIC_SYMBOLS
        .iter()
        .position(|&s| s == normalize_symbol(symbol))?
        + 1;
    *atom_count += 1;
    Some(Atom::new(*atom_count, atomic_number as u8, x, y, z))
}

/// Parses a single line of an MOL file and returns a `Bond` object.
/// The line should contain the atoms ids and the bond order where 4 is aromatic bond.
/// Example line: "  1  2  2  0  0  0  0"
fn parse_bond_line(line: &str) -> Option<Bond> {
    let mut iter = line.split_whitespace();

    let atom1 = iter.next()?.parse().ok()?;
    let atom2 = iter.next()?.parse().ok()?;
    let mut order = iter.next()?.parse().ok()?;
    let is_aromatic = order == 4;
    if is_aromatic {
        order = 1;
    }

    Some(Bond {
        atom1,
        atom2,
        order,
        is_aromatic,
    })
}

/// Parses a mol file and returns a vector of `Atom` and a vector of `Bond` objects.
/// # Examples
/// ```
/// use chelate::mol;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let file = File::open("data/corrole.mol").unwrap();
/// let reader = BufReader::new(file);
/// let (atoms, bonds) = mol::parse(reader).unwrap();
///
/// assert_eq!(atoms.len(), 37);
/// assert_eq!(bonds.len(), 41);
/// ```
pub fn parse<P: Read>(reader: BufReader<P>) -> io::Result<(Vec<Atom>, Vec<Bond>)> {
    let mut atom_count = 0;

    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    for line in reader.lines().skip(4) {
        let line = line?;
        if line.contains("M  END") {
            break;
        }

        if let Some(atom) = parse_atom_line(&line, &mut atom_count) {
            atoms.push(atom);
        }
        if let Some(bond) = parse_bond_line(&line) {
            bonds.push(bond);
        }
    }
    Ok((atoms, bonds))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::fs::File;

    #[rstest]
    #[case("data/benzene_3d.mol", 12, 12)]
    #[case("data/benzene_arom.mol", 12, 12)]
    #[case("data/benzene.mol", 6, 6)]
    #[case("data/tep.mol", 46, 50)]
    #[case("data/corrole.mol", 37, 41)]
    fn test_mol_files(#[case] filename: &str, #[case] atom_len: usize, #[case] bond_len: usize) {
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let (atoms, bonds) = parse(reader).unwrap();

        assert_eq!(atoms.len(), atom_len);
        assert_eq!(bonds.len(), bond_len);
    }
}

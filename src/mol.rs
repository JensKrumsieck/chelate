//! Functions for parsing MDL MOL files (chemical/x-mdl-molfile)
//! Documentation can be found here: <https://en.wikipedia.org/wiki/Chemical_table_file#Molfile>
use super::normalize_symbol;
use crate::atom::{ATOMIC_SYMBOLS, Atom, Bond};
use std::io::{self, BufRead, BufReader, Read};

/// Parses a single line of an MOL file and returns an `Atom` object.
/// The line should contain the x, y, and z coordinates followed by the atomic symbol.
/// Example line: `    1.3194   -1.2220   -0.8506 N   0  0  0  0  0  0  0  0  0  0  0  0`
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

    let mut atom = Atom::new(*atom_count, atomic_number as u8, x, y, z);
    atom.name = symbol.to_string();

    Some(atom)
}

/// Parses the counts line of a MOL file and returns the number of atoms and bonds.
/// The fields are 3 characters wide and run together for counts above 99.
/// Example line: `101100  0  0  0  0  0  0  0  0999 V2000` (101 atoms, 100 bonds)
fn parse_counts_line(line: &str) -> Option<(usize, usize)> {
    let atoms = line.get(0..3)?.trim().parse().ok()?;
    let bonds = line.get(3..6)?.trim().parse().ok()?;
    Some((atoms, bonds))
}

/// Parses a single line of an MOL file and returns a `Bond` object.
/// The line should contain the atoms ids and the bond order where 4 is aromatic bond.
/// The fields are 3 characters wide and run together for atom ids above 99.
/// Example lines: `  1  2  2  0  0  0  0`, `100101  1  0  0  0  0`
fn parse_bond_line(line: &str) -> Option<Bond> {
    let atom1 = line.get(0..3)?.trim().parse().ok()?;
    let atom2 = line.get(3..6)?.trim().parse().ok()?;
    let mut order = line.get(6..9)?.trim().parse().ok()?;
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
/// use nalgebra::Point3;
/// use std::io::BufReader;
///
/// let file = File::open("data/corrole.mol").unwrap();
/// let reader = BufReader::new(file);
/// let (atoms, bonds) = mol::parse(reader).unwrap();
///
/// assert_eq!(atoms.len(), 37);
/// assert_eq!(bonds.len(), 41);
/// assert_eq!(atoms[0].atomic_number, 7);
/// assert_eq!(atoms[0].coord, Point3::new(1.3194, -1.2220, -0.8506));
/// assert_eq!(atoms[0].resname, "UNK");
/// assert_eq!(atoms[0].resid, 0);
/// assert_eq!(atoms[0].chain, char::default());
/// assert_eq!(atoms[0].occupancy, 1.0);
/// assert_eq!(atoms[0].name, "N");
/// ```
pub fn parse<P: Read>(reader: BufReader<P>) -> io::Result<(Vec<Atom>, Vec<Bond>)> {
    let mut atom_count = 0;

    // the counts line follows the 3 header lines and tells how many atom and bond lines follow
    let mut lines = reader.lines().skip(3);
    let counts_line = lines.next().transpose()?.unwrap_or_default();
    let (atom_len, bond_len) = parse_counts_line(&counts_line)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid counts line"))?;

    let mut atoms = Vec::with_capacity(atom_len);
    for line in lines.by_ref().take(atom_len) {
        if let Some(atom) = parse_atom_line(&line?, &mut atom_count) {
            atoms.push(atom);
        }
    }

    let mut bonds = Vec::with_capacity(bond_len);
    for line in lines.take(bond_len) {
        if let Some(bond) = parse_bond_line(&line?) {
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

    #[test]
    fn test_mol_more_than_99_atoms() {
        // carbon chain with 101 atoms, the last bond line reads `100101  1  0  0  0  0`
        let mut mol = String::from("chain\n\n\n101100  0  0  0  0  0  0  0  0999 V2000\n");
        for i in 0..101 {
            mol += &format!(
                "{:10.4}{:10.4}{:10.4} C   0  0  0  0  0  0  0  0  0  0  0  0\n",
                i as f32 * 1.5,
                0.0,
                0.0
            );
        }
        for i in 1..101 {
            mol += &format!("{:3}{:3}  1  0  0  0  0\n", i, i + 1);
        }
        mol += "M  END\n";

        let (atoms, bonds) = parse(BufReader::new(mol.as_bytes())).unwrap();

        assert_eq!(atoms.len(), 101);
        assert_eq!(bonds.len(), 100);
        let last = bonds.last().unwrap();
        assert_eq!((last.atom1, last.atom2, last.order), (100, 101, 1));
    }
}

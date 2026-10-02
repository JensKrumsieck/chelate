//! Functions for parsing MDL MOL files (chemical/x-mdl-molfile)
//! Documentation can be found here: <https://en.wikipedia.org/wiki/Chemical_table_file#Molfile>
use super::{atomic_number, column};
use crate::atom::{Atom, Bond};
use crate::error::{FileError, ParseError};
use std::io::{BufRead, BufReader, Read};

/// Parses a single line of an MOL file and returns an `Atom` object,
/// `None` for symbols that are no elements like `R#` or `*`.
/// The line should contain the x, y, and z coordinates (10 characters each) followed by the atomic symbol (columns 32 - 34).
/// Example line: `    1.3194   -1.2220   -0.8506 N   0  0  0  0  0  0  0  0  0  0  0  0`
fn parse_atom_line(line: &str, atom_count: &mut usize) -> Result<Option<Atom>, ParseError> {
    let coord = |start| {
        column(line, start, start + 10)
            .parse::<f32>()
            .map_err(|_| ParseError::new("Invalid coordinates"))
    };
    let (x, y, z) = (coord(0)?, coord(10)?, coord(20)?);

    let symbol = column(line, 31, 34);
    let Some(atomic_number) = atomic_number(symbol) else {
        return Ok(None);
    };
    *atom_count += 1;

    let mut atom = Atom::new(*atom_count, atomic_number, x, y, z);
    atom.name = symbol.into();

    Ok(Some(atom))
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
fn parse_bond_line(line: &str) -> Result<Bond, ParseError> {
    let invalid = |_| ParseError::new("Invalid bond");
    let atom1 = column(line, 0, 3).parse().map_err(invalid)?;
    let atom2 = column(line, 3, 6).parse().map_err(invalid)?;
    let mut order = column(line, 6, 9).parse().map_err(invalid)?;
    let is_aromatic = order == 4;
    if is_aromatic {
        order = 1;
    }

    Ok(Bond {
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
pub fn parse<P: Read>(reader: BufReader<P>) -> Result<(Vec<Atom>, Vec<Bond>), FileError> {
    let mut lines = reader.lines().zip(1..).skip(3);
    // the file must not end before all atom and bond lines are read
    let mut next_line = || -> Result<(String, usize), FileError> {
        match lines.next() {
            Some((line, line_number)) => Ok((line?, line_number)),
            None => Err(ParseError::new("Unexpected end of file").into()),
        }
    };

    // the counts line follows the 3 header lines and tells how many atom and bond lines follow
    let (counts_line, line_number) = next_line()?;
    let (atom_len, bond_len) = parse_counts_line(&counts_line)
        .ok_or_else(|| ParseError::new("Invalid counts line").with_line(line_number))?;

    let mut atom_count = 0;
    let mut atoms = Vec::with_capacity(atom_len);
    // ids of the parsed atoms by their position in the file, `None` for skipped atoms
    let mut ids = Vec::with_capacity(atom_len);
    for _ in 0..atom_len {
        let (line, line_number) = next_line()?;
        let atom = parse_atom_line(&line, &mut atom_count).map_err(|e| e.with_line(line_number))?;
        ids.push(atom.as_ref().map(|a| a.id));
        atoms.extend(atom);
    }

    let mut bonds = Vec::with_capacity(bond_len);
    for _ in 0..bond_len {
        let (line, line_number) = next_line()?;
        let bond = parse_bond_line(&line).map_err(|e| e.with_line(line_number))?;
        let id = |atom: usize| {
            atom.checked_sub(1)
                .and_then(|index| ids.get(index).copied())
                .ok_or_else(|| ParseError::new("Bond refers to missing atom").with_line(line_number))
        };
        // bonds to skipped atoms are skipped as well
        if let (Some(atom1), Some(atom2)) = (id(bond.atom1)?, id(bond.atom2)?) {
            bonds.push(Bond {
                atom1,
                atom2,
                ..bond
            });
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

    #[test]
    fn test_mol_invalid_counts_line() {
        let mol = "name\n\n\nnot a counts line\n";

        let error = parse(BufReader::new(mol.as_bytes())).unwrap_err();

        assert!(matches!(
            error,
            FileError::Parse(ParseError {
                line_number: Some(4),
                ..
            })
        ));
    }
}

//! Functions for parsing PDB files (chemical/x-pdb)
//! The PDB file format is documented here <https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html>
use super::normalize_symbol;
use crate::atom::{ATOMIC_SYMBOLS, Atom};
use std::io::{self, BufRead, BufReader, Read};

/// Parses a single line of a PDB file and returns an `Atom` object.
/// Accoding to the PDB format specification, the line should contain the atomic symbol followed by the x, y, and z coordinates.
/// <http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM>
///
/// | Columns  | Data Type   | Field        | Definition                                  |
/// |----------|-------------|--------------|---------------------------------------------|
/// | 1 - 6    | Record name | "ATOM  "     | (or HETATM)                                 |
/// | 7 - 11   | Integer     | serial       | Atom serial number.                         |
/// | 13 - 16  | Atom        | name         | Atom name.                                  |
/// | 17       | Character   | altLoc       | Alternate location indicator.               |
/// | 18 - 20  | Residue name| resName      | Residue name.                               |
/// | 22       | Character   | chainID      | Chain identifier.                           |
/// | 23 - 26  | Integer     | resSeq       | Residue sequence number.                    |
/// | 27       | AChar       | iCode        | Code for insertion of residues.             |
/// | 31 - 38  | Real(8.3)   | x            | Orthogonal coordinates for X in Angstroms.  |
/// | 39 - 46  | Real(8.3)   | y            | Orthogonal coordinates for Y in Angstroms.  |
/// | 47 - 54  | Real(8.3)   | z            | Orthogonal coordinates for Z in Angstroms.  |
/// | 55 - 60  | Real(6.2)   | occupancy    | Occupancy.                                  |
/// | 61 - 66  | Real(6.2)   | tempFactor   | Temperature factor.                         |
/// | 77 - 78  | LString(2)  | element      | Element symbol, right-justified.            |
/// | 79 - 80  | LString(2)  | charge       | Charge on the atom.                         |
fn parse_atom_line(line: &str) -> Option<Atom> {
    let symbol = line[76..78].trim();
    let atomic_number = ATOMIC_SYMBOLS
        .iter()
        .position(|&s| s == normalize_symbol(symbol))?
        + 1;

    let id = line[8..11].trim().parse().ok()?;

    let x = line[30..38].trim().parse().ok()?;
    let y = line[38..46].trim().parse().ok()?;
    let z = line[46..54].trim().parse().ok()?;

    Some(Atom::new(id, atomic_number as u8, x, y, z))
}

/// Parses an XYZ file and returns a vector of `Atom` objects.
/// # Examples
/// ```
/// use chelate::pdb;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let file = File::open("data/0001.pdb").unwrap();
/// let reader = BufReader::new(file);
/// let atoms = pdb::parse(reader).unwrap();
///
/// assert_eq!(atoms.len(), 15450);
/// ```
pub fn parse<P: Read>(reader: BufReader<P>) -> io::Result<Vec<Atom>> {
    let mut atoms = Vec::new();
    for line in reader.lines() {
        let line = line?;
        // Skip lines that are not ATOM or HETATM records
        if !line.starts_with("ATOM") && !line.starts_with("HETATM") {
            continue;
        }
        if let Some(atom) = parse_atom_line(&line) {
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
    #[case("data/oriluy.pdb", 130)]
    #[case("data/2spl.pdb", 1437)]
    #[case("data/1hv4.pdb", 9288)]
    #[case("data/0001.pdb", 15450)]
    fn test_pdb_files(#[case] filename: &str, #[case] len: usize) {
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let atoms = parse(reader).unwrap();

        assert_eq!(atoms.len(), len);
    }
}

//! Functions for parsing PDB files (chemical/x-pdb)
//! The PDB file format is documented here <https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html>
use super::{atomic_number, column};
use crate::types::Atom;
use crate::error::{FileError, ParseError};
use std::io::{BufRead, BufReader, Read};

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
fn parse_atom_line(line: &str, atom_count: &mut usize) -> Result<Option<Atom>, ParseError> {
    let coord = |start, end| {
        column(line, start, end)
            .parse::<f32>()
            .map_err(|_| ParseError::new("Invalid coordinates"))
    };
    let (x, y, z) = (coord(30, 38)?, coord(38, 46)?, coord(46, 54)?);

    // older files often end after the coordinates or the temperature factor
    let symbol = match column(line, 76, 78) {
        "" => line
            .get(12..16)
            .and_then(element_from_atom_name)
            .unwrap_or_default(),
        symbol => symbol,
    };
    let Some(atomic_number) = atomic_number(symbol) else {
        return Ok(None);
    };

    let chain = column(line, 21, 22).parse().unwrap_or_default();
    let resname = column(line, 17, 20).to_string();
    let resid = column(line, 22, 26).parse().unwrap_or_default();
    let occ = column(line, 54, 60).parse().unwrap_or(1.0);

    *atom_count += 1;

    let mut atom = Atom::new(*atom_count, atomic_number, x, y, z);
    atom.data.chain = chain;
    atom.data.resname = resname.into();
    atom.data.resid = resid;
    atom.data.occupancy = occ;
    atom.data.name = symbol.into();

    Ok(Some(atom))
}

/// Guesses the element symbol from the atom name (columns 13 - 16) for lines without element columns.
/// The element symbol is right-justified in columns 13 - 14, so ` CA ` is a carbon and `CA  ` is calcium.
/// Hydrogen names may start with a digit (`1HB `) or fill all four columns (`HG21`).
fn element_from_atom_name(name: &str) -> Option<&str> {
    if name.starts_with(|c: char| c == ' ' || c.is_ascii_digit()) {
        name.get(1..2)
    } else if name.get(2..4)?.chars().all(|c| c.is_ascii_digit()) {
        name.get(0..1)
    } else {
        // non-standard left-justified names like `C1  `
        Some(
            name.get(0..2)?
                .trim_end_matches(|c: char| !c.is_ascii_alphabetic()),
        )
    }
}

/// Parses an PDB file and returns a vector of `Atom` objects.
/// # Examples
/// ```
/// use chelate::pdb;
/// use std::fs::File;
/// use nalgebra::Point3;
/// use std::io::BufReader;
///
/// let file = File::open("data/0001.pdb").unwrap();
/// let reader = BufReader::new(file);
/// let atoms = pdb::parse(reader).unwrap();
///
/// assert_eq!(atoms.len(), 15450);
/// assert_eq!(atoms[0].atomic_number, 7);
/// assert_eq!(atoms[0].coord, Point3::new(58.667, 69.671, 7.056));
/// assert_eq!(atoms[0].resname, "ASP");
/// assert_eq!(atoms[0].resid, 25);
/// assert_eq!(atoms[0].occupancy, 1.0);
/// assert_eq!(atoms[0].name, "N");
/// ```
pub fn parse<P: Read>(reader: BufReader<P>) -> Result<Vec<Atom>, FileError> {
    let mut atom_count = 0;

    let mut atoms = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        // Skip lines that are not ATOM or HETATM records
        if !line.starts_with("ATOM") && !line.starts_with("HETATM") {
            continue;
        }
        let atom = parse_atom_line(&line, &mut atom_count).map_err(|e| e.with_line(i + 1))?;
        atoms.extend(atom);
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

    #[test]
    fn test_pdb_without_element_columns() {
        let pdb = "\
ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00
ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00  0.00
HETATM    3 ZN    ZN A 101      10.000   5.000  -5.000
";
        let atoms = parse(BufReader::new(pdb.as_bytes())).unwrap();

        assert_eq!(atoms.len(), 3);
        assert_eq!(atoms[0].atomic_number, 7);
        assert_eq!(atoms[1].atomic_number, 6);
        assert_eq!(atoms[2].atomic_number, 30);
        assert_eq!(atoms[2].data.resname, "ZN");
        assert_eq!(atoms[2].data.occupancy, 1.0);
    }

    #[rstest]
    #[case(" CA ", Some("C"))]
    #[case("CA  ", Some("CA"))]
    #[case("1HB ", Some("H"))]
    #[case("HG21", Some("H"))]
    #[case("FE  ", Some("FE"))]
    #[case("C1  ", Some("C"))]
    fn test_element_from_atom_name(#[case] name: &str, #[case] element: Option<&str>) {
        assert_eq!(element_from_atom_name(name), element);
    }
}

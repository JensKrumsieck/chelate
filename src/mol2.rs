//! Functions for parsing TRIPOS MOL2 files (TRIPOS now is Certara) (chemical/x-mol2)
//! Native format of the SYBYL cheminformatics application.
//! The MOL2 Format is documented here <https://paulbourke.net/dataformats/mol2/>
use super::normalize_symbol;
use crate::atom::{ATOMIC_SYMBOLS, Atom, Bond};
use std::io::{self, BufRead, BufReader, Read};

/// Parses a single line of an TRIPOS MOL2 file and returns an `Atom` object.
/// The line should contain the x, y, and z coordinates followed by the atomic symbol.
/// Example line: `     1 N       58.6644  69.6736   7.0558   N.3       1 ASP25  32.7500`
fn parse_atom_line(line: &str) -> Option<Atom> {
    let mut iter = line.split_whitespace();

    let id = iter.next()?.parse().ok()?;
    let name = iter.next()?;
    let x = iter.next()?.parse().ok()?;
    let y = iter.next()?.parse().ok()?;
    let z = iter.next()?.parse().ok()?;

    let type_ = iter.next()?;
    let mut symbol = type_.split('.').next()?;

    let chain_id = if let Some(next) = iter.next() {
        next.parse().unwrap_or_default()
    } else {
        Default::default()
    };

    let residue_raw = iter.next().unwrap_or_default();
    let id_pos = residue_raw.chars().position(|c| c.is_numeric());
    let residue = if let Some(pos) = id_pos {
        &residue_raw[..pos]
    } else {
        residue_raw
    };

    let res_id = if let Some(pos) = id_pos {
        residue_raw[pos..].parse().unwrap_or(0)
    } else {
        0
    };

    if !ATOMIC_SYMBOLS.contains(&symbol) {
        symbol = get_symbol_from_name(name)
    }

    let atomic_number = ATOMIC_SYMBOLS
        .iter()
        .position(|&s| s == normalize_symbol(symbol))?
        + 1;

    let mut atom = Atom::new(id, atomic_number as u8, x, y, z);
    atom.chain = chain_id;
    atom.resname = residue.to_string();
    atom.resid = res_id;
    atom.name = name.to_string();

    Some(atom)
}

fn get_symbol_from_name(s: &str) -> &str {
    s.char_indices()
        .find(|(_, c)| c.is_ascii_digit())
        .map(|(i, _)| &s[..i])
        .unwrap_or(s)
}

/// Parses a single line of an MO2L file and returns a `Bond` object.
/// The line should contain the bond id, the atoms ids by the bond order where "ar" is aromatic bond.
/// Example line: `     1     1     2   un`
fn parse_bond_line(line: &str) -> Option<Bond> {
    let mut iter = line.split_whitespace();

    iter.next()?; // Skip the bond id
    let atom1 = iter.next()?.parse().ok()?;
    let atom2 = iter.next()?.parse().ok()?;
    let raw_order = iter.next()?;
    let order = raw_order.parse().unwrap_or(1);
    let is_aromatic = raw_order.starts_with("ar");

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
/// use chelate::mol2;
/// use nalgebra::Vector3;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let file = File::open("data/ptcor.mol2").unwrap();
/// let reader = BufReader::new(file);
/// let (atoms, bonds) = mol2::parse(reader).unwrap();
///
/// assert_eq!(atoms.len(), 129);
/// assert_eq!(bonds.len(), 127);
/// 
/// assert_eq!(atoms[0].atomic_number, 78);
/// assert_eq!(atoms[0].coord, Vector3::new(8.7088, 6.0412,5.1685));
/// assert_eq!(atoms[0].resname, "RES");
/// assert_eq!(atoms[0].resid, 1);
/// assert_eq!(atoms[0].chain, '1');
/// assert_eq!(atoms[0].occupancy, 1.0);
/// assert_eq!(atoms[0].name, "Pt1");
/// ```
pub fn parse<P: Read>(reader: BufReader<P>) -> io::Result<(Vec<Atom>, Vec<Bond>)> {
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();

    let mut pick_atoms = false;
    let mut pick_bonds = false;

    for line in reader.lines().skip_while(|s| {
        s.as_ref()
            .unwrap_or(&String::from(""))
            .contains("@<TRIPOS>ATOM")
    }) {
        let line = line?;
        let line_trimmed = line.trim();
        match line_trimmed {
            "@<TRIPOS>ATOM" => {
                pick_atoms = true;
                pick_bonds = false;
                continue;
            }
            "@<TRIPOS>BOND" => {
                pick_atoms = false;
                pick_bonds = true;
                continue;
            }
            _ => {
                if line.starts_with("@<TRIPOS>") {
                    pick_atoms = false;
                    pick_bonds = false;
                    continue;
                }
            }
        }

        if pick_atoms {
            if let Some(atom) = parse_atom_line(&line) {
                atoms.push(atom);
            }
        } else if pick_bonds {
            if let Some(bond) = parse_bond_line(&line) {
                bonds.push(bond);
            }
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
    #[case("data/0001.mol2", 15450, 14898)]
    #[case("data/benzene.mol2", 12, 12)]
    #[case("data/myo.mol2", 1437, 1312)]
    #[case("data/ptcor.mol2", 129, 127)]
    #[case("data/tep.mol2", 46, 50)]
    #[case("data/VATTOC.mol2", 130, 146)]
    fn test_mol2_files(#[case] filename: &str, #[case] atom_len: usize, #[case] bond_len: usize) {
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let (atoms, bonds) = parse(reader).unwrap();

        assert_eq!(atoms.len(), atom_len);
        assert_eq!(bonds.len(), bond_len);
    }
}

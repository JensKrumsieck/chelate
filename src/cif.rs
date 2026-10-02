//! Functions for parsing CIF files
//! Multiple types of CIF files are supported:
//! - CCDC/IUCr CIF: <https://www.iucr.org/resources/cif> (chemical/x-cif)
//! - PDBx/mmCIF: <https://mmcif.wwpdb.org/docs/user-guide/guide.html> (chemical/x-mmcif)
//!
//! See also: <https://en.wikipedia.org/wiki/Crystallographic_Information_File>
use super::normalize_symbol;
use crate::error::{FileError, ParseError};
use crate::types::{ATOMIC_SYMBOLS, Atom, Bond};
use nalgebra::Matrix4;
use std::{
    collections::HashMap,
    f32::consts::PI,
    io::{BufRead, BufReader, Read},
};

/// Represents the type of CIF file (dialect) being parsed.
#[derive(Default, PartialEq)]
enum CIFDialect {
    #[default]
    Undefined,
    #[allow(clippy::upper_case_acronyms)]
    CCDC,
    #[allow(non_camel_case_types)]
    mmCIF,
    #[allow(non_camel_case_types)]
    compCIF,
}

/// Stores position of columns based on column headers, `None` if the column is missing
#[derive(Default)]
struct CIFAtomHeader {
    symbol: Option<usize>,
    x: Option<usize>,
    y: Option<usize>,
    z: Option<usize>,
    id: Option<usize>,
    disorder: Option<usize>,
    residue: Option<usize>,
    chain: Option<usize>,
    seq_id: Option<usize>,
    occupancy: Option<usize>,
    /// coordinates are fractional and need the cell parameters to be converted
    fractional: bool,
}

/// Parses a single line of a CCDC CIF (Crystallographic information file) file and returns an `Atom` object.
/// In comparison to other file formats, cif files can have 3 different types: CCDC, mmCIF and compCIF
///
/// CCDC:       `N1 N 0.0662(3) 0.55056(14) 0.1420(3) 0.066(2) Uani 1 1 d . . . . .`
///
/// mmCIF:      `ATOM   2    C  CA  . MET A 1 13  ? -16.763 -22.990 22.365  1.00 30.45  ? 1   MET A CA`
///
/// compCIF:    `H20 CA   CA   C  0 1 N N N 1.747  -36.297 22.990 -1.853 -0.230 0.046  CA   H20 1`
fn parse_atom_line(
    line: &str,
    header: &CIFAtomHeader,
    atom_count: &mut usize,
    label_map: &mut HashMap<String, usize>,
) -> Option<Atom> {
    let vec = line.split_whitespace().collect::<Vec<_>>();
    let column = |index: Option<usize>| vec.get(index?).copied();

    let id = column(header.id)?;
    //the type symbol is optional as labels start with the element symbol
    let symbol = column(header.symbol).unwrap_or_else(|| element_from_label(id));
    let x = column(header.x)?.split('(').next()?.parse().ok()?;
    let y = column(header.y)?.split('(').next()?.parse().ok()?;
    let z = column(header.z)?.split('(').next()?.parse().ok()?;
    let disorder_group = column(header.disorder)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0);
    let residue = column(header.residue).unwrap_or("UNK");
    let chain_id = column(header.chain).unwrap_or_default();
    let seq_id = column(header.seq_id)
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or_default();
    let occ = column(header.occupancy)
        .and_then(|s| s.parse::<f32>().ok())
        .unwrap_or(1.0);

    let atomic_number = ATOMIC_SYMBOLS
        .iter()
        .position(|&s| s == normalize_symbol(symbol))?
        + 1;
    *atom_count += 1;
    label_map.insert(id.to_owned(), *atom_count);
    let mut atom = Atom::new(*atom_count, atomic_number as u8, x, y, z);
    //additional info
    atom.data.disorder_group = disorder_group;
    atom.data.name = id.into();
    atom.data.resname = residue.into();
    atom.data.chain = chain_id.into();
    atom.data.resid = seq_id;
    atom.data.occupancy = occ;
    Some(atom)
}

/// Derives the element symbol from an atom label for files without `_atom_site_type_symbol`.
/// Labels start with the element symbol in its usual case, so `Cu1` is copper while `C1A` and `CA1` are carbon.
fn element_from_label(label: &str) -> &str {
    match label.as_bytes() {
        [_, second, ..] if second.is_ascii_lowercase() && ATOMIC_SYMBOLS.contains(&&label[..2]) => {
            &label[..2]
        }
        _ => label.get(..1).unwrap_or_default(),
    }
}

/// Parses a single line of an CIF file and returns a `Bond` object.
/// The line should contain the atoms names which need to be mapped to ids
/// Example line: `C4A N21A 1.370(3) . ?`
fn parse_bond_line(line: &str, map: &HashMap<String, usize>, dialect: &CIFDialect) -> Option<Bond> {
    let mut iter = line.split_whitespace();
    if *dialect == CIFDialect::compCIF {
        iter.next()?;
    }
    let atom1 = iter.next()?;
    let atom2 = iter.next()?;
    Some(Bond {
        atom1: *map.get(atom1)?,
        atom2: *map.get(atom2)?,
        order: 1,
        is_aromatic: false,
    })
}

/// Parses a CIF file and returns a vector of `Atom` and a vector of `Bond` objects.
/// # Examples
/// ```
/// use chelate::cif;
/// use std::fs::File;
/// use nalgebra::Point3;
/// use approx::relative_eq;
/// use std::io::BufReader;
///
/// let file = File::open("data/147288.cif").unwrap();
/// let reader = BufReader::new(file);
/// let (atoms, bonds) = cif::parse(reader).unwrap();
///
/// assert_eq!(atoms.len(), 206);
/// assert_eq!(bonds.len(), 230);
/// assert_eq!(atoms[0].atomic_number, 31);
/// assert!(relative_eq!(atoms[0].coord, Point3::new(11.377683611607571, 1.637743396392762, 3.827447754962335), epsilon = 1.0e-5));
/// assert_eq!(atoms[0].resname, "UNK");
/// assert_eq!(atoms[0].resid, 0);
/// assert_eq!(atoms[0].chain, char::default());
/// assert_eq!(atoms[0].occupancy, 1.0);
/// assert_eq!(atoms[0].name, "Ga1A");
/// ```
pub fn parse<P: Read>(reader: BufReader<P>) -> Result<(Vec<Atom>, Vec<Bond>), FileError> {
    let mut dialect = CIFDialect::default();

    let mut atoms = Vec::new();
    let mut bonds = Vec::new();

    let mut pick_atoms = false;
    let mut pick_bonds = false;

    let mut headers = CIFAtomHeader::default();
    let mut header_idx = 0;
    let mut atom_count = 0;

    //index of the first atom of the current data block
    let mut block_start = 0;
    let mut cell_params: [Option<f32>; 6] = Default::default();

    let mut map: HashMap<String, usize> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let line_trimmed = line.trim();
        if line_trimmed.is_empty() {
            continue;
        }

        if line_trimmed.starts_with("data_") {
            if headers.fractional {
                fractional_to_cartesian(&mut atoms[block_start..], &cell_params)?;
            }
            //each data block is a separate structure with its own cell, columns and labels
            dialect = CIFDialect::default();
            pick_atoms = false;
            pick_bonds = false;
            headers = CIFAtomHeader::default();
            header_idx = 0;
            block_start = atoms.len();
            cell_params = Default::default();
            map.clear();
            continue;
        }

        if dialect == CIFDialect::Undefined {
            dialect = set_dialect(line_trimmed);
        }
        if line_trimmed.starts_with("_cell_") {
            parse_cell_param(line_trimmed, &mut cell_params);
        }

        if line_trimmed.starts_with("loop_") {
            pick_atoms = false;
            pick_bonds = false;
        }
        if line_trimmed.starts_with("_atom_site_label")
            || line_trimmed.starts_with("_atom_site.")
            || line_trimmed.starts_with("_chem_comp_atom.")
        {
            pick_atoms = true;
            pick_bonds = false;
        }
        if line_trimmed.starts_with("_geom_bond") || line_trimmed.starts_with("_chem_comp_bond") {
            pick_atoms = false;
            pick_bonds = true;
        }

        if pick_atoms {
            if line_trimmed.starts_with("_") {
                set_header_indices(line_trimmed, header_idx, &mut headers);
                header_idx += 1;
            } else if let Some(atom) =
                parse_atom_line(line_trimmed, &headers, &mut atom_count, &mut map)
            {
                atoms.push(atom);
            }
        } else if pick_bonds && let Some(bond) = parse_bond_line(line_trimmed, &map, &dialect) {
            bonds.push(bond);
        }
    }

    if headers.fractional {
        fractional_to_cartesian(&mut atoms[block_start..], &cell_params)?;
    }
    Ok((atoms, bonds))
}

fn set_dialect(line: &str) -> CIFDialect {
    if line.starts_with("_chem_comp") {
        CIFDialect::compCIF
    } else if line.starts_with("_atom_type_symbol") || line.starts_with("_symmetry") {
        CIFDialect::CCDC
    } else if line.starts_with("_pdbx") {
        CIFDialect::mmCIF
    } else {
        CIFDialect::Undefined
    }
}

fn set_header_indices(header: &str, index: usize, headers: &mut CIFAtomHeader) {
    match header {
        h if h.contains("symbol") => headers.symbol = Some(index),
        h if h.contains("fract_x") || h.contains("Cartn_x") => {
            headers.x = Some(index);
            headers.fractional = h.contains("fract_x");
        }
        h if h.contains("fract_y") || h.contains("Cartn_y") => headers.y = Some(index),
        h if h.contains("fract_z") || h.contains("Cartn_z") => headers.z = Some(index),
        h if h.contains("label_atom_id")
            || h.contains("atom.atom_id")
            || h.contains("_site_label") =>
        {
            headers.id = Some(index)
        }
        h if h.contains("disorder_group") => headers.disorder = Some(index),
        h if h.contains("auth_comp_id") || h.contains("comp_id") => headers.residue = Some(index),
        h if h.contains("auth_asym_id") => headers.chain = Some(index),
        h if h.contains("auth_seq_id") => headers.seq_id = Some(index),
        h if h.contains("occupancy") => headers.occupancy = Some(index),
        _ => {}
    }
}

/// Stores the value of a cell parameter line like `_cell_length_a 8.1707(5)`.
/// The order is a, b, c, alpha, beta, gamma, independent of the order in the file.
fn parse_cell_param(line: &str, cell_params: &mut [Option<f32>; 6]) {
    let mut iter = line.split_whitespace();
    let index = match iter.next() {
        Some("_cell_length_a") => 0,
        Some("_cell_length_b") => 1,
        Some("_cell_length_c") => 2,
        Some("_cell_angle_alpha") => 3,
        Some("_cell_angle_beta") => 4,
        Some("_cell_angle_gamma") => 5,
        _ => return,
    };
    cell_params[index] = iter.next().and_then(get_value_from_uncertainity);
}

fn get_value_from_uncertainity(input: &str) -> Option<f32> {
    input.split('(').next()?.parse().ok()
}

/// Converts fractional coordinates to cartesian coordinates.
/// Runs once a data block is read completely, as the cell parameters may follow the atoms.
fn fractional_to_cartesian(
    atoms: &mut [Atom],
    cell_params: &[Option<f32>; 6],
) -> Result<(), ParseError> {
    if atoms.is_empty() {
        return Ok(());
    }
    let [
        Some(a),
        Some(b),
        Some(c),
        Some(alpha),
        Some(beta),
        Some(gamma),
    ] = *cell_params
    else {
        return Err(ParseError::new(
            "Missing cell parameters for fractional coordinates",
        ));
    };
    let matrix = conversion_matrix(a, b, c, alpha, beta, gamma);
    for atom in atoms {
        atom.coord = matrix.transform_vector(&atom.coord.coords).into();
    }
    Ok(())
}

fn conversion_matrix(a: f32, b: f32, c: f32, alpha: f32, beta: f32, gamma: f32) -> Matrix4<f32> {
    let cos_alpha = (alpha * PI / 180.0).cos();
    let cos_beta = (beta * PI / 180.0).cos();
    let cos_gamma = (gamma * PI / 180.0).cos();
    let sin_gamma = (gamma * PI / 180.0).sin();

    Matrix4::new(
        a,
        b * cos_gamma,
        c * cos_beta,
        0.0,
        0.0,
        b * sin_gamma,
        c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma,
        0.0,
        0.0,
        0.0,
        c * ((1.0 - cos_alpha.powi(2) - cos_beta.powi(2) - cos_gamma.powi(2)
            + 2.0 * cos_alpha * cos_beta * cos_gamma)
            .sqrt())
            / sin_gamma,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::relative_eq;
    use nalgebra::Point3;
    use rstest::rstest;
    use std::fs::File;

    #[rstest]
    #[case("data/4n4n.cif", 15450, 0)] //mmcif = no bonds
    #[case("data/4r21.cif", 6752, 0)] //mmcif = no bonds
    #[case("data/147288.cif", 206, 230)]
    #[case("data/1484829.cif", 466, 528)]
    #[case("data/cif_noTrim.cif", 79, 89)]
    #[case("data/cif.cif", 79, 89)]
    #[case("data/CuHETMP.cif", 85, 0)] //no bonds in file
    #[case("data/ligand.cif", 44, 46)]
    #[case("data/mmcif.cif", 1291, 0)] //mmcif = no bonds
    fn test_cif_files(#[case] filename: &str, #[case] atom_len: usize, #[case] bond_len: usize) {
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let (atoms, bonds) = parse(reader).unwrap();
        assert_eq!(
            atoms.iter().filter(|a| a.data.disorder_group != 2).count(),
            atom_len
        );
        assert_eq!(
            bonds
                .iter()
                .filter(|b| atoms[b.atom1 - 1].data.disorder_group != 2
                    && atoms[b.atom2 - 1].data.disorder_group != 2)
                .count(),
            bond_len
        );
    }

    #[test]
    fn test_cif_multiple_data_blocks() {
        let parse_file =
            |filename: &str| parse(BufReader::new(File::open(filename).unwrap())).unwrap();
        let (atoms1, bonds1) = parse_file("data/147288.cif");
        let (atoms2, bonds2) = parse_file("data/cif.cif");

        let file = File::open("data/147288.cif")
            .unwrap()
            .chain(File::open("data/cif.cif").unwrap());
        let (atoms, bonds) = parse(BufReader::new(file)).unwrap();

        assert_eq!(atoms.len(), atoms1.len() + atoms2.len());
        assert_eq!(bonds.len(), bonds1.len() + bonds2.len());
        //second block uses its own cell parameters
        for (atom, expected) in atoms[atoms1.len()..].iter().zip(&atoms2) {
            assert_eq!(atom.data.name, expected.data.name);
            assert_eq!(atom.coord, expected.coord);
        }
        //bonds of the second block refer to atoms of the second block
        for (bond, expected) in bonds[bonds1.len()..].iter().zip(&bonds2) {
            assert_eq!(bond.atom1, expected.atom1 + atoms1.len());
            assert_eq!(bond.atom2, expected.atom2 + atoms1.len());
        }
    }

    #[test]
    fn test_cif_without_type_symbol_and_markers() {
        //no _atom_site_type_symbol, no _symmetry_* or _atom_type_symbol and the cell follows the atoms
        let cif = "\
data_test
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Cu1 0.1 0.2 0.3
Cl1 0.2 0.2 0.3
C1A 0.3 0.2 0.3
H1A 0.4 0.2 0.3
_cell_length_a 10
_cell_length_b 20
_cell_length_c 30
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
";
        let (atoms, _) = parse(BufReader::new(cif.as_bytes())).unwrap();

        let atomic_numbers: Vec<_> = atoms.iter().map(|a| a.symbol.atomic_number()).collect();
        assert_eq!(atomic_numbers, [29, 17, 6, 1]);
        assert!(relative_eq!(
            atoms[0].coord,
            Point3::new(1.0, 4.0, 9.0),
            epsilon = 1.0e-5
        ));
    }

    #[test]
    fn test_cif_without_dialect_markers() {
        //cif.cif uses _space_group_* instead of _symmetry_*, so only _atom_type_symbol marks it as CCDC
        let reference = parse(BufReader::new(File::open("data/cif.cif").unwrap()))
            .unwrap()
            .0;
        let cif = std::fs::read_to_string("data/cif.cif")
            .unwrap()
            .replace("_atom_type_symbol", "_atom_type_renamed");
        let (atoms, _) = parse(BufReader::new(cif.as_bytes())).unwrap();

        assert_eq!(atoms, reference);
    }

    #[test]
    fn test_cif_fractional_without_cell() {
        let cif = "data_test\nloop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nC1 0.1 0.2 0.3\n";

        let error = parse(BufReader::new(cif.as_bytes())).unwrap_err();

        assert!(matches!(error, FileError::Parse(_)));
    }

    #[rstest]
    #[case("Cu1", "Cu")]
    #[case("Cl2", "Cl")]
    #[case("C1A", "C")]
    #[case("CA1", "C")]
    #[case("H1A", "H")]
    #[case("Ow1", "O")]
    fn test_element_from_label(#[case] label: &str, #[case] element: &str) {
        assert_eq!(element_from_label(label), element);
    }

    #[rstest]
    #[case("data/147288.cif")]
    #[case("data/1484829.cif")]
    #[case("data/cif.cif")]
    #[case("data/CuHETMP.cif")]
    fn test_element_from_label_matches_type_symbol(#[case] filename: &str) {
        let (atoms, _) = parse(BufReader::new(File::open(filename).unwrap())).unwrap();

        for atom in atoms {
            let symbol = normalize_symbol(element_from_label(&atom.data.name));
            let atomic_number = ATOMIC_SYMBOLS.iter().position(|&s| s == symbol).unwrap() + 1;
            assert_eq!(
                atomic_number as u8, atom.symbol.atomic_number(),
                "{}",
                atom.data.name
            );
        }
    }
}

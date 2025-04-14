//! Functions for parsing CIF files
//! Multiple types of CIF files are supported:
//! - CCDC/IUCr CIF: <https://www.iucr.org/resources/cif> (chemical/x-cif)
//! - PDBx/mmCIF: <https://mmcif.wwpdb.org/docs/user-guide/guide.html> (chemical/x-mmcif)
//!
//! See also: <https://en.wikipedia.org/wiki/Crystallographic_Information_File>
use super::normalize_symbol;
use crate::atom::{ATOMIC_SYMBOLS, Atom, Bond};
use nalgebra::{Matrix4, Vector3};
use std::{
    collections::HashMap,
    f32::consts::PI,
    io::{self, BufRead, BufReader, Read},
};

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
    header: &[usize; 6],
    atom_count: &mut usize,
    label_map: &mut HashMap<String, usize>,
) -> Option<Atom> {
    let vec = line.split_whitespace().collect::<Vec<_>>();

    if vec.len() <= header[3] {
        return None;
    }

    let symbol = vec[header[0]];
    let x = vec[header[1]].split('(').next()?.parse().ok()?;
    let y = vec[header[2]].split('(').next()?.parse().ok()?;
    let z = vec[header[3]].split('(').next()?.parse().ok()?;
    let id = vec[header[4]];

    let is_disordered = if header[5] != 0 && header[5] < vec.len() {
        if let Ok(disorder_group) = vec[header[5]].parse::<usize>() {
            disorder_group == 2
        } else {
            false
        }
    } else {
        false
    };

    let atomic_number = ATOMIC_SYMBOLS
        .iter()
        .position(|&s| s == normalize_symbol(symbol))?
        + 1;
    *atom_count += 1;
    label_map.insert(id.to_owned(), *atom_count);
    Some(Atom {
        id: *atom_count,
        atomic_number: atomic_number as u8,
        coord: [x, y, z],
        is_disordered,
    })
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
        atom1: map[atom1],
        atom2: map[atom2],
        order: 1,
        is_aromatic: false,
    })
}

/// Parses a CIF file and returns a vector of `Atom` and a vector of `Bond` objects.
/// # Examples
/// ```
/// use chelate::cif;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let file = File::open("data/147288.cif").unwrap();
/// let reader = BufReader::new(file);
/// let (atoms, bonds) = cif::parse(reader).unwrap();
///
/// assert_eq!(atoms.len(), 206);
/// assert_eq!(bonds.len(), 230);
/// ```
pub fn parse<P: Read>(reader: BufReader<P>) -> io::Result<(Vec<Atom>, Vec<Bond>)> {
    let mut dialect = CIFDialect::default();

    let mut atoms = Vec::new();
    let mut bonds = Vec::new();

    let mut pick_atoms = false;
    let mut pick_bonds = false;

    let mut headers = [0; 6];
    let mut header_idx = 0;
    let mut atom_count = 0;

    let mut param_idx = 0;
    let mut cell_params: [f32; 6] = Default::default();
    let mut fract_mtrx: Matrix4<f32> = Default::default();

    let mut map: HashMap<String, usize> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let line_trimmed = line.trim();
        if line_trimmed.is_empty() {
            continue;
        }

        if dialect == CIFDialect::Undefined {
            dialect = set_dialect(line_trimmed);
        } else if dialect == CIFDialect::CCDC {
            //collect cell parameters
            if line_trimmed.starts_with("_cell_length_") || line_trimmed.starts_with("_cell_angle_")
            {
                let mut iter = line_trimmed.split_whitespace();
                iter.next(); //discard name

                cell_params[param_idx] =
                    get_value_from_uncertainity(iter.next().unwrap_or_default())
                        .unwrap_or_default();
                param_idx += 1;

                if param_idx == 6 {
                    fract_mtrx = conversion_matrix_arr(cell_params);
                }
            }
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
                if line_trimmed.contains("symbol") {
                    headers[0] = header_idx;
                } else if line_trimmed.contains("fract_x") || line_trimmed.contains("Cartn_x") {
                    headers[1] = header_idx
                } else if line_trimmed.contains("fract_y") || line_trimmed.contains("Cartn_y") {
                    headers[2] = header_idx
                } else if line_trimmed.contains("fract_z") || line_trimmed.contains("Cartn_z") {
                    headers[3] = header_idx
                } else if line.starts_with("_atom_site.label_atom_id")
                    || line_trimmed.contains("atom.atom_id")
                    || line_trimmed.starts_with("_atom_site_label")
                {
                    headers[4] = header_idx
                } else if line_trimmed.contains("disorder_group") {
                    headers[5] = header_idx
                }
                header_idx += 1;
            } else if let Some(mut atom) =
                parse_atom_line(line_trimmed, &headers, &mut atom_count, &mut map)
            {
                if dialect == CIFDialect::CCDC {
                    atom.coord = fractional_to_cartesian(&atom.coord, fract_mtrx);
                }
                atoms.push(atom);
            }
        } else if pick_bonds {
            if let Some(bond) = parse_bond_line(line_trimmed, &map, &dialect) {
                bonds.push(bond);
            }
        }
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

fn get_value_from_uncertainity(input: &str) -> Option<f32> {
    input.split('(').next()?.parse().ok()
}

fn conversion_matrix_arr(array: [f32; 6]) -> Matrix4<f32> {
    conversion_matrix(array[0], array[1], array[2], array[3], array[4], array[5])
}

fn conversion_matrix(a: f32, b: f32, c: f32, alpha: f32, beta: f32, gamma: f32) -> Matrix4<f32> {
    let cos_alpha = (alpha * PI / 180.0).cos();
    let cos_beta = (beta * PI / 180.0).cos();
    let cos_gamma = (gamma * PI / 180.0).cos();
    let sin_gamma = (gamma * PI / 180.0).sin();

    let m = Matrix4::new(
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
    );

    m.transpose()
}

fn fractional_to_cartesian(fractional: &[f32; 3], matrix: Matrix4<f32>) -> [f32; 3] {
    let pos = Vector3::from_column_slice(fractional);
    let cartesian = matrix.transform_vector(&pos);
    [cartesian.x, cartesian.y, cartesian.z]
}

#[cfg(test)]
mod tests {
    use super::*;
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
        assert_eq!(atoms.iter().filter(|a| !a.is_disordered).count(), atom_len);
        assert_eq!(
            bonds
                .iter()
                .filter(|b| !atoms[b.atom1 - 1].is_disordered && !atoms[b.atom2 - 1].is_disordered)
                .count(),
            bond_len
        );
    }
}

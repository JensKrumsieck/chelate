use super::{get_next_id, normalize_symbol, reset_counter};
use crate::{ATOMIC_SYMBOLS, Atom, Bond};
use nalgebra::{Matrix4, Vector3};
use std::{
    collections::HashMap,
    f32::consts::PI,
    hash::Hash,
    io::{self, BufRead, BufReader, Read},
};

#[derive(Hash, Eq, PartialEq, Debug)]
enum CIFDialect {
    #[allow(clippy::upper_case_acronyms)]
    CCDC,
    #[allow(non_camel_case_types)]
    mmCIF,
    #[allow(non_camel_case_types)]
    compCIF,
}

#[derive(Hash, Eq, PartialEq, Debug)]
enum CIFAtomField {
    Id,
    Symbol,
    XCoord,
    YCoord,
    ZCoord,
}

fn get_field_name(dialect: &CIFDialect, field: CIFAtomField) -> &'static str {
    match (dialect, field) {
        //id
        (CIFDialect::CCDC, CIFAtomField::Id) => "_atom_site_label",
        (CIFDialect::mmCIF, CIFAtomField::Id) => "_atom_site.label_atom_id ",
        (CIFDialect::compCIF, CIFAtomField::Id) => "_chem_comp_atom.atom_id",
        //symbol
        (CIFDialect::CCDC, CIFAtomField::Symbol) => "_atom_site_type_symbol",
        (CIFDialect::mmCIF, CIFAtomField::Symbol) => "_atom_site.type_symbol",
        (CIFDialect::compCIF, CIFAtomField::Symbol) => "_chem_comp_atom.type_symbol",
        //x coord
        (CIFDialect::CCDC, CIFAtomField::XCoord) => "_atom_site_fract_x",
        (CIFDialect::mmCIF, CIFAtomField::XCoord) => "_atom_site.Cartn_x",
        (CIFDialect::compCIF, CIFAtomField::XCoord) => "_chem_comp_atom.model_Cartn_x",
        //y coord
        (CIFDialect::CCDC, CIFAtomField::YCoord) => "_atom_site_fract_y",
        (CIFDialect::mmCIF, CIFAtomField::YCoord) => "_atom_site.Cartn_y",
        (CIFDialect::compCIF, CIFAtomField::YCoord) => "_chem_comp_atom.model_Cartn_y",
        //z coord
        (CIFDialect::CCDC, CIFAtomField::ZCoord) => "_atom_site_fract_z",
        (CIFDialect::mmCIF, CIFAtomField::ZCoord) => "_atom_site.Cartn_z",
        (CIFDialect::compCIF, CIFAtomField::ZCoord) => "_chem_comp_atom.model_Cartn_z",
    }
}

/// Parses a single line of a CIF (Crystallographic information file) file and returns an `Atom` object.
/// In comparison to other file formats, cif files can have 3 different types: CCDC, mmCIF and compCIF
/// CCDC:       N1 N 0.0662(3) 0.55056(14) 0.1420(3) 0.066(2) Uani 1 1 d . . . . .
/// mmCIF:      ATOM   2    C  CA  . MET A 1 13  ? -16.763 -22.990 22.365  1.00 30.45  ? 1   MET A CA
/// compCIF:    H20 CA   CA   C  0 1 N N N 1.747  -36.297 22.990 -1.853 -0.230 0.046  CA   H20 1 
/// As columns are not fixed, this function takes a slice as header positions of [symbol, x, y, z]
fn parse_atom_line(line: &str, header: &[usize; 4]) -> Option<Atom> {
    let vec = line.split_whitespace().collect::<Vec<_>>();

    if vec.len() <= header[3] {
        return None;
    }

    let symbol = vec[header[0]];
    let x = vec[header[1]].split('(').next()?.parse().ok()?;
    let y = vec[header[2]].split('(').next()?.parse().ok()?;
    let z = vec[header[3]].split('(').next()?.parse().ok()?;

    let atomic_number = ATOMIC_SYMBOLS
        .iter()
        .position(|&s| s == normalize_symbol(symbol))?
        + 1;

    Some(Atom::new(get_next_id(), atomic_number as u8, x, y, z))
}


/// Parses a single line of an CIF file and returns a `Bond` object.
/// The line should contain the atoms names which need to be mapped to ids
/// Example line: "C4A N21A 1.370(3) . ?"
fn parse_bond_line(line: &str, map: &HashMap<String, usize>) -> Option<Bond> {
    let mut iter = line.split_whitespace();
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
/// use chelate::format::cif;
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
    reset_counter();
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();

    let mut dialect = CIFDialect::CCDC;
    let mut dialect_set = false;
    let mut headers: [usize; 4] = Default::default();
    let mut header_idx = 0;

    let mut ccdc_matrix_set = false;
    let mut param_idx = 0;
    let mut cell_params: [f32; 6] = Default::default();
    let mut fract_mtrx: Matrix4<f32> = Default::default();

    let mut map: HashMap<String, usize> = HashMap::new();

    let mut pick_atoms = false;
    let mut pick_bonds = false;

    for line in reader.lines() {
        let line = line?;
        let line_trimmed = line.trim();
        if line_trimmed.is_empty() {
            continue;
        }

        if !dialect_set {
            if let Some(d) = set_dialect(line_trimmed) {
                dialect = d;
                dialect_set = true;
            }
        } else if dialect == CIFDialect::CCDC && !ccdc_matrix_set {
            //collect cell parameters
            if line.starts_with("_cell_length_") || line.starts_with("_cell_angle_") {
                let mut iter = line.split_whitespace();
                iter.next(); //discard name

                cell_params[param_idx] =
                    get_value_from_uncertainity(iter.next().unwrap_or_default())
                        .unwrap_or_default();
                param_idx += 1;

                if param_idx == 6 {
                    ccdc_matrix_set = true;
                    fract_mtrx = conversion_matrix_arr(cell_params);
                }
            }
        }

        if line.starts_with("loop_") {
            pick_atoms = false;
            pick_bonds = false;
        }
        if line.starts_with(get_field_name(&dialect, CIFAtomField::Id)) {
            pick_atoms = true;
            pick_bonds = false;
        }
        if line.starts_with("_geom_bond") || line.starts_with("_chem_comp_bond") {
            pick_atoms = false;
            pick_bonds = true;
        }

        if pick_atoms {
            if line.starts_with("_") {
                //pick headers
                match &line {
                    l if l.starts_with(get_field_name(&dialect, CIFAtomField::Symbol)) => {
                        headers[0] = header_idx
                    }
                    l if l.starts_with(get_field_name(&dialect, CIFAtomField::XCoord)) => {
                        headers[1] = header_idx
                    }
                    l if l.starts_with(get_field_name(&dialect, CIFAtomField::YCoord)) => {
                        headers[2] = header_idx
                    }
                    l if l.starts_with(get_field_name(&dialect, CIFAtomField::ZCoord)) => {
                        headers[3] = header_idx
                    }
                    _ => {}
                }
                header_idx += 1;
            } else if let Some(mut atom) = parse_atom_line(&line, &headers) {
                if dialect == CIFDialect::CCDC {
                    atom.coord = fractional_to_cartesian(&atom.coord, fract_mtrx);
                }
                let id = line.split_whitespace().next().unwrap_or_default();
                map.insert(id.to_string(), atom.id);
                atoms.push(atom);
            }
        }

        if pick_bonds {
            if let Some(bond) = parse_bond_line(&line, &map) {
                bonds.push(bond);
            }
        }
    }

    Ok((atoms, bonds))
}

fn get_value_from_uncertainity(input: &str) -> Option<f32> {
    input.split('(').next()?.parse().ok()
}

fn set_dialect(line: &str) -> Option<CIFDialect> {
    match line {
        l if l.starts_with("_chem_comp") => Some(CIFDialect::compCIF),
        l if l.starts_with("_atom_type_symbol") || l.starts_with("_symmetry") => {
            Some(CIFDialect::CCDC)
        }
        l if l.starts_with("_pdbx") => Some(CIFDialect::mmCIF),
        _ => None,
    }
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
    #[case("data/4r21.cif", 6752, 0)]//mmcif = no bonds
    #[case("data/147288.cif", 206, 230)]
    #[case("data/1484829.cif", 466, 528)]
    #[case("data/cif_noTrim.cif", 79, 89)]
    #[case("data/cif.cif", 79, 89)]
    #[case("data/CuHETMP.cif", 85, 92)]
    #[case("data/ligand.cif", 44, 46)]
    #[case("data/mmcif.cif", 1291, 1256)]
    fn test_mol_files(#[case] filename: &str, #[case] atom_len: usize, #[case] bond_len: usize) {
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let (atoms, bonds) = parse(reader).unwrap();
        assert_eq!(atoms.len(), atom_len);
        assert_eq!(bonds.len(), bond_len);
    }
}

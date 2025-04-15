# Chelate
[![Crates.io Version](https://img.shields.io/crates/v/chelate)
![Crates.io Total Downloads](https://img.shields.io/crates/d/chelate)](https://crates.io/crates/chelate)
[![🦀 Continuous Integration](https://github.com/JensKrumsieck/chelate/actions/workflows/CI.yaml/badge.svg)](https://github.com/JensKrumsieck/chelate/actions/workflows/CI.yaml)
[![Docs.rs](https://img.shields.io/docsrs/chelate/latest)](https://docs.rs/chelate)
[![License](https://img.shields.io/badge/license-MIT%2FApache-blue.svg)](https://github.com/jenskrumsieck/chelate#license)

Chelate is a simple parser for a bunch of molecular file formats.

| Format |  MIME | Import | Export |
|--------|-------|--------|--------|
|[IUCr CIF](https://www.iucr.org/resources/cif) | `chemical/x-cif` | ✅  | ❌ |
|[PDBx/mmCIF](https://mmcif.wwpdb.org/docs/user-guide/guide.html)| `chemical/x-mmcif`|✅  | ❌ |
|[PDB](https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)| `chemical/x-pdb`|✅  | 🏗️ |
|[MOL](https://en.wikipedia.org/wiki/Chemical_table_file#Molfile)|`chemical/x-mdl-molfile` |✅  | 🏗️ |
|[MOL2](https://paulbourke.net/dataformats/mol2/) (TRIPOS)| `chemical/x-mol2`|✅  | 🏗️ |
|[XYZ](https://en.wikipedia.org/wiki/XYZ_file_format)|`chemical/x-xyz` |✅  | 🏗️ |

✅ = implemented 🏗️ = planned ❌ = not available

## Example
### Simple Example returning `(Vec<Atom>, Vec<Bond>)`
```rust
use chelate;

let (atoms, bonds) = chelate::from_file("data/147288.cif").unwrap();

assert_eq!(atoms.len(), 206);
assert_eq!(bonds.len(), 230);
```
### Example using `petgraph` returning `Graph<Atom, Edge, Undirected>`
Formats that do not have bond information like xyz are able to generate `Bond` objects when `petgraph` feature is active (default).
```rust
use chelate;
//its a Molecule (type alias for Graph<Atom, Edge, Undirected>)
let mol = chelate::molecule_from_file("data/oriluy.pdb").unwrap(); 

assert_eq!(mol.node_count(), 130);
assert_eq!(mol.edge_count(), 151);
```

## Installation
Run the following Cargo command in your project directory:
```
cargo add chelate
```
Or add the following line to your Cargo.toml:
```toml
chelate = "0.2.0"
```

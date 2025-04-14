# Chelate
[![Crates.io Version](https://img.shields.io/crates/v/chelate)
![Crates.io Total Downloads](https://img.shields.io/crates/d/chelate)](https://crates.io/crates/chelate)
[![🦀 Continuous Integration](https://github.com/JensKrumsieck/chelate/actions/workflows/CI.yaml/badge.svg)](https://github.com/JensKrumsieck/chelate/actions/workflows/CI.yaml)
[![Docs.rs](https://img.shields.io/docsrs/chelate/latest)](https://docs.rs/chelate)
[![License](https://img.shields.io/badge/license-MIT%2FApache-blue.svg)](https://github.com/jenskrumsieck/chelate#license)

Chelate is a simple parser for a bunch of molecular file formats.

Supported formats:
* PDB
* MOL
* MOL2 (TRIPOS)
* CIF (CCDC, mmCIF)
* XYZ

## Example
```rust
use chelate;

let (atoms, bonds) = chelate::from_file("data/147288.cif").unwrap();

assert_eq!(atoms.len(), 206);
assert_eq!(bonds.len(), 230);
```

## Installation
Run the following Cargo command in your project directory:
```
cargo add chelate
```
Or add the following line to your Cargo.toml:
```toml
chelate = "0.1.0"
```
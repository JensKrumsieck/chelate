# Fuzzing

Fuzz targets for the parsers, using [cargo-fuzz](https://github.com/rust-fuzz/cargo-fuzz) (libFuzzer, needs a nightly toolchain).

| Target  | Input |
|---------|-------|
| `cif`, `mol`, `mol2`, `pdb`, `xyz` | the file content, parsed with `chelate::parse` |
| `bonds` | first byte selects the format (`byte % 5`: 0 CIF, 1 MOL, 2 MOL2, 3 PDB, 4 XYZ), the rest is parsed and converted with `to_molecule` |

## Running
```
cargo install cargo-fuzz
cargo +nightly fuzz run cif fuzz/corpus/cif fuzz/seed/cif -- -dict=fuzz/dict/cif.dict -timeout=10
```
- `fuzz/corpus/<target>` collects generated inputs and is not committed.
- `fuzz/seed/<target>` holds small example files to start from.
- `fuzz/dict/<target>.dict` lists the keywords of the format, which helps libFuzzer reach the parsing code behind them.
- Add `-max_total_time=300` to stop after 5 minutes.

## Crashes
Crashing inputs are written to `fuzz/artifacts/<target>/`. To reproduce and shrink one:
```
cargo +nightly fuzz run mol fuzz/artifacts/mol/crash-<hash>
cargo +nightly fuzz tmin mol fuzz/artifacts/mol/crash-<hash>
```
Add a unit test with the shrunk input next to the fix.

## [Unreleased]

## [0.1.0] - 2026-10-09
### Added
- Project skeleton, packaging setup (`pyproject.toml`, `pixi.toml`), and CI configuration.
- Profile fetch utility for encapsulin shell HMMs (`PF04454` / Linocin_M18 family).
- Genomic I/O layer with FASTA, GFF3 parser and Pyrodigal gene calling.
- HMMER-based shell detection via `pyhmmer`.
- Genomic neighborhood window extraction around detected shells.
- Terminal peptide sequence extraction and physicochemical/motif feature computation.
- MMseqs2 clustering wrapper for cross-system candidate recurrence.
- Transparent, weighted additive scoring engine.
- Candidate table export (TSV and JSON) and neighborhood visualization (matplotlib).
- Unified CLI (`cargofinder run` and `cargofinder fetch-profiles`).
- Validation dataset with 6 experimentally characterized encapsulin systems and benchmark analysis.

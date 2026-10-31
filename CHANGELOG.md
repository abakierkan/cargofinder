# Changelog

## [Unreleased]

## [0.1.0]
### Added
- Genome FASTA (Pyrodigal) and protein FASTA + GFF3 input, including several genomes in one run.
- Shell detection with pyhmmer; `fetch-profiles` downloads PF04454 and optional extra profiles and records their versions.
- Neighborhood extraction with gene and bp distance, strand and operon likelihood.
- Terminal peptide features and user-supplied motif regexes.
- MMseqs2 clustering for recurrence across shells.
- Additive scoring with all components in the output, and ranks both across the run and within each shell.
- TSV and JSON tables, shell table, neighborhood figures and run summary.
- Validation on whole NCBI replicons, with baselines and a false-flag list.

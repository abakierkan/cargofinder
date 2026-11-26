# Changelog

## [Unreleased]
### Added
- `score_cargo_class` and `cargo.py`: candidates are scanned against 15 Pfam
  families marking the encapsulin cargo classes named by Andreas & Giessen
  (2021), who assigned cargo by identity this way. Two tiers: PF22277
  (EncFtn-like) and PF24309 (IMEF cargo Flp) are defined by Pfam as encapsulin
  cargo and score 1.0; the 13 broad enzyme families score 0.6. Accessory
  components (rhodanese, serine O-acetyltransferase) are excluded, because that
  paper states they are not cargo proteins. `fetch-profiles` downloads these
  into `profiles/cargo_classes/`. Each row reports the class and Pfam matched.
- `ScoringWeights.without_cargo_class()`: when no cargo class profiles are
  present, the term is dropped and the remaining weights rescaled to sum to 1,
  so totals stay comparable instead of being uniformly depressed.
- Validation gained a cargo-class-only baseline and a cargo class recognition
  section (`validation/cargo_class_calls.tsv`).
- `score_tp_consensus`: BLOSUM62 similarity of each terminus to the published
  Family 1 targeting-peptide consensus `SDGSLGIGSLKRS` (Giessen & Silver 2017),
  reproducing the comparison Andreas & Giessen (2021) made with `blastp`. New
  `substitution.py` holds BLOSUM62 and the ungapped window scan; no new
  dependency.
- `data/published_tp_motifs.json`: targeting-peptide anchor motifs from
  Altenburg et al. (2021), each quoting its source sentence. Opt-in via
  `--motifs`.
- Prophage screening (`prophage.py`): the genes within +/-20 of each shell are
  scanned for four phage marker families, and a locus with >=2 distinct families
  is flagged in `shells.tsv` and on every candidate row. New
  `--drop-prophage-loci` excludes flagged loci. `fetch-profiles` now also
  downloads the marker profiles into `profiles/phage_markers/`.
- Validation gained a TP-consensus-only baseline and a prophage screening
  control (`validation/prophage_control.tsv`).

### Changed
- Default weights, now: `w_distance` 0.25, `w_operon` 0.15, `w_tp_consensus`
  0.15, `w_cargo_class` 0.30, `w_tail_features` 0.05, `w_recurrence` 0.10
  (v0.1 was 0.35/0.25/-/-/0.25/0.15). The invented composition heuristic was
  down-weighted hard in favour of the two literature-grounded signals, after the
  benchmark showed it ranked only 2 of 4 adjacent cargos first against 4 of 4
  for both cargo class and consensus similarity. Scores from v0.1 are not
  comparable with these.
- README and `docs/methods.md` now state plainly that the method is not new and
  credit Andreas & Giessen (2021) as its source, and document that Family 2B and
  Family 3 have no Pfam profile and so cannot be detected at all.


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

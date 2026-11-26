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
- Two systems added to `validation/known_systems.tsv`, taking the benchmark from
  4 to 6 adjacent cargos: *Rhodospirillum rubrum* ATCC 11170 Rru_A0974/Rru_A0973
  (encapsulated ferritin, He et al. 2016, doi:10.7554/eLife.18972) and
  *Mycolicibacterium smegmatis* MC2 155 MSMEG_5830/MSMEG_5829 (MsDyP, Tang
  et al. 2021, doi:10.1073/pnas.2025658118). Both locus-tag pairs are stated in
  those papers and independently confirmed by the current RefSeq annotation.
- `score_targeting`, a family-aware targeting signal replacing v0.1's
  composition-only tail score. At a Family 1 shell (PF04454) it is BLOSUM62
  similarity of the terminus to the published targeting-peptide consensus
  `SDGSLGIGSLKRS` (Giessen & Silver 2017), reproducing the comparison Andreas &
  Giessen (2021) made with `blastp`; new `substitution.py` holds BLOSUM62 and
  the ungapped window scan, with no new dependency. At a Family 2A shell
  (PF19307) there is no short targeting peptide at all, so it is instead the
  length of the unannotated N-terminal domain, bounded by the CyD truncation
  series in Nichols et al. (2021), doi:10.7554/eLife.59288: 0 at 30 residues,
  1.0 at 225. Unknown families fall back to C-terminal mode, and every row
  records which was used in `targeting_mode`. This moved the S. elongatus
  Family 2A cargo from rank 3 to rank 1.
- `data/published_tp_motifs.json`: targeting-peptide anchor motifs from
  Altenburg et al. (2021), each quoting its source sentence. Opt-in via
  `--motifs`.
- Prophage screening (`prophage.py`): the genes within +/-20 of each shell are
  scanned for four phage marker families, and a locus with >=2 distinct families
  is flagged in `shells.tsv` and on every candidate row. New
  `--drop-prophage-loci` excludes flagged loci. `fetch-profiles` now also
  downloads the marker profiles into `profiles/phage_markers/`.
- Validation gained a targeting-signal-only baseline and a prophage screening
  control (`validation/prophage_control.tsv`).

- `fetch-profiles --extra` now accepts any syntactically valid Pfam accession,
  not only the three curated ones, with the description stating that cargofinder
  makes no claim about it. Custom non-Pfam HMMs are scanned by dropping them in
  the `--profiles` directory, which is the only way to reach Family 2B and
  Family 3.
- Added the Family 2A CyD N-terminal motif AASPYYFLDG to published_tp_motifs.json,
  with the paper's own finding that the motif alone is not sufficient for
  encapsulation.

### Changed
- **PF19307 (Family 2) is now a default shell profile** alongside PF04454.
  Family 2 is the larger group in the Andreas & Giessen (2021) census (3,523
  systems against Family 1's 2,383), so scanning Family 1 alone missed most
  encapsulin systems. On the benchmark genomes it adds three further shells,
  every one annotated by RefSeq as an encapsulin shell protein, with no false
  positives. Despite being named for Family 2A it also detects Family 2B: the
  M. xanthus MXAN_6248/MXAN_6249 pair is found at E ~ 1e-136. Benchmark recovery
  rises from 6/11 to 7/11, since the S. elongatus cargo is now reached without
  a special run.
- Recurrence counts **independent genes** rather than distinct shells. Tandem
  shells such as the Family 2B pair share their windows, so counting shells let
  a single gene appear to recur on its own; the pooled run was reporting
  clusters like "MXAN_6245, MXAN_6245". Each gene now counts once however many
  shells it neighbours. This also surfaced a real cross-genome cluster, the two
  Family 2A desulfurase cargos of M. smegmatis and S. elongatus.
- Default weights, now: `w_distance` 0.25, `w_operon` 0.15, `w_targeting`
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

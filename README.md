# cargofinder

A transparent, rule-based tool that ranks genes near encapsulin shell genes as candidate cargo.

This is a reproducible ranking tool, not a discovery claim. Every candidate is a hypothesis to test.

## What it does

`cargofinder` finds encapsulin shell genes in bacterial and archaeal genomes with HMMER, then scores the genes around each shell using:

1. **Position**: distance in genes, strand, and intergenic gaps.
2. **Terminal peptide features**: hydrophobic and Gly/Ser content, plus optional motif regexes you supply.
3. **Recurrence**: whether a gene family sits next to shells in several genomes (MMseqs2).

It writes ranked tables (TSV and JSON) with every score component, a shell table, a gene-neighborhood figure per shell, and a run summary with parameters, tool versions and profile checksums.

## Why it exists

Encapsulins enclose cargo enzymes that are often loaded through short targeting peptides. Those peptides are weakly conserved, so sequence searches alone miss cargo, and a score you cannot break down is hard to trust. `cargofinder` keeps the evidence for each candidate visible so you can judge it yourself.

## Install

### Pixi (recommended)

[Pixi](https://pixi.sh) installs everything, including MMseqs2:

```bash
git clone https://github.com/abakierkan/cargofinder.git
cd cargofinder
pixi install
pixi run pytest
```

Prefix the commands below with `pixi run`.

### pip

`cargofinder` is not on PyPI, so install from a clone:

```bash
pip install .
```

Runs with two or more shells also need [MMseqs2](https://github.com/soedinglab/MMseqs2) on your `PATH`.

## Quick start

```bash
# Download the shell profile (Pfam PF04454, with its version and checksum recorded)
cargofinder fetch-profiles --out profiles/

# Genome FASTA: genes are called with Pyrodigal
cargofinder run genome.fasta --out results/ --window 5 --tail 20

# Protein FASTA with annotation
cargofinder run proteins.faa --gff genes.gff3 --out results/

# Several genomes, so that recurrence can be measured
cargofinder run genome1.fasta genome2.fasta genome3.fasta --out results/

# Add the Family 2A profile and score both termini
cargofinder fetch-profiles --out profiles/ --extra PF19307
cargofinder run genome.fasta --profiles profiles/ --termini both --out results/
```

Outputs in `results/`:

- `candidates.tsv`, `candidates.json`: ranked candidates with all score components.
- `shells.tsv`: detected shell genes and their E-values.
- `run_summary.json`, `run_summary.txt`: counts, parameters, tool versions and profiles.
- `plots/neighborhood_<shell>.png`: one figure per shell.

*Thermotoga maritima* MSB8, whole-replicon run. The shell (`gene_0780`, locus `TM0785`) is blue and the top-ranked gene (`gene_0781`, locus `TM0786`, the known ferritin-like cargo) is orange. Figures show Pyrodigal gene numbers, not locus tags:

![Thermotoga maritima encapsulin neighborhood](docs/images/neighborhood_example.png)

## How scoring works

The score is a weighted sum. Weights are in [default_weights.json](src/cargofinder/data/default_weights.json) and can be overridden with `--weights`.

| Component | Weight | Value |
| :--- | :---: | :--- |
| `score_distance` | 0.35 | 1.0 for the adjacent gene, falling by 1/W per gene (W = window) |
| `score_operon` | 0.25 | 1.0 likely operon, 0.5 same strand, 0.0 opposite strand |
| `score_tail_features` | 0.25 | 0.6 × hydrophobic fraction + 0.4 × min(1, Gly/Ser fraction / 0.35); motif matches, if used, take 20% |
| `score_recurrence` | 0.15 | 0.2 for one shell, 0.6 for two, 1.0 for three or more |

The weights are judgment calls, not fitted values. An adjacent same-strand operon gene scores 0.63 before any tail signal in a single-shell run, so a high score mostly reflects position. Read the components, not just the total. See [docs/methods.md](docs/methods.md).

## Validation results

Run on one complete NCBI replicon per organism, with default parameters, the PF04454 profile and no motifs. Truth is [known_systems.tsv](validation/known_systems.tsv), matched by locus tag. The full generated output is in [validation/report.md](validation/report.md).

| Cargo recovered at rank 1 within its shell | Result |
| :--- | :---: |
| All cited, experimentally shown, separately encoded cargos | **4 / 9** |
| Those encoded next to the shell | 4 / 4 |
| Adjacent cargos, position alone (ties count against) | 3 / 4 |
| Adjacent cargos, tail features alone (ties count against) | 2 / 4 |

| Organism | Known cargo | Outcome |
| :--- | :--- | :--- |
| *T. maritima* MSB8 | `TM0786` Flp | Rank 1 of 10 |
| *M. tuberculosis* H37Rv | `Rv0799c` DyP | Rank 1 of 10 |
| *M. tuberculosis* H37Rv | `Rv3841` BfrB | **Missed**: outside window |
| *M. tuberculosis* H37Rv | `Rv3607c` FolB | **Missed**: outside window |
| *M. xanthus* DK 1622 | `MXAN_3557` EncB | Rank 1 of 10 |
| *M. xanthus* DK 1622 | `MXAN_4464` EncC | **Missed**: outside window |
| *M. xanthus* DK 1622 | `MXAN_2410` EncD | **Missed**: outside window |
| *R. jostii* RHA1 | `RHA1_ro02407` DypB | Rank 1 of 10 |
| *S. elongatus* PCC 7942 | `Synpcc7942_B2661` | **Missed**: shell not detected with the default profile |

Failures and caveats:

- *S. elongatus* with `PF19307` added: the shell is detected, but the cargo ranks only 3rd of 10, with or without `--termini both`.
- *P. furiosus* is a **false positive of the benchmark design**, not a recovery. `PF1192` ranks first, but it and `PF1191` are one 345-aa protein in PDB 2E0Z, split in two by the annotation.
- *R. erythropolis* PR4 `RER_59910` ranks first, but it is not counted: the experiments we found used strain N771, so this row is homology-inferred.
- **False flags**: in *M. xanthus* a sigma-54 regulator scores 0.71, close to the true cargo at 0.76. Non-cargo genes elsewhere reach 0.69. All top-3 lists are in [false_flags.tsv](validation/false_flags.tsv).
- **Recurrence**: pooling the genomes clusters only the three DyP homologues. The ferritin-like cargos do not cluster at the default thresholds.
- **Not blind**: parameters were not fitted to this set, but the weights and features were chosen with these systems in mind, and low complexity was removed from the score after reviewing an earlier benchmark.
- The set is small and related systems count more than once, so these numbers do not predict performance on new genomes. Position alone does nearly as well here.
- Only the listed replicon was scanned per organism, not whole assemblies.

Reproduce it (requires network access):

```bash
pixi run python validation/run_validation.py
```

## Assumptions and limitations

- **Targeting peptides are weakly conserved**, so motif signals are noisy.
- **Adjacency is not proof.** Some cargo is encoded elsewhere or fused to the shell, and fusions annotated as two genes can look like adjacent cargo.
- **Cargo without a recognizable terminal peptide will be missed** or get no tail credit. By default only the C-terminus is scored.
- **Only PF04454 is scanned by default.** Family 2A requires `--extra PF19307`. PF05065 is a general phage capsid family, so expect prophage hits. Other encapsulin families are not covered in v0.1.
- Contigs are treated as linear, so windows do not wrap around circular replicons.
- Scores from single-genome and multi-genome runs are not comparable, because recurrence is constant with one shell.

## Roadmap

- More encapsulin family profiles.
- Structure prediction and peptide-shell docking.
- Machine-learning scoring, if enough validated examples exist.
- Web UI.
- Metagenome-scale parallelization.

## How to cite

```bibtex
@software{cargofinder,
  author = {Erkan, Abdulbaki},
  title = {cargofinder: candidate cargo discovery and ranking for encapsulin nanocompartments},
  year = {2026},
  url = {https://github.com/abakierkan/cargofinder}
}
```

See [CITATION.cff](CITATION.cff).

This project was developed with AI coding assistance, and the author reviewed the methods.

## License

[MIT](LICENSE)

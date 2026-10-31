# Methods: `cargofinder` v0.1

This document describes what the code does, why, and where it falls short. The formulas match `src/cargofinder/scoring.py`.

## 1. Biological background

Encapsulins are prokaryotic protein shells with an HK97-like fold that enclose cargo enzymes (Sutter 2008; McHugh 2014). Cargo is usually loaded through a short targeting peptide, often a C-terminal extension (Sutter 2008; Cassidy-Amstutz 2016), that binds the inner surface of the shell. Other routes exist: the *S. elongatus* cysteine desulfurase uses an N-terminal targeting sequence (Nichols 2021), and the *P. furiosus* shell protein has a ferritin-like segment in the same polypeptide (Akita 2007, PDB 2E0Z).

Cargo genes are often next to the shell gene, but not always. *M. tuberculosis* BfrB and FolB (Contreras 2014) and *M. xanthus* EncC and EncD (McHugh 2014) are encoded elsewhere.

Targeting peptides are weakly conserved, so the output is a ranking of candidates, never a definitive call.

## 2. Pipeline

1. **Input** (`io.py`, `genecalling.py`). One or more nucleotide FASTA files, with genes called by Pyrodigal, or one protein FASTA with a GFF3.
   - Inputs under 100 kb use Pyrodigal's meta mode, because Pyrodigal warns that training needs at least 100 kb.
   - A protein FASTA passed without `--gff` is an error, not an empty result.
   - GFF3 features with no matching protein are skipped with a warning and counted in the run summary. They also shorten gene distances, so treat that warning seriously.
   - With several genomes, IDs become `<genome>|<id>`.
2. **Shell detection** (`shells.py`). `hmmsearch` through pyhmmer. A hit needs a full-sequence E-value of at most 1e-5 and a protein of at least 150 aa. A gene hitting several profiles keeps its best hit. Results are in `shells.tsv`.
3. **Neighborhood** (`neighborhood.py`). Up to ±N genes (default 5) on the same contig. Contigs are treated as linear, so a window does not wrap around the origin of a circular replicon. Other detected shells are not candidates.
   - A neighbor is a likely operon member when every gene from the shell to the neighbor shares its strand and every gap is at most 150 bp. The gap threshold is a heuristic.
4. **Clustering** (`cluster.py`). With two or more shells in a run, all neighbors are clustered with MMseqs2 `easy-cluster` (30% identity, 80% coverage, one thread for reproducibility). With one shell there is no recurrence to measure, so MMseqs2 is skipped.
5. **Terminal peptide features** (`features.py`). The last and first K residues (default 20).
6. **Scoring** (`scoring.py`). See below.
7. **Output** (`reporting.py`, `plotting.py`). Tables, figures and a run summary with parameters, tool versions and profile checksums.

## 3. Scoring

$$\text{Score} = 0.35\,S_\text{distance} + 0.25\,S_\text{operon} + 0.25\,S_\text{tail} + 0.15\,S_\text{recurrence}$$

Weights must be non-negative and sum to 1. Defaults are in `src/cargofinder/data/default_weights.json`; override them with `--weights`. Each row reports all four components.

| Component | Definition |
| :--- | :--- |
| $S_\text{distance}$ | $(W - d + 1) / W$ for gene distance $d \le W$. With $W=5$: 1.0, 0.8, 0.6, 0.4, 0.2 |
| $S_\text{operon}$ | 1.0 for a likely operon member, 0.5 for same strand otherwise, 0.0 for opposite strand |
| $S_\text{tail}$ | $0.6\,h + 0.4\,\min(1, g / 0.35)$, where $h$ is the fraction of A, V, I, L, M, F, W, P and $g$ is the fraction of G, S |
| $S_\text{tail}$ with motifs | $0.8 \times$ the above $+\ 0.2\,m$, where $m$ is 1 if any motif matches the C-terminus |
| $S_\text{recurrence}$ | 0.2 for one shell, 0.6 for two, 1.0 for three or more shells in the cluster |

**Why these weights.** They are judgment calls, not fitted values. Adjacency gets the most weight because it is the best documented signal. The tail gets less because it is weakly conserved. Recurrence gets the least because it depends on what else is in the run.

**Consequences to keep in mind:**

- In a single-shell run the recurrence term is a constant 0.03.
- An adjacent same-strand operon gene therefore scores 0.63 before any tail signal. A high score mostly reflects position.
- Scores from single-genome and multi-genome runs are not comparable.
- Sequence complexity is reported but not scored. We found no evidence that it should count in either direction.
- `--termini both` uses the better terminus. It is off by default because signal peptides also make N-termini hydrophobic.
- `rank` is across the run. `rank_in_shell` is within one shell's neighborhood, and it is the rank the figures show and validation uses.

## 4. Validation

`validation/run_validation.py` runs one complete replicon per organism and matches genes to `known_systems.tsv` by locus tag, using strand and stop coordinate. `validation/report.md` is generated from the results. It compares against adjacency-only and tail-only rankings and lists the top non-cargo genes.

## 5. Limitations

1. **Adjacency is not proof of encapsulation.** Neighbors may be unrelated genes.
2. **Distal cargo is missed**, and no window size fixes that.
3. **Fusions can look like adjacent cargo.** *P. furiosus* PF1191 and PF1192 are annotated as two ORFs, but PDB 2E0Z shows one protein. PF1192 is ranked first as a "neighbor", which is an artifact. Whether the split is biological or a sequencing error was not determined.
4. **Cargo without a recognizable terminal peptide gets no tail credit.** Targeting peptides are weakly conserved, so motif signals are noisy, and any motif match is only a hint.
5. **Profile coverage.** Only PF04454 is scanned by default. Family 2A needs PF19307. PF05065 is a general phage capsid family, so expect prophage hits if you add it. Other encapsulin families may have no suitable profile here. See Andreas & Giessen (2021) for the wider diversity.
6. **E-value depends on the sequences searched**, including the length filter, so values are not comparable across differently sized inputs.
7. **Recurrence can be inflated** by tandem shells whose windows overlap, because it counts shells rather than independent loci.

## 6. References

- Akita F et al. (2007) The crystal structure of a virus-like particle from the hyperthermophilic archaeon *Pyrococcus furiosus* provides insight into the evolution of viruses. *J Mol Biol* 368:1469-1483. doi:10.1016/j.jmb.2007.02.075
- Andreas MP, Giessen TW (2021) Large-scale computational discovery and analysis of virus-derived microbial nanocompartments. *Nat Commun* 12. doi:10.1038/s41467-021-25071-y
- Cassidy-Amstutz C et al. (2016) Identification of a minimal peptide tag for in vivo and in vitro loading of encapsulin. *Biochemistry* 55:3461-3468. doi:10.1021/acs.biochem.6b00294
- Contreras H et al. (2014) Characterization of a *Mycobacterium tuberculosis* nanocompartment and its potential cargo proteins. *J Biol Chem* 289:18279-18289. doi:10.1074/jbc.M114.570119
- Giessen TW, Silver PA (2017) Widespread distribution of encapsulin nanocompartments reveals functional diversity. *Nat Microbiol* 2:17029. doi:10.1038/nmicrobiol.2017.29
- Kyte J, Doolittle RF (1982) A simple method for displaying the hydropathic character of a protein. *J Mol Biol* 157:105-132.
- McHugh CA et al. (2014) A virus capsid-like nanocompartment that stores iron and protects bacteria from oxidative stress. *EMBO J* 33:1896-1911. doi:10.15252/embj.201488566
- Nichols RJ et al. (2021) Discovery and characterization of a novel family of prokaryotic nanocompartments involved in sulfur metabolism. *eLife* 10:e59288. doi:10.7554/eLife.59288
- Rahmanpour R, Bugg TDH (2013) Assembly in vitro of *Rhodococcus jostii* RHA1 encapsulin and peroxidase DypB to form a nanocompartment. *FEBS J* 280:2097-2104. doi:10.1111/febs.12234
- Sutter M et al. (2008) Structural basis of enzyme encapsulation into a bacterial nanocompartment. *Nat Struct Mol Biol* 15:939-947. doi:10.1038/nsmb.1473
- Tamura A et al. (2015) Packaging guest proteins into the encapsulin nanocompartment from *Rhodococcus erythropolis* N771. *Biotechnol Bioeng* 112:13-20. doi:10.1002/bit.25322

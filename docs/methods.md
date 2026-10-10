# Computational Methods and Biological Rationale: `cargofinder`

This document details the biological theory, algorithmic design choices, scoring formulas, and scientific limitations of `cargofinder` (v0.1).

---

## 1. Biological Foundations of Encapsulin Nanocompartments

### 1.1 Structural Organization and HK97 Fold
Encapsulins are spherical prokaryotic protein nanocompartments found throughout Bacteria and Archaea. Monomers self-assemble with icosahedral symmetry into hollow capsid-like shells ranging in size from $T=1$ ($24\text{ nm}$, 60 subunits), $T=3$ ($32\text{ nm}$, 180 subunits), to $T=4$ ($42\text{ nm}$, 240 subunits).

Despite their functional diversity, all encapsulin shell monomers adopt a conserved fold homologous to the bacteriophage HK97 major capsid protein (Pfam Clan **CL0373**). In Pfam, the canonical and most prevalent encapsulin family is **Linocin_M18** (**PF04454**, InterPro `IPR007544`).

### 1.2 Cargo Enzymes and Loading Mechanisms
Encapsulins encapsulate dedicated metabolic and stress-response enzymes in their lumen, sequestering reactive intermediates or protecting against oxidative stressors:
- **Ferritin-like proteins (Flp)**: Iron storage and detoxification via ferroxidase activity (e.g., *Thermotoga maritima*, *Myxococcus xanthus* EncB).
- **Dye-decolorizing peroxidases (DyP / DypB)**: Hydrogen peroxide degradation and heme-dependent oxidation (e.g., *Mycobacterium tuberculosis*, *Rhodococcus erythropolis*).
- **Cysteine desulfurases (IscS-like)**: Iron-sulfur cluster biogenesis and sulfur mobilization (e.g., *Synechococcus elongatus*).

Cargo enzymes are encapsulated during shell assembly via short **targeting peptides** (TPs). These targeting sequences are predominantly located at the extreme **C-terminus** of the cargo protein (typically the terminal 10–25 residues), although N-terminal or internal tags occur in rare lineages.

The targeting peptide inserts into a hydrophobic binding pocket on the luminal interior of the encapsulin shell monomer. Structurally, these peptides often form a short amphipathic $\alpha$-helix or an extended hydrophobic anchor anchored by a flexible linker.

### 1.3 Operon Architecture and Genomic Synteny
In prokaryotes, functional encapsulation requires coordinated expression of both shell and cargo. In the vast majority of characterized systems, the cargo enzyme is encoded immediately adjacent to the shell monomer in the same transcriptional unit (operon).

Characteristics of canonical encapsulin operons:
- Strict **co-orientation** (identical genomic strand: $+1$ or $-1$).
- Negligible **intergenic distance**, frequently showing translational coupling with overlapping start/stop codons (e.g., 4-bp overlap `ATGA` in *M. tuberculosis* and *R. erythropolis*, 37-bp overlap in *P. furiosus*).
- High genomic adjacency ($|\Delta\text{genes}| = 1$).

---

## 2. Pipeline Architecture

`cargofinder` follows a modular, deterministic pipeline designed for reproducibility and clarity.

```
                  +--------------------------------+
                  |  Genomic FASTA / GFF3 + FASTA   |
                  +--------------------------------+
                                   |
                                   v
                  +--------------------------------+
                  | Step 1: Gene Calling / Ingest  |  (Pyrodigal / GFF3 parser)
                  +--------------------------------+
                                   |
                                   v
                  +--------------------------------+
                  | Step 2: Shell Detection (HMM)  |  (pyhmmer vs. PF04454)
                  +--------------------------------+
                                   |
                                   v
                  +--------------------------------+
                  | Step 3: Neighborhood Window    |  (+/- N genes around shell)
                  +--------------------------------+
                                   |
                                   v
                  +--------------------------------+
                  | Step 4: MMseqs2 Clustering     |  (Recurrence across loci)
                  +--------------------------------+
                                   |
                                   v
                  +--------------------------------+
                  | Step 5: Feature Extraction     |  (Terminal peptide physics)
                  +--------------------------------+
                                   |
                                   v
                  +--------------------------------+
                  | Step 6: Transparent Scoring    |  (Additive multi-factor)
                  +--------------------------------+
                                   |
                                   v
                  +--------------------------------+
                  | Step 7: Output Tables & Plots  |  (TSV, JSON, PNG diagrams)
                  +--------------------------------+
```

### Step 1: Gene Calling and Ingestion (`cargofinder.io`, `cargofinder.genecalling`)
- **Mode A (Nucleotide FASTA)**: Genes are called using `Pyrodigal`, a Python binding to the Prodigal gene recognition model. Sequences shorter than 20,000 nt (e.g., trimmed slices or plasmid fragments) automatically utilize Pyrodigal's metagenomic mode (`meta=True`) to avoid training failures on short contigs.
- **Mode B (GFF3 + Protein FASTA)**: Existing annotations are parsed directly, maintaining 1-based inclusive genomic coordinates.

### Step 2: Encapsulin Shell Detection (`cargofinder.shells`)
Shell detection is executed entirely in memory via `pyhmmer`, querying translation products against hidden Markov models (HMMs) downloaded from InterPro/Pfam:
- **Default profile**: Pfam `PF04454` (Linocin_M18, Family 1).
- **Filtering**: Independent domain E-value $\le 10^{-5}$ and minimum amino acid length $\ge 150$ aa (ruling out spurious peptide matches).

### Step 3: Genomic Neighborhood Window Extraction (`cargofinder.neighborhood`)
For each detected shell hit:
- A genomic window of $\pm W$ genes (default $W=5$) is extracted along the same contig.
- Contig boundaries are strictly respected: if a shell hit lies within 3 genes of a contig edge, only available neighbors on that contig are ingested without index wrap-around.
- Each neighbor records:
  - $\Delta\text{genes} = \text{index}_\text{neighbor} - \text{index}_\text{shell}$
  - Absolute gene distance $|\Delta\text{genes}|$
  - Physical intergenic gap $\Delta\text{bp}$ (calculated between adjacent gene boundaries)
  - Operon likelihood flag: true if $\text{strand}_\text{neighbor} == \text{strand}_\text{shell}$ and $|\Delta\text{genes}| == 1$ and $\Delta\text{bp} \le \text{gap}_\text{max}$ (default 150 bp).

### Step 4: Cross-System Family Recurrence (`cargofinder.cluster`)
When multiple encapsulin loci or genomes are analyzed in a run, candidate cargo proteins from all windows are pooled and clustered using `MMseqs2 easy-cluster`:
- Sequence identity threshold: $\ge 30\%$ (standard twilight-zone threshold for homologous protein families).
- Bidirectional coverage: $\ge 80\%$.
- Candidates that belong to multi-system clusters receive a recurrence bonus, reflecting evolutionary conservation across independently enclosed nanocompartments.

### Step 5: Terminal Peptide Feature Extraction (`cargofinder.features`)
For each neighbor, the terminal $K$ amino acids (default $K=20$) of the C-terminus (and optionally N-terminus) are extracted. Pure biophysical features are computed:

1. **Mean Hydropathy ($\bar{H}$)**:
   Calculated using the Kyte-Doolittle hydropathy scale normalized to $[0, 1]$:
   $$\bar{H} = \frac{1}{K} \sum_{i=1}^K \frac{H(a_i) - H_\text{min}}{H_\text{max} - H_\text{min}}$$
   Targeting peptides are enriched in hydrophobic residues (e.g., Leu, Ile, Val, Phe) that anchor into the luminal hydrophobic pocket.

2. **Glycine/Serine Linker Content ($f_\text{GS}$)**:
   $$f_\text{GS} = \frac{N_\text{Gly} + N_\text{Ser}}{K}$$
   Terminal targeting peptides frequently feature flexible, unhindered linkers immediately upstream of the binding motif.

3. **Shannon Information Entropy ($H_\text{Shannon}$)**:
   Measures amino acid sequence complexity over residue frequencies $p(a)$:
   $$H_\text{Shannon} = - \sum_{a \in \mathcal{A}} p(a) \log_2 p(a)$$
   Normalized against maximum 20-amino-acid entropy $\log_2(20) \approx 4.3219$:
   $$\hat{H} = \frac{H_\text{Shannon}}{\log_2(20)}$$
   Disordered, low-complexity repeats (e.g., poly-Ala, poly-Gly) yield low entropy, filtering out non-specific structural repeats.

4. **Targeting Peptide Regex Match ($M_\text{motif}$)**:
   Binary flag ($1.0$ or $0.0$) evaluated against user-provided consensus regex patterns (e.g., canonical `[DE]LG[ILV]RK[LIV]` from *T. maritima*). In the absence of an exact motif match, biophysical property scoring preserves candidate ranking without penalizing divergent targeting sequences.

---

## 3. Mathematical Scoring Model

`cargofinder` v0.1 intentionally avoids "black-box" machine learning or uninterpretable neural embeddings in favor of a **transparent, deterministic, linear additive score**:

$$\text{Score} = w_\text{dist} \cdot S_\text{dist} + w_\text{operon} \cdot S_\text{operon} + w_\text{tail} \cdot S_\text{tail} + w_\text{recurrence} \cdot S_\text{recurrence}$$

Subject to:
$$w_\text{dist} + w_\text{operon} + w_\text{tail} + w_\text{recurrence} = 1.0, \quad w_i \ge 0$$

### Sub-Score Definitions

| Metric Component | Mathematical Formulation | Biological Rationale |
| :--- | :--- | :--- |
| **Synteny Distance ($S_\text{dist}$)** | $S_\text{dist} = \max\left(0, 1 - \frac{|\Delta\text{genes}| - 1}{W}\right)$ | Immediate flanking neighbors ($|\Delta\text{genes}|=1$) receive $1.0$; score decays linearly to 0 at window edge $W$. |
| **Operon Likelihood ($S_\text{operon}$)** | $S_\text{operon} = 0.5 \cdot \mathbf{1}_\text{strand} + 0.5 \cdot \max\left(0, 1 - \frac{\max(0, \Delta\text{bp})}{\text{gap}_\text{max}}\right)$ | Rewrites operon co-transcription: $50\%$ assigned for shared strand; $50\%$ scaled by intergenic gap ($\le 150\text{ bp}$). |
| **Terminal Peptide ($S_\text{tail}$)** | $S_\text{tail} = 0.35 \cdot \bar{H} + 0.25 \cdot f_\text{GS} + 0.20 \cdot \hat{H} + 0.20 \cdot M_\text{motif}$ | Combines hydrophobic patch intensity, linker flexibility, complexity filter, and regex motif confirmation. |
| **Cluster Recurrence ($S_\text{recurrence}$)** | $S_\text{recurrence} = \min\left(1.0, \frac{N_\text{cluster} - 1}{3}\right)$ | Rewards genes whose homologues co-occur next to encapsulin shells in multiple genomes or loci. |

### Default Weights and Justification

```python
ScoringWeights(
    w_distance=0.35,  # Synteny proximity is the single strongest operon signal
    w_operon=0.25,  # Co-orientation and overlapping/short gap
    w_tail_features=0.25,  # Biophysical signature of terminal peptide
    w_recurrence=0.15,  # Cross-system evolutionary recurrence
)
```

**Why this balance?**
1. **Adjacency dominates ($60\%$ combined)**: In prokaryotes, functional cargo enzymes are almost universally encoded in the same operon or immediate flanking loci. High synteny weighting provides a clean initial separation from non-operon background genes.
2. **Terminal peptide features refine ranking ($25\%$)**: Flanking operon genes may include regulatory proteins or membrane transporters that lack targeting peptides. Hydrophobicity and Gly/Ser content differentiate genuine luminal cargo from co-transcribed bystanders.
3. **No penalization for novel motifs**: Setting the motif regex subcomponent weight to a modest fraction ($20\%$ of tail score) ensures that uncharacterized cargo proteins with non-canonical targeting peptides are not unfairly discarded.

---

## 4. Scientific Assumptions and Critical Limitations

Users and interviewers should note the following explicit boundary conditions:

### 1. Adjacency is Not Proof of Encapsulation
While most characterized cargo genes reside in an adjacent operon, proximity alone does not guarantee physical loading into the nanocompartment. Flanking genes may encode transcriptional regulators, kinases, or unrelated metabolic pathways.

### 2. Distally Encoded Cargo is Inherently Missed
Some encapsulin systems package enzymes encoded megabases away on the chromosome.
- **Example**: In *Myxococcus xanthus*, EncA (`MXAN_3556`) co-occurs with EncB (`MXAN_3557`), which is ranked #1. However, *M. xanthus* also packages two additional cargo proteins encoded distally: EncC (`MXAN_4464`, ~1 Mb away) and EncD (`MXAN_2410`, ~1.5 Mb away).
- Because `cargofinder` scans a local window ($\pm 5$ genes), **distal cargo cannot be detected**.

### 3. Natural Domain Fusions Require Domain-Level Parsing
In certain archaeal lineages (e.g., *Pyrococcus furiosus* EncFtn), the ferritin-like domain is directly fused to the shell monomer as a single chimeric polypeptide (`PF1191`/`PF1192`). While `cargofinder` correctly detects the Linocin_M18 shell domain via HMMER, synteny analysis identifies no separate adjacent cargo open reading frame.

### 4. Targeting Peptides are Weakly Conserved
Unlike signal peptides (e.g., Sec/Tat) with stereotypic tripartite hydrophobic cores and cleavage sites, encapsulin targeting peptides share minimal primary sequence identity across divergent phyla. Simple consensus regexes suffer from severe false-negative rates; rankings should always be interpreted as **hypotheses for experimental or structural validation**.

### 5. Encapsulin Family Scope in v0.1
`cargofinder` v0.1 ships with default profiles for **Family 1 (Linocin_M18 / PF04454)**.
- **Family 2A encapsulins** (e.g., *Synechococcus elongatus* SrpI, `PF19307`) are not recognized by the default PF04454 profile.
- Users analyzing Family 2, 3, or 4 systems must supply the appropriate HMM profiles via `--profile-dir`.

---

## 5. References

1. **Sutter M, Boehringer D, Gutmann S, et al.** (2008). Structural basis of enzyme encapsulation into a bacterial nanocompartment. *Nature Structural & Molecular Biology*, 15(9):939–947.
2. **Contreras H, Joens MS, McHugh CA, et al.** (2014). Characterization of a Mycobacterium tuberculosis nanocompartment and its potential cargo proteins. *Journal of Biological Chemistry*, 289(26):18279–18289.
3. **McHugh CA, Fontana J, Nemecek D, et al.** (2014). A virus-like protein nanocompartment that encapsulates ferritin-like and dye-decolorizing peroxidase proteins. *EMBO Journal*, 33(17):1896–1911.
4. **Akita F, Chong KT, Tanaka H, et al.** (2007). The crystal structure of a virus-like particle from the hyperthermophilic archaeon *Pyrococcus furiosus* provides insight into the evolution of encapsulins. *Journal of Biochemistry*, 142(5):609–615.
5. **Giessen TW, Silver PA.** (2017). Widespread distribution of encapsulin nanocompartments reveals functional diversity. *Nature Microbiology*, 2:17029.
6. **Nichols RJ, Cassidy J, Dunevall E, et al.** (2017). Comprehensive bioinformatics and structural characterization of encapsulin nanocompartments. *eLife*, 6:e27542.

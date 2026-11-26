# Methods: `cargofinder` v0.2

This document describes what the code does, why, and where it falls short. The formulas match `src/cargofinder/scoring.py`.

## 1. Biological background

Encapsulins are prokaryotic protein shells with an HK97-like fold that enclose cargo enzymes (Sutter 2008; McHugh 2014). Cargo is usually loaded through a short targeting peptide, often a C-terminal extension (Sutter 2008; Cassidy-Amstutz 2016), that binds the inner surface of the shell. Other routes exist: the *S. elongatus* cysteine desulfurase uses an N-terminal targeting sequence (Nichols 2021), and the *P. furiosus* shell protein has a ferritin-like segment in the same polypeptide (Akita 2007, PDB 2E0Z).

Cargo genes are often next to the shell gene, but not always. *M. tuberculosis* BfrB and FolB (Contreras 2014) and *M. xanthus* EncC and EncD (McHugh 2014) are encoded elsewhere.

Targeting peptides are weakly conserved, so the output is a ranking of candidates, never a definitive call.

## 1a. Relation to prior work

The method here is not new. The genome-mining approach it implements was
established by **Andreas & Giessen (2021)**, who searched prokaryotic genomes
with HMM profiles, computed genome neighborhoods extending 20 ORFs either side
of each shell gene, and assigned Family 1 cargo by comparing proteins encoded
within 5000 bp of the shell against the targeting-peptide consensus
`SDGSLGIGSLKRS` from Giessen & Silver (2017). They curated over 6,000
encapsulin-like systems that way and defined the family numbering used here.
Trujillo & Komeili (2025) apply comparable comparative genomics across
prokaryotic organelles more broadly.

cargofinder deliberately reimplements that logic rather than extending it. What
it adds is engineering and transparency, not biological insight:

1. An installable, deterministic tool with recorded tool versions and profile
   checksums, rather than an analysis performed once.
2. A score broken into named components, all of which appear in every output
   row, so a ranking can be inspected instead of trusted.
3. A benchmark against literature-curated systems that reports its own failures,
   including the components that do not help.

Anyone needing a curated catalogue of encapsulin systems should use the
published datasets. cargofinder is for ranking candidates in a genome that is
not in those datasets, and for seeing why each candidate ranks where it does.

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

$$\text{Score} = 0.25\,S_\text{distance} + 0.15\,S_\text{operon} + 0.15\,S_\text{target} + 0.30\,S_\text{cargo} + 0.05\,S_\text{tail} + 0.10\,S_\text{recurrence}$$

Weights must be non-negative and sum to 1. Defaults are in `src/cargofinder/data/default_weights.json`; override them with `--weights`. Each row reports all four components.

| Component | Definition |
| :--- | :--- |
| $S_\text{distance}$ | $(W - d + 1) / W$ for gene distance $d \le W$. With $W=5$: 1.0, 0.8, 0.6, 0.4, 0.2 |
| $S_\text{operon}$ | 1.0 for a likely operon member, 0.5 for same strand otherwise, 0.0 for opposite strand |
| $S_\text{cargo}$ | 1.0 if the candidate matches an encapsulin-cargo-specific Pfam family (PF22277, PF24309), 0.6 for a broad cargo-class family, 0.0 for no match |
| $S_\text{target}$ | **Family 1 shell:** best normalised ungapped BLOSUM62 similarity of the terminal slice to the published Family 1 TP consensus `SDGSLGIGSLKRS`. **Family 2A shell:** length of the unannotated N-terminal extension, ramped from 0 at 30 residues to 1 at 225. See section 3c |
| $S_\text{tail}$ | $0.6\,h + 0.4\,\min(1, g / 0.35)$, where $h$ is the fraction of A, V, I, L, M, F, W, P and $g$ is the fraction of G, S |
| $S_\text{tail}$ with motifs | $0.8 \times$ the above $+\ 0.2\,m$, where $m$ is 1 if any motif matches the C-terminus |
| $S_\text{recurrence}$ | 0.2 for one gene, 0.6 for two, 1.0 for three or more independent genes in the cluster |

**Why these weights.** They are judgment calls, not fitted values. Adjacency gets
the most weight because it is the best documented signal. $S_\text{target}$ gets more
than $S_\text{tail}$ because it rests on published, family-specific evidence, whereas
$S_\text{tail}$ is a cargofinder heuristic with no literature basis. Recurrence
gets the least because it depends on what else is in the run. $S_\text{cargo}$
gets the largest single weight because cargo identity is the evidence the
literature actually assigns cargo on, but it is held to 0.30 rather than higher
for two reasons given in section 3b.

If no cargo class profiles are available, the term is dropped and the remaining
weights are rescaled to sum to 1, so totals stay on $[0,1]$ and remain
comparable with runs that had the profiles.

**On $S_\text{target}$ in C-terminal mode.** This reproduces, in a small deterministic form, the
comparison Andreas & Giessen (2021) made with `blastp` against the same
consensus. It is not BLAST: there is no E-value and no significance claim, only a
bounded similarity. Two caveats matter when reading the number:

- The consensus is itself Gly/Ser-rich, so generic Gly/Ser linkers pick up around
  0.3 as background. The `GGSGGSGGSGG` control peptide, which Altenburg et al.
  (2021) showed does not bind the shell pocket, scores about as high as some
  native targeting peptides. Similarity alone does not separate a real anchor
  from a flexible linker.
- The consensus is a **Family 1** consensus, so it is only applied at Family 1
  loci. Family 2 cargo has no reason to match it: before targeting became
  family-aware, the *S. elongatus* cysteine desulfurase scored 0.24 here and
  ranked third. Family 4 has no documented mechanism and falls back to this
  mode, which is a default rather than a claim.

The anchor motifs in `src/cargofinder/data/published_tp_motifs.json` address the
first caveat in principle, since Altenburg et al. resolved which residues
actually contact the shell. They are opt-in via `--motifs` and were not used to
produce the benchmark numbers. Each one is three to six residues long and will
match many non-cargo termini by chance.

**What the benchmark says about these components.** On the six adjacent cargos,
$S_\text{target}$ alone ranks all 6 first, as does $S_\text{cargo}$; adjacency alone
ranks 5 of 6 first; the composition heuristic $S_\text{tail}$ alone ranks 3 of 6
first. Six cases is far too few to call that a result, and several of the systems
are related to one another, but it is the reason $S_\text{tail}$ was cut from
0.25 to 0.05 in favour of the two signals with a literature basis. See
`validation/report.md`.

## 3c. Family-aware targeting

Families do not recruit cargo the same way, so one targeting formula cannot
serve them all. Andreas & Giessen (2021): "In classical encapsulins (Family 1),
encapsulation is mediated by short C-terminal peptide sequences referred to as
targeting peptides (TPs) ... while for Family 2 systems, larger N-terminal
protein domains are proposed to mediate encapsulation." For Family 2 cysteine
desulfurases specifically, "no specific targeting peptide has been identified",
and "the unannotated N-terminal domain has been shown to be responsible for
mediating encapsulation."

Scoring a Family 2 candidate against a Family 1 C-terminal consensus therefore
measures nothing at all. $S_\text{target}$ picks its mode from the Pfam family
of the shell that found the locus:

| Shell profile | Family | Targeting mode |
| :--- | :--- | :--- |
| PF04454 | 1 | C-terminal TP consensus similarity |
| PF19307 | 2A, and 2B in practice | N-terminal extension length |
| anything else | unknown | C-terminal, as a documented fallback |

The mode used is recorded in `targeting_mode` on every row, so the fallback is
never silent.

**Measuring the extension.** It is the number of residues before the *earliest*
domain found by *any* cargo class profile. Taking the best-scoring domain
instead would be wrong: a DyP peroxidase matches `Dyp_perox_C` around residue
140, so the best hit alone implies a spurious ~140-residue "extension" for an
ordinary two-domain protein. Measured from the earliest domain, the benchmark's
Family 1 cargos have extensions of 3, 7, 8 and 11 residues, while the Family 2A
cargo has 242.

**The bounds come from an experiment, not from taste.** Nichols et al. (2021)
fused truncations of the *S. elongatus* CyD N-terminal domain to sfGFP: the full
225-residue domain "yielded the highest loading signal", the first 100 residues
"also functioned, albeit not as efficiently", a 29-residue fragment carrying a
conserved motif "was not sufficient to sequester cargo", and deleting the domain
abolished loading. So the score is 0 at and below 30 residues and 1.0 at 225,
ramping linearly between; 100 residues lands at about 0.36. The linear shape
between the anchors is a judgment call, the anchors are not.

**Limits.** The extension can only be measured where a cargo class domain is
recognised, so a Family 2 candidate of an unknown class scores 0 — the same as
a candidate with no extension. Length is also a crude proxy for what is really
a disordered-domain interaction; the paper used DISOPRED3 to show the domain is
disordered, which cargofinder does not attempt. And Family 2B has no Pfam shell
profile, so it never reaches this code path at all.

## 3b. Cargo class recognition

v0.1 had no notion of *what a gene is*. That is not how cargo is assigned in the
literature: Andreas & Giessen (2021) built "HMM templates representative of each
cargo class" and classified hits "based on their Pfam or Interpro annotation".

`cargo.py` scans candidates against the Pfam families marking the cargo classes
that paper names. Family 1: DyP peroxidases ("the most abundant cargo class in
Family 1", 1,505 systems) and the ferritin-like superfamily, which they report
corresponds to "Flp, rubrerythrin, CoQ7, Mn catalase, IMEF, hemerythrin, and
bacterioferritin". Family 2: cysteine desulfurases (they give the domain as
PF00266 directly), terpene cyclases, polyprenyl transferases and xylulose
kinases.

Serine *O*-acetyltransferases and rhodaneses are **deliberately excluded**
despite being the most enriched Family 2A operon components, because that paper
states such accessory components "are not cargo proteins".

Two tiers, because the families differ enormously in specificity:

| Tier | Families | Score |
| :--- | :--- | :---: |
| Defined by Pfam as encapsulin cargo | PF22277 (EncFtn-like), PF24309 (IMEF cargo Flp) | 1.0 |
| Broad enzyme families containing cargo among much else | the other 13 | 0.6 |

Where a protein matches both tiers, the encapsulin-specific call wins regardless
of E-value, since it is the more informative one.

**Two limits that cap how much this can be trusted.**

1. *It cannot separate cargo from a co-located paralogue.* In the benchmark,
   every non-cargo class match was in the *P. furiosus* locus: a further
   ferritin-like protein and three rubrerythrins. Andreas & Giessen note that
   "Flp-fusion encapsulins are often located in operons containing other
   ferritin-like proteins or rubrerythrins", so those genes genuinely are that
   class of enzyme. Position and targeting peptide are what distinguish them.
2. *It can only recognise cargo classes that are already known.* Novel cargo
   scores zero, so a candidate with no class match is not excluded. This
   component must never be read as a filter. The benchmark cannot measure this
   failure mode at all, because every system in it is already characterised.

## 3a. Prophage screening

Encapsulin shell proteins and bacteriophage major capsid proteins share the HK97
fold and sit in the same Pfam clan, CL0373. A profile search for a shell
therefore also returns phage capsids, and no E-value threshold separates them
because the similarity is real. Andreas & Giessen (2021) report that their Family
1 and Family 2 datasets "were significantly contaminated with bacteriophage
capsid proteins", and removed those hits using genome neighborhood context.

`prophage.py` does the same locally. It scans the genes within ±20 of each shell
for four phage marker families and flags a locus carrying at least two distinct
ones. The screening radius is wider than the ±5 cargo window on purpose: a
prophage spans far more genes than could plausibly be cargo, and within ±5 even
an obvious prophage capsid rarely shows two marker families.

| Accession | Name | Role |
| :--- | :--- | :--- |
| PF03354 | TerL_ATPase | Terminase large subunit |
| PF04586 | Peptidase_S78 | Caudovirus prohead protease |
| PF00589 | Phage_integrase | Integrase |
| PF05135 | Phage_connect_1 | Head-tail connector |

The flag never enters the score. It is evidence about whether the **locus** is an
encapsulin, not about whether a given neighbour is cargo, and mixing the two
would conflate them. It appears in `shells.tsv` and on every candidate row;
`--drop-prophage-loci` excludes flagged loci.

**Known gap in our source.** Andreas & Giessen list the domains they filtered on
as "PF0860, PF03354, PF04586, PF00589, PF05135". `PF0860` is not a well-formed
Pfam accession, which have five digits, and the nearest valid accession, PF00860,
is `Xan_ur_permease`, a xanthine/uracil permease family with no phage
association. We could not determine what was intended and omitted it rather than
guess. **Please check this against the paper.**

**Measured behaviour, not assumed.** In the benchmark control, 0 of 8 true
encapsulin shells were wrongly flagged and 2 of 9 phage capsid hits were caught.
Sensitivity is therefore weak: the misses carry one marker family or none within
the window. An unflagged locus means "not evidently prophage", never "confirmed
encapsulin".

**Consequences to keep in mind:**

- In a single-shell run the recurrence term is a constant 0.02.
- An adjacent same-strand operon gene scores 0.42 with no tail or cargo signal.
  Recognised cargo in that position reaches roughly 0.72, so position alone no
  longer dominates the way it did in v0.1.
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
5. **Profile coverage, and a remaining gap.** PF04454 and PF19307 are scanned by
   default, covering Family 1 and Family 2; Family 4B needs `--extra PF08967`.
   Checking Pfam clan CL0373 against the InterPro API on 2026-11-26, those three
   are the **only** encapsulin-specific families in Pfam; every other member of
   the clan is a phage capsid family.

   PF19307 was promoted to a default in v0.2. Family 2 is the larger group in the
   Andreas & Giessen (2021) census (3,523 systems against Family 1's 2,383), so
   scanning Family 1 alone missed most encapsulin systems. On the benchmark
   genomes it added three further shells, all annotated by RefSeq as encapsulin
   shell proteins, with no false positives. Although Pfam names it for Family 2A,
   it also detects Family 2B shells: the *M. xanthus* MXAN_6248/MXAN_6249 pair is
   annotated Family 2B and is found at E ~ 1e-136.

   **Family 3 still has no Pfam profile**, so cargofinder cannot detect it out of
   the box. Andreas & Giessen built custom HMMs for Family 2 and Family 3 for
   exactly this reason. Any `.hmm` placed in the directory given to `--profiles`
   is scanned as a shell profile, which is how to cover it; `fetch-profiles
   --extra` also accepts any Pfam accession, though cargofinder makes no claim
   that an arbitrary family detects encapsulins. PF05065 is a general phage
   capsid family, so expect prophage hits if you add it.
6. **E-value depends on the sequences searched**, including the length filter, so values are not comparable across differently sized inputs.
7. **Recurrence counts independent genes, not shells.** Tandem shells, such as the two adjacent Family 2B shell genes of *M. xanthus*, have overlapping windows and share neighbours, so counting shells let a single gene appear to recur on its own. Each gene now counts once however many shells it neighbours. Recurrence can still be inflated by closely related genomes in one run.

## 6. References

- Altenburg WJ, Rollins N, Silver PA, Giessen TW (2021) Exploring targeting peptide-shell interactions in encapsulin nanocompartments. *Sci Rep* 11:4951. doi:10.1038/s41598-021-84329-z
- Bibb LA, Hatfull GF (2002) Integration and excision of the *Mycobacterium tuberculosis* prophage-like element, phiRv1. *Mol Microbiol* 45:1515-1526. doi:10.1046/j.1365-2958.2002.03130.x
- Fan X, Abd Alla AA, Xie J (2016) Distribution and function of prophage phiRv1 and phiRv2 among *Mycobacterium tuberculosis* complex. *J Biomol Struct Dyn* 34:233-238. doi:10.1080/07391102.2015.1022602
- Henikoff S, Henikoff JG (1992) Amino acid substitution matrices from protein blocks. *PNAS* 89:10915-10919. doi:10.1073/pnas.89.22.10915
- Trujillo HA, Komeili A (2025) Revealing the diversity of bacterial and archaeal organelles via comparative genomics. *Mol Biol Cell* 36:pe4. doi:10.1091/mbc.e20-08-0564
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

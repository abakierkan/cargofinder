"""Input/Output layer for genomic sequences, GFF3 annotations, and data models."""

from __future__ import annotations

import urllib.parse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple


@dataclass
class Gene:
    """Represents a protein-coding gene in a genome.

    Coordinates follow standard 1-based inclusive genomic convention (GFF3 style).
    Strand is +1 (forward) or -1 (reverse).
    """

    gene_id: str
    contig_id: str
    start: int
    end: int
    strand: int
    gene_index: int
    protein_seq: str
    product: Optional[str] = None

    def __post_init__(self) -> None:
        if self.start > self.end:
            raise ValueError(
                f"Gene start ({self.start}) cannot exceed end ({self.end}) "
                f"for gene {self.gene_id}"
            )
        if self.strand not in (+1, -1):
            raise ValueError(
                f"Gene strand must be +1 or -1 (got {self.strand}) "
                f"for gene {self.gene_id}"
            )
        # Strip trailing stop codon character if present
        if self.protein_seq.endswith("*"):
            self.protein_seq = self.protein_seq[:-1]

    @property
    def length_nt(self) -> int:
        """Length in nucleotides."""
        return self.end - self.start + 1

    @property
    def length_aa(self) -> int:
        """Length in amino acids."""
        return len(self.protein_seq)


@dataclass
class GenomeData:
    """Container holding ordered genes for a single genome."""

    genome_id: str
    genes: List[Gene]
    contigs: Dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self._gene_map: Dict[str, Gene] = {g.gene_id: g for g in self.genes}
        self._by_contig: Dict[str, List[Gene]] = {}
        for g in self.genes:
            self._by_contig.setdefault(g.contig_id, []).append(g)
        # Ensure genes per contig are sorted by start coordinate
        for ctg in self._by_contig:
            self._by_contig[ctg].sort(key=lambda x: (x.start, x.end))
            for idx, g in enumerate(self._by_contig[ctg]):
                g.gene_index = idx

    def get_gene(self, gene_id: str) -> Optional[Gene]:
        """Lookup a gene by its unique identifier."""
        return self._gene_map.get(gene_id)

    def genes_on_contig(self, contig_id: str) -> List[Gene]:
        """Return ordered list of genes on the specified contig."""
        return self._by_contig.get(contig_id, [])

    @property
    def contig_ids(self) -> List[str]:
        """List of contig identifiers present in this genome."""
        return list(self._by_contig.keys())

    def __len__(self) -> int:
        return len(self.genes)


def read_fasta_sequences(path: Path | str) -> Dict[str, str]:
    """Read a multi-record FASTA file into a mapping of ID -> sequence.

    Lines are stripped and concatenated. Header comments (after whitespace)
    are ignored; the identifier is taken as the first token after '>'.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    records: Dict[str, List[str]] = {}
    current_id: Optional[str] = None

    with p.open("r", encoding="utf-8", errors="replace") as fh:
        for line_num, line in enumerate(fh, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip()
                if not header:
                    raise ValueError(
                        f"Malformed FASTA at {path}:{line_num}: empty header"
                    )
                # ID is the first whitespace-delimited token
                current_id = header.split()[0]
                if current_id in records:
                    raise ValueError(
                        f"Duplicate FASTA ID '{current_id}' at {path}:{line_num}"
                    )
                records[current_id] = []
            else:
                if current_id is None:
                    raise ValueError(
                        f"Malformed FASTA at {path}:{line_num}: sequence before header"
                    )
                records[current_id].append(line.replace(" ", "").upper())

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    return {k: "".join(v) for k, v in records.items()}


def write_fasta_sequences(
    sequences: Dict[str, str] | Iterable[Tuple[str, str]],
    path: Path | str,
    line_width: int = 60,
) -> None:
    """Write sequences to a FASTA file."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)

    items: Iterator[Tuple[str, str]]
    if isinstance(sequences, dict):
        items = iter(sequences.items())
    else:
        items = iter(sequences)

    with p.open("w", encoding="utf-8") as fh:
        for header, seq in items:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                fh.write(f"{seq[i : i + line_width]}\n")


def parse_gff3(gff_path: Path | str, feature_type: str = "CDS") -> List[Dict]:
    """Parse GFF3 features into structured dictionaries.

    Raises ValueError if GFF3 rows are malformed (wrong columns or bad coordinates).
    """
    p = Path(gff_path)
    if not p.is_file():
        raise FileNotFoundError(f"GFF3 file not found: {gff_path}")

    features: List[Dict] = []
    with p.open("r", encoding="utf-8", errors="replace") as fh:
        for line_num, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) != 9:
                raise ValueError(
                    f"Malformed GFF3 line {line_num} in {gff_path}: "
                    f"expected 9 tab-separated fields, got {len(parts)}"
                )

            seqid, source, ftype, start_s, end_s, score, strand_s, phase, attrs_s = (
                parts
            )

            if feature_type and ftype != feature_type:
                continue

            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid integer coordinates at line {line_num} in {gff_path}: "
                    f"start='{start_s}', end='{end_s}'"
                ) from exc

            strand = -1 if strand_s == "-" else +1

            attrs: Dict[str, str] = {}
            for item in attrs_s.split(";"):
                item = item.strip()
                if not item:
                    continue
                if "=" in item:
                    key, val = item.split("=", 1)
                    attrs[key.strip()] = urllib.parse.unquote(val.strip())

            # Determine gene identifier: check ID, locus_tag, Name, Parent
            gid = (
                attrs.get("ID")
                or attrs.get("locus_tag")
                or attrs.get("Name")
                or attrs.get("Parent")
            )
            if not gid:
                gid = f"{seqid}_{start}_{end}"

            product = (
                attrs.get("product") or attrs.get("Note") or attrs.get("description")
            )

            features.append(
                {
                    "gene_id": gid,
                    "contig_id": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "product": product,
                    "attributes": attrs,
                }
            )

    if not features:
        raise ValueError(
            f"No features of type '{feature_type}' found in GFF3 file {gff_path}"
        )

    return features


def load_protein_and_gff3(
    protein_fasta: Path | str,
    gff_path: Path | str,
    genome_id: Optional[str] = None,
    feature_type: str = "CDS",
) -> GenomeData:
    """Load genome data by combining a protein FASTA and a GFF3 annotation file.

    Features and proteins are matched by gene_id or matching attribute tokens.
    """
    p_fasta = Path(protein_fasta)
    p_gff = Path(gff_path)

    if genome_id is None:
        genome_id = p_fasta.stem

    proteins = read_fasta_sequences(p_fasta)
    gff_records = parse_gff3(p_gff, feature_type=feature_type)

    # Map candidate IDs to protein sequences
    genes: List[Gene] = []
    missing_proteins: List[str] = []

    for rec in gff_records:
        gid = rec["gene_id"]
        attrs = rec["attributes"]
        seq = None

        # Try direct match
        if gid in proteins:
            seq = proteins[gid]
        else:
            # Check other common attribute keys: locus_tag, Name, ID, protein_id
            for alt_key in ("locus_tag", "Name", "ID", "protein_id"):
                alt_val = attrs.get(alt_key)
                if alt_val and alt_val in proteins:
                    seq = proteins[alt_val]
                    gid = alt_val
                    break

        if seq is None:
            missing_proteins.append(gid)
            continue

        genes.append(
            Gene(
                gene_id=gid,
                contig_id=rec["contig_id"],
                start=rec["start"],
                end=rec["end"],
                strand=rec["strand"],
                gene_index=0,  # will be assigned by GenomeData post_init
                protein_seq=seq,
                product=rec.get("product"),
            )
        )

    if not genes:
        sample_gff = [r["gene_id"] for r in gff_records[:5]]
        sample_faa = list(proteins.keys())[:5]
        raise ValueError(
            f"Could not match any GFF3 features with protein FASTA entries.\n"
            f"Sample GFF3 IDs: {sample_gff}\n"
            f"Sample FASTA IDs: {sample_faa}"
        )

    return GenomeData(genome_id=genome_id, genes=genes)

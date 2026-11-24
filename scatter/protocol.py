"""Protocol-level read selection and transcript-orientation helpers.

Alignment layout, evidence mate, transcript strand, and TSO source are kept
independent deliberately.  A paired BAM does not imply that Read 1 is the
assay endpoint, and alignment strand does not necessarily equal RNA strand.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import pysam


ALIGNMENT_LAYOUTS = {"paired", "single"}
PAIRED_R1_MIN_BP = 100
EVIDENCE_MATES = {"auto", "read1", "read2", "single"}
ASSAY_TYPES = {"3prime", "5prime", "full_length"}
LIBRARY_STRANDS = {"gene", "forward", "reverse"}
TSO_SOURCES = {"evidence", "tag", "none"}
ENDPOINT_RULES = {"auto", "physical-5prime", "transcript-5prime"}
LEGACY_LIBRARY_LAYOUT = {"10x_5prime_pe": "paired", "10x_5prime_se": "single"}


@dataclass(frozen=True)
class GeneRecord:
    gene_id: str
    chrom: str
    strand: str


def forced_alignment_layout_for_r1(r1_length_bp: int) -> str:
    """Return the mandatory alignment layout for one raw Read 1 length.

    SCATTER's project-wide policy is deliberately binary: datasets whose R1
    is shorter than 100 bp must use an unpaired/single alignment, while an R1
    of at least 100 bp is eligible for paired alignment.  The threshold is
    applied to the complete sequenced R1, including any CB/UMI prefix.
    """
    if r1_length_bp < 1:
        raise ValueError("R1 length must be positive")
    return "single" if r1_length_bp < PAIRED_R1_MIN_BP else "paired"


def clean_gene(value: object) -> str | None:
    """Return one versionless gene ID, rejecting STAR ambiguous assignments."""
    token = str(value) if value is not None else ""
    if not token or token == "-" or "," in token or ";" in token:
        return None
    return token.split(".")[0]


def resolve_alignment_layout(library_type: str | None, alignment_layout: str | None) -> str:
    """Resolve the explicit layout while retaining the legacy library alias."""
    if alignment_layout is not None and alignment_layout not in ALIGNMENT_LAYOUTS:
        raise ValueError(f"Unknown alignment layout: {alignment_layout}")
    if library_type is not None and library_type not in LEGACY_LIBRARY_LAYOUT:
        raise ValueError(f"Unknown legacy library type: {library_type}")
    legacy_layout = LEGACY_LIBRARY_LAYOUT.get(library_type) if library_type else None
    if alignment_layout and legacy_layout and alignment_layout != legacy_layout:
        raise ValueError(
            f"Conflicting layout options: --alignment-layout {alignment_layout} "
            f"versus --library-type {library_type}"
        )
    resolved = alignment_layout or legacy_layout
    if resolved is None:
        raise ValueError("Specify --alignment-layout or the legacy --library-type")
    return resolved


def resolve_evidence_mate(alignment_layout: str, evidence_mate: str | None) -> str:
    """Resolve and validate the BAM mate used as one molecular endpoint."""
    requested = evidence_mate or "auto"
    if requested not in EVIDENCE_MATES:
        raise ValueError(f"Unknown evidence mate: {requested}")
    resolved = ("read1" if alignment_layout == "paired" else "single") if requested == "auto" else requested
    if alignment_layout == "paired" and resolved not in {"read1", "read2"}:
        raise ValueError("Paired layout requires --evidence-mate read1 or read2")
    if alignment_layout == "single" and resolved != "single":
        raise ValueError("Single layout requires --evidence-mate single")
    return resolved


def record_matches_mate(read: pysam.AlignedSegment, evidence_mate: str) -> bool:
    if evidence_mate == "read1":
        return bool(read.is_read1)
    if evidence_mate == "read2":
        return bool(read.is_read2)
    if evidence_mate == "single":
        return not read.is_read1 and not read.is_read2
    raise ValueError(f"Unresolved evidence mate: {evidence_mate}")


def load_gene_catalog(reference_tss_path: str | Path) -> dict[str, GeneRecord]:
    """Load the versioned reference's sibling genes.tsv as a stable-ID map."""
    genes_path = Path(reference_tss_path).expanduser().resolve().parent / "genes.tsv"
    if not genes_path.is_file():
        raise FileNotFoundError(f"Expected sibling genes.tsv beside reference TSS catalog: {genes_path}")
    result: dict[str, GeneRecord] = {}
    with genes_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = clean_gene(row.get("gene_id"))
            strand = row.get("strand")
            if gene is None or strand not in {"+", "-"}:
                continue
            record = GeneRecord(gene, row["chrom"], strand)
            previous = result.get(gene)
            if previous is not None and previous != record:
                raise ValueError(f"Reference gene ID has conflicting loci or strands: {gene}")
            result[gene] = record
    if not result:
        raise ValueError(f"Reference gene catalog contains no usable records: {genes_path}")
    return result


def alignment_strand(read: pysam.AlignedSegment) -> str:
    return "-" if read.is_reverse else "+"


def transcript_strand(read: pysam.AlignedSegment, gene: GeneRecord, library_strand: str) -> str:
    """Resolve RNA strand from GX annotation or a forward/reverse library rule."""
    if library_strand == "gene":
        return gene.strand
    observed = alignment_strand(read)
    if library_strand == "forward":
        return observed
    if library_strand == "reverse":
        return "+" if observed == "-" else "-"
    raise ValueError(f"Unknown library strand: {library_strand}")


def evidence_five_position(read: pysam.AlignedSegment) -> int:
    """Return the selected sequencing record's physical 5′ genomic endpoint.

    Library-strand conversion changes the RNA strand assigned to the evidence,
    not which end of the selected read was sequenced first.  Keeping these two
    operations independent is essential for reverse-stranded libraries.
    """
    if not read.is_reverse:
        return int(read.reference_start) + 1
    return int(read.reference_end)


def resolve_endpoint_rule(requested: str, assay_type: str, alignment_layout: str) -> str:
    """Use the transcript-facing boundary for single-read 5-prime evidence."""
    if requested not in ENDPOINT_RULES:
        raise ValueError(f"Unknown endpoint rule: {requested}")
    if requested == "auto":
        return "transcript-5prime" if assay_type == "5prime" and alignment_layout == "single" else "physical-5prime"
    return requested


def transcript_five_position(read: pysam.AlignedSegment, strand: str) -> int:
    """Return the alignment boundary toward RNA 5-prime, in 1-based coordinates.

    For reverse-stranded R2 this is its physical 3-prime aligned boundary.
    This selects a boundary of the observed alignment; it does not infer an
    unsequenced R1 start using fragment length or extrapolate soft clips.
    """
    if strand == "+":
        return int(read.reference_start) + 1
    if strand == "-":
        return int(read.reference_end)
    raise ValueError(f"Invalid transcript strand: {strand}")


def evidence_five_softclip(read: pysam.AlignedSegment) -> str:
    """Return the soft clip at the selected sequencing record's physical 5′ end."""
    cigar = read.cigartuples or []
    sequence = read.query_sequence or ""
    if not cigar or not sequence:
        return ""
    if not read.is_reverse:
        operation, length = cigar[0]
        return sequence[:length] if operation == 4 else ""
    operation, length = cigar[-1]
    return sequence[-length:] if operation == 4 else ""


def biological_five_softclip(read: pysam.AlignedSegment, technical_prefix_bp: int = 0) -> str:
    """Return the physical 5′ soft clip after removing a sequenced prefix.

    STAR retains bases removed with ``--clip5pNbases`` in the terminal soft
    clip. On a reverse alignment those original 5′ bases occur at the right
    side of the stored soft clip, whereas on a forward alignment they occur at
    the left side.
    """
    if technical_prefix_bp < 0:
        raise ValueError("technical_prefix_bp must be non-negative")
    clip = evidence_five_softclip(read)
    if technical_prefix_bp == 0:
        return clip
    if len(clip) <= technical_prefix_bp:
        return ""
    return clip[:-technical_prefix_bp] if read.is_reverse else clip[technical_prefix_bp:]


def has_five_softclip_alignment_anchor(
    read: pysam.AlignedSegment, min_aligned_bases: int = 6
) -> bool:
    """Require a terminal soft clip adjacent to a substantive alignment."""
    if min_aligned_bases < 1:
        raise ValueError("min_aligned_bases must be positive")
    cigar = read.cigartuples or []
    if len(cigar) < 2:
        return False
    adjacent = cigar[-2] if read.is_reverse else cigar[1]
    terminal = cigar[-1] if read.is_reverse else cigar[0]
    # M, =, and X all consume aligned query/reference bases.
    return terminal[0] == 4 and adjacent[0] in {0, 7, 8} and adjacent[1] >= min_aligned_bases


def boolean_tag(read: pysam.AlignedSegment, tag: str) -> bool:
    """Interpret an explicitly transferred BAM tag as a boolean TSO result."""
    if not read.has_tag(tag):
        raise KeyError(tag)
    value = read.get_tag(tag)
    if isinstance(value, (bool, int, float)):
        return bool(value)
    token = str(value).strip().lower()
    if token in {"1", "true", "t", "yes", "y", "positive", "pos"}:
        return True
    if token in {"0", "false", "f", "no", "n", "negative", "neg"}:
        return False
    raise ValueError(f"BAM tag {tag} is not boolean-like: {value!r}")

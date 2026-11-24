"""Strand-aware span annotation separated from evidence-based GX assignment."""
from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass, field
import csv
from pathlib import Path


@dataclass
class GeneIntervalIndex:
    intervals: list[tuple[int, int, str, str]]
    starts: list[int] = field(init=False)
    prefix_max_end: list[int] = field(init=False)
    genes: dict[str, str] = field(init=False)

    def __post_init__(self) -> None:
        self.intervals.sort()
        self.starts = [row[0] for row in self.intervals]
        self.prefix_max_end = []
        self.genes = {}
        maximum = 0
        for _, end, gene, strand in self.intervals:
            maximum = max(maximum, end)
            self.prefix_max_end.append(maximum)
            if gene in self.genes and self.genes[gene] != strand:
                raise ValueError(f"Conflicting reference strands for {gene}")
            self.genes[gene] = strand

    def overlaps(self, start: int, end: int) -> list[tuple[int, int, str, str]]:
        if start > end:
            raise ValueError("Interval start must not exceed its end")
        right = bisect_right(self.starts, end)
        left = bisect_right(self.prefix_max_end, start - 1, 0, right)
        return [row for row in self.intervals[left:right] if row[1] >= start]


def load_genes(reference_tss_path: str | Path) -> dict[str, GeneIntervalIndex]:
    path = Path(reference_tss_path).resolve().parent / "genes.tsv"
    if not path.is_file():
        raise FileNotFoundError(f"Expected sibling genes.tsv beside reference TSS catalog: {path}")
    grouped: dict[str, list[tuple[int, int, str, str]]] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            grouped.setdefault(row["chrom"], []).append((
                int(row["start_1"]), int(row["end_1"]), row["gene_id"].split(".")[0], row["strand"],
            ))
    return {chrom: GeneIntervalIndex(rows) for chrom, rows in grouped.items()}


def annotate_gene_overlap(rows: list[dict[str, object]], genes: dict[str, GeneIntervalIndex]) -> None:
    """Annotate spans without invalidating an already unique GX assignment.

    The caller's evidence contract requires unique GX, matching chromosome,
    and an explicitly resolved gene strand before records enter its database.
    Antisense, nested, and same-strand span overlaps are annotations, not new
    ambiguous read assignments. Unknown GX/chromosome/strand remains invalid.
    Endpoints outside a gene span are flagged, not silently reassigned.
    """
    for row in rows:
        index = genes.get(str(row["chrom"]))
        strand = str(row["strand"])
        gene = str(row["gene_id"]).split(".")[0]
        overlaps = index.overlaps(int(row["start_1"]), int(row["end_1"])) if index else []
        all_ids = sorted({item[2] for item in overlaps})
        same_ids = sorted({item[2] for item in overlaps if item[3] == strand})
        opposite_ids = sorted({item[2] for item in overlaps if item[3] != strand})
        valid = index is not None and index.genes.get(gene) == strand
        row["overlapping_gene_ids"] = ";".join(all_ids)
        row["same_strand_overlapping_gene_ids"] = ";".join(same_ids)
        row["opposite_strand_overlapping_gene_ids"] = ";".join(opposite_ids)
        row["gene_assignment_evidence"] = "unique_GX" if valid else "invalid_GX_reference"
        row["gene_assignment_status"] = "unique" if valid else "ambiguous"
        row["assigned_gene_span_overlap"] = gene in all_ids
        if not valid and row.get("filter_status") == "pass":
            row["filter_status"] = "fail_ambiguous_gene"

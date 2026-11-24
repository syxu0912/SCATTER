"""Strand-aware GTF reference construction for SCATTER."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import csv
import re
from .io import open_text


_ATTRIBUTE = re.compile(r'\s*([^\s;]+)\s+"([^"]*)"\s*')


@dataclass(frozen=True)
class TranscriptTSS:
    gene_id: str
    transcript_id: str
    chrom: str
    strand: str
    tss_0: int


def parse_attributes(value: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for token in value.rstrip(";\n").split(";"):
        match = _ATTRIBUTE.fullmatch(token)
        if match:
            attributes[match.group(1)] = match.group(2)
    return attributes


def build_reference(gtf_path: str | Path, output_dir: str | Path, merge_distance: int = 50) -> dict[str, int]:
    """Build gene and reference-TSS catalogs from a GTF.

    Input GTF coordinates are one-based inclusive.  CSV coordinates written by
    this function are one-based biological positions, while BED output follows
    the conventional zero-based half-open scheme.
    """
    source = Path(gtf_path).expanduser().resolve()
    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    genes: dict[str, dict[str, object]] = {}
    transcripts: list[TranscriptTSS] = []
    transcript_spans = {}

    if merge_distance < 0:
        raise ValueError('Reference TSS merge distance must be non-negative')
    with open_text(source) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[6] not in {"+", "-"}:
                continue
            chrom, _, feature, start, end, _, strand, _, raw_attributes = fields
            attributes = parse_attributes(raw_attributes)
            gene_id = attributes.get("gene_id")
            if not gene_id:
                continue
            start_0 = int(start) - 1
            end_0 = int(end)
            if feature == "gene":
                genes[gene_id] = {
                    "gene_id": gene_id,
                    "gene_name": attributes.get("gene_name", ""),
                    "gene_type": attributes.get("gene_type", attributes.get("gene_biotype", "")),
                    "chrom": chrom,
                    "start_0": start_0,
                    "end_0": end_0,
                    "strand": strand,
                }
            elif feature == "transcript":
                transcript_id = attributes.get("transcript_id")
                if not transcript_id:
                    continue
                tss_0 = start_0 if strand == "+" else end_0 - 1
                transcripts.append(TranscriptTSS(gene_id, transcript_id, chrom, strand, tss_0))
                key = (gene_id, chrom, strand)
                a, b = transcript_spans.get(key, (start_0, end_0))
                transcript_spans[key] = (min(a, start_0), max(b, end_0))

    # Fall back to transcript spans for a GTF that omits an explicit gene line.
    for transcript in transcripts:
        record = genes.get(transcript.gene_id)
        if record is None:
            genes[transcript.gene_id] = {
                "gene_id": transcript.gene_id,
                "gene_name": "",
                "gene_type": "",
                "chrom": transcript.chrom,
                "start_0": transcript_spans[(transcript.gene_id, transcript.chrom, transcript.strand)][0],
                "end_0": transcript_spans[(transcript.gene_id, transcript.chrom, transcript.strand)][1],
                "strand": transcript.strand,
            }

    gene_path = destination / "genes.tsv"
    with gene_path.open("wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["gene_id", "gene_name", "gene_type", "chrom", "start_1", "end_1", "strand"],
            delimiter="\t",
        )
        writer.writeheader()
        for record in sorted(genes.values(), key=lambda item: (str(item["chrom"]), int(item["start_0"]))):
            writer.writerow(
                {
                    "gene_id": record["gene_id"],
                    "gene_name": record["gene_name"],
                    "gene_type": record["gene_type"],
                    "chrom": record["chrom"],
                    "start_1": int(record["start_0"]) + 1,
                    "end_1": int(record["end_0"]),
                    "strand": record["strand"],
                }
            )

    grouped: dict[tuple[str, str, str], list[TranscriptTSS]] = defaultdict(list)
    for transcript in transcripts:
        grouped[(transcript.gene_id, transcript.chrom, transcript.strand)].append(transcript)

    groups: list[dict[str, object]] = []
    for (gene_id, chrom, strand), members in grouped.items():
        members.sort(key=lambda item: item.tss_0)
        current: list[TranscriptTSS] = []
        for member in members:
            if current and member.tss_0 - current[-1].tss_0 > merge_distance:
                groups.append(_make_tss_group(current, gene_id, chrom, strand))
                current = []
            current.append(member)
        if current:
            groups.append(_make_tss_group(current, gene_id, chrom, strand))

    tss_path = destination / "reference_tss.tsv"
    bed_path = destination / "reference_tss.bed"
    fields = ["reference_tss_id", "gene_id", "chrom", "tss_1", "strand", "transcript_ids"]
    with tss_path.open("wt", newline="", encoding="utf-8") as tsv, bed_path.open("wt", encoding="utf-8") as bed:
        writer = csv.DictWriter(tsv, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for group in sorted(groups, key=lambda item: (str(item["chrom"]), int(item["tss_0"]))):
            writer.writerow(
                {
                    "reference_tss_id": group["reference_tss_id"],
                    "gene_id": group["gene_id"],
                    "chrom": group["chrom"],
                    "tss_1": int(group["tss_0"]) + 1,
                    "strand": group["strand"],
                    "transcript_ids": ";".join(group["transcript_ids"]),
                }
            )
            bed.write(
                f"{group['chrom']}\t{group['tss_0']}\t{int(group['tss_0']) + 1}\t"
                f"{group['reference_tss_id']}\t0\t{group['strand']}\n"
            )
    return {"genes": len(genes), "transcripts": len(transcripts), "reference_tss_groups": len(groups)}


def _make_tss_group(members: list[TranscriptTSS], gene_id: str, chrom: str, strand: str) -> dict[str, object]:
    representative = members[len(members) // 2].tss_0
    return {
        "reference_tss_id": f"{gene_id}|Known|{chrom}_{representative + 1}_{strand}",
        "gene_id": gene_id,
        "chrom": chrom,
        "strand": strand,
        "tss_0": representative,
        "transcript_ids": sorted(member.transcript_id for member in members),
    }

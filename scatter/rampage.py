"""Parsing utilities for strand-aware RAMPAGE peak benchmarks."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


class RampageParseError(ValueError):
    """Raised when a RAMPAGE BED record cannot be interpreted."""


_SUMMIT_IN_NAME = re.compile(r"^[^|]+\|[+-]\|(\d+)(?:\||$)")


@dataclass(frozen=True)
class RampagePeak:
    chrom: str
    start_0: int
    end_0: int
    summit_0: int
    strand: str
    name: str

    @property
    def summit_1(self) -> int:
        """One-based biological coordinate used by SCATTER TSV/CSV outputs."""
        return self.summit_0 + 1


def parse_high_peak_bed(path: str | Path) -> list[RampagePeak]:
    """Parse ENCODE-like high-confidence RAMPAGE BED files.

    The K562 ENCSR000AER file encodes its zero-based summit coordinate as the
    third pipe-delimited field in BED column 4.  We retain both BED and
    SCATTER coordinate conventions to make the conversion auditable.
    """
    peak_path = Path(path).expanduser().resolve()
    if not peak_path.is_file():
        raise RampageParseError(f"RAMPAGE BED does not exist: {peak_path}")
    peaks: list[RampagePeak] = []
    with peak_path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise RampageParseError(f"Line {line_number} has fewer than six BED columns.")
            chrom, start, end, name, _, strand = fields[:6]
            if strand not in {"+", "-"}:
                raise RampageParseError(f"Line {line_number} has invalid strand {strand!r}.")
            match = _SUMMIT_IN_NAME.match(name)
            if match is None:
                raise RampageParseError(
                    f"Line {line_number} does not contain an encoded summit in column 4."
                )
            start_0, end_0, summit_0 = int(start), int(end), int(match.group(1))
            if not start_0 <= summit_0 < end_0:
                raise RampageParseError(
                    f"Line {line_number} summit {summit_0} lies outside [{start_0}, {end_0})."
                )
            peaks.append(RampagePeak(chrom, start_0, end_0, summit_0, strand, name))
    if not peaks:
        raise RampageParseError("RAMPAGE BED contains no peaks.")
    return peaks

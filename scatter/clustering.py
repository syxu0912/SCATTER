"""Valley-aware subdivision of adaptive-gap components.

The width limit is an explicit engineering constraint, not a biological
definition of a promoter. Every input site is assigned to exactly one child;
the caller subsequently applies the existing read threshold to each child.
"""
from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks


@dataclass(frozen=True)
class SplitConfig:
    max_width_bp: int = 200
    valley_ratio: float = 0.25
    smoothing_sigma_bp: float = 5.0
    min_peak_distance_bp: int = 50

    def __post_init__(self) -> None:
        if self.max_width_bp < 0:
            raise ValueError("max_width_bp must be non-negative; zero selects legacy clustering")
        if not 0 <= self.valley_ratio < 1:
            raise ValueError("valley_ratio must be in [0, 1)")
        if self.smoothing_sigma_bp <= 0 or self.min_peak_distance_bp < 1:
            raise ValueError("Smoothing sigma and peak distance must be positive")


def split_sites(
    sites: list[tuple[int, int]],
    config: SplitConfig,
    min_cluster_reads: int,
) -> list[tuple[list[tuple[int, int]], str]]:
    """Partition one sorted component without losing or duplicating sites.

    Supported modes separated by a deep valley are split first. Components
    still wider than the cap are bisected at a low-density boundary in their
    central half. Equal-density boundaries prefer balanced coordinate spans.
    Recursion is iterative, and even unsupported fragments are returned so
    that the caller can account for their discarded reads explicitly.
    """
    if not sites:
        return []
    if any(c <= 0 for _, c in sites) or any(a[0] >= b[0] for a, b in zip(sites, sites[1:])):
        raise ValueError("Sites must have positive counts and strictly increasing coordinates")
    if config.max_width_bp == 0:
        return [(sites, "legacy")]
    origin = sites[0][0]
    signal = np.zeros(sites[-1][0] - origin + 1, dtype=float)
    for position, count in sites:
        signal[position - origin] = count
    smooth = gaussian_filter1d(signal, config.smoothing_sigma_bp, mode="constant")
    positions = [p for p, _ in sites]
    prefix = np.concatenate(([0], np.cumsum([c for _, c in sites], dtype=np.int64)))
    result = []
    pending = [(0, len(sites), ())]
    while pending:
        left, right, reasons = pending.pop()
        start, end = positions[left], positions[right - 1]
        width = end - start + 1
        cut = None
        reason = None
        if width >= config.min_peak_distance_bp and right - left >= 2:
            profile = smooth[start - origin:end - origin + 1]
            peaks, _ = find_peaks(np.pad(profile, (1, 1)), distance=config.min_peak_distance_bp)
            peaks = peaks - 1
            # Tiny ripples on the bridge between two supported modes must
            # not hide their intervening valley by becoming adjacent peaks.
            radius = max(1, config.min_peak_distance_bp // 2)
            peaks = [p for p in peaks if signal[
                max(start - origin, start - origin + int(p) - radius):
                min(end - origin + 1, start - origin + int(p) + radius + 1)
            ].sum() >= min_cluster_reads]
            options = []
            for p, q in zip(peaks, peaks[1:]):
                trough = profile[p:q + 1]
                minimum = float(trough.min())
                peak_floor = min(float(profile[p]), float(profile[q]))
                if peak_floor <= 0 or minimum > config.valley_ratio * peak_floor:
                    continue
                minima = np.flatnonzero(trough == minimum)
                coordinate = start + int(p + minima[len(minima) // 2])
                boundary = bisect_right(positions, coordinate, left, right)
                if not left < boundary < right:
                    continue
                if (prefix[boundary] - prefix[left] < min_cluster_reads
                        or prefix[right] - prefix[boundary] < min_cluster_reads):
                    continue
                options.append((minimum / peak_floor, abs(coordinate - (start + end) / 2), boundary))
            if options:
                cut = min(options)[2]
                reason = "valley"
        if cut is None and width > config.max_width_bp and right - left > 1:
            middle = (start + end) / 2
            central_left = start + width * 0.25
            central_right = end - width * 0.25
            candidates = []
            fallback = []
            for boundary in range(left + 1, right):
                midpoint = (positions[boundary - 1] + positions[boundary]) // 2
                density = float(smooth[midpoint - origin])
                score = (density, abs(midpoint - middle), boundary)
                fallback.append((abs(midpoint - middle), density, boundary))
                if central_left <= midpoint <= central_right:
                    candidates.append(score)
            cut = min(candidates)[2] if candidates else min(fallback)[2]
            reason = "max_width"
        if cut is None:
            result.append((sites[left:right], "+".join(sorted(set(reasons))) or "unchanged"))
        else:
            child_reasons = reasons + (reason,)
            pending.append((cut, right, child_reasons))
            pending.append((left, cut, child_reasons))
    return result

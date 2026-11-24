"""Versioned, packaged default model and length-aware batched inference."""
from __future__ import annotations
import hashlib
from importlib.resources import files
import json
from pathlib import Path

import numpy as np
import pysam
import torch

from .deeptssfilter import load_model, reverse_complement, score_sequences

DEFAULT_THRESHOLD = .9


def bundled_model() -> tuple[Path, dict]:
    root = files('scatter').joinpath('models')
    path = Path(str(root.joinpath('shift_hg38_1000bp_sigma10.pth')))
    metadata = json.loads(root.joinpath('shift_hg38_1000bp_sigma10.json').read_text())
    if hashlib.sha256(path.read_bytes()).hexdigest() != metadata['sha256']:
        raise ValueError('Packaged Shift model checksum mismatch; reinstall scatter-tss')
    return path, metadata


def resolve_device(requested: str) -> str:
    if requested == 'auto':
        return 'cuda' if torch.cuda.is_available() else 'cpu'
    if requested.startswith('cuda') and not torch.cuda.is_available():
        raise ValueError('CUDA requested but unavailable; choose --device cpu or auto')
    if requested != 'cpu' and not requested.startswith('cuda'):
        raise ValueError('Device must be cpu, cuda, cuda:N, or auto')
    return requested


def score_loci(rows, fasta_path, model_path, *, device='auto', batch_size=256, threads=4):
    if batch_size < 1 or threads < 1:
        raise ValueError('Batch size and threads must be positive')
    device = resolve_device(device)
    torch.set_num_threads(threads)
    model = load_model(model_path, device=device)
    length = int(model.sequence_length_bp)
    if length % 2 or length < 2:
        raise ValueError('Model requires a positive even sequence length')
    unique = list(dict.fromkeys((str(r['chrom']), int(r['summit_1']), str(r['strand'])) for r in rows))
    scores, statuses, valid_loci, sequences = {}, {}, [], []
    def flush():
        if sequences:
            values = score_sequences(model, sequences, device=device, batch_size=batch_size)
            for loc, value in zip(valid_loci, values):
                scores[loc] = float(value)
            valid_loci.clear(); sequences.clear()
    with pysam.FastaFile(str(fasta_path)) as fasta:
        sizes = dict(zip(fasta.references, fasta.lengths))
        for loc in unique:
            chrom, summit, strand = loc
            start, end = summit - 1 - length // 2, summit - 1 + length // 2
            status = 'valid'
            if strand not in ('+', '-'):
                raise ValueError(f'Invalid RNA strand: {strand}')
            if chrom not in sizes:
                status = 'missing_contig'
            elif start < 0 or end > sizes[chrom]:
                status = 'out_of_bounds'
            else:
                sequence = fasta.fetch(chrom, start, end).upper()
                if len(sequence) != length:
                    status = 'invalid_length'
                elif set(sequence) - set('ATCG'):
                    status = 'ambiguous_bases'
                else:
                    sequences.append(reverse_complement(sequence) if strand == '-' else sequence)
                    valid_loci.append(loc)
            statuses[loc] = status
            if status != 'valid': scores[loc] = None
            if len(sequences) >= batch_size: flush()
        flush()
    return scores, statuses, dict(sequence_length_bp=length, device=device, batch_size=batch_size,
        threads=threads, unique_loci=len(unique), valid_loci=sum(v=='valid' for v in statuses.values()),
        unscorable_loci=sum(v!='valid' for v in statuses.values()),
        model_sha256=hashlib.sha256(Path(model_path).read_bytes()).hexdigest(),
        sequence_rule='[summit_1-1-L/2,summit_1-1+L/2), reverse complement on RNA minus strand, ATCG channels',
        inference='one centered window; no shift averaging')

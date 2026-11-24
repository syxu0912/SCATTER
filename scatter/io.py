"""Small, shared input readers for plain and gzip-compressed text."""
import gzip
from pathlib import Path


def open_text(path):
    path = Path(path)
    return gzip.open(path, 'rt', encoding='utf-8') if path.suffix == '.gz' else path.open(encoding='utf-8')


def read_barcodes(path):
    cells, seen = [], set()
    with open_text(path) as handle:
        for line in handle:
            token = line.rstrip('\r\n').split('\t', 1)[0].strip()
            if not token or token.startswith('#'):
                continue
            if not cells and token.lower() in {'barcode', 'cell_barcode', 'cell_id', 'cb'}:
                continue
            if token in seen:
                raise ValueError(f'Duplicate cell barcode in {path}: {token}')
            cells.append(token); seen.add(token)
    if not cells:
        raise ValueError(f'No cell barcodes in {path}')
    return cells

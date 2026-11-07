#!/usr/bin/env python3
"""Reconstitute `results/trials.csv` from its GitHub-friendly CSV partitions."""
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parts = sorted((root/'results'/'trial_parts').glob('trials_*.csv'))
if not parts:
    raise SystemExit('No trials_*.csv parts found')
out = root/'results'/'trials.csv'
with out.open('w', encoding='utf-8', newline='') as target:
    header = None
    for part in parts:
        with part.open(encoding='utf-8', newline='') as source:
            first = source.readline()
            if header is None:
                header = first
                target.write(first)
            elif first != header:
                raise ValueError(f'Header mismatch in {part}')
            for line in source:
                target.write(line)
print(f'Rebuilt {out} from {len(parts)} parts ({out.stat().st_size} bytes)')
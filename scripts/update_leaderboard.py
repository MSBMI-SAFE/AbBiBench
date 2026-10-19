"""Validate the manuscript snapshot and refresh its README table (stdlib only)."""

import argparse
import csv
import json
import math
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
START = '<!-- manuscript-leaderboard:start -->'
END = '<!-- manuscript-leaderboard:end -->'


def render():
    metadata = json.loads((ROOT / 'leaderboard/provenance.json').read_text())
    with (ROOT / 'leaderboard/manuscript_spearman.csv').open(newline='') as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        rows = list(reader)
    models = metadata['rank_order']
    if len(models) != 23 or len(set(models)) != len(models):
        raise ValueError('Expected 23 unique model implementations.')
    if fields is None or fields[0] != 'assay' or set(fields[1:]) != set(models):
        raise ValueError('CSV model columns do not match the manuscript registry.')
    if len(rows) != 16 or len({row['assay'] for row in rows}) != 16:
        raise ValueError('Expected 16 unique assays.')
    for row in rows:
        if None in row or any(value is None for value in row.values()):
            raise ValueError('Malformed CSV row.')
        for model in models:
            value = float(row[model])
            if not math.isfinite(value) or not -1 <= value <= 1:
                raise ValueError(f'Invalid correlation for {model}: {value}')
    lines = [START,
             '| Rank¹ | Model Type | Model / score | ' + ' | '.join(row['assay'] for row in rows) + ' | Approx. mean ρ² ↑ |',
             '|' + '|'.join(['---'] * (len(rows) + 4)) + '|']
    for rank, model in enumerate(models, 1):
        label = {'AF3': 'AF3 (ipSAE)', 'Boltz-2': 'Boltz-2 (ipSAE)',
                 'ESMFold2': 'ESMFold2 (ipSAE)', 'Delta G': 'FoldX (−ΔΔG)',
                 'SASA': 'FreeSASA (−ΔSASA)', 'Antifold': 'AntiFold'}.get(model, model)
        values = [row[model] for row in rows]
        mean = sum(Decimal(value) for value in values) / len(rows)
        formatted = [f'{Decimal(value):.2f}' if Decimal(value) else '0.00' for value in values]
        lines.append('| ' + ' | '.join([str(rank), metadata['families'][model], label,
                                         *formatted, f'{mean:.2f}']) + ' |')
    lines.extend([END, '',
                  '¹ Rank follows the manuscript’s Figure 3b ordering. Rounded displayed means can tie.',
                  '² Approximate arithmetic mean across the 16 assays, recalculated from the two-decimal Figure 3a labels. '
                  'These means can differ from estimates calculated from unrounded results; no confidence intervals are reconstructed.'])
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Validate without modifying files.')
    args = parser.parse_args()
    path = ROOT / 'README.md'
    original = path.read_text()
    start = original.index(START)
    end = original.index('\n# Installation', start)
    updated = original[:start] + render() + '\n\n' + original[end:]
    if args.check:
        if updated != original:
            raise SystemExit('README leaderboard is stale; run scripts/update_leaderboard.py.')
        print('Validated 23 models × 16 assays and README synchronization.')
    else:
        path.write_text(updated)
        print('Updated README from the validated manuscript snapshot.')


if __name__ == '__main__':
    main()

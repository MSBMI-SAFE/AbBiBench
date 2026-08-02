#!/usr/bin/env python3
"""Replay frozen D4-PAR1/IF10 label-blind AbBiBench scores.

The frozen score files in this directory were generated before AbBiBench
labels were opened for evaluation. This script only copies the locked
candidate scores into the CSV format expected by AbBiBench.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


SCORE_COLUMN = "d4_par1_if10_strict_internal_guard_score"
PREFIX = {'aayl50': 'aayl50_LC', 'aayl52': 'aayl52_LC'}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", required=True)
    parser.add_argument("--score-root", default=None)
    parser.add_argument("--output-dir", default="./notebooks/scoring_outputs")
    args = parser.parse_args()

    score_root = Path(args.score_root) if args.score_root else Path(__file__).resolve().parent / "frozen_scores"
    score_path = score_root / f"{args.name}.frozen_scores.csv"
    if not score_path.is_file():
        raise FileNotFoundError(score_path)
    frame = pd.read_csv(score_path)
    out = pd.DataFrame({
        "row_id": frame["row_id"],
        "sequence_sha256": frame["sequence_sha256"],
        "log-likelihood": frame["log-likelihood"],
        "d4_par1_if10_policy": frame["d4_par1_if10_policy"],
    })
    prefix = PREFIX.get(args.name, args.name)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{prefix}_benchmarking_data_D4-PAR1_IF10_scores.csv"
    out.to_csv(out_path, index=False)
    print(out_path)


if __name__ == "__main__":
    main()

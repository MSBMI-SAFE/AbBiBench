# D4-PAR1/IF10 AbBiBench Submission Audit

- Model: `D4-PAR1_IF10`
- Model type: `Geometry/Physics-Aware Sequence Ranker`
- Official-style preview rank: **1**
- Official-style preview macro Spearman: **0.310716**
- Public leaderboard source: `/root/mutdock_ranker/outputs/d4_par1_top3_20260801/abbibench_official_leaderboard_audit/remote_public_leaderboard.csv`
- Label policy: frozen scores only; no `binding_score` column is copied into model inputs or score outputs.
- Entrypoint: `models/D4-PAR1_IF10/run.py --name <dataset>`
- Output directory: `notebooks/scoring_outputs/`

The packaged runner is intentionally portable: it reads the frozen label-blind score table
from `models/D4-PAR1_IF10/frozen_scores/` and writes the standard AbBiBench score CSV.

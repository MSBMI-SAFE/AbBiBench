# D4-PAR1 IF10 AbBiBench Submission Preview

This runner replays frozen, label-blind D4-PAR1/IF10 scores.
It does not read `binding_score`; labels are used only by the external evaluator after scores are frozen.
The default score source is `models/D4-PAR1_IF10/frozen_scores/`, so the command is portable inside the submission package.

Example:

```bash
python models/D4-PAR1_IF10/run.py --name 3gbn_h1
```

# D4-PAR1/IF10 Portable Runner Validation

Status: **portable_runner_validation_passed**

This validation replays every packaged frozen score file through
`models/D4-PAR1_IF10/run.py` and compares the generated CSV with the
prebuilt AbBiBench scoring output. It does not read AbBiBench labels.

## Summary

- Datasets checked: 19
- Total rows replayed: 216650
- Failures: 0

## Dataset Checks

| Dataset | Rows | Same Values | Forbidden Columns |
|---|---:|---:|---|
| 1mhp | 39 | True |  |
| 1mhp_LC | 69 | True |  |
| 1mlc | 1229 | True |  |
| 1mlc_LC | 912 | True |  |
| 1n8z | 419 | True |  |
| 2fjg | 2223 | True |  |
| 3gbn_h1 | 1887 | True |  |
| 3gbn_h9 | 1842 | True |  |
| 4d5_her2 | 2080 | True |  |
| 4fqi_h1 | 65094 | True |  |
| 4fqi_h3 | 65535 | True |  |
| 5a12_ang2 | 944 | True |  |
| 5a12_vegf | 29981 | True |  |
| aayl49 | 4312 | True |  |
| aayl49_ML | 8953 | True |  |
| aayl50 | 11473 | True |  |
| aayl51 | 4320 | True |  |
| aayl52 | 13324 | True |  |
| g6_LC | 2014 | True |  |

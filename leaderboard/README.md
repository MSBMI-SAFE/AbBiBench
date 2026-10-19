# Manuscript leaderboard snapshot

The main README leaderboard reflects the manuscript supplied on 15 October 2026, updated here on 19 October 2026. It covers **20 base methods / 23 implementations** across **16 assays**. Fixed-backbone variants of DiffAb, MEAN and dyMEAN count as additional implementations, not additional base methods.

## Provenance and precision

- `manuscript_spearman.csv`: 368 two-decimal Spearman correlations extracted from the numerical annotations in manuscript Figure 3a, stored as `figures/Figure2_complete_nature.pdf` in the LaTeX package.
- `provenance.json`: source-figure SHA-256, model families and rank order transcribed from Figure 3b.
- README means are approximate equal-assay arithmetic averages of rounded correlations. They are not the unrounded analysis results or new model evaluations. Rank follows Figure 3b; rounded averages can tie.
- The underlying new per-variant results and model runners are not included in this snapshot. Existing `notebooks/scoring_outputs/` files predate this update and should not be assumed to reproduce it.
- The supplied precision-at-k appendix omits ESMFold2, although the main figure includes it. No precision-at-k leaderboard or confidence intervals have been inferred from that inconsistent appendix.

## Score definitions

AF3, Boltz-2 and ESMFold2 use **ipSAE** as their primary structure-prediction score. The previous AF3/Boltz-2 entries used older confidence-proxy results and are superseded here. ipTM and interface pLDDT remain auxiliary metrics; none of these confidence scores is a calibrated affinity estimate. Boltz-2's small-molecule affinity output is not used as an antibody-affinity measurement.

ESMDance is the manuscript's biophysics-informed dynamics score. FoldX and FreeSASA use the negative mutation-induced changes in binding energy and epitope solvent accessibility. Other methods use their manuscript-defined sequence/structure scores; the leaderboard is not uniformly a comparison of log-likelihoods.

Assay correlations compare scores with experimental binding within each assay. Direct Kd and enrichment readouts are not pooled onto a single thermodynamic affinity scale. Light-chain components are merged into their biological assays (including `1mlc`/`1mlc_LC` and `2fjg`/`g6_LC`). AF3, Boltz-2 and ESMFold2 were evaluated on the manuscript's common-case population; this does not establish identical variant coverage for every other model.

Structure-conditioned methods use experimental complex templates where available and AF3-predicted templates for the five AAYL assays. The structure predictors' ipSAE entries are confidence scores from their predictions. The proposed comparison of experimental templates versus predicted templates/alternative poses is **not** a completed benchmark result and is not represented by these rows.

## Updating

```bash
python3 scripts/update_leaderboard.py
python3 scripts/update_leaderboard.py --check
```

When the original analysis outputs become available, replace this rounded snapshot with versioned per-assay estimates and coverage, update the provenance and ranking rule, and regenerate the table. Do not silently mix iPTM and ipSAE results, or infer missing model results from an older score file.

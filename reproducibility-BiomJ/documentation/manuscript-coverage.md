# Manuscript Coverage Check

Checked on 19 November 2026 against the author-selected
[arXiv:2606.16080v1](https://arxiv.org/abs/2606.16080v1), including its supplement.
The versioned PDF contains 31 pages: the main article occupies PDF pages 1-18;
the supplement occupies PDF pages 19-31 and restarts its printed numbering at 1.
The June 11 local drafts were not used as the final numbering reference.

Use [manuscript_map.csv](../manuscript_map.csv) for exact numbered locations,
generated filenames, archived references, commands and verification status.
Generated filenames such as `table2_1_...` are internal identifiers, not the
published table numbers. Blank paths/commands mean that no numerical reproduction
is supplied for that item; they are not instructions to run a missing analysis.

## Verified Simulation Results

- Main Table 4 and Supplement Tables S1-S8: all nine generated LaTeX tables have
  the same displayed data cells as their counterparts in the arXiv source.
  Whitespace and the unavailable-value markers `NA` and `--` were normalised for
  this comparison. Font sizes, captions, model labels and table scaling differ;
  byte-identical manuscript layout is not claimed.
- Main Figure 1 and Supplement Figures S1-S2: each generated PDF renders
  pixel-identically to the corresponding arXiv source figure using the same
  Poppler renderer at a maximum dimension of 1600 pixels. PDF bytes differ.
- Main Table 2: the five coefficient scenarios and fixed truths correspond to
  `simulation_parameters()` in `production/Bunya/ENZAMET/Data-Gen-ENZAMET.R`.
  The configuration grid and generated design CSV support Section 5.1; this
  design CSV is not the real-data Table 1.
- Main Table 3: prior specifications are traceable to the shared Gamma5 Stan
  source and the R driver's prior-data preparation. The longitudinal intercept
  prior is specified on the centered-response scale; the induced original-scale
  alpha prior SD is `log(2)/1.96`. Tables 2 and 3 are manually typeset, not emitted
  by the manuscript-output generator.
- Section 5.3 runtimes: archived records reproduce the printed mean, median and
  95th percentile after rounding. In minutes these are 27.2/21.4/67.2 for sAFT-JM,
  0.005/0.004/0.008 for LMM and 1.2/1.1/1.7 for Cox-JM. Use the main runtime
  summary's `OVERALL` row for sAFT-JM; use `total_model_time/60` grouped by `model`
  for LMM and Cox in the comparator per-simulation runtime CSV. Each uses 40,000
  recorded datasets. These are historical Bunya runtimes, not local-run targets.

This covers every numbered simulation table and figure in the specified version.
It complements the saved reaggregation checks in `validation-report.md`; those
checks compare against archived outputs, whereas this check compares against the
manuscript source. No model fitting or scientific-code changes were needed.

## Explicit Exceptions And Extra Outputs

- Main Table 1 and Section 4's real-data estimates cannot be numerically
  reproduced from the public mock data. Workflow 1 reproduces the fitting
  procedure with synthetic data and compares against mock references only. Its
  additional five BP-weight rows are not printed in the current Table 1.
- Section 4's PH diagnostic p-values (0.013 and 0.025; PDF page 6) also require
  restricted ENZAMET data. The supplied mock workflow does not reproduce these
  tests. Authors should disclose this exception explicitly in the audit response
  and agree any further verification with the editor; no agreement is assumed.
- Supplement Section S1 gives methodological guidance on model checks and
  information criteria, not reported numerical predictive checks or DIC results.
  Saved conditional LOO/WAIC scores do not establish new-patient prediction.
- `fig9_longitudinal_LL_coverage.pdf` and
  `fig11_alpha_prior_bias_coverage_heatmap.pdf` are additional outputs, not
  numbered figures in this version. Other exploratory outputs may also be
  generated. They are not additional manuscript requirements.

If the journal revision differs from this arXiv version, repeat the numbering
and coverage comparison before submission. Data/code permissions, independent
machine testing and archive delivery remain separate author sign-off items.

## Source Identity

The PDF and TeX source were downloaded from the versioned arXiv links. Local
review downloads and rendered pages are under the ignored
`output/manuscript-map-review/`, not required package inputs.

| Download | SHA-256 |
| --- | --- |
| [PDF](https://arxiv.org/pdf/2606.16080v1) | `f6c831d3e3554a9efa9013f56061cc28531506bb5c356a4be2eab39a8561fd04` |
| [TeX source](https://arxiv.org/src/2606.16080v1) | `bc4d1fe1511315ab9c5d411b34856c8cd1f1e95f5d223bbd0aa89ec44f2e15a9` |

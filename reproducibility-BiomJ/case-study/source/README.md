# Mock Case-Study Inputs

These synthetic records mimic the structure of the restricted ENZAMET data,
not the manuscript's real-data results.

Run the example from the package root using `run_checks.R` with
`task <- "case-study"`. See the root README for dependencies and output paths.
Do not use a separate case-study project.

- `mock-data.rds`: 17,623 assessments and 1,105 synthetic survival records.
- `fit-case-study-model.R`: local fitting/table procedure, seed 10235.
- Shared Stan model: `production/Stan/Bernstein-Polynomials-JM-Hist-alt-DM-Gamma5.stan`.
- Supplied references: `case-study/reference/`.
- Generated estimates, inputs, draws and diagnostics: `output/case-study/`.

The script's input/initial-value equivalence was checked against the saved
completed mock fit. Review sampling diagnostics before interpreting estimates.
Redistribution permissions and licensing require author confirmation.

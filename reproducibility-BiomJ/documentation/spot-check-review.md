# Reviewing Simulation Spot Checks

Use this guide after `task <- "records"` in `run_checks.R`. This task combines
existing local outputs; it does not fit models. Paths below are relative to the
package folder. Mock case-study results are separate: see README Workflow 1.
These local records are the endpoint of Workflow 2. Workflow 3's `"aggregate"`
and `"results"` tasks read archived production results/summaries, not these fits.

## 1. Confirm Which Runs Are Included

Open `output/local-simulations/editor-records/progress_counts.csv` for counts,
then `progress.csv` for individual configuration/replicate IDs and errors.
Records can include earlier selections: filter by your chosen `collection`,
`full_config_id` and `replicate_id`. A `completed` checkpoint means execution
finished, not that every model or diagnostic is acceptable.

## 2. Compare Matching Estimates

Open **`comparison_with_historical_production.csv`** in the same folder. The
local and archived results are already paired; no manual RDS extraction is needed.

| Column | Meaning |
|---|---|
| `collection`, `full_config_id`, `replicate_id`, `unique_seed` | Identity of the selected dataset and analysis |
| `variable` | Parameter being compared |
| `model` | Comparator model, when applicable; may be blank for main/sensitivity rows |
| `mean_fresh` | Local point estimate (posterior mean for joint-model parameters) |
| `mean_production` | Archived point estimate for the same configuration, replicate and model |
| `mean_difference` | Local minus archived estimate |

For example, from the package's R project:

```r
folder <- file.path("output", "local-simulations", "editor-records")
comparison <- read.csv(file.path(folder, "comparison_with_historical_production.csv"))
selected <- subset(comparison, collection == "main" & full_config_id == 1 &
                   replicate_id %in% 1:5 &
                   variable %in% c("beta_long[3]", "beta_surv[1]", "alpha"))
selected[c("replicate_id", "variable", "mean_fresh", "mean_production", "mean_difference")]
```

Change the selection to match your spot check. `beta_long[3]` is the longitudinal
treatment-by-time coefficient; `beta_surv[1]` is the survival treatment coefficient
(manuscript gamma); Stan `gamma[1:5]` are the Bernstein weights (manuscript theta).
Intercept summaries in these records are already on the original response scale.

Compare **the same replicate**, not a five-replicate average against the published
1,000-replicate average or each estimate against the generating truth. Exact MCMC
equality across machines is not required; differences need interpretation alongside
sampling diagnostics and Monte Carlo uncertainty, not a universal decimal tolerance.

## 3. Review Reliability Before Judging Agreement

| File in `editor-records/` | What to inspect |
|---|---|
| `stan_chain_diagnostics.csv` | Divergences, maximum-treedepth hits and E-BFMI for each Stan chain |
| `replicate_parameter_records.csv` | Estimates, interval endpoints (`q2.5`, `q97.5`), R-hat, ESS and available MCSEs |
| `comparator_model_status.csv` | Status/error for each comparator model; inspect even when its dataset checkpoint completed |

Flags such as divergences, persistent treedepth hits or poor mixing require review
even when point estimates look close. Do not silently exclude flagged replicates.
`uncertainty_sd` describes posterior uncertainty or an LMM standard error when
available; it is **not** the Monte Carlo standard error `mcse_mean`.
Historical results and some comparator records lack MCSEs, so a formal
combined-MCSE agreement calculation is not always possible. Interval overlap alone
does not establish convergence or numerical equivalence.

For details, open the matching folder under
`output/local-simulations/fits/<collection>/configNNN/repNNN/`: it contains the
production comparison, saved inputs, initial values, diagnostics and warnings.
`progress.csv` also records the archived reference filename and within-batch position.

## Optional Files And What Not To Infer

- `replicate_runtime_records.csv`: local timings, not an agreement criterion.
- `reduced_run_summary.csv`: descriptive means/SDs across available local runs;
  a small spot check cannot recover precise full-study bias, RMSE or coverage.
- `saft_information_criteria.csv`: conditional joint LOO/WAIC, not required for
  checking parameter reproducibility. Retain score warnings; unreliable scores
  do not change the fitted posterior, but should not be used to certify agreement.
- `local-audit-records.rds`: the combined records in R format; CSVs suffice for review.

Record the selected IDs, environment, estimate differences and diagnostic concerns.
Author-run examples are in [local-fit-notes.md](local-fit-notes.md); software and
cross-platform qualifications are in [software-environment.md](software-environment.md)
and the [technical guide](technical-guide.md#what-agreement-means).

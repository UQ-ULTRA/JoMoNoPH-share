# Optional Command-Line Workflow

The recommended RStudio route is `run_checks.R` in the package root; see the
root `README.md`. No terminal commands are needed for that route. The commands
below are an optional equivalent interface to the same fitting functions.

The audit workflow is to generate a small, clearly identified reference set on
the authors' machine, while allowing the editor to choose other configurations
and original replicate IDs. It is **not necessary to refit every configuration**
locally. Full saved production results remain available for additional checks
and regeneration of the published simulation summaries.

No unattended all-settings launcher is provided. Nothing starts automatically
when the project is opened.

## Choose Settings Before Fitting

From the parent `JoMoNoPH-share` directory:

```sh
# List all candidates; planning does not fit anything:
Rscript reproducibility-BiomJ/main.R --mode=plan --collection=all --n-reps=5
```

`output/simulation-plans/configuration_candidates.csv` lists all 40 original
configurations, including baseline distribution, coefficient scenario,
censoring mechanism/rate and follow-up. The collections are:

| Collection | Available configurations | Model/prior |
|---|---|---|
| main | IDs 1 to 40 | Production sAFT-JM; alpha SD log(2)/1.96 |
| alpha125 | IDs 1 to 10 | Log-logistic sAFT-JM; alpha SD log(1.25)/1.96 |
| alpha150 | IDs 1 to 10 | Log-logistic sAFT-JM; alpha SD log(1.5)/1.96 |
| alpha300 | IDs 1 to 10 | Log-logistic sAFT-JM; alpha SD log(3)/1.96 |
| comparator | IDs 1 to 40 | Cox-, Weibull-, log-logistic JM and LMM |

`--collection=alpha` selects all three extra alpha-prior settings for the chosen
log-logistic configuration IDs. Configuration IDs are never renumbered after
selection. Choose settings for their scientific coverage, not because their
previous estimates or diagnostics happened to look favourable.

## Run Five Replicates Of Selected Settings

These are syntax examples, not a prescribed set of audit settings:

```sh
# First inspect five replicates of main configuration 1:
Rscript reproducibility-BiomJ/main.R --mode=plan --collection=main --config-ids=1 --n-reps=5
# Optional data/initial-value check without MCMC:
Rscript reproducibility-BiomJ/main.R --mode=prepare --collection=main --config-ids=1 --n-reps=5 --no-pch=true
# Run those five replicates:
Rscript reproducibility-BiomJ/main.R --mode=spot-check --collection=main --config-ids=1 --n-reps=5 --no-pch=true
# Same selected setting under another alpha prior or comparator framework:
Rscript reproducibility-BiomJ/main.R --mode=spot-check --collection=alpha150 --config-ids=1 --n-reps=5 --no-pch=true
Rscript reproducibility-BiomJ/main.R --mode=comparator --config-ids=1 --n-reps=5
```

Every fitting/preparation command requires an explicit `--config-ids`. Multiple
IDs can be comma-separated. Five replicates means original replicate IDs 1 to 5
**for each selected configuration**, not five datasets overall. Sampling settings,
initial-value procedures, priors and seed rules remain those of production.

For a different replication count or a random spot check chosen by the editor:

```sh
Rscript reproducibility-BiomJ/main.R --mode=plan --collection=main --config-ids=1,2 --n-reps=10
Rscript reproducibility-BiomJ/main.R --mode=spot-check --collection=main --config-ids=7 --replicates=42 --no-pch=true
```

Use either `--n-reps` or `--replicates`, not both. Replicate IDs up to 1,000 have
historical production references. Explicitly selecting all configurations remains
possible with `--config-ids=all`, but is not required for the agreed spot-check
strategy. Even five sAFT fits can take hours; plan a small check first.

## Saved Records And Comparison

For file-by-file review instructions, see the [spot-check review guide](spot-check-review.md).

Each run saves its RDS result and provenance under
`output/local-simulations/fits/<collection>/configNNN/repNNN/`. This supports any
selected replicate count; actual IDs/counts are recorded explicitly in the
plan and progress table. Existing compatible completed
checkpoints are skipped. Errors are recorded and stop the selected run for review.

```sh
Rscript reproducibility-BiomJ/main.R --mode=local-status
```

Editor-facing CSV files and a combined RDS are in
`output/local-simulations/editor-records/`. They contain per-replicate estimates,
intervals, uncertainty, diagnostics, runtimes, model failures, conditional joint
LOO/WAIC where available, and comparisons with the matching production records.
The reduced-run summary is descriptive; five replicates cannot recreate precise
1,000-replicate bias, RMSE or coverage estimates.

For comparison, use the same configuration ID, replicate ID and seed, not the
mean over all 1,000 production replicates. A different computer may produce
different MCMC draws; inspect diagnostics and Monte Carlo uncertainty rather
than demanding identical rounded estimates. Historical summaries do not always
contain exact MCSEs, so formal combined-MCSE agreement cannot always be assessed.

Exact Stan data/initializations, chain IDs, selected draws, raw CSVs, warnings,
build and software/hardware metadata are retained. Comparator generated data,
preliminary-model inputs and seeds are also saved. This follows the categories
in the [Stan reproducibility guide](https://mc-stan.org/docs/reference-manual/reproducibility.html),
without claiming bitwise equivalence across different environments.

The mock case-study example and full-results manuscript regeneration workflow
are unchanged. All new output stays inside `reproducibility-BiomJ`; previous
figures and original production source/data are not overwritten.

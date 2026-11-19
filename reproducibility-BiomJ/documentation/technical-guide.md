# Technical Reference

Start with the package-root `README.md` and `run_checks.R` for the recommended
RStudio workflow. This document retains the detailed implementation, provenance
and optional command-line instructions. All commands below assume the package
root as the working directory, not this documentation folder.

Manuscript: **Bayesian joint modelling using semiparametric accelerated failure
time approaches**, bimj.202500502.

This package separates three tasks. They are not interchangeable:

1. Recreate the manuscript's simulation tables and figures from the supplied
   **full, post-rescue simulation results** (1,000 datasets per configuration).
2. Fit **five replicates per selected configuration** locally as computational
   spot checks. The author/editor chooses a small subset; fitting all settings
   is not required. Five
   replicates do not reproduce 1,000-replicate bias, RMSE or coverage estimates.
3. Run the case-study code on supplied **synthetic mock data**, not ENZAMET data.
   The mock results are not expected to equal the manuscript's real-data results.

## Quick Start

Unzip the entire folder. No RStudio project, SLURM, shell script or parent
`JoMoNoPH-share` checkout is required. Run the following commands from a terminal
in this folder; alternatively, give the full path to `main.R`. On Windows use an
R installation with `Rscript` on PATH, or its full executable path.

```sh
Rscript main.R --mode=check
Rscript main.R --mode=results
Rscript main.R --mode=plan --collection=all --n-reps=5
# A small selected check; fitting requires explicit configuration IDs:
Rscript main.R --mode=prepare --collection=main --config-ids=1 --n-reps=5
Rscript main.R --mode=spot-check --collection=main --config-ids=1 --n-reps=5 --no-pch=true
Rscript main.R --mode=case-study
```

`check` inventories the five saved collections, mock inputs, dependencies and
CmdStan/basic toolchain. Missing components and failed stages are reported without
stopping the other checks. The final PASS/WARNING table is returned and saved as
`output/check/check-summary.csv` where writable, with one concluding warning if
anything remains unresolved. No model is compiled; passing does not guarantee
successful model linking or convergence. `results`
regenerates manuscript-facing CSV, LaTeX and PDF/PNG files from the supplied
summary tables. `prepare` regenerates the selected datasets, fits the preliminary
models and saves all Stan inputs/initial values, **without sampling**. It checks
the response-centering constant and exponential-AFT starting estimate against
corresponding historical replicates. `spot-check` runs the requested full
four-chain Stan fits sequentially, using production MCMC settings.

All newly computed outputs go under `output/`. Supplied data and reference
results are not overwritten. Completed local checkpoints are skipped by default.
Preparation ignores the `no_pch` build option because it never compiles Stan;
actual Stan fit checkpoints still check it. Other identity differences or missing
files stop resumption and are named in the error. Inspect them before deliberately
setting `resume <- FALSE` in `run_checks.R` (CLI: `--resume=false`).
Use `--resume=false` only when intentionally repeating and overwriting a local
check. See [command-line.md](command-line.md) for selected-run commands and exports.

## Software And Installation

Use a current R installation with the following packages installed before
running the commands. This script does not install packages or modify an
existing CmdStan installation automatically.

```r
install.packages(c("tidyverse", "survival", "MASS", "nlme", "loo",
                  "posterior", "patchwork", "JMbayes2"))
install.packages("cmdstanr",
                 repos = c("https://stan-dev.r-universe.dev",
                           getOption("repos")))
cmdstanr::check_cmdstan_toolchain()
# Install CmdStan only if needed. This is an explicit user action:
# cmdstanr::install_cmdstan()
```

The Stan fits require CmdStan and a supported C++ toolchain (RTools on Windows;
Xcode command-line tools on macOS; a compatible C++ compiler and make on Linux).
`--mode=results` needs R packages but does not compile or sample a model.

CmdStanR's configured installation is used, rather than the hard-coded Bunya
path in the historical source. To choose an installed version explicitly:

```sh
Rscript main.R --mode=spot-check --config-ids=1 --cmdstan=/path/to/cmdstan-2.37.0
```

Historical production used CmdStan 2.37.0. The supplied mock reference README
records its own, different historical environment. New-run software versions,
`sessionInfo()`, installed-package versions, compiler information, CmdStan
library versions and hardware are written below each output's `environment/`
folder. The locally tested environment is recorded in `documentation/validation/`.
Completed local fit timings and diagnostic caveats are in
`documentation/local-fit-notes.md`.
Do not assume this is the original cluster's complete environment specification.

If a compiler update makes an existing CmdStan precompiled header incompatible,
this package can compile its model without that header, without rebuilding the
user's CmdStan installation:

```sh
Rscript main.R --mode=spot-check --config-ids=1 --no-pch=true
```

## Simulation Collections

| Collection | Configurations | Replicates/config | Files | Contents |
|---|---:|---:|---:|---|
| `main` | 40 | 1,000 | 10,000 | Proposed sAFT-JM; four replicates/batch |
| `alpha125` | 10 log-logistic | 1,000 | 2,500 | Alpha prior SD `log(1.25)/1.96` |
| `alpha150` | 10 log-logistic | 1,000 | 2,500 | Alpha prior SD `log(1.5)/1.96` |
| `alpha300` | 10 log-logistic | 1,000 | 2,500 | Alpha prior SD `log(3)/1.96` |
| `comparator` | 40 | 1,000 | 4,000 | JMbayes2 Cox/Weibull/log-logistic JMs and LMM; ten datasets/batch |

The default `log(2)/1.96` sensitivity comparator comes from the main collection,
not a sixth collection. Alpha sensitivity is restricted to the **ten original
log-logistic configurations**; it was not run for all 40 configurations.

The 40 configurations cross four event-time distributions, five coefficient
scenarios and two censoring mechanisms. Each dataset has 1,100 patients and
maximum follow-up 120 months. The baseline distributions are log-logistic
(shape 1.2, scale 23) and Weibull (shapes 0.9, 1.3, 2.1; scale 38). `check`
writes the exact configuration grid, including censoring rates, to
`output/check/configurations.csv` using the production grid expressions.

The candidate grid is in `output/simulation-plans/configuration_candidates.csv`.
The default reduced run uses original replicate IDs 1 to 5 for each selected
configuration; IDs are never renumbered after filtering. The original seed
is `1e6 * full_config_id + replicate_id`. Production used `furrr` workers with
R's `L'Ecuyer-CMRG` RNG; retaining that RNG kind is necessary when regenerating
the datasets locally. Stan uses the same numerical seed with chain IDs 1 to 4.
Current completion is recorded in `output/local-simulations/editor-records/progress.csv`.

### Model Settings Preserved

- Observed-mean response centering, with intercept summaries returned to the
  original outcome scale. No true simulation coefficient is passed to the fitter.
- Five Bernstein baseline weights, each positive with Normal(0, 5 squared) prior;
  no shape roughness penalty.
- Zero-centered Normal priors with SD 10 for longitudinal and survival regression
  coefficients (the longitudinal intercept prior is on the centered-response scale).
- Raw-scale association prior SD `log(2)/1.96`, implemented through the unchanged
  scaled `alpha_tilde` parameterization. Sensitivity overrides alter only this SD.
- Half-Cauchy(0, 5) scale priors and LKJ(2) random-effect correlation prior, as
  coded in the copied Stan file; random-effect SDs have the coded lower bound .001.
- Initial values from `nlme::lme()` and exponential `survival::survreg()`, with
  production jitter and non-finite survival-start fallback. Baseline weights are
  initialized near .1; association near zero; initial random-effect correlation .2.
- Four chains, 1,000 warmup and 1,000 retained iterations per chain,
  `max_treedepth=10`, one thread/chain, default Stan `adapt_delta` (0.8).

### Spot-Check Options

```sh
# Inspect candidates and planned counts without sampling:
Rscript main.R --mode=plan --collection=all --n-reps=5
# Five replicates for configuration 1; still four chains, two at a time:
Rscript main.R --mode=spot-check --collection=main --config-ids=1 --n-reps=5 --parallel-chains=2 --no-pch=true
# Other selected checks (configuration 1 is an example, not a prescribed set):
Rscript main.R --mode=spot-check --collection=alpha125 --config-ids=1 --n-reps=5 --no-pch=true
Rscript main.R --mode=spot-check --collection=alpha150 --config-ids=1 --n-reps=5 --no-pch=true
Rscript main.R --mode=spot-check --collection=alpha300 --config-ids=1 --n-reps=5 --no-pch=true
Rscript main.R --mode=comparator --config-ids=1 --n-reps=5
# Increase replication or select exact original replicate IDs:
Rscript main.R --mode=plan --collection=main --config-ids=1,2 --n-reps=10
Rscript main.R --mode=plan --collection=main --config-ids=1 --replicates=7,42
```

`--collection=alpha` selects all three additional prior widths on the ten
log-logistic configurations; `--collection=all` selects all five collections.
Alpha collections retain full-grid configuration IDs 1 to 10. The comparator uses all three original
JMbayes2 survival models, 3,500 iterations, 500 burn-in, thinning 1 and four
chains. Each selected comparator dataset therefore includes three JM fits plus an LMM.
Fitting/preparation commands require `--config-ids`; planning may list all
candidates without initiating MCMC. No unattended all-settings launcher is supplied.
Local-only fitting definitions are under `production/`. SLURM dispatch,
workspace clearing and fixed workstation paths were removed, and one shared
generator replaces five scenario copies. Use `run_checks.R` or `main.R` rather
than sourcing individual modules. Numerical equivalence checks and source
provenance are described in `documentation/code-equivalence.md`.

Each new Stan check saves exact data and all four initialization lists (including
random-effect starts), the call configuration, raw Stan CSVs, coefficient draws,
diagnostics, mean and quantile MCSEs, build metadata, and original-scale summaries.
The correct historical batch/replicate is read from `archived-production-results/` for every
collection.
Each new run saves its reference provenance. The regenerated-data check must pass
before comparing posterior estimates. Historical summaries do not contain exact
MCSEs, so comparisons with those summaries do not claim an exact MCSE-based pass.

## Rebuild Summaries From All Saved Batch Results

For a full check of the aggregation stage, rather than simply reading the
supplied summary CSVs, use the following chain. Inputs are the individual saved
production results in `archived-production-results/output_*/`, not local spot
checks or previously aggregated summary CSVs:

```sh
Rscript main.R --mode=aggregate --collection=main
Rscript main.R --mode=aggregate --collection=alpha125
Rscript main.R --mode=aggregate --collection=alpha150
Rscript main.R --mode=aggregate --collection=alpha300
Rscript main.R --mode=aggregate --collection=comparator
Rscript main.R --mode=results --summaries=regenerated
```

`--collection=all` runs all five aggregation commands. This reads tens of GB and
may take many minutes or longer, depending on disk speed and R/package versions.
Wait for aggregation to finish before running `results --summaries=regenerated`.
That final command reads `output/aggregated/` and writes `output/manuscript/`;
no further aggregation follows it. The shorter `results --summaries=saved` route
instead reads the supplied summaries under `archived-production-results/auxiliary/`.
It does not refit models. The archive retains original model failures and
diagnostics; no problematic replicate is silently excluded by the audit adapter.
The original aggregation code determines which estimates are available.

Only each collection's **top-level merged result directory** is aggregated.
Never point aggregation recursively at all of `archived-production-results/`: auxiliary rescue
folders contain duplicate historical batches.

The manuscript-output generator writes descriptive filenames rather than relying
on figure numbers that may change during revision. See `manuscript_map.csv`.
The `fig_slide_*` outputs and other additional diagnostics are auxiliary outputs,
not claims that every exploratory graph is included in the article. Supplied
reference figures remain in `archived-production-results/auxiliary/manuscript-outputs/`.

## Case Study: Synthetic Data Only

The mock data and fitting expressions originate from Andrew's supplied
`audit-mock-case-study 2.zip`. The R script now uses package-relative paths and
does not clear the workspace. Both simulation and case-study fitting use the
single unchanged Stan model in `production/Stan/`. `case-study/reference/`
contains the supplied CSV and three screenshots.

`--mode=case-study` runs those same expressions without clearing the workspace or
changing the working directory. The adapter redirects reads/writes and saves
portable numeric results and Stan CSVs instead of an R6 object pointing to
temporary files. It preserves the example's seed 10235, centering, priors,
initialization draw order and four-chain MCMC settings. The numerical model is
not changed to make its output agree with a reference table.

The supplied example removes missing HRQoL measurements when fitting and adds
1e-6 months to survival times to allow zero-time events. These are retained mock
example conventions, not newly imposed changes to the production generators.
Data dictionary and scope are in `documentation/data_dictionary.md`.

## Runtime And Parallelization

Production ran on The University of Queensland's Bunya HPC, using EPYC3/EPYC4
CPU nodes. Proposed-model batches requested eight cores and 8 GB RAM for four
replicates, generally two concurrent four-chain fits. Historical runtime records
report a mean of **27.16 minutes**, median **21.38 minutes**, and maximum
**153.57 minutes** per main-run replicate (40,000 replicates). This recorded
quantity comprises LMM fit, model loading/compilation and Stan sampling; it does
not include all data-generation, LOO/WAIC and file-saving overhead. Do not treat
the sum of these times as elapsed calendar time for the parallel cluster study.

Allow several hours for five local fits plus compilation; difficult warmup can
take considerably longer. Running all configurations would take many days,
but is not required for the selected-setting audit. A single configuration
first is recommended. The mock case
study is another full fit, not an instantaneous table-rendering operation.

The local runner uses one dataset at a time and up to four cores for its four
chains. `--parallel-chains=1` or `2` reduces simultaneous CPU/RAM demand while
retaining all four chains; fewer simultaneous chains generally increases elapsed
time. No automatic cluster or multisession workers are launched for Stan checks.
The optional JMbayes2 checks use that package's own parallel-chain mechanism.

## What Agreement Means

Saved-result aggregation should recover the same numerical tables, subject to
rounding and software-dependent rendering. PDF creation timestamps and plotting
package changes can change file hashes without changing the underlying results.
Check CSV values rather than insisting that PDF binaries have identical hashes.

Fresh MCMC results need not be bitwise identical across machines. Stan's
[reproducibility guide](https://mc-stan.org/docs/reference-manual/reproducibility.html)
requires matching versions, interface/language, included libraries, OS, hardware,
compiler/build flags and libraries, data, seed, chain IDs and initial values for
exact reproducibility. This package records what is obtainable for new runs;
unrecorded historical motherboard/compiler details cannot be reconstructed.

Inspect divergences, treedepth hits, R-hat, ESS and MCSE before interpreting
differences. Small rounded-table differences can be ordinary Monte Carlo
variation; large differences or poor diagnostics need investigation. For two
independent good fits, differences in posterior means can be assessed against
the square root of the sum of their squared mean MCSEs. Interval endpoints have
their own MCSEs. Neither closeness of means nor interval overlap proves convergence.

LOO/WAIC in the proposed-model outputs use the **patient-level conditional joint
log likelihood**, conditional on fitted patient random effects. They are not a
strict new-patient cross-validation analysis. A warning replacing missing `r_eff`
values by 1 does not affect fitting, but does warrant checking LOO diagnostics;
it is not suppressed or treated as an automatic pass here.

## Folder Contents And Provenance

| Path | Purpose |
|---|---|
| `run_checks.R` | Recommended RStudio entry script with editable settings |
| `main.R`, `R/` | Optional CLI entry point and shared portability/recording adapters |
| `production/` | Verified local production definitions and unchanged Stan model |
| `analysis/` | Copies of aggregation/presentation scripts with path/workspace edits only |
| `output/simulation-plans/` | Candidate configuration grid and explicit selected run plans |
| `documentation/command-line.md` | Optional selected-setting commands and CSV/RDS exports |
| `case-study/source/` | Path-adapted mock-example script and unchanged data from the uploaded ZIP |
| `case-study/reference/` | Author-supplied mock reference table/screenshots |
| `archived-production-results/` | Moved, consolidated post-rescue collections and historical auxiliary material |
| `documentation/` | Checklists, data dictionary, change log, validation evidence |
| `output/` | Newly generated results; separate from all supplied references |
| `files_manifest.csv` | Optional file-size/MD5 inventory; generated with `--mode=manifest` |

The source is the authors' existing research implementation; no change in model,
likelihood, priors, censoring calibration or production estimands was intended.
See `documentation/changes.md` for the exact adapter scope. Data provenance and
permissions still require author confirmation before redistribution. No actual
ENZAMET participant data are included.

## Before Sending To The Editor

This is a prepared audit package, not an author certification that every
checklist item has already been completed. Consult
`documentation/essential-checklist-response.md` and the validation report.
Have a coauthor test the instructions on a second machine. Confirm that the
manuscript-map labels correspond to the final submitted PDF and confirm sharing
permissions. Do not describe mock case-study estimates as recovered real-data
estimates.

The full saved archive is approximately 37 GB before auxiliary material. It is
not suitable for email. Either include it in the ZIP if the delivery mechanism
permits, or host it in an editor-accessible repository (for example Zenodo, OSF
or Figshare) and put the actual download URL and checksum in this README before
submission. A OneDrive link requiring an institutional account is not suitable.
No public upload or email submission has been performed by this script.

To make a complete inventory before archiving:

```sh
Rscript main.R --mode=manifest
```

Create the ZIP only after validation and the final author checks. Do not send
unnecessary compiled executables, failed compilation logs, `.DS_Store`, RStudio
caches, or large temporary rerun CSVs unless they are needed for the agreed audit.
Retain numeric local reference outputs and their environment information.

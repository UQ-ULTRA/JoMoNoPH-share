# Biometrical Journal Reproducibility Package

**Bayesian joint modelling using semiparametric accelerated failure time
approaches**, manuscript bimj.202500502.

## Get The Files And Reserve Disk Space

**Recommended: clone the whole code repository, then work only inside
`reproducibility-BiomJ/`.** First create a local parent folder using Finder or
File Explorer, for example:

| System | Example parent folder (replace the username) |
|---|---|
| Windows | `C:/Users/YourName/Documents/BiomJ-audit/` |
| macOS | `/Users/yourname/Documents/BiomJ-audit/` |

Open a terminal in that parent folder and run:

```sh
git clone --depth 1 https://github.com/UQ-ULTRA/JoMoNoPH-share.git
cd JoMoNoPH-share
```

Git creates **`JoMoNoPH-share/` inside the folder where the clone command runs**.
`--depth 1` includes all current project folders but limits downloaded history;
see the [Git documentation](https://git-scm.com/docs/git-clone). On Windows, for
example, open this project in RStudio:

```text
C:/Users/YourName/Documents/BiomJ-audit/JoMoNoPH-share/reproducibility-BiomJ/reproducibility-BiomJ.Rproj
```

On macOS the corresponding location begins
`/Users/yourname/Documents/BiomJ-audit/JoMoNoPH-share/`. Open the **nested
`reproducibility-BiomJ.Rproj`**, not the parent research project. A supplied code
ZIP can instead be extracted into your chosen parent folder; the same nested
project is used.

The current tracked/Git-eligible working files total approximately **7.2 MB**;
allow additional space for Git metadata. This is a small code checkout, not the
large results archive. The mock dataset and mock reference outputs are included;
ignored simulation results, generated outputs and local `.git` caches are not.

Download **`archived-production-results`** using the separate results link supplied
with the audit submission. Extract/place the complete folder, including all five
collections and `auxiliary/`, directly inside the cloned audit folder:

```text
your-chosen-parent/
  JoMoNoPH-share/
    reproducibility-BiomJ/
      run_checks.R
      reproducibility-BiomJ.Rproj
      archived-production-results/
        output_.../             (five production collections)
        auxiliary/              (summaries and reference tables/figures)
```

The extracted archive contains approximately **40.6 GB (37.8 GiB)**. Allow at
least **100 GB free disk space before downloading** for a typical selected audit,
including room for both a downloaded ZIP and its extracted files, R/CmdStan and
generated outputs. More may be needed for many additional fits/raw Stan CSVs.
Keep the archive locally available, not cloud-only placeholders. A mock-case-study
check alone does not need this archive; in that case skip the archive inventory
task `"check"` and go directly to Workflow 1 after installing dependencies.

## Start Here

1. After obtaining the files above, open `reproducibility-BiomJ.Rproj`.
2. Open **`run_checks.R`** and edit its task/settings.
3. Run sections 1-4 in order, or click **Source**. Repeat after changing tasks.

Edit settings in the script, not just the Console. Opening the project starts
nothing automatically. Running the R workflows requires no terminal, SLURM or parent repository;
without RStudio, work in this folder and use `source("run_checks.R")`.

Install missing dependencies using [Software Setup](#software-setup), then run
**`task <- "check"` once** to inventory files/packages and record the environment
in `output/check/`. This setup check does not fit models or certify agreement.

## Three Main Workflows

These are **independent workflows, not three mandatory stages**. All use `run_checks.R`.

| Purpose | Task sequence | What to compare |
|---|---|---|
| [1. Mock case study](#workflow-1-mock-case-study) | `"case-study"` | Fresh mock-data estimates against supplied mock references |
| [2. Simulation spot checks](#workflow-2-selected-simulation-spot-checks) | `"plan"` -> `"fit"` -> `"records"` | Selected original replicates against saved production results |
| [3. Manuscript tables/figures](#workflow-3-recreate-manuscript-tables-and-figures) | Route A: `"results"`; Route B: `"aggregate"` -> `"results"` | Archived production results, not fresh spot checks; no MCMC |

## Workflow 1: Mock Case Study

```r
task <- "case-study"
```
Source the script to fit the supplied synthetic data. Inspect
`output/case-study/mock_table1_update.csv` and the saved diagnostics. Compare
with `case-study/reference/`, **not** the restricted real-data manuscript table.

No other task is needed. Simulation selection settings do not affect this workflow.

## Workflow 2: Selected Simulation Spot Checks

**`"plan"` -> `"fit"` -> `"records"` -> inspect comparisons and diagnostics.**

Example settings in `run_checks.R`:
```r
task <- "plan"
collection <- "main"
config_ids <- 1L          # Example, not a prescribed audit selection.
replicates <- 1:5         # Five datasets PER selected configuration.
```

1. Run `"plan"` to inspect the workload and candidate configuration table.
2. Change only `task <- "fit"` and source the script to fit those datasets.
3. Change only `task <- "records"` and source it again to combine estimates,
   diagnostics, runtimes and comparisons with matching production replicates.
4. Review `output/local-simulations/editor-records/` and detailed files under
   `output/local-simulations/fits/`.

**What to open and compare:** follow the [spot-check review guide](documentation/spot-check-review.md)
for the file review order, comparison columns, example R commands and diagnostics.

**Optional:** `task <- "prepare"` inspects data/preliminary fits/initial values
without MCMC. `"fit"` regenerates its inputs; it does not load preparation files.

Keep the same selection between tasks; inspect a new plan if it changes.
All 40 main configurations with five replicates request **200 fits** without
confirmation. Five replicates cannot recover precise 1,000-replicate performance estimates.

### Available Settings

| `collection` | Configuration IDs | Model / association prior SD |
|---|---|---|
| `"main"` | 1-40 | `log(2)/1.96` |
| `"alpha125"` | 1-10, log-logistic | `log(1.25)/1.96` |
| `"alpha150"` | 1-10, log-logistic | `log(1.5)/1.96` |
| `"alpha300"` | 1-10, log-logistic | `log(3)/1.96` |
| `"comparator"` | 1-40 | JMbayes2 Cox/Weibull/log-logistic JMs and LMM |

Archives contain 1,000 datasets per configuration. `"alpha"` selects all three
sensitivity collections; `"all"` selects all five collections. View `configurations`
or `output/simulation-plans/` for the full grid. You need not refit every setting.

For example, `config_ids <- c(1, 2)` selects two configurations;
`replicates <- 56:60` selects five original replicate IDs. The seed rule is
`1000000 * config_id + replicate_id` (IDs 1-1,000 have archived references).
`config_ids <- NULL` is planning-only; `config_ids <- "all"` selects all for fitting.

`parallel_chains` controls how many of the four chains run simultaneously.
Datasets run sequentially. `resume <- TRUE` skips compatible completed
simulation checkpoints. Set `resume <- FALSE` to rerun supplied completed local
checks; this overwrites local outputs only, not archived production results.

## Workflow 3: Recreate Manuscript Tables And Figures

This workflow uses the **old production simulations in
`archived-production-results/`**, not freshly fitted spot checks. Choose either
Route A or Route B below. Both recreate the same manuscript outputs without MCMC.

### Route A: Use The Supplied Production Summaries

```text
archived-production-results/auxiliary/.../summary_tables_*/
    -> results (saved)
    -> output/manuscript/
```

Set these values in `run_checks.R` and source it once:
```r
task <- "results"
summaries <- "saved"
```

This reads the already-calculated production summary CSVs and creates manuscript
tables, figures and derived CSVs. **This route is complete; no aggregation is needed.**

### Route B: Recalculate The Production Summaries First

**Route B is an alternative to Route A, not a required additional step.** Choose
it only if you also want to check the calculations that turn the old production
replicate results into summary tables; otherwise Route A is sufficient.

```text
archived-production-results/output_*/ (per-replicate results in batch RDS files)
    -> aggregate
    -> output/aggregated/ (recalculated summaries and diagnostic plots)
    -> results (regenerated)
    -> output/manuscript/
```

**Step 1.** Set these values and source the script:
```r
task <- "aggregate"
collection <- "all"
```

This rereads the saved batch `.rds` files for all five production collections
and recalculates performance, diagnostic and runtime summaries. It does **not**
summarise the old summary CSVs again or refit models. Allow **tens of minutes or
longer** to read/process the full archive (tested machine: Apple M2 Pro, 32 GB RAM;
an exact aggregation elapsed time was not recorded). **Wait for this task to finish.**

Aggregation is sequential in **one R process, with no parallel workers**.
`parallel_chains <- 4L` controls model-fitting chains only and does not speed up
aggregation; underlying libraries may use internal threads.

**Step 2.** After Step 1 finishes, change these values and source the script again:
```r
task <- "results"
summaries <- "regenerated"
```

This reads the recalculated summaries in `output/aggregated/` and creates the
manuscript outputs. **The route is now complete. Do not run another aggregation
after this step.** `summaries` chooses the input summaries for `"results"`; it
does not control `"aggregate"`.

The verified parameter-summary tables for all five collections agree with the
supplied production summaries within numerical tolerance `1e-10`. The regenerated
16 manuscript CSVs and 16 LaTeX tables were byte-identical to the references.
This is numerical agreement, not a promise that every file (especially figures
with timestamps or software-dependent rendering) has an identical file checksum.

Both routes write to `output/manuscript/`, so running one after the other replaces
those generated presentation files, not the archived references. For `"results"`,
`collection` does not restrict the manuscript outputs; for both tasks,
`config_ids` and `replicates` do not limit the archived production data used.

### Fresh Spot Checks Are A Separate Chain

Use **Workflow 2**: `"plan"` -> `"fit"` -> `"records"` -> review
`output/local-simulations/editor-records/`. This includes local per-replicate
comparisons and `reduced_run_summary.csv`. Neither `"aggregate"` nor `"results"`
uses these fresh local fits. A five-replicate spot check is not an alternative
input for recreating the full production-study manuscript tables.

**Optional:** `task <- "validate"` checks sources, saved summaries/draws and output
agreement, writing `documentation/validation-report.md`. Review its `REVIEW` and
`PENDING` entries; execution alone does not certify acceptable diagnostics.

## Software Setup

Install dependencies once if needed; scripts do not install them automatically.
```r
install.packages(c("tidyverse", "survival", "MASS", "nlme", "loo",
                  "posterior", "patchwork", "JMbayes2"))
install.packages("cmdstanr",
                 repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::check_cmdstan_toolchain()
# Only if needed:
# cmdstanr::install_cmdstan()
```

Stan fits need CmdStan and a working C++ toolchain. Leave `cmdstan <- ""` to
use the configured installation, or supply its path. Bunya production used
CmdStan 2.37.0. The editable `no_pch <- TRUE` bypasses cached-header compiler
incompatibility; it does not change the model. Local tests used R 4.5.2,
CmdStanR 0.9.0 and CmdStan 2.40.0. Full package/compiler versions and loaded-package
`sessionInfo()` are in [software environment](documentation/software-environment.md).

## Outputs And References

**`archived-production-results/` contains the supplied historical results and
is input-only.** Fresh outputs are kept separately:

| Location | Contents |
|---|---|
| `output/case-study/` | Mock fit, table and diagnostics |
| `output/local-simulations/fits/<collection>/configNNN/repNNN/` | Actual local fits and per-replicate records |
| `output/local-simulations/prepared/` | Optional preparation-only inputs |
| `output/local-simulations/editor-records/` | Combined estimates, diagnostics, runtimes and comparisons |
| `output/simulation-plans/` | Candidates and selected plans |
| `output/aggregated/` | Recalculated full-result summaries/plots |
| `output/manuscript/` | Manuscript `tables/`, `figures/` and `derived_csv/` |

Local results are organised by IDs, not by invocation or replicate count.
Combined records can include earlier selections. Rerunning tasks may replace
generated files, but not the archived references.

## Runtime And Agreement

On an Apple M2 Pro (32 GB RAM), the mock fit took about 15 minutes, five main
configuration-1 fits about 76 sampling minutes combined, and one comparator
dataset about 2 minutes. Difficult fits can take hours. Presentation generation
takes minutes; full aggregation reads tens of GB and can take tens of minutes or longer.

Bunya production used UQ EPYC3/EPYC4 nodes; main-replicate times averaged
27.16 minutes (median 21.38; maximum 153.57), excluding some overhead. See
[local checks](documentation/local-fit-notes.md) for timings and diagnostics.

Fresh MCMC draws need not match across platforms; see the
[Stan reproducibility guide](https://mc-stan.org/docs/reference-manual/reproducibility.html).
Review divergences, treedepth, R-hat, ESS and Monte Carlo uncertainty. Historical
summaries lack some MCSEs; rounded agreement or interval overlap does not prove convergence.

### LOO/WAIC Warnings

Scores are computed after fitting and do not change posterior estimates. They
use patient-level joint likelihoods conditional on fitted random effects, not
strict new-patient prediction. High Pareto-k or many `p_waic > 0.4` values warn
against relying on the scores. Retain warnings and do not use flagged scores to
judge parameter agreement or select models without further validation. See
[LOO diagnostics](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html).

## Adaptations From The Production Workflow

`run_checks.R` replaces SLURM dispatch, using a shared scenario generator and
package-relative paths. Seeds, generated data, centering, priors, initial values
and sampling settings are preserved; spot checks reduce dataset counts, not the
model. Temporary task settings are restored and run records saved. The case
study substitutes mock data for restricted ENZAMET data. See
[verification](documentation/code-equivalence.md) and the retained
[development history](documentation/changes.md).

## Supporting Files

Only `run_checks.R` needs editing for routine use. `R/` contains execution
helpers; `production/` contains local fitting/generation code; `analysis/`
contains aggregation/presentation code; `case-study/` contains mock inputs and
references. [manuscript_map.csv](manuscript_map.csv) links outputs to manuscript
results; `files_manifest.csv` records checksums.

Data provenance/definitions are in the [data dictionary](documentation/data_dictionary.md).
Optional details: [technical guide](documentation/technical-guide.md),
[function interfaces](documentation/function-reference.md) and
[terminal commands](documentation/command-line.md).

The [essential](documentation/essential-checklist-response.md) and
[detailed](documentation/detailed-checklist-response.md) responses identify pending
author sign-off: manuscript mapping, permissions/licensing, independent-machine
testing and delivery. The approximately 40.6 GB archive requires an editor-accessible
repository or transfer link without an institutional login, not email attachment.

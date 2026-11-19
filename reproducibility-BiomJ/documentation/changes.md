# Changes Made For The Audit Package

## Check Wording And Independent Mock Check (19 November 2026)

The concluding `check` message now describes a setup-only check, with detailed
result validation and MCMC diagnostics assessed separately. Check behaviour is
unchanged. Checklist responses record PM's independently run mock-case-study
table screenshot without implying independent testing of all workflows or a
sampling-diagnostic assessment. No scientific calculations or results changed.

## Overall Preflight Check (19 November 2026)

`check` now continues through missing archive collections, mock inputs, packages,
CmdStan/toolchain components and failed recording/export stages. It returns and
saves a PASS/WARNING table and emits one final summary warning, rather than
stopping at the first missing archive. Invalid explicit CmdStan paths and a
missing configuration preview are reported within this task; other tasks retain
their existing strict input checks. Tests cover incomplete and complete fixtures
without model compilation or MCMC. No fitting, data-generation or aggregation
calculations changed, and no installed software is modified automatically.

## onBunya Linux Linker Workaround (19 November 2026)

README now documents the reported PIE/non-PIE mismatch involving CmdStan's
`main_threads.o`. The author confirmed that appending `LDFLAGS += -no-pie` to
their user-owned CmdStan 2.37.0 configuration allowed compilation and sampling
to start on the onBunya desktop. Completed Linux results have not yet been
verified. This is an opt-in installation-level workaround, not a change to the
default build, Stan source, data generation, priors or sampling settings.

## Mock Reference Filenames (19 November 2026)

The four author-supplied reference files now identify the Windows-AM, Windows-DM
and MacBook-M2Pro-DM runs. Documentation, the manuscript map, manifest paths and
the CSV Git-ignore exception were updated. All four files retain their original
checksums. Fresh case-study outputs remain named `mock_table1_update.*`; no
fitting code, data or result calculations changed.

## Manuscript Output Mapping (19 November 2026)

The map now uses author-selected arXiv:2606.16080v1, including its supplement.
It identifies all numbered tables/figures, adds four previously omitted
supplementary tables, separates extra plots from published figures, and records
restricted real-data exceptions. Nine simulation table bodies match the arXiv
source; all three published figures render identically under the same renderer.
See manuscript-coverage.md. Only documentation/metadata changed; no models,
generated data, saved results or calculation code were altered or rerun.

## Case-Study Table Namespace Fix (18 November 2026)

The final table explicitly calls `dplyr::select()` so it also works after a
simulation has loaded MASS. A non-sampling regression test covers that conflict.
No data, model, initial values or numerical summary calculations are changed.

## Detailed-Checklist Compliance Pass (18 November 2026)

All changes in this pass are confined to this reproducibility folder. Local-only
main/comparator definitions replace the full cluster drivers, five scenario
copies are consolidated into a shared simulator/coefficient table, and unused
source copies are removed. The original project outside this folder is not
edited. One unchanged Stan source is shared by both fitting workflows. Temporary
task settings are restored, function interfaces/software versions are documented,
and the project no longer automatically restores or saves a workspace.

Before/after checks cover 80 generated datasets, 43 complete Stan input/start
cases, all 40 comparator preparation cases and the saved mock inputs. Numerical
values and RNG states match exactly on the recorded environment. Recognised
source-only checkpoint mappings retain old fit provenance without refitting.
See code-equivalence.md and source-provenance.csv; earlier raw source hashes are
retained as original-source-checksums.csv. Independent-machine testing, final
manuscript mapping, licensing and archive delivery remain author actions.

## Final Layout And Development Cleanup

Historical results are now under `archived-production-results/`; `production/`
contains source code only. Current local outputs remain under
`output/local-simulations/` and `output/case-study/`. The superseded five-case
development outputs, extracted development references/selection list, abandoned
launcher records, unused binaries and generated test fixtures were removed with
permission. Unused legacy dispatch/extraction helpers and their validation paths
were removed as well. Tests now create disposable fixtures in R's temporary
directory. Validation inspects the current saved fits without rewriting them.
The README is workflow-focused; current execution notes replace the earlier
development narrative. No current model fit or archived numerical result was
deleted or refitted.

## Clearer Directory Names

`output/five-per-config/` was renamed to `output/local-simulations/`, since
selected runs can contain any number of replicates. `Sim-Results/` was renamed
to `production-results/` to distinguish supplied historical production results
from fresh local outputs. The two directories were moved in place, not copied
or regenerated. Original dated collection and batch names inside the archive
are retained for traceability. Code, tests, documentation, current plan/index
CSVs and the file manifest use the new paths. Raw sampler files, original RDS
run snapshots and historical logs retain the paths recorded at execution time;
their numerical contents and provenance are not rewritten. Original production
and mock-example source copies remain unchanged, including historical comments.

## Local Compiler Header Workaround

The editable `no_pch` setting in `run_checks.R` is enabled after the local
case-study compilation reported a cached Apple Clang 2100.3.34.2 header versus
the installed compiler 2100.4.1.1. The existing wrapper passes
`PRECOMPILED_HEADERS = FALSE` to CmdStan compilation for simulations and the
mock case study. This does not rebuild the CmdStan installation or change the
statistical model. The setting remains editable for other machines.

## Preparation Checkpoint Compatibility

Preparation-only checkpoints no longer reject a difference in `no_pch`:
preparation does not compile or run Stan. The build choice remains checked for
actual Stan fits, and all other identity fields remain checked. JMbayes2 also
continues to ignore this Stan-only flag. Saved checkpoint provenance is not
modified. Errors now name differing identity fields or missing files and give
the RStudio `resume <- FALSE` instruction as well as its CLI equivalent.

## RStudio Entry Point

`run_checks.R` is the recommended entry point, with numbered sections and an
editable settings block. Its default only plans five replicates of an example
configuration; it does not start MCMC. `R/entry.R` is shared with the optional
`main.R` CLI, so both interfaces call the existing fitting/recording functions.
Task aliases `fit` and `records` mean `spot-check` and `export-local`.
The no-PCH environment option is restored after dispatch in an interactive
session. All original source copies, statistical settings, reference results
and output paths are retained. No automatic model fits were added.

The root README now prioritizes RStudio. Its former detailed contents and the
former root `RUN_LOCAL.md` are retained as `documentation/technical-guide.md`
and `documentation/command-line.md`. Supporting folders were not renamed or
merged, avoiding broken paths and unnecessary copies of the large archive.

Validation: `tests/test-rstudio-entry.R` checks `source()` from both supported
project directories, task routing, ID validation and build-option restoration
without MCMC. The existing dispatch/export tests and package validation also
pass. All 11 recorded production source checksums remain unchanged. This entry
point change was not validated by starting another posterior sampling run.

## Selected-Setting Audit

The unattended all-settings coordinator and Mac launcher were removed on
17 November 2026 at the user's request. Fitting/preparation now requires an
explicit configuration selection. The existing result exports, metadata,
checkpoints and arbitrary original-replicate selection are retained so the
author can supply a few local reference settings and the editor can choose
additional checks against production. Existing outputs and production data
were not deleted. No new fits were started as part of this reversion.

The original project R and Stan files were not edited. The existing cleaned
`Bunya/ENZAMET/Sim-Results` directory was moved into this package with permission,
instead of making a second approximately 37 GB copy.

## Unmodified Copies

`production/` contains the primary simulation driver, comparator driver, five
scenario generators, selected Gamma5 Stan model, two aggregation scripts and
the manuscript-output generator. `case-study/source/` contains the relevant
unmodified source/data/README extracted from the uploaded ZIP. Original platform
binaries, R history and RStudio caches were not extracted.

## Small Adapters

`R/audit.R` loads production expressions only up to the SLURM dispatch block.
It skips workspace clearing and fixed CmdStan-path selection, redirects generator
source paths to the included copies, and substitutes a recording wrapper around
`cmdstan_model()`. It does not rewrite the model or generator calculations.
The configuration grid is evaluated from the production expressions without
the array-job mapping or environment filters. Configuration values are passed
explicitly and retain the production function's argument-forcing protections.

The wrapper records complete inputs/initial values, uses persistent Stan CSV
paths, records chain IDs explicitly as 1 to 4, and optionally reduces the number
of simultaneously running chains without changing the total four chains. It
does not change warmup, retained iterations, treedepth, adapt_delta, priors or
likelihood. The optional `--no-pch=true` build flag bypasses incompatible cached
compiler headers; it is recorded as a build choice, not a statistical-model change.

The mock wrapper redirects only file I/O, skips workspace clearing and replaces
serialization of the R6 fit with portable numeric summaries plus persistent CSVs.
It retains the supplied script's initialization random-number draw order, which
is not identical to the production script's order. The initialization strategy
and statistical target are the same; exact random starts across the two distinct
datasets/scripts are not claimed to match.

`analysis/` contains working copies of the production aggregation and presentation
scripts. Edits remove workspace clearing and supply explicit package-relative
input/output paths. Calculations and plotting expressions are unchanged. Some
historical exploratory slide figures no longer produced by the current generator
remain in the archived reference folder; they are not silently substituted into
newly generated results.

## Deliberately Not Changed

No refactoring of duplicated historical generator functions, no adjustment of
Stan priors or initial-value distributions, no change to empirical centering,
no fit exclusion rule, no numerical likelihood rewrite, no new data imputation,
no attempt to force cross-platform bitwise equality, and no changes to the
restricted real-data analysis. These would be scientific or implementation
changes beyond the requested minimal-edit audit preparation.

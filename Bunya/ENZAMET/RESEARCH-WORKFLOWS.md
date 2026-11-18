# ENZAMET Research Workflow Index

The editor-facing entry point is `reproducibility-BiomJ/README.md`. This index
organises the older research files by experimental stage, broadly from earlier
tests to production and publication. Paths below are relative to the repository.

This is a retrospective grouping of current files, not reconstructed historical
code versions. Some launchers were updated for later reruns, including joint
LOO/WAIC output. Their shared R/Stan files also evolved. Old filenames and output
labels therefore do not guarantee that running them now recreates the old fit.
Do not backdate commits, silently replace shared code, or submit old jobs merely
to test these archives. Files remain at their existing paths to preserve calls.

## Completed Commit And Push Record

The nine approved commits below were created in this order and pushed to
`origin/main` (`https://github.com/UQ-ULTRA/JoMoNoPH-share.git`). They record the
repository organisation, not the original dates when the experiments ran.

| Scope | Commit | Contents |
| --- | --- | --- |
| Housekeeping | `ee70c5b` | Exclude generated results; retain small audit inputs |
| Group 1 | `edafd84` | Early test/rescue wrappers, this index and generated-HTML ignore rule |
| Group 2 | `f172d14` | Standard BP prior and shape-penalty experiments |
| Group 3 | `bce9825` | Alternative BP implementation experiments |
| Group 4 | `f9e1f89` | Zero-prior pilot and final sAFT rescue workflows |
| Group 5 | `60bd707` | JMbayes2 comparator workflows and aggregation |
| Group 6A | `9d02c1c` | Original application, aggregation, presentation and HPC setup |
| Group 6B | `7959dfa` | Historical reproducibility review notes |
| Group 7 | `fbb8523` | Self-contained reproducibility package: 45 files |

Group 7 includes editor-facing defaults (`task = "check"`, `collection = "main"`,
`summaries = "saved"`) and refreshed manifest entries. Both non-sampling regression
tests passed in a temporary package copy, and all 19 R files parsed. This commit
preparation did not rerun MCMC or change scientific settings or saved results.
Large result archives and generated outputs were not uploaded. The six deferred
files listed below remain locally modified or untracked and were not part of
these commits or the push.

## Group 1: Earlier Test And Rescue Wrappers

- `Bunya/ENZAMET/ENZAMET-DM.slurm`
- `Bunya/ENZAMET/ENZAMET-DM-test-rescue.slurm`

The committed change to the first wrapper reduces memory from 20 GB to 16 GB. These
40-config/10-replicate launchers use the historical td12-labelled output folder
and current `ENZAMET-DM.R`; they do not themselves pin the old model, treedepth
or prior settings. Preserve them as historical wrappers, not current production.

## Group 2: Standard BP Prior/Penalty Tests And Scaling Pilot

- `Stan/Bernstein-Polynomials-JM-Hist-Gamma05.stan`
- `Stan/Bernstein-Polynomials-JM-Hist-Gamma5.stan`
- `Stan/Bernstein-Polynomials-JM-Hist-Penalty-DM.stan`
- `Bunya/ENZAMET/ENZAMET-DM-gamma05-td10-10rep.slurm`
- `Bunya/ENZAMET/ENZAMET-DM-gamma1-td10-10rep.slurm`
- `Bunya/ENZAMET/ENZAMET-DM-gamma5-td10-10rep.slurm`
- `Bunya/ENZAMET/ENZAMET-DM-gamma1-shapepenalty-td10-10rep.slurm`
- `Bunya/ENZAMET/ENZAMET-DM-gamma05-td10-100rep.slurm`

Tests cover positive-normal weight-prior SDs 0.5, 1 and 5, a shape penalty,
and a larger SD-0.5 pilot. The SD-1 standard model and main R driver are already
tracked. Some 10-replicate wrappers now target later joint-IC rerun outputs.

## Group 3: Alternative BP Implementation Tests

- `Stan/Bernstein-Polynomials-JM-Hist-alt-DM-Gamma05.stan`
- `Bunya/ENZAMET/ENZAMET-DM-altDM-gamma05-td10-10rep.slurm`
- `Bunya/ENZAMET/ENZAMET-DM-altDM-gamma1-td10-10rep.slurm`
- `Bunya/ENZAMET/ENZAMET-DM-altDM-gamma5-td10-10rep.slurm`

The alternative implementation constructs BP bases directly. Its SD-1 and
SD-5 Stan sources are already tracked. These are implementation/prior tests,
not additional published production collections.

## Group 4: Zero-Mean-Prior Pilot And Final sAFT Rescues

- `Bunya/ENZAMET/ENZAMET-DM-altDM-gamma5-100rep-zeroprior-initinfo-1.slurm`
- `Bunya/ENZAMET/ENZAMET-DM-altDM-gamma5-100rep-zeroprior-initinfo-2.slurm`
- `Bunya/ENZAMET/rescue_slurms_20260527_timeout_20rep_1000rep_zeroprior_initinfo/ENZAMET-DM-altDM-gamma5-1000rep-zeroprior-initinfo-timeout-rescue-standard.slurm`
- `Bunya/ENZAMET/rescue_slurms_20260527_timeout_20rep_1000rep_zeroprior_initinfo/ENZAMET-DM-altDM-gamma5-1000rep-zeroprior-initinfo-timeout-rescue-W210-scen3-batch171.slurm`
- `Bunya/ENZAMET/rescue_slurms_20260527_timeout_20rep_1000rep_zeroprior_initinfo/README.md`
- `Bunya/ENZAMET/rescue_slurms_20260530_alphaSens_LL_13batch/ENZAMET-DM-altDM-gamma5-alphaSens-LL-1000rep-rescue-standard.slurm`
- `Bunya/ENZAMET/rescue_slurms_20260530_alphaSens_LL_13batch/ENZAMET-DM-altDM-gamma5-alphaSens-LL-1000rep-rescue-highTime.slurm`
- `Bunya/ENZAMET/rescue_slurms_20260530_alphaSens_LL_13batch/README.md`

The main 40-config/1000-replicate production wrappers, 10-config log-logistic
alpha-sensitivity wrappers, shared launchers, generators and Gamma5 model are
already tracked. Keep these pilot and rescue records alongside that provenance.
The named rescue arrays are historical targets, not fresh production submissions.

## Group 5: JMbayes2 Comparator Pilots, Production And Rescues

- `Bunya/ENZAMET/ENZAMET-JMbayes2-DM.R`
- `Bunya/ENZAMET/run-enzamet-jmbayes2-dm-common.sh`
- `Bunya/ENZAMET/ENZAMET-JMbayes2-DM-LL-pilot.slurm`
- `Bunya/ENZAMET/ENZAMET-JMbayes2-DM-all40-pilot.slurm`
- `Bunya/ENZAMET/run-enzamet-jmbayes2-dm-all40-1000rep-common.sh`
- `Bunya/ENZAMET/ENZAMET-JMbayes2-DM-all40-1000rep-01.slurm`
- `Bunya/ENZAMET/ENZAMET-JMbayes2-DM-all40-1000rep-02.slurm`
- `Bunya/ENZAMET/ENZAMET-JMbayes2-DM-all40-1000rep-03.slurm`
- `Bunya/ENZAMET/ENZAMET-JMbayes2-DM-all40-1000rep-04.slurm`
- `Bunya/ENZAMET/ENZAMET-JMbayes2-DM-all40-1000rep-05.slurm`
- `Bunya/ENZAMET/ENZAMET-JMbayes2-DM-all40-1000rep-rescue.slurm`
- `Bunya/ENZAMET/ENZAMET-JMbayes2-DM-all40-1000rep-rescue2.slurm`
- `Bunya/ENZAMET/ENZAMET-JMbayes2-DM-all40-1000rep-randSDsummary-rescue.slurm`
- `Bunya/SumSimRes-ENZAMET-JMbayes2-DM.R`

The pilot launcher defaults differ from production (including chain counts),
so they must not be described as interchangeable. The all-40 production launcher
sets four chains and targets the final random-effect-SD-summary collection.
The first two rescue files target earlier comparator outputs; the randSDsummary
rescue targets the final collection. All use the already-tracked generators.

## Group 6A: Original Application, Presentation And HPC Setup

- `Bunya/single-stan-run-ENZAMET-DM.R`
- `Bunya/SumSimRes-ENZAMET-DM.R`
- `Bunya/Manuscript/generate_ENZAMET_manuscript_outputs.R`
- `Bunya/Bunya-simulation-setup.md`
- `Stan/compile-stan-models.slurm`

The application script requires restricted ENZAMET data and is not the public
mock example. The committed change to the original sAFT aggregation script
selects alpha300 instead of the main collection as the default input, without
altering its aggregation calculations.
The presentation generator reads summaries in original project locations. The
committed compilation-wrapper change redirects logs. These files preserve original
workflows; they do not replace the self-contained audit copies.

## Group 6B: Historical Review Notes (Separate Commit)

- `Bunya/Manuscript/audit_mock_case_study_review_20261111.md`
- `Bunya/Manuscript/sim_results_archive_audit_20261116.md`

Preserve these dated reviews as provenance, not current running instructions.
Historical banners point readers to the current audit README. Findings, paths
and pending actions describe the state at the original review date and may have
since been resolved or superseded; the historical text is retained.

## Group 7: Editor-Facing Reproducibility Package

Git-eligible files under `reproducibility-BiomJ/` were committed separately after
the earlier groups. Its local sources are the versions tested for the audit;
do not replace them with the older drivers above. Large production results and
generated outputs remain excluded, while the synthetic mock input is included.

## Deferred: Separate Review, Not Part Of The Above Commits

- `Bunya/ENZAMET/ENZAMET-single-aft-only-DM.R`: substantial exploratory changes,
  including log-time BP support and additional AFT comparators.
- `Stan/Bernstein-Polynomials-JM-Hist-alt.stan`: a different model/interface change,
  including renamed survival data/parameters; not the production alt-DM-Gamma5 model.
- `Stan/Bernstein-Polynomials-JM-Hist-DM.stan`: separate shape/amplitude penalty model,
  not the fixed shape-penalty variant in Group 2.
- `Bunya/ENZAMET/ENZAMET-single.R`: older multi-model local prototype with workstation
  paths and optional model dependencies; not needed by published production.
- `Stan/Gaussian-Basis-JM-Hist-Quantile-PM.stan`: separate exploratory basis model.
- `Stan/general_aft-PM.stan`: includes `aft_functions.stan`, absent from this checkout;
  do not present this as a complete runnable model without resolving that dependency.

Deferred means retained locally for an explicit later decision, not deleted,
reverted or silently added to an unrelated commit. Future changes require
approval of their exact staged scope before committing; pushes require separate approval.

## Generated Side Product: Excluded From Git

`Bunya/Bunya-simulation-setup.html` is the rendered counterpart of the Markdown
guide. It is ignored and retained locally; the source `.md` was committed in
Group 6A. The generated-HTML ignore rule was committed in Group 1, without
amending the earlier housekeeping commit.

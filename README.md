# Bayesian Joint Modelling With Semiparametric AFT Models

Research code for **Bayesian joint modelling using semiparametric accelerated
failure time approaches** (Biometrical Journal manuscript bimj.202500502),
including ENZAMET-based simulations, model comparisons and prior sensitivity
analyses.

## Reproducibility Audit: Start Here

Use the self-contained [reproducibility package](reproducibility-BiomJ/README.md).
Open `reproducibility-BiomJ/reproducibility-BiomJ.Rproj` in RStudio, then edit
and run [run_checks.R](reproducibility-BiomJ/run_checks.R). Its README covers
software setup, expected runtimes and three independent workflows:

- Fit the supplied mock case-study data and compare with mock reference results.
- Run selected simulation replicates locally and compare with matching production results.
- Recreate manuscript tables and figures from archived production results.

The mock data and mock references are included. Restricted real ENZAMET data
are not shared. Large simulation results, generated outputs and compiled Stan
binaries are excluded from Git. Obtain `archived-production-results/` separately
using the link supplied with the audit submission and place it inside
`reproducibility-BiomJ/`. The archive is approximately 40.6 GB; see the package
README for storage requirements. The mock case-study workflow does not require it.

## Research Files

| Location | Purpose |
| --- | --- |
| `Bunya/ENZAMET/` | Original simulation drivers, SLURM launchers, pilots and rescue workflows |
| `Bunya/SumSimRes-ENZAMET-DM.R` | Original sAFT-JM simulation aggregation |
| `Bunya/SumSimRes-ENZAMET-JMbayes2-DM.R` | Original comparator aggregation |
| `Bunya/Manuscript/` | Manuscript-output generation and historical review notes |
| `Stan/` | Stan models, experimental variants and compilation helpers |
| `reproducibility-BiomJ/` | Portable audit workflow with its own source copies and documentation |

The [research workflow index](Bunya/ENZAMET/RESEARCH-WORKFLOWS.md) groups the
historical experiments and records the archive commits. Older launchers share
code that evolved over time; their filenames alone do not reproduce historical
settings. For audit runs, follow the reproducibility-package README rather than
submitting historical SLURM jobs.

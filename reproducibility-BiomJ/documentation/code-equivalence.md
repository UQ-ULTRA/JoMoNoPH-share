# Numerical Equivalence Of The Local Cleanup

Verified on the authors' recorded local environment on 18 November 2026.
No new posterior sampling was used for these source-cleanup checks.

## What Changed

- One simulator and five-row coefficient table replace the duplicated scenario scripts.
- Local main/comparator definitions replace the historical SLURM drivers.
- Unused main-model branches, fixed workstation paths and workspace clearing were removed.
- The case-study and simulation share one Stan source. Its checksum is unchanged.
- Local recording, session restoration and function documentation were added.

## Checks Performed

| Check | Coverage | Result |
|---|---|---|
| Generated observations and final R RNG state | All 40 configurations, replicate IDs 1 and 56 (80 datasets) | Exact match |
| Stan data, all four initialization lists, response center and pre-sampling RNG state | All 40 configurations at replicate 1; configuration 1/replicate 56 at alpha factors 1.25, 1.5 and 3 (43 cases) | Exact match |
| Comparator pre-MCMC inputs | All 40 configurations at replicate 1: observations, preliminary LMM coefficients/covariance/random effects/residuals, controls and RNG | Exact match |
| Mock example inputs | Data, four initialization lists, centering and RNG versus the saved completed fit | Exact match |
| Main/comparator checkpoint compatibility | Existing completed local checks | Reused without refitting |
| Saved posterior summaries | Current local draws versus saved summaries | Checked by task validate |
| Summary/presentation outputs | Full saved-result aggregation and manuscript tables | See validation-report.md |

Hashes and per-case outcomes are retained under `output/compliance-checks/`.
`tests/check-code-equivalence.R` is an optional exact-hash regression check on
the recorded author environment, not the recommended cross-platform audit
entry point. Different numerical libraries can change such hashes; use the
README's workflows and diagnostic/Monte Carlo comparisons on other machines.
`source-provenance.csv` records original/current source hashes. Old raw source
hashes can be accepted for existing local checkpoints only through the explicit
`verified-source-cleanup.csv` mapping; data IDs, seeds, priors, software-version
checks and the unchanged Stan checksum are not bypassed. New checkpoints record
the current sources, including the extracted configuration grid.

The exact before/after checks concern this controlled cleanup on one environment.
They are not a promise of bitwise-identical MCMC draws across platforms. Full
new posterior runs for all configurations were not performed. Future source or
software changes require their own checks and must not be silently added to the
compatibility mapping.

The existing generation mechanism is intentionally preserved, including
Normal(0, 1)-month visit jitter (baseline reset to zero). No timing windows,
priors, likelihoods, stopping rules or diagnostics were altered to improve results.

# Data And Saved-Result Dictionary

## Synthetic Case Study

Source: Andrew's uploaded `audit-mock-case-study 2.zip`. These are simulated
mock records, not the ENZAMET clinical-trial data. The authors must confirm
permission to redistribute the mock data and their code before submission.

`case-study/source/mock-data.rds` is an R list:

| Object/field | Meaning |
|---|---|
| `QoL` | One row per synthetic assessment |
| `QoL$id` | Synthetic patient identifier |
| `QoL$time` | Assessment time in months |
| `QoL$D6` | Mock HRQoL score; may be missing |
| `QoL$arm` | Treatment indicator 0/1 |
| `cPFS` | One row per synthetic patient |
| `cPFS$id` | Matching synthetic identifier |
| `cPFS$time` | Observed follow-up time in months |
| `cPFS$status` | Event indicator 1/0 |
| `cPFS$arm` | Treatment indicator 0/1 |

The supplied file has 1,105 survival rows and 17,623 assessment rows, including
176 missing D6 values. The supplied script handles missing D6 through its LMM
model frame and shifts follow-up times by 1e-6 months. Some mock assessments
extend beyond recorded follow-up; this example is supplied to exercise the
fitting code and is not a claim about the chronology or results of real ENZAMET.
No records were silently removed or altered by the audit adapter.

## Simulation Batches

The shared generator is `production/Bunya/ENZAMET/Data-Gen-ENZAMET.R`; its
five-row coefficient table defines the scenarios. The original numerical
algorithm and RNG order are retained. In particular, scheduled assessment times
use independent Normal(0, 1)-month jitter, with baseline reset to zero; this is
not fixed-window uniform jitter. Visits at or beyond observed follow-up are
excluded. These details document existing generation, not new changes.

Each top-level simulation output `.rds` is a list of replicate results. The
filename contains model, basis count, sample size, censoring setting, baseline
distribution, coefficient scenario, follow-up and batch number. The main and
alpha-sensitivity batches have four entries; comparator batches have ten.

Main/sensitivity fields include `unique_seed`, `censor_prop`,
`longitudinal_response_center`, `joint_summary`, `joint_diagnostics`, `lme_coefs`,
`prior_info`, `init_info`, `runtimes`, `loo_joint`, and `waic_joint`. They are
compact fitted-result summaries, not raw generated datasets or full posterior
draws. `init_info` deliberately omits the large subject-level random-effect
starting matrix in historical production. The local audit saves complete
initializations separately.

Comparator entries include `comparison_tidy`, `jmb2_results`, `model_status`,
`runtimes` and preliminary-fit summaries. Missing model-specific estimates must
be interpreted with their recorded fit status, not as zero estimates.

The saved aggregate CSVs define their quantities through descriptive column
names, including bias, empirical SD, mean posterior SD, RMSE, coverage,
Monte Carlo SE where recorded, fit diagnostics and runtimes. The aggregation
code is the definitive record of the denominators and missing-result handling.

Parameter names in Stan differ from manuscript notation:

| Stan name | Manuscript meaning |
|---|---|
| `beta_long[1]` | Longitudinal fixed intercept beta_0 |
| `beta_long[2]` | Time fixed effect beta_1 |
| `beta_long[3]` | Treatment-by-time effect beta_2 |
| `beta_surv[1]` | Survival treatment coefficient gamma |
| `gamma[1]` to `gamma[5]` | Positive Bernstein weights theta_1 to theta_5 |
| `sd_1_long[1]`, `[2]` | Random intercept and random slope SDs |
| `sigma_long` | Residual SD |
| `alpha` | Original-scale longitudinal-survival association |
| `alpha_tilde` | Computationally scaled association parameter |

Pure response centering leaves slopes and SDs on their existing scales. Saved
posterior intercept means and 95% intervals are back-transformed. Raw Stan draws
and CSVs remain on the centered-response parameterization and are labelled so.

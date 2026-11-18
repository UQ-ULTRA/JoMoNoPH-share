# Mock case-study reproducibility review

> Historical review of the original supplied example, not current instructions.
> Several recommendations below have since been addressed. For the current
> package and verified status, see the [audit README](../../reproducibility-BiomJ/README.md).

Reviewed 2026-11-11 against the supplied `audit-mock-case-study 2.zip`.

## Scope and verification

- Compared the supplied Stan source against
  `Stan/Bernstein-Polynomials-JM-Hist-alt-DM-Gamma5.stan`, used by the
  40-configuration, 1000-replicate production launcher.
- Read `fit-case-study-model.R`, `README.md`, the reference CSV, and the mock RDS.
- Read both editor-provided reproducibility checklists.
- Executed the mock data preparation, staging LMM/exponential AFT fits, and
  initialization code. No MCMC was run.
- Compared the resulting data with the production `convert_data_from_models()`
  function on the same complete observations. All 28 shared fields were
  numerically equal, ignoring non-numerical attributes.
- Stan syntax passed using the locally installed CmdStan 2.40.0. This is not a
  new posterior fit or verification under the supplied reference environment.
- Checked prior metadata in one saved production batch, confirming zero beta
  prior means, SD 10, and raw alpha prior SD approximately 0.3536465.

The original OneDrive directory was inaccessible to this session. All mock
file inspection and checks used a temporary extraction of the supplied ZIP.
The supplied files and production source files were not changed.

## Settings comparison

The mock Stan source is identical to the production source after normalizing
Windows versus Unix line endings. Both include the longitudinal and survival
likelihoods in the posterior and generate patient-level conditional joint
log likelihoods.

| Component | Production and mock example |
|---|---|
| Longitudinal fixed effects | Intercept, time, time by treatment |
| Random effects | Patient intercept and time slope; correlated Gaussian effects |
| Response centring | Subtract the mean of non-missing observed responses before staging/fitting |
| Intercept reporting | Add the observed response mean back to posterior/LMM intercept summaries |
| Regression priors | Independent Normal(0, 10^2), on the centred response scale where applicable |
| BP weights | Five positive weights; Normal(0, 5^2) restricted to positive values |
| Residual SD prior | Cauchy(0, 5) restricted to positive values |
| Random-effect SD priors | Cauchy(0, 5), with the shared Stan lower bound of 0.001 |
| Correlation prior | LKJ correlation prior with shape 2 |
| Alpha prior | Raw alpha SD log(2)/1.96; implemented using inverse scaling of alpha_tilde SD |
| Longitudinal staging fit | nlme::lme(), REML; same model and optimizer controls |
| Survival staging fit | survival::survreg(), exponential AFT, used for starting values |
| Chains | Four, one thread each |
| Iterations | 1000 warm-up plus 1000 retained iterations per chain |
| Sampling controls | max_treedepth 10; adapt_delta not supplied (CmdStan default) |
| Roughness penalty | None; r and tau_h are accepted but unused |

Initial-value distributions also agree: fixed effects from the LMM with SD 0.1
jitter; survival treatment coefficient from exponential survreg with SD 0.05
jitter; residual SD drawn around the LMM SD with SD 1 and made positive;
random-effect SDs augmented by absolute Normal(0, 0.05) jitter; correlation
initialized to 0.2; standard-normal latent effects; alpha_tilde drawn from
Normal(0, 0.1); BP weights drawn as abs(Normal(0.1, 0.01)).

The mock script draws the BP initial values before alpha_tilde, whereas the
production code draws alpha_tilde first. This changes the realized initial
values under a common R seed but not their distributions or the posterior model.
The single-example seed 10235 is intentionally different from the per-replicate
production seeds. The mock example uses 1105 patients rather than the
simulation's 1100 patients.

## Mock-data checks and reproduced staging values

- 1105 survival records and 17623 longitudinal rows, including 176 missing D6
  values; 17447 measurements enter the fit.
- 668 events; observed censoring proportion approximately 39.55%.
- No duplicate survival IDs; every survival patient has a non-missing D6 value.
- Treatment assignment agrees between the two datasets.
- Survival rows match the random-effect level order expected by the Stan
  likelihood. Preserve this ordering or explicitly reorder before fitting;
  part of the likelihood uses b_long directly rather than applying a mapping.
- D6 range among non-missing values: 0.9420715 to 100.
- Response mean subtracted: 72.56573.
- Reproduced raw-scale LMM fixed effects: 73.28292199, -0.03083351,
  -0.04625706. These round to the values shown on all three machines.
- LMM residual SD: approximately 11.52495; random-effect SDs approximately
  14.20437 and 0.2077552.
- Exponential staging treatment coefficient: approximately 0.8319389.
- s_surv = 1; alpha_tilde prior SD approximately 4.075757; implied raw alpha
  prior SD approximately 0.3536465.
- All generated initial values were finite.

Two data/preprocessing differences need clarification:

1. Two patients have an event recorded at time zero. The mock script adds
   1e-6 to every survival time after replacing negative times by zero. This
   avoids a zero input but should be documented as a mock-data convention.
   The production converter does not perform this time shift.
2. There are 172 non-missing HRQoL observations from 172 patients after their
   recorded follow-up time. The excess ranges from 0.110 to 2.999 months, with
   median 1.374 months. The production generator retains measurements only
   before observed follow-up. Confirm whether the mock-data construction
   intentionally permits this. If the intended case-study dataset stops at
   progression/follow-up, make the mock chronology consistent and regenerate
   all mock reference results after any change.

## Cross-machine screenshots

The following are ranges of printed posterior means, not MCSE calculations:

| Parameter | Windows Ding | Windows Andrew | Mac M2 Pro | Maximum minus minimum |
|---|---:|---:|---:|---:|
| Longitudinal intercept | 73.276 | 73.303 | 73.265 | 0.038 |
| Longitudinal time slope | -0.040 | -0.041 | -0.041 | 0.001 |
| Treatment by time | -0.041 | -0.041 | -0.041 | 0.000 |
| Random intercept SD | 14.200 | 14.228 | 14.205 | 0.028 |
| Survival treatment coefficient | 0.996 | 0.999 | 0.996 | 0.003 |
| Alpha | 0.015 | 0.015 | 0.015 | 0.000 |

The intervals and substantive interpretations look similar. These differences
are compatible with ordinary Monte Carlo variation; screenshots cannot establish
that they are within MCSE or rule out convergence problems. Matching rounded
alpha values does not imply identical unrounded estimates.

Exact Stan draws are not guaranteed across CPUs, operating systems, compilers,
or Stan versions, even with the same seed. The README reports CmdStan 2.36.0
for the reference example; the production script selects 2.37.0. Document each
machine's actual environment and verify data/code/settings before attributing
all variation to hardware.

For each machine, retain unrounded means, interval endpoints, R-hat, bulk/tail
ESS, mean MCSE, and chain diagnostics. For an approximate cross-run comparison
of means use

    abs(mean_A - mean_B) / sqrt(mcse_A^2 + mcse_B^2)

as a diagnostic, not an absolute pass/fail guarantee. Quantile agreement needs
quantile MCSE; mean MCSE does not measure interval-endpoint error.

## Items to address before sending the audit package

1. **Check the treedepth claim.** README line 64 says one chain reached maximum
   treedepth on every iteration and that this did not affect results. This is
   not demonstrated by the ZIP, which contains no saved fitted diagnostics.
   Treedepth saturation concerns efficiency and warrants checking ESS, MCSE,
   chain agreement, R-hat, and divergences. Similar printed means alone do not
   establish a harmless warning. Replace the categorical claim with the actual
   diagnostic results and a qualified interpretation.
2. **Supply numerical reference diagnostics.** README line 60 refers to MCSE
   "below", but none are supplied, and fit$summary() does not include mean MCSE
   by default. Export unrounded reference summaries including mean MCSE and
   sampler diagnostics. Keep reference output separate from rerun output to
   avoid overwriting it.
3. **Correct the alpha-prior comment.** R lines 121-124 describe a 1.96-unit
   D6 contrast as giving a factor of two at the 95% prior scale. With raw SD
   log(2)/1.96, 1.96 is the Normal 97.5% quantile, not an extra D6 contrast.
   Algebraically the approximate central prior interval for alpha is
   [-log(2), log(2)], corresponding to a factor-of-two range for a one-unit
   fixed contrast. For a contrast Delta, the corresponding calibration is
   SD(alpha) = log(F)/(1.96 * abs(Delta)). This is a comment correction; the
   implemented raw prior matches production.
4. **Report runtime and environments.** Add total runtime, compilation versus
   sampling time, CPU, core use, and sessionInfo(). The existing package lists
   software versions but gives no runtime or CPU specification. Explain how
   to reduce parallel_chains on a smaller machine.
5. **Clean the ZIP.** Remove rm(list = ls()) from the distributed script.
   Exclude .Rhistory, .Rproj.user, __MACOSX files, and stale platform-specific
   binaries from the normal source workflow. Compile locally from the supplied
   source once; provide explicit instructions if optional prebuilt binaries
   are retained. The current ZIP has 148 entries, many from editor metadata.
6. **Preserve portable posterior output.** The script saves an RDS containing
   the CmdStan fit but leaves CmdStan CSV output in temporary storage. Use
   fit$save_object() or preserve the needed draws/diagnostics explicitly, and
   use a permanent output directory if distributing CSVs. A plain saveRDS()
   of an R6 fit is not the strongest guarantee that every lazily loaded
   component survives after temporary files disappear.
7. **Clarify manuscript mapping.** README and script comments say BP weights
   are not in the published table. Check against the final manuscript; the
   earlier Table 1 shared in this conversation did include theta_1 to theta_5.
8. **Keep the overall audit scope explicit.** This package addresses the mock
   case study only. The editor also requires manuscript-linked simulation code,
   intermediate production results and summarization scripts, and a documented
   small spot-check workflow with supplied reference results. The mock example
   alone cannot cover the simulation tables/figures or sensitivity/comparator
   analyses. Full cluster rescue infrastructure is unnecessary if a portable
   workflow reproduces the relevant computations.

## Sources

- Editor-supplied `Reproducible-Research-Author-Checklist-2026.pdf`, especially
  Sections 3-5, and `Checklist_for_authors_essential_V1.0.pdf`.
- Stan reproducibility guidance:
  https://mc-stan.org/docs/reference-manual/reproducibility.html
- Stan convergence/efficiency diagnostics:
  https://mc-stan.org/learn-stan/diagnostics-warnings.html
- posterior mean MCSE:
  https://mc-stan.org/posterior/reference/mcse_mean.html
- CmdStanR persistent model saving:
  https://mc-stan.org/cmdstanr/reference/fit-method-save_object.html

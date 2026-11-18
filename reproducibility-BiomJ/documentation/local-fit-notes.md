# Completed Local Checks

Checked on 17 November 2026. Results are under `output/local-simulations/`
and `output/case-study/`. Successful execution is not a declaration that
all predictive scores or sampling diagnostics are acceptable.

## Main Simulation: Configuration 1

Five original replicates completed on an Apple M2 Pro with 32 GB RAM,
R 4.5.2, CmdStanR 0.9.0 and CmdStan 2.40.0. Each fit retained 1,000 draws
from each of four chains after 1,000 warmup iterations per chain.

| Replicate | Seed | Sampling minutes | Divergences | Treedepth hits | Maximum coefficient R-hat |
|---|---|---|---|---|---|
| 1 | 1000001 | 17.30 | 0 | 0 | 1.00742 |
| 2 | 1000002 | 8.35 | 0 | 0 | 1.00950 |
| 3 | 1000003 | 28.74 | 0 | 1000 | 1.00544 |
| 4 | 1000004 | 6.51 | 0 | 0 | 1.00328 |
| 5 | 1000005 | 14.78 | 0 | 0 | 1.00368 |

Replicate 3 reached depth 10 in every retained transition of chain 4.
Its minimum coefficient bulk/tail ESS were about 571/1,155; the warning is
retained for review, not hidden or treated as an automatic failure.
All five fits warned about high Pareto-k and many pointwise WAIC variances
above 0.4. These predictive scores need separate validation and are not the
basis for declaring agreement of parameter estimates.

Matching archived-production comparisons, inputs, initial values, draws and
diagnostics are saved per replicate. Historical summaries lack some MCSEs;
no formal combined-MCSE equivalence claim is made.

## Comparator Check

One configuration-1, replicate-1 comparator check completed the Cox, Weibull
and log-logistic JM fits successfully, in approximately 107 seconds combined.
Its statuses and estimates remain in the local editor records.

## Mock Case Study

The mock example completed with seed 10235 in 14.83 minutes overall,
including 13.79 minutes of sampling. No divergences or treedepth hits occurred.
Maximum reported-coefficient R-hat was 1.00593; E-BFMI ranged from .586 to .736.
The intercept, survival coefficient and alpha means were 73.2654, .996488
and .0145264. The LMM table matches the supplied reference; the joint-model
estimates are close to the Windows references, not identical.

Inputs, four initializations, raw chain CSVs, numeric results and
software/build/hardware records are retained. This is a synthetic-data check,
not recovery of the manuscript's real ENZAMET estimates.

## Saved-Result Checks

The regenerated 16 manuscript CSVs and 16 LaTeX tables match the archived
reference files. Full regenerated parameter-summary tables for all five
collections agree with the supplied summaries within tolerance 1e-10.
See `validation-report.md` for the latest non-sampling checks.

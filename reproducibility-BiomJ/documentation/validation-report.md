# Validation Report

Recorded: 2026-11-18 12:40:00 UTC

These checks distinguish completed execution from acceptable sampling diagnostics.
EXECUTED is not an automatic declaration of convergence or statistical equivalence.

- **PASS** R syntax: 19 R scripts parsed
- **PASS** Stan source consistency: One shared Stan source; original model checksum unchanged
- **PASS** Configuration grid: 40 unique original configuration IDs
- **PASS** Merged archives: 10000, 2500, 2500, 2500, 4000
- **PASS** Available selection space: All configuration/replicate mappings verified; running all 550 candidate combinations is not required
- **EXECUTED** Recorded local checks comparator: 1 of 1 recorded dataset checkpoints completed; inspect diagnostics separately
- **EXECUTED** Recorded local checks main: 5 of 5 recorded dataset checkpoints completed; inspect diagnostics separately
- **PASS** Saved summary main/config001/rep001: Compared saved coefficient summaries with saved posterior draws
- **CHECKED** Sampling diagnostics main/config001/rep001: 0 divergences; 0 treedepth hits; max coefficient R-hat 1.0074; 17.3 sampling min
- **PASS** Saved summary main/config001/rep002: Compared saved coefficient summaries with saved posterior draws
- **CHECKED** Sampling diagnostics main/config001/rep002: 0 divergences; 0 treedepth hits; max coefficient R-hat 1.0095; 8.3 sampling min
- **PASS** Saved summary main/config001/rep003: Compared saved coefficient summaries with saved posterior draws
- **REVIEW** Sampling diagnostics main/config001/rep003: 0 divergences; 1000 treedepth hits; max coefficient R-hat 1.0054; 28.7 sampling min
- **PASS** Saved summary main/config001/rep004: Compared saved coefficient summaries with saved posterior draws
- **CHECKED** Sampling diagnostics main/config001/rep004: 0 divergences; 0 treedepth hits; max coefficient R-hat 1.0033; 6.5 sampling min
- **PASS** Saved summary main/config001/rep005: Compared saved coefficient summaries with saved posterior draws
- **CHECKED** Sampling diagnostics main/config001/rep005: 0 divergences; 0 treedepth hits; max coefficient R-hat 1.0037; 14.8 sampling min
- **PASS** Manuscript derived_csv: 16 of 16 reference files byte-identical
- **PASS** Manuscript tables: 16 of 16 reference files byte-identical
- **PASS** Reaggregate main: Full regenerated parameter table agrees within 1e-10
- **PASS** Reaggregate alpha125: Full regenerated parameter table agrees within 1e-10
- **PASS** Reaggregate alpha150: Full regenerated parameter table agrees within 1e-10
- **PASS** Reaggregate alpha300: Full regenerated parameter table agrees within 1e-10
- **PASS** Reaggregate comparator: Full regenerated parameter table agrees within 1e-10
- **PASS** Mock saved summary: Compared saved coefficient summaries with saved posterior draws
- **EXECUTED** Mock case-study fit: Inspect output/case-study diagnostics and original supplied mock reference table

## Still Requires Author Sign-Off

- Confirm final manuscript/supplement numbering and coverage of every reported result.
- Have a coauthor test the final package on another machine.
- Confirm data/code redistribution permissions and the large-archive delivery URL.
- Review all new-fit diagnostics before claiming successful numerical reproduction.

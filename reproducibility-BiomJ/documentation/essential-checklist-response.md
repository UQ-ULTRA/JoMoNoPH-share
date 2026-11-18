# Essential Checklist: Draft Author Response

Status is deliberately qualified. The authors should review and sign off before
submitting the supplied essential PDF checklist; this file is not a false claim
that every final-submission requirement has already been met.

| Essential item | Status and evidence |
|---|---|
| README and software versions | Implemented. README, `--mode=check`, per-run environment records. Historical complete compiler/hardware versions were not all retained. |
| Exact running instructions / one main file | Start with `run_checks.R` and the root README in RStudio; `main.R` is an optional terminal equivalent. |
| All code/data for results | Full simulation intermediates and relevant production/aggregation/presentation code included. Mock case-study substitute supplied; actual ENZAMET data not included. Author must confirm that the final manuscript has no additional results outside the map. |
| Outputs linked to manuscript | `manuscript_map.csv` links descriptive filenames to analysis/results. Final figure/table numbering requires checking against the final submitted PDF. |
| Seeded simulations | Implemented with original configuration IDs, per-replicate seeds, R RNG kind and chain IDs. Exact cross-platform Stan equality cannot be certified. |
| Runtime and hardware | Historical runtime summaries included; README describes Bunya EPYC3/EPYC4 and new-run metadata. New local test timings recorded when completed. |
| Parallelization explained | README documents four chains, sequential datasets and `--parallel-chains`. |
| Avoid absolute paths | Runnable local sources use package-relative paths. Original hard-coded driver paths and SLURM dispatch were removed. Recorded run metadata may retain the actual machine paths used during execution. |
| Organized files | Separate local source, archived references and regenerated output. One shared scenario generator replaces five duplicated scripts. |
| Helpful comments | Local runner and scientific stages are commented; function-reference.md describes inputs, outputs and side effects. Numerical cleanup checks are documented in code-equivalence.md. |
| Intermediate results / reduced checks | Full post-rescue collections plus a local runner supporting five replicates per selected configuration, all three alpha settings and all comparator models. Current local completion is reported explicitly in `output/local-simulations/editor-records/progress.csv`. |
| Single ZIP, no unnecessary files | PENDING author packaging/delivery. Full archive is too large for normal email; upload to an accessible repository and add URL/checksum if supplied separately. No cloud link has been invented. |

## Before Checking The Remaining Boxes

1. Inspect `validation-report.md` and complete any pending full fits or
   aggregation checks before claiming they were executed.
2. Have a coauthor run the final archive independently on another machine.
3. Confirm final manuscript/supplement output mapping and data/code permissions.
4. Resolve the final delivery location for the large archive and check that the
   editor can access it without an institutional login.
5. Fill/sign the essential PDF checklist using these verified statuses. It is
   intentionally not auto-signed or universally checked by this preparation step.

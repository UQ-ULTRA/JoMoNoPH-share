# Essential Checklist: Responses And Mark Review

Reviewed on 19 November 2026 against
`Checklist_for_authors_essential_V1.0_completed.pdf`. All twelve boxes are now
checked. The final two marks depend on supplying the final ZIP and an accessible
results download with the submission; this review does not verify delivery.
The PDF was not changed. These notes accompany, rather than replace, the
completed checklist.

| Essential item (PDF order) | Observed mark | Response / qualification |
|---|---|---|
| 1. Include a README file | Checked | README and software-environment.md document dependencies and the recorded local environment. Historical Bunya OS/compiler/package details were not all retained; complete historical environment identity is not claimed. |
| 2. Describe how to run your code | Checked | Open run_checks.R in the supplied R project. The README gives three independent workflows; main.R is an optional terminal interface to the same implementation. |
| 3. Include all code/data | Checked; mock-data qualification | All numbered simulation results in author-selected arXiv:2606.16080v1 are mapped. Restricted ENZAMET Table 1 and the PH diagnostic p-values cannot be reproduced from the public inputs. As permitted by the checklist, synthetic data and mock references exercise the case-study fitting procedure instead. See manuscript-coverage.md. |
| 4. Link outputs to manuscript results | Checked | manuscript_map.csv gives exact table/figure numbers and pages. Nine simulation table bodies match the source and three figures render identically with the same renderer. Table typography, captions and unavailable-value markers can differ. |
| 5. Ensure simulations are reproducible | Checked; Stan qualification | Original configuration/replicate seeds, R RNG kinds and chain IDs are retained. Controlled generator/input checks match on the recorded environment. Numerical agreement is assessed with diagnostics and Monte Carlo uncertainty; bitwise-identical MCMC across different hardware/toolchains is not promised. |
| 6. Report runtime and hardware | Checked | README, software-environment.md and recorded runtime summaries distinguish historical Bunya EPYC3/EPYC4 runs from local audit checks. |
| 7. Explain parallelisation | Checked | README explains parallel_chains, four total chains, sequential datasets and single-process aggregation. |
| 8. Avoid absolute paths | Checked | Executable input/output paths are package-relative. Recorded provenance and environment files may contain the actual historical machine paths; these are not paths the editor must edit. |
| 9. Organise code files | Checked | Sources, archived results, mock references and newly generated outputs are separate. A shared generator/coefficient table serves all five scenarios. |
| 10. Add helpful comments | Checked | Entry scripts and scientific stages are commented; manuscript_map.csv, function-reference.md and manuscript-coverage.md link calculations and outputs to their purpose. |
| 11. Provide intermediate results | Checked; delivery qualification | Complete post-rescue production collections and selected local reference results are prepared for separate sharing. Supply the editor-accessible download link with the submission; link access was not verified in this review. |
| 12. Single ZIP without unnecessary files | Checked; packaging qualification | Supply the final clean code ZIP together with the separately hosted large archive. This mark records the intended submission contents; ZIP creation and delivery were not verified in this review. |

## Before Submission

- Finalise the ZIP, external archive access and selected local reference records.
  The detailed checklist's page 3 footnote specifically excludes OneDrive;
  use a suitable repository or agree a different transfer method with the editor.
- A coauthor (PM) independently ran Workflow 1; the supplied
  `case-study/reference/mock_table1_update_Windows_PM.png` records the mock
  estimates. This is case-study evidence, not an independent check of the full
  simulation workflows or their sampling diagnostics. The author's simulation
  and manuscript-regeneration checks are documented separately.
- Confirm redistribution permissions and whether the journal revision differs
  from the mapped arXiv version. Retain the restricted-data and Stan qualifications
  in the submission.

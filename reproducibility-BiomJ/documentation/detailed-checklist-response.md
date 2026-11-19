# Detailed Checklist: Responses And Mark Review

Reviewed on 19 November 2026 against all six rendered pages of
`Reproducible-Research-Author-Checklist-2026_completed.pdf`. The PDF marks were
not changed. There are 34 checked and four unchecked boxes. Locations below
use the printed page numbers. These notes record observed marks, evidence and
qualifications; delivery is not verified merely by a checkmark.

## Unchecked Boxes

| Page / item opening | Current mark | Response / recommended action |
|---|---|---|
| 2, RNG: "By re-running any code... multiple times with different seeds..." | Unchecked | Appropriate. Production uses 1,000 independently seeded datasets per configuration, but an additional full-study rerun demonstrating unchanged rankings/conclusions has not been established. No such result is claimed. |
| 3, restricted data: "The original data has been made available..." | Unchecked | Appropriate: restricted ENZAMET data are not supplied. The preceding synthetic/pseudo-data option is checked and is the chosen alternative. |
| 4, documentation: "contains full documentation for all data sets including their provenance..." | Unchecked | Data provenance/dictionary are supplied, but copyright, redistribution permissions and final licensing need author/coauthor confirmation. |
| 6, compiled code: "The supplement includes both the documented source code..." | Unchecked | Explain the exception: Stan source and local CmdStan/toolchain instructions are supplied, not machine-specific precompiled executables. Do not claim the executable requirement is fulfilled. |

## Checked Boxes Requiring Qualification

| Page / item opening | Observed mark | Qualification / recommendation |
|---|---|---|
| 1, Section 2: "We have verified that re-running..." | Checked | Author-run simulation reaggregation and manuscript-output checks passed; PM independently ran Workflow 1, recorded in mock_table1_update_Windows_PM.png. Independent checking covers the mock case study, not all workflows. Restricted real-data results have a mock-data substitute, and cross-platform MCMC identity is not claimed. |
| 1, Section 2: "The code and data supplement has been uploaded..." | Checked; delivery qualification | Supply the final clean code ZIP and external results link with the submission. ZIP creation, upload and recipient access were not verified in this review. |
| 1, Section 3: "The code and data supplement contains all code and data..." | Checked | All numbered simulation results in arXiv:2606.16080v1 are mapped. Real-data Table 1 and Section 4's PH p-values require restricted ENZAMET data; mock data reproduce the fitting procedure and mock references only. See manuscript-coverage.md. |
| 2: "They create dataframes or tables with identical structure..." | Checked | All nine simulation table bodies match at displayed precision. Presentation differs in font size, captions/model labels and NA versus dash markers; identical PDF layout is not claimed. Exact numbers/pages are in manuscript_map.csv. |
| 2, RNG: "Any code whose output relies on... RNG..." | Checked | Seeds are retained, but the box also asserts exact identity to manuscript results. Add an explicit Stan cross-platform qualification or leave it unchecked with this explanation; do not present the current tick as unconditional exact numerical identity. |
| 3, long computations: "The code supplement provides intermediate results..." | Checked; delivery qualification | Production results, local reference results and regeneration code are prepared for sharing. Supply an accessible download link; access was not verified here. The page 3 footnote excludes OneDrive, so use a suitable repository or obtain the editor's agreement to another transfer route. |
| 3, spot checks: "We have verified that the code and intermediate results... enable spot checks..." | Checked | Selected original configurations/replicates are independently runnable and paired with archived references. The sentence's exact-identity claim needs the same Stan qualification. Use uncertainty and diagnostics, not bitwise equality, to assess different environments. |
| 3, software: "contains version information for all the software..." | Checked | The local session and per-fit build/environment records are supplied. Complete historical Bunya OS/compiler/package versions were not retained. No Docker/renv environment or recovery of missing historical versions is claimed. |
| 4, style: "All code files use a consistent, sensible format and style..." | Checked | The stated maximum of 80-100 characters is not met throughout: several R/Stan files retain longer lines. Explicitly qualify this criterion or leave it unchecked; strict line-length compliance cannot be certified. No formatting changes were made in this review. |
| 5: "Code files do not contain extraneous or outdated code..." | Checked | The manuscript generator retains additional unpublished plots/tables and the comparator includes extra fitted models. These are documented research/audit extras, not additional published results. Qualify the literal claim that every function/output is only for published results, or leave it unchecked. |
| 5: "Code files do not contain commands that affect users' global workspace..." | Checked | No automatic package installation or whole-workspace clearing is used. Tasks do temporarily set process environment variables/options/RNG controls and restore the recorded task settings on exit; the literal assertion that no variables are modified needs this explanation. |

## Supporting Evidence And Remaining Sign-Off

README describes the three workflows, editable settings, package-relative paths,
runtimes and hardware. The parallelisation item on page 4 is now checked,
consistent with the documented settings. The shared generator avoids near-duplicate scenario
scripts. Function interfaces and data fields are documented in
function-reference.md and data_dictionary.md. The RStudio project is optional;
the same scripts can run through R directly. Platform-specific environment
queries are guarded, not model inputs.

Source/input cleanup checks and existing local fits are documented in
code-equivalence.md and validation-report.md. Nine simulation table bodies and
all three published figures were compared with the arXiv source; see
manuscript-coverage.md. These checks do not establish exact cross-platform
sampling or independent testing of the final submission ZIP. A coauthor (PM)
has independently run Workflow 1; its mock estimate screenshot is included as
case-study evidence, not a full-workflow or convergence assessment.

Retain the exact-identity, line-length and published-code-only qualifications
above when submitting these checked forms. Permissions and archive delivery
remain to be confirmed; independent testing currently covers Workflow 1 only.
The copied PDFs were deliberately left unchanged.

# Detailed 2026 Checklist: Verification And Remaining Author Actions

This is evidence for completing the journal checklist, not a signed declaration
of universal compliance. Section numbers refer to the supplied six-page
`Reproducible-Research-Author-Checklist-2026.pdf`.

| Item | Evidence / status |
|---|---|
| Sections 2-3: code/data and intermediate results | Five complete saved simulation collections; source and three workflows supplied. Mock data substitute is explicitly identified. Final manuscript completeness still needs author confirmation. |
| Repeated settings without copied scenario scripts | Shared simulator/coefficient table and local task runner; arbitrary configuration and replicate selection. |
| Output/result association | manuscript_map.csv links generated outputs and references. Exact final manuscript/supplement numbering remains to be confirmed. |
| RNG and selected-replicate reproduction | Original configuration/replicate seed rule and RNG kinds retained. Controlled cleanup checks passed; see code-equivalence.md. Cross-platform MCMC equality is not claimed. |
| Monte Carlo stability across seeds | Archived studies use 1,000 independently seeded datasets per configuration. No new full-study reruns or comprehensive method-ranking stability study were performed during cleanup; author must assess/report this item honestly. |
| Restricted case-study data | Mock data, code, reference results and a completed local mock run supplied. No actual ENZAMET participant data included. |
| Long computations / spot checks | Archived intermediates regenerate results without sampling; selected IDs run independently, normally five replicates per configuration. |
| Section 4: software documentation | README links software-environment.md, including loaded-package sessionInfo; per-fit compiler/build/library records supplied. Unrecorded historical versions are explicitly identified, not invented. |
| Running instructions and manual edits | One R settings block selects tasks and datasets; no repeated edits to numerical algorithms. Optional terminal interface uses the same engine. |
| Data provenance / licensing | Data dictionary supplied. Copyright holders, redistribution permissions and final license terms require author/coauthor confirmation. |
| Parallelization / runtimes / hardware | README and local-fit-notes.md; four chains, configurable concurrency, sequential datasets. |
| Section 5: organization / function documentation | Clean local sources, shared generator and function-reference.md. Source formatting changes do not alter numerical expressions. |
| Workspace / software installation | No workspace-clearing command in the runnable sources; no automatic installation. Temporary task options/environment/RNG/parallel state restored. |
| Paths / platform dependence | Package-relative inputs and separate outputs; no RStudio API requirement. Platform-specific metadata queries are guarded and do not control inference. Independent second-platform testing is still required. |
| Compiled code | Stan source and local CmdStan/toolchain instructions supplied. Confirm with the editor whether locally compiled Stan models satisfy the checklist's precompiled-executable wording; platform-specific binaries are not supplied. |
| Single ZIP / large files | Pending final packaging, external archive URL/checksum and unauthenticated editor access. |
| Independent reproduction | Pending a coauthor's test of the final package on another machine and final comparison/reporting of diagnostics. |

## Do Not Certify Yet

Final manuscript coverage, licensing, large-archive access and independent
reproduction are author-controlled steps. Missing historical environment details
and Stan's cross-platform numerical limitations should be disclosed to the RR
editor rather than represented as exact reproducibility guarantees.

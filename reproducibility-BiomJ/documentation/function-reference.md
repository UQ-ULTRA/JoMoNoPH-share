# Function Interfaces

This reference lists the actual formal arguments and outputs of top-level
functions. Defaults are visible in their source signatures. Analysis helpers
run inside their script's private environment and use its prepared tables,
labels and output directories; they are not separate entry points.

## Argument Conventions

- Paths: root is the package folder; out/output_dir/fig_dir/folder are output or
  input folders as specified by the caller; file/path identify a file; stem is a filename prefix.
- Selection: collection is a named simulation family; config_ids/config_hash identify
  original grid rows; i/replicates/sim_id/idx identify original replicates/rows;
  row is one plan row; scenario is the integer coefficient setting 1:5.
- Fitting: fit/fit_long/fit_surv are Stan/JMbayes2, nlme::lme and survreg objects;
  draws/chain_draws/D_mat contain posterior draws; center/longitudinal_response_center
  is the observed response mean; surv_data is the patient-level survival data frame.
- Generation: D is the 2x2 random-effect covariance; beta_0/beta_1/beta_2 are
  longitudinal fixed effects; sigma_e is residual SD; alpha_AFT is association;
  log_AF is treatment log time ratio. log_HR/alpha_PH are unused legacy call
  arguments. Baseline shapes/scales and cfg_* parameters set the chosen family;
  visit/max_FU/fu are month-scale times, n_patients/n sample size, lambda_c/cens
  censoring controls. aft_mode/dist/baseline_label name the baseline family;
  link_type is value. seed/unique_seed give the original RNG seed. k_bases is BP count.
- Execution: task/mode select an operation; parallel_chains controls simultaneous
  chains; cmdstan is an installed path; no_pch/prepare/resume/allow_all/quiet/
  continue_on_error/comparator/preparation_only are logical switches. summaries
  chooses saved or regenerated tables. saved/current/reference/reference_result/
  result/fresh/sim are checkpoint or replicate lists; source_aliases maps verified source hashes.
- Tables: x/y/df/per_sim/all_per_sim/runtime_df/prior_info_df are the scalar,
  vector or data-frame inputs appropriate to the named helper; parsed/config_info
  contain configuration metadata. component/model/model_name/source/parameter/
  variable name model components and parameters. pattern/prefix select columns;
  diag_id selects a covariance diagonal; n/upper bound lengths or IDs; value/name/
  default specify environment or validation inputs; args/packages are command arguments/package names.
- Formatting: config/configs/config_str/cfg/lbl/label/code/cens/config_label/cm_label/
  censoring_label/censoring_value describe configurations/censoring. digits/column_digits
  control precision; caption/note/super_header/na_text are display text; *_labels
  set display group order/text; align/table_env/font_size control LaTeX layout;
  plot is a ggplot object and width/height are inches; group_vars define grouping;
  maxF identifies alpha-prior width; file_name is an input basename; aft_only
  requests historical AFT-only truth labels; log_lik_var names the Stan score variable.

## R/audit.R

| Function | Inputs (formal names) | Output / side effects |
|---|---|---|
| `audit_collections` | root | Collection/archive/summary-path data frame; no writes |
| `audit_parameter_summary` | draws, center | Posterior summary with original-scale intercept and MCSEs |
| `audit_require` | packages | Stops for missing packages; no installation |
| `audit_command` | command, args | Captured command output/error text; read-only system query |
| `audit_environment` | out | Writes session, package, compiler and hardware records under out |
| `audit_load_production` | root, comparator | Private environment of local fitting functions/settings |
| `audit_grid` | root | Fixed 40-row configuration data frame |
| `audit_model` | root, out, parallel_chains, prepare | Recording model factory; writes inputs/build/draws/diagnostics under out |
| `audit_compare_reference` | root, out, reference_result | Writes regeneration/initial-value/estimate comparisons; stops on staging mismatch |
| `audit_case_study` | root, parallel_chains | Fits mock data; writes case-study table/draws/diagnostics |
| `audit_analysis_environment` | None | Private script environment with explicit dplyr aliases; no writes |
| `audit_results` | root, summaries | Writes manuscript tables/figures from summary files |
| `audit_aggregate` | root, collection | Writes collection summary tables and diagnostic plots |
| `audit_check` | root, cmdstan | Returns PASS/WARNING table; writes check-summary.csv and available inventories/environment/grid; continues after missing components and emits one concluding warning |
| `audit_manifest` | root | Writes file-size/checksum CSV |

## R/entry.R

| Function | Inputs (formal names) | Output / side effects |
|---|---|---|
| `audit_run` | root, task, collection, config_ids, replicates, parallel_chains, cmdstan, no_pch, resume, summaries | Dispatches one task; returns its result invisibly where provided |

## R/export-local.R

| Function | Inputs (formal names) | Output / side effects |
|---|---|---|
| `audit_write_csv` | x, file | Atomically writes x to file |
| `audit_record_one` | root, row, out, result | Writes annotated parameter/runtime/status/diagnostic records for one fit |
| `audit_export_local` | root, quiet | Writes merged CSV/RDS records; returns progress counts |

## R/session.R

| Function | Inputs (formal names) | Output / side effects |
|---|---|---|
| `audit_session` | None; reads scoped settings/script tables where noted | Restoration callback for temporary environment/options/RNG/parallel settings |

## R/simulation.R

| Function | Inputs (formal names) | Output / side effects |
|---|---|---|
| `audit_parse_ids` | value, name, upper, allow_all | Validated positive integer IDs, or NULL for all |
| `audit_simulation_plan` | root, collection, config_ids, replicates | One row per selected collection/configuration/replicate |
| `audit_plan_reference` | root, row | Matching archived replicate result list |
| `audit_comparator_check` | out, reference, fresh | Writes comparator preparation/estimate comparisons; stops on center mismatch |
| `audit_checkpoint_differences` | saved, current, preparation_only, comparator, source_aliases | Character vector of incompatible identity field names |
| `audit_simulations` | root, collection, config_ids, replicates, mode, parallel_chains, resume, continue_on_error | Writes plans and selected fit/preparation records; returns plan |

## R/validate.R

| Function | Inputs (formal names) | Output / side effects |
|---|---|---|
| `audit_validate` | root | Writes consistency report; stops on failed checks without changing fits |

## analysis/generate_manuscript_outputs.R

| Function | Inputs (formal names) | Output / side effects |
|---|---|---|
| `%||%` | x, y | x unless missing/empty, otherwise y |
| `read_summary_csv` | path | Data frame loaded from a saved summary CSV |
| `safe_mean` | x | Mean of nonmissing values, or NA |
| `safe_median` | x | Median of nonmissing values, or NA |
| `safe_min` | x | Minimum nonmissing value, or NA |
| `safe_max` | x | Maximum nonmissing value, or NA |
| `latex_escape` | x | LaTeX-escaped character vector |
| `format_number` | x, digits | Rounded display strings with missing-value handling |
| `format_table` | df, digits, column_digits | Formatted display table |
| `write_latex_table` | df, path, caption, label, digits, column_digits, align, table_env, font_size | Writes a labelled LaTeX table |
| `parameter_math_label` | variable | Plotmath/LaTeX parameter-label string |
| `write_alpha_prior_estimate_table` | df, path, caption, label, prior_labels, cm_labels, prior_display_labels, cm_display_labels, super_header, na_text, note | Writes grouped alpha-prior estimate LaTeX table |
| `write_distribution_model_estimate_table` | df, path, caption, label, distribution_labels, cm_labels, model_labels, distribution_display_labels, cm_display_labels, model_display_labels, na_text, note | Writes grouped distribution/model estimate LaTeX table |
| `write_alpha_prior_detail_table` | df, path, caption, label, prior_labels, metric_labels, censoring_label, note, prior_display_labels | Writes grouped alpha-prior performance LaTeX table |
| `write_csv` | df, path | Writes a data frame to CSV |
| `clean_config` | config | Cleaned display configuration string |
| `parse_dist_code` | code | Baseline distribution metadata |
| `parse_config` | configs | Parsed sample-size/censoring/baseline/scenario/follow-up metadata |
| `parse_short_config` | configs | Parsed abbreviated configuration metadata |
| `parameter_label` | variable | Human-readable parameter-label string |
| `parameter_plot_label` | variable | Plotmath parameter-label string |
| `summarise_config_rows` | df, group_vars | Grouped manuscript performance/estimate table |
| `saft_bp_distribution_detail_note` | cm_label | LaTeX table-note text for censoring mechanism |
| `make_saft_bp_distribution_detail_rows` | censoring_value | Distribution-specific manuscript performance rows |
| `bp_lmm_parameter_math_label` | variable | BP/LMM parameter-label string |
| `longitudinal_ll_detail_note` | cm_label | LaTeX note for longitudinal log-logistic comparison |
| `make_longitudinal_ll_detail_rows` | censoring_value | Log-logistic longitudinal performance rows |
| `read_alpha_sensitivity` | folder, maxF, label, file_name | Alpha-prior summary data frame with prior labels |
| `alpha_prior_detail_note` | cm_label | LaTeX note for alpha-prior comparison |
| `make_alpha_prior_detail_rows` | censoring_value | Alpha-prior performance rows for selected censoring mechanism |
| `save_figure` | plot, stem, width, height | Writes the manuscript figure at the requested dimensions |

## analysis/SumSimRes-ENZAMET-DM.R

| Function | Inputs (formal names) | Output / side effects |
|---|---|---|
| `format_cens_short` | cens | Abbreviated censoring label |
| `format_cens_label` | cens | Display censoring label |
| `format_cens_table` | cens | Manuscript censoring label |
| `parse_baseline_label` | label | Baseline family/shape/scale metadata |
| `truth_values_for_scenario` | scenario, aft_only | Named data-generating coefficients/SDs for the scenario |
| `shorten_config` | cfg | Abbreviated configuration label |
| `parse_config` | config_str | Parsed sample-size/censoring/baseline/scenario/follow-up metadata |
| `parse_batch_number` | path | Batch number parsed from filename |
| `summarise_runtime_table` | runtime_df, config_label | Replicate/model runtime summary data frame |
| `empty_prior_info_per_sim` | None; reads scoped settings/script tables where noted | Zero-row prior-information table with expected columns |
| `safe_prior_scalar` | x | First finite numeric scalar, or NA |
| `prior_vector_cols` | x, prefix, n | Named list of indexed prior-vector values |
| `extract_prior_info_row` | sim, idx | One replicate's prior-information row |
| `summarise_prior_info_table` | prior_info_df, config_label | Configuration-level prior-information summary table |
| `parse_short_config` | lbl | Parsed abbreviated configuration metadata |
| `process_simulation_results` | bf, k_bases, n, cens, baseline_label, dist, scenario, fu, output_dir | Main-model parameter/diagnostic/runtime/prior/failed-replicate tables |
| `shorten_config_1` | config | Abbreviated display configuration label |

## analysis/SumSimRes-ENZAMET-JMbayes2-DM.R

| Function | Inputs (formal names) | Output / side effects |
|---|---|---|
| `env_optional_int` | name | Integer environment setting or NA |
| `env_bool` | name, default | Logical environment setting |
| `parse_batch_number` | path | Batch number parsed from filename |
| `parse_baseline_label` | label | Baseline family/shape/scale metadata |
| `parse_config` | config_str | Parsed sample-size/censoring/baseline/scenario/follow-up metadata |
| `truth_values_for_scenario` | scenario | Named data-generating coefficients/SDs for the scenario |
| `truth_for_jmbayes2_summary` | variable, source, parsed, scenario | Relevant longitudinal truth or NA for noncomparable estimands |
| `format_cens_short` | cens | Abbreviated censoring label |
| `shorten_config` | cfg | Abbreviated configuration label |
| `empty_tidy` | None; reads scoped settings/script tables where noted | Zero-row common-format comparator table |
| `empty_initial_survreg_tidy` | None; reads scoped settings/script tables where noted | Zero-row preliminary-survival estimate table |
| `empty_initial_survreg_summary` | None; reads scoped settings/script tables where noted | Zero-row preliminary-survival performance table |
| `standardise_summary_table` | x, component, model, sim_id | Long parameter-summary data frame |
| `map_parameter_to_common` | component, parameter | Common manuscript-compatible parameter name |
| `normalise_saved_jmb2_tidy` | x, model_name, sim_id | Common-format saved comparator summary table |
| `tidy_jmb2_from_raw_summary` | sim, sim_id | Common-format comparator table reconstructed from raw summaries |
| `repair_comparison_tidy` | sim, sim_id | Common-format table with missing available summaries repaired |
| `tidy_initial_survreg_from_sim` | sim, sim_id | Common-format stored preliminary survival estimates |
| `truth_for_initial_survreg` | variable, model, parsed, scenario | Preliminary survival estimand truth where defined, otherwise NA |
| `summarise_initial_survreg` | per_sim, parsed | Preliminary survival-model bias/RMSE/coverage summary |
| `process_jmbayes2_config` | config, output_dir | Comparator estimate/performance/status/runtime tables for a configuration |
| `write_outputs` | x, stem | Writes comparator summary CSV/RDS files |
| `generate_jmbayes2_figures` | all_per_sim, config_info, fig_dir | Writes comparator diagnostic/performance figures |

## production/Bunya/ENZAMET/Data-Gen-ENZAMET.R

| Function | Inputs (formal names) | Output / side effects |
|---|---|---|
| `simulation_parameters` | scenario | Coefficient/design list for scenario 1:5 |
| `simulate_joint_dataset` | D, beta_0, beta_1, beta_2, sigma_e, log_HR, alpha_PH, log_AF, alpha_AFT, weibull_shape, weibull_scale, loglogistic_shape, loglogistic_scale, visit, seed, n_patients, max_FU, lambda_c, aft_mode, link_type | List of longitudinal and survival data frames; sets the supplied seed |

## production/Bunya/ENZAMET/ENZAMET-DM.R

| Function | Inputs (formal names) | Output / side effects |
|---|---|---|
| `compute_waic_from_log_lik` | fit, log_lik_var | loo::waic object from saved pointwise log likelihoods |
| `alpha_prior_calibration` | None; reads scoped settings/script tables where noted | Scoped alpha contrast/max-factor/SD list |
| `convert_data_from_models` | fit_long, fit_surv, surv_data, k_bases | Stan data plus fixed prior and initial-value inputs |
| `run_simulation` | i, n_patients, lambda_c, aft_mode, max_FU, config_hash, joint_model_type, k_bases, scenario, baseline_label, cfg_weibull_shape, cfg_weibull_scale, cfg_loglogistic_shape, cfg_loglogistic_scale | One dataset fit: coefficient summaries, diagnostics/statuses, priors/starts and runtimes |

## production/Bunya/ENZAMET/ENZAMET-JMbayes2-DM.R

| Function | Inputs (formal names) | Output / side effects |
|---|---|---|
| `env_int` | name, default | Validated positive integer environment setting |
| `env_nonnegative_int` | name, default | Validated nonnegative integer environment setting |
| `parse_jmb2_models` | value | Canonical JMbayes2 model-name vector |
| `standardise_summary_table` | x, component, model, sim_id | Long parameter-summary data frame |
| `as_mcmc_chain_matrices` | x | List of numeric draw matrices, one per chain |
| `find_D_diag_column` | D_mat, diag_id | Column identifier for a random-effect covariance diagonal |
| `summarise_draws` | draws, chain_draws | Mean/SD/interval/ESS/R-hat draw summary |
| `jmb2_random_sd_summary` | fit, model, sim_id | Random-effect SD posterior summary derived from covariance draws |
| `map_parameter_to_common` | component, parameter | Common manuscript-compatible parameter name |
| `tidy_jmb2_summary` | fit, model, sim_id, longitudinal_response_center | Long JMbayes2 parameter-summary table on original response scale |
| `tidy_lme_summary` | fit_long, sim_id, longitudinal_response_center | Long LMM coefficient/SD summary on original response scale |
| `empty_initial_survreg_tidy` | None; reads scoped settings/script tables where noted | Zero-row preliminary-survival estimate table |
| `tidy_initial_survreg` | fit_surv, model, sim_id | Long preliminary survival-fit coefficient/shape summary |
| `compact_jmb2_draws` | fit, pattern | Selected posterior draws only |
| `save_full_fit_if_requested` | fit, model, unique_seed, sim_id | Optional saved-fit path, or NULL; only writes when requested |
| `fit_one_jmb2_model` | model, fit_long, surv_data, unique_seed, sim_id, longitudinal_response_center | Model fit/status/summary/draw/runtime list, retaining failures |
| `run_simulation` | i, n_patients, lambda_c, aft_mode, max_FU, config_hash, k_bases, scenario, baseline_label, cfg_weibull_shape, cfg_weibull_scale, cfg_loglogistic_shape, cfg_loglogistic_scale | One dataset fit: coefficient summaries, diagnostics/statuses, priors/starts and runtimes |


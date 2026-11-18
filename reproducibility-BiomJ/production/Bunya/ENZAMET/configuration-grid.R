# Fixed production grid: four baselines x five coefficient scenarios x two censoring settings.

baseline_grid <- data.frame(
  baseline_label = c("LL120_230", "W090_380", "W130_380", "W210_380"),
  aft_mode = c("loglogistic", "Weibull", "Weibull", "Weibull"), weibull_shape = c(
    NA_real_,
    0.9, 1.3, 2.1
  ), weibull_scale = c(NA_real_, 38, 38, 38), loglogistic_shape = c(
    1.2,
    NA_real_, NA_real_, NA_real_
  ), loglogistic_scale = c(23, NA_real_, NA_real_, NA_real_),
  stringsAsFactors = FALSE
)

config_base <- merge(
  expand.grid(n_patients = c(1100), k_bases = c(5), scenario = c(
    1, 2, 3,
    4, 5
  ), max_FU = 120, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE), baseline_grid,
  by = NULL
)

lambda_c_table <- data.frame(
  baseline_label = rep(baseline_grid$baseline_label, each = 5),
  scenario = rep(1:5, times = nrow(baseline_grid)), lambda_c = c(
    0.00994873, 0.00390625,
    0.004272461, 0.017028809, 0.017211914, 0.009643555, 0.003479004, 0.00378418, 0.017211914,
    0.017456055, 0.008666992, 0.002227783, 0.002593994, 0.015563965, 0.015716553, 0.008056641,
    0.001098633, 0.001464844, 0.014099121, 0.014160156
  ), stringsAsFactors = FALSE
)

cfg_default <- merge(config_base, lambda_c_table,
  by = c("baseline_label", "scenario"), all.x = TRUE,
  sort = FALSE
)

cfg_default$cens_label <- "50"

cfg_admin <- transform(config_base, lambda_c = -1, cens_label = "admin")

config_grid <- rbind(cfg_default, cfg_admin)

config_grid <- config_grid[order(
  match(config_grid$baseline_label, baseline_grid$baseline_label),
  config_grid$scenario, match(config_grid$cens_label, c("50", "admin"))
), ]

row.names(config_grid) <- NULL

config_grid$full_config_id <- seq_len(nrow(config_grid))

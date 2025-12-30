# ENZAMET Joint Modelling Simulation – Short Guide

## Purpose

This repository contains code to simulate and fit joint longitudinal–survival models for ENZAMET-based scenarios. The primary aim is to compare Bayesian joint models that differ only in their **baseline hazard representation**, under a common data-generating mechanism.

The code is designed for batch execution on HPC systems using SLURM.

---

## Main Components

### Data Generation

**Files**
- `Data-Gen-ENZAMET-01.R` to `Data-Gen-ENZAMET-07.R`

These scripts simulate individual-level trial data:

- Continuous longitudinal outcome (linear mixed model)

- Time-to-event outcome linked to the longitudinal process

Event times are generated using AFT-style models (Weibull or log-logistic). A PH formulation exists for completeness but is not used in standard runs.

---

### Simulation and Inference Driver

**File**
- `I2a-ENZAMET.R`

This script:

1. Selects an ENZAMET scenario

2. Simulates joint longitudinal–survival data

3. Fits classical models (`lme`, `survreg`)

4. Converts fitted models to Stan data

5. Fits a Bayesian joint model (Stan)

6. Saves results per batch as `.rds`

The script is intended to be run as a SLURM array job.

---

### Bayesian Joint Models (Stan)

**Files**

- `Bernstein-Polynomials-JM-Hist.stan`

- `Gaussian-Basis-JM-Hist.stan`

- (optional) Gaussian Process variants

All Stan models share the same structure and differ only in how the baseline hazard is represented.

---

## How to Run

### Local testing

For development or debugging, `I2a-ENZAMET.R` can be run locally by manually setting:

- `joint_model_type`

- `scenario`

- `k_bases`

- `lambda_c`

Reduce sample size and number of simulations for local runs.

### HPC execution
Production runs are launched as SLURM array jobs. Each array task runs one batch of simulations for a specific configuration. Output files are written to the configured output directory and skipped if already present.

---

## Common Configuration Changes

Most changes can be made in `I2a-ENZAMET.R`:

- Baseline hazard model: `joint_model_type`

- Number of basis functions: `k_bases`

- ENZAMET scenario: `scenario`

- Censoring rate: `lambda_c`

- MCMC settings: `iter_warmup`, `iter_sampling`, `chains`

---

## Outputs

Each batch produces one `.rds` file containing:

- Joint model posterior summaries

- Convergence diagnostics

- Classical LME coefficients

- Runtime information

- LOO and WAIC estimates (when available)

These files are intended for downstream aggregation and analysis.

---

## Notes

- Data preparation and prior calibration are handled in R; Stan models assume fully prepared inputs.

- The PH formulation is retained for flexibility but is not used in standard ENZAMET runs.

- The code prioritizes reproducibility and consistency across configurations.




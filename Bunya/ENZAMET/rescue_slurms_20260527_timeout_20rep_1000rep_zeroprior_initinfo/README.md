# Timeout Rescue Slurms

These slurms rescue the 20 null/time-limited replicates from the previous
40-config x 1000-replicate ENZAMET alt-DM gamma5 zero-prior/init-info run.

The rescue jobs rerun the original full 4-replicate batch files rather than
writing single-replicate `.rds` files. This is intentional: it preserves the
original batch numbering and list length expected by `SumSimRes-ENZAMET-DM.R`.
Successful rescue jobs overwrite the incomplete `.rds` files in:

```text
$HOME/JoMoNoPH-share/Bunya/ENZAMET/output_20260525_ENZAMET-altDM-gamma5-1000rep-zeroprior-initinfo
```

Submit both files:

```bash
sbatch ENZAMET-DM-altDM-gamma5-1000rep-zeroprior-initinfo-timeout-rescue-standard.slurm
sbatch ENZAMET-DM-altDM-gamma5-1000rep-zeroprior-initinfo-timeout-rescue-W210-scen3-batch171.slurm
```

Settings:

- Standard rescue: 18 original batches, `--time=10:00:00`.
- W210 scen3 batch171 rescue: 1 original batch with two timed-out reps, `--time=15:00:00`.
- Both use `SIM_SAFE_ELAPSED=24600` seconds (`6h50m`) per replicate.
- Both use `SIM_BATCH_SIZE=4`, `SIMS_PER_CONFIG=1000`, `8` CPUs, and `8GB` RAM.

# Alpha-Sensitivity LL Rescue Slurms

This folder contains the rescue SLURMs for the remaining problematic batches
from the LL120_230 alpha-prior sensitivity simulation.

Submit both:

```bash
sbatch ENZAMET-DM-altDM-gamma5-alphaSens-LL-1000rep-rescue-standard.slurm
sbatch ENZAMET-DM-altDM-gamma5-alphaSens-LL-1000rep-rescue-highTime.slurm
```

Each array task reruns one original 4-replicate batch and writes to the original
alpha-specific output folder:

```text
$HOME/JoMoNoPH-share/Bunya/ENZAMET/output_20260528_ENZAMET-alphaSens-LL-${ENZ_ALPHA_LABEL}-altDM-gamma5-1000rep-zeroprior-initinfo
```

Standard-time rescue targets (`4h50m` per replicate, `10h` walltime):

| Array task | Alpha | Original global sim_id | Config/batch |
|---:|---|---:|---|
| 1 | maxF125 | 522 | scen2, cens50, batch022 |
| 2 | maxF125 | 1261 | scen3, admin, batch011 |
| 3 | maxF125 | 1264 | scen3, admin, batch014 |
| 4 | maxF125 | 1398 | scen3, admin, batch148 |
| 5 | maxF300 | 742 | scen2, cens50, batch242 |
| 6 | maxF300 | 861 | scen2, admin, batch111 |
| 7 | maxF300 | 902 | scen2, admin, batch152 |
| 8 | maxF300 | 1110 | scen3, cens50, batch110 |
| 9 | maxF300 | 1123 | scen3, cens50, batch123 |

High-time rescue targets (`5h50m` per replicate, `12h30m` walltime):

| Array task | Alpha | Original global sim_id | Config/batch |
|---:|---|---:|---|
| 1 | maxF125 | 1262 | scen3, admin, batch012 |
| 2 | maxF125 | 1315 | scen3, admin, batch065 |
| 3 | maxF150 | 1047 | scen3, cens50, batch047 |
| 4 | maxF300 | 1482 | scen3, admin, batch232 |

Settings:

- Standard rescue: `--time=10:00:00`, `SIM_SAFE_ELAPSED_OVERRIDE=17400`
  seconds (`4h50m`) per replicate.
- High-time rescue: `--time=12:30:00`, `SIM_SAFE_ELAPSED_OVERRIDE=21000`
  seconds (`5h50m`) per replicate.
- Both use `--cpus-per-task=8` and `--mem=8GB`.
- `SIM_BATCH_SIZE=4`
- `SIMS_PER_CONFIG=1000`
- `SIM_BASELINE_FILTER=LL120_230` inherited from the common launcher

The `maxF125` scen2/cens50 batch022 target rescues the missing file whose
rerun log showed a `serverSocket` failure. The other targets overwrite saved
files that contain one or two `NULL` replicates.

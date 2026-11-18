# Bunya Simulation Setup Guide

Last updated: 2026-05-08

## File format recommendation

Use Markdown (`.md`) for this guide. The document is mostly setup notes, shell commands, R commands, and Slurm reminders, so `.md` is easier to read in GitHub, onBunya, VS Code, and plain terminals. Use `.Rmd` only if the guide needs to execute R chunks and generate a report.

If copying commands into Apple Notes, use `Cmd + Shift + M` to create code-style blocks.

## Scope

This guide records the known-working Bunya setup for running the ENZAMET simulation scripts in this repository. It combines project-specific notes with the official UQ RCC Bunya user guide.

Examples below use `uqdma8` as the Bunya username. Replace it with your own Bunya username if adapting the guide for another user.

Official references:

- UQ RCC Bunya user guide: <https://github.com/UQ-RCC/hpc-docs/blob/main/guides/Bunya-User-Guide.md>
- CmdStanR installation notes: <https://github.com/stan-dev/cmdstanr>

## Access options

### Option 1: onBunya web portal

Open the onBunya dashboard:

<https://onbunya.rcc.uq.edu.au/pun/sys/dashboard>

For a graphical desktop session:

1. Go to `Interactive Apps`.
2. Choose `Desktop`.
3. Set the maximum running time, for example `3` hours.
4. Tick `Constrained to epyc3 Host`.
5. Click `Launch`.
6. Open a terminal inside the desktop session.

<!-- For shell-only access from the onBunya dashboard, use `Clusters` then `Bunya Shell Access`. -->

### Option 2: SSH from macOS or Linux

Open Terminal and run:

```bash
ssh uqdma8@bunya.rcc.uq.edu.au
```

Enter the Bunya/UQ password when prompted, then enter the Duo passcode.

To leave the SSH session:

```bash
exit
```

### Option 3: SSH from a Windows laptop

The simplest Windows setup is Windows Terminal or PowerShell with the built-in OpenSSH client.

Open PowerShell or Windows Terminal and run:

```powershell
ssh uqdma8@bunya.rcc.uq.edu.au
```

Connection details:

- Host name: `bunya.rcc.uq.edu.au`
- Port: `22`
- Username: `uqdma8`
- Authentication: UQ/Bunya password, then Duo MFA passcode

On the first connection, SSH may ask whether to trust the host key. Continue only if the host is `bunya.rcc.uq.edu.au`.

Use the Duo passcode for Bunya MFA unless UQ RCC announces a change. The official Bunya guide says not to use Okta MFA for Bunya until advised.

PuTTY is also fine:

1. Open PuTTY.
2. Set `Host Name` to `bunya.rcc.uq.edu.au`.
3. Set `Port` to `22`.
4. Set `Connection type` to `SSH`.
5. Click `Open`.
6. Log in with username, password, and Duo MFA.

If using MobaXterm, disable `Remote-monitoring (Experimental)` in SSH settings. The official Bunya guide warns that this feature should not be used.

Do not share passwords, SSH keys, or MFA credentials.

## File transfer

### Recommended for Windows: WinSCP

WinSCP is a convenient graphical SFTP client for Windows and works with MFA.

Use these settings:

- File protocol: `SFTP`
- Host name: `bunya.rcc.uq.edu.au`
- Port: `22`
- Username: `uqdma8`

Log in with password and Duo MFA. WinSCP is useful because it can transfer multiple files without asking for the MFA passcode for every individual file.

### Command-line transfer from Windows, macOS, or Linux

Upload one file:

```bash
scp local_file.txt uqdma8@bunya.rcc.uq.edu.au:/home/uqdma8/
```

Upload a folder recursively:

```bash
scp -r local_folder uqdma8@bunya.rcc.uq.edu.au:/home/uqdma8/
```

Start an interactive SFTP session:

```bash
sftp uqdma8@bunya.rcc.uq.edu.au
```

Useful SFTP commands:

```text
pwd      # Bunya working directory
ls       # list Bunya files
cd       # change Bunya directory
lpwd     # local laptop working directory
lls      # list local laptop files
lcd      # change local laptop directory
put x    # upload x to Bunya
get x    # download x from Bunya
exit     # leave sftp
```

### onBunya file upload

In the onBunya dashboard, go to `Files` then `Home Directory`. Files can be dragged and dropped into the current folder.

### macOS graphical transfer

Cyberduck can be used on macOS, but it may ask for a passcode each time files are uploaded or downloaded.

### Avoid FileZilla

The official Bunya guide no longer recommends FileZilla.

## Repository location

The current Slurm and R scripts use paths under:

```bash
~/JoMoNoPH-share
```

On Bunya, confirm that the repository exists there:

```bash
cd ~/JoMoNoPH-share
pwd
ls
```

The official Bunya guide recommends running jobs from `/scratch` rather than `/home`, especially for larger workloads. This repository currently uses `$HOME` paths in several scripts. If moving the project to `/scratch/user/$USER`, update the hard-coded paths in the Slurm and R scripts before submitting jobs.

Check Bunya quota and usage:

```bash
rquota
```

## Set up the R and CmdStan environment

Use the standalone R terminal on Bunya. For this project, do not use RStudio for package installation or Stan compilation, because the R version in RStudio may differ from the standalone R module used by the Slurm scripts.

Known-working module:

```bash
module purge
module load r/4.4.0-gfbf-2023a
R
```

Inside R, install CmdStanR using the recommended R-universe route:

```r
# Run in a fresh R session if possible.
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
```

Install CmdStan:

```r
library(cmdstanr)
install_cmdstan(version = "2.37.0") # instead of install_cmdstan()
cmdstan_path()
cmdstan_version()
```

The project file `Stan/compile-stan-models.R` currently contains:

```r
set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")
```

If `install_cmdstan()` installs CmdStan somewhere else, either install CmdStan 2.37.0 or update `set_cmdstan_path()` before compiling models. The exact CmdStan version appears less important than keeping the installed path and the project script consistent.

To install the version currently expected by the project script:

```r
install_cmdstan(version = "2.37.0")
```

Avoid installing the development version unless there is a specific reason:

```r
# Not recommended for the stable Bunya setup.
# install.packages("remotes")
# remotes::install_github("stan-dev/cmdstanr")
```

## Install required R packages

Start a fresh standalone R session after loading `r/4.4.0-gfbf-2023a`, then run:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c(
  "tidyverse",
  "ggplot2",
  "nlme",
  "survival",
  "MASS",
  "simsurv",
  "future",
  "furrr",
  "R.utils"
)

lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(lib)) {
  dir.create(lib, recursive = TRUE)
}

.libPaths(c(lib, .libPaths()))
install.packages(
  setdiff(pkgs, rownames(installed.packages())),
  dependencies = TRUE
)
```

After installation, quit R:

```r
q()
```

## Compile Stan models

Compile Stan executable files once before submitting the main simulations.

Recommended interactive method:

```bash
module purge
module load r/4.4.0-gfbf-2023a
R
```

Inside R:

```r
source("~/JoMoNoPH-share/Stan/compile-stan-models.R")
```

Then quit R:

```r
q()
```

Optional Slurm method:

```bash
cd ~/JoMoNoPH-share/Stan
sbatch compile-stan-models.slurm
```

If compilation fails after changing CmdStan versions, remove the old CmdStan installation and reinstall:

```bash
rm -rf ~/.cmdstan
```

Then repeat the CmdStanR and CmdStan installation steps in standalone R.

## Submit ENZAMET simulation jobs

If a `.slurm` file was created or edited on Windows, convert it to Unix line endings before submitting it on Bunya. Otherwise `sbatch` may fail with errors such as `bad interpreter`, `^M`, or `Batch script contains DOS line breaks`.

Convert one Slurm file:

```bash
dos2unix ~/JoMoNoPH-share/Bunya/ENZAMET/ENZAMET-DM-2.slurm
```

Convert all ENZAMET Slurm launchers:

```bash
dos2unix ~/JoMoNoPH-share/Bunya/ENZAMET/*.slurm
```

From any Bunya shell:

```bash
sbatch ~/JoMoNoPH-share/Bunya/ENZAMET/ENZAMET-DM.slurm
```

For the split full run, submit the numbered launchers as needed (Min's is -2.slurm):

```bash
sbatch ~/JoMoNoPH-share/Bunya/ENZAMET/ENZAMET-DM-1.slurm
sbatch ~/JoMoNoPH-share/Bunya/ENZAMET/ENZAMET-DM-2.slurm
sbatch ~/JoMoNoPH-share/Bunya/ENZAMET/ENZAMET-DM-3.slurm
sbatch ~/JoMoNoPH-share/Bunya/ENZAMET/ENZAMET-DM-4.slurm
```

## Monitor and manage jobs

Check queued and running jobs:

```bash
squeue -u $USER
```

Cancel a job:

```bash
scancel JOB_ID
```

Check accounting details after a job starts or finishes:

```bash
sacct -j JOB_ID
```

This repository redirects Slurm output to folders like:

```bash
~/out_YYYYMMDD_JOBID/
```

Simulation outputs are written by the current ENZAMET Slurm scripts to folders like:

```bash
~/JoMoNoPH-share/Bunya/ENZAMET/output_YYYYMMDD_JOBID/
```

Inspect output:

```bash
ls ~/out_YYYYMMDD_JOBID
tail -f ~/out_YYYYMMDD_JOBID/ENZAMET_JOBID_TASKID.out
```

## Useful shell commands

```bash
cd ~/folderX          # change directory to folderX
cd ..                # go to the parent directory
pwd                  # print current directory
ls                   # list files
mkdir -p folderX     # create folderX if needed
sbatch xx.slurm      # submit a Slurm job
dos2unix xx.slurm    # convert Windows line endings before sbatch
squeue -u $USER      # check your jobs
scancel JOB_ID       # terminate a job
rm -rf ~/.cmdstan    # remove CmdStan installation
exit                 # leave ssh, sftp, or shell session
```

## Slurm syntax notes

In Slurm scripts, lines beginning with uppercase `#SBATCH` are special Slurm directives. They look like shell comments, but `sbatch` reads them before the shell runs the script.

Active directive:

```bash
#SBATCH --cpus-per-task=16
```

Disabled directive:

```bash
##SBATCH --cpus-per-task=40
```

A normal shell command is commented out with one `#`:

```bash
# sbatch other_script.slurm
```

The special behavior applies to `#SBATCH` directives, not to every line that starts with `#`.

## Stability checklist

Before a production run, check:

- Use onBunya with `Constrained to epyc3 Host` when using the desktop.
- Load `r/4.4.0-gfbf-2023a` before installing packages, compiling Stan, or running Slurm jobs.
- Install CmdStanR from R-universe, not the GitHub development version.
- Keep `cmdstan_path()` consistent with `Stan/compile-stan-models.R`.
- Compile Stan models in standalone R or Slurm Rscript, not RStudio.
- Avoid running simulations directly on login nodes; submit jobs with `sbatch`.
- Confirm quotas with `rquota` before large runs.

## Troubleshooting

If login fails, confirm that Bunya is not down for scheduled maintenance. The official guide says Bunya is taken offline four times per year, on the second Tuesday of February, May, August, and November.

If Duo MFA fails repeatedly, wait before retrying. If the problem persists, contact UQ ITS for MFA help.

If Bunya login still fails and it is not a maintenance day, email `rcc-support@uq.edu.au` with:

- Bunya username
- Where you are connecting from, for example campus, home, or VPN
- Tool used, for example Windows Terminal, PuTTY, onBunya, or WinSCP
- Any error messages

If package loading fails inside jobs, reinstall packages from standalone R after loading the same R module used in the Slurm scripts:

```bash
module purge
module load r/4.4.0-gfbf-2023a
R
```

If Stan compilation fails after a new CmdStan installation, re-check:

```r
library(cmdstanr)
cmdstan_path()
cmdstan_version()
```

Then update `Stan/compile-stan-models.R` or reinstall the expected CmdStan version.

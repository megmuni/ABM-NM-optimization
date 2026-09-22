# Overview
Contains the Nelder-Mead optimization workflow for VUA Lab ABMs. This is based off of [https://github.com/mintary/abm-opt](https://github.com/mintary/abm-opt), but updated to work with the latest ABM input/output structure (for example, the [IVDBM-ABM](https://github.com/megmuni/IVDBM-ABM)).

## Setup
1. Update `simulation_config.template.json` with your specific ABM's config
1. Update `parameters.xlsx` so its rows match the tagged template lines in the same order, and set `Vary?` for the parameters you want optimized
1. Put `experimental_config.csv` in place with your measurements
1. Copy your ABM's `bin/` and `configFiles/` into this directory, then chmod +x bin/testRun
1. Edit OUTPUT_METRICS in ABM.py if you're fitting different timepoints or biomarkers
1. Review the #SBATCH headers in the job scripts (account, walltime, resources)
1. Install the Python dependencies for your account (see below); one-time step

### Python requirements
`requirements.txt` has everything you need to install for the optimization scripts to import. Do this once for your account on a given cluster:
```bash
module load StdEnv/2020 gcc/9.3.0 python/3.10
pip install --no-index --user -r requirements.txt
```

Then, at the beginning of every session, you only need to import modules:
```bash
module load StdEnv/2023 gcc/12.3 cuda/12.2 python/3.11
```
# Workflow
There are 3 parts to the optimization workflow, which have to be run in order.

## Part 1: Generate samples 

### Step 1 - Generate parameter sets

To generate parameter sets for the ABM, run this from the main optimization directory:

```bash
python ABM_generate_samples.py write-configs --n-samples 10
```
This draws the parameter sets, writes one ready-to-run JSON per (sample, condition) into `configFiles/samples/`, and records the sampled values in `sample_parameters.csv`. Parameters are drawn once per sample and shared across its conditions, so high/low differ only in the scaffold fields (world_init in the JSON).

### Step 2 - Run the ABM to create output
```bash
mkdir -p logs   # one-time
NTASKS=$(python ABM_generate_samples.py count --n-samples 10)
export EMAIL="you@mail.com"
sbatch --array=0-$((NTASKS-1))%20 --mail-user $EMAIL ABM_generate_samples_job.sh
```

Setting `NTASKS` is important for ensuring that the right number of array jobs are submitted by the `sbatch` statement. Make sure you are consistent with what you selected for `--n-samples` in Part 1.

### Step 3 - Collect metrics
```bash
python ABM_generate_samples.py collect          # skips missing tasks with a warning
python ABM_generate_samples.py collect --strict # fail if any task is incomplete
```
Joins `sample_parameters.csv` with each task's metrics into `generated_samples_with_outputs.csv`. To re-run individual failed tasks before collecting:
```bash
python ABM_generate_samples.py run --task 7
```

## Part 2: Running Verification

To verify the ABM model outputs against the previously generated results, run:

```bash
python ABM_verify.py
```

This will re-execute the ABM for each sample, compare the simulated outputs to the expected values, and write a summary of the results (per-sample errors and summary metrics) to `verification_results.json`. 
Useful flags:
- --reps N: averages several runs per sample before scoring
- --limit N: for a very quick test
- --samples / --out / --template: to override paths if necessary

## Part 3: Parameter optimization

### Overview

**Goal**: To determine the values of parameters that minimize the error of ABM outputs

**Logistics**: Uses scipy.optimize package from Python (https://docs.scipy.org/doc/scipy/reference/optimize.html); runs on the cluster using a job script

### Execution

1. Export your email address for notifications like so:

```bash
export EMAIL="youremail@mail.com"
```

2. Submit the script

```bash
sbatch --mail-user $EMAIL ABM_optimize_job.sh [args]
```
| Flag | Default | Meaning |
| --- | --- | --- |
| `--mode` | `joint` | `joint` = one fit, objective summed over all conditions; `separate` = independent fit and parameter set per condition |
| `--conditions` | `high low` | Which scaffold conditions to fit |
| `--iters` | `3` | ABM runs per condition per function evaluation, averaged (the ABM is stochastic) |
| `--maxiter` | `50` | Maximum Nelder-Mead **iterations** (not evaluations) |
| `--maxfev` | none | Hard cap on objective evaluations |
| `--max-params` | none | **Testing only.** Optimize just the first N varying parameters |
| `--numticks` | from `OUTPUT_METRICS` | **Testing only.** Shorten each execution |
| `--tol` | `1e-4` | Convergence tolerance |
| `--parallel` | all | Concurrent ABM executions. Within one evaluation the (conditions x iters) executions are independent; `1` for serial. |
| `--devices` | `$SLURM_GPUS_ON_NODE` | GPUs to spread executions across, round-robin. |
| `--cache` | `output/objective_cache.json` | Checkpoint of evaluated parameter vectors, enables resuming |
| `--no-cache` | off | Disable the cache, re-run every evaluation |
| `--workdir-root` | `$SLURM_TMPDIR/abm_opt` | Where per-execution working directories go |
| `--snapshots` | off | Save periodic biomarker snapshots |
| `--out` | `output/optimization_results.json` | Results path |

3. Results will be outputed in `output`. You should also see a tarball archive of the entire directory.

#### Notes
**Parallelism.** Nelder-Mead is sequential; that is, each iteration needs the
previous result. However, (conditions x `--iters`) can be parallel and
now run concurrently by default. Each runs in
its own working directory under `$SLURM_TMPDIR`, since the ABM resolves
`./bin/testRun` and its output paths relative to the working directory.

**Resuming after a timeout.** Every completed evaluation is appended to
`output/objective_cache.json`. Nelder-Mead is deterministic given the
objective values, so **resubmitting the same command** replays the
previous run's path from the cache in minutes and continues from where it
stopped:

```bash
sbatch --mail-user $EMAIL ABM_optimize_job.sh --mode joint   # timed out
sbatch --mail-user $EMAIL ABM_optimize_job.sh --mode joint   # picks up where it left off
```

Other arguments must match, since the cache key includes the conditions and
`--iters`. Delete
the cache to start clean. Parameters optimized are every row of `parameters.xlsx` with a blank
`Vary?` column. Initial
simplex is built from the parameter bounds (all minimums, then one maximum
at a time, then all maximums).

### Analysis

- Output includes all the variables you included for each ABM execution after the “print” statement, which you can use to monitor how they change with optimization
- Output ends with results of optimization, i.e. optimal parameter values, minimum error, and stopping criteria met
- We are typically most interested in evaluating how much error decreased (absolute and % decrease in error)

## Testing

TO BE UPDATED

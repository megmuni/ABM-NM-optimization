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
module load StdEnv/2020 gcc/9.3.0 python/3.10
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

# Part 2: Running Verification

To verify the ABM model outputs against the previously generated results, run:

```bash
python ABM_verify.py
```

This will re-execute the ABM for each sample, compare the simulated outputs to the expected values, and write a summary of the results (per-sample errors and summary metrics) to `verification_results.json`. 
Useful flags:
- --reps N: averages several runs per sample before scoring
- --limit N: for a very quick test
- --samples / --out / --template: to override paths if necessary

# Part 3: Parameter optimization

## Overview

**Goal**: To determine the values of parameters that minimize the error of ABM outputs

**Logistics**: Uses scipy.optimize package from Python (https://docs.scipy.org/doc/scipy/reference/optimize.html)

## Setup

1. Update the Sensitivity Analysis.xslx file with your parameters
2. Update the `ABM_optimize.py` file as needed.
   a. sam: vector of parameters you are optimizing
   b. Your ABM specifications in the subprocess.call
   c. Y: stores each term of your objective function (i.e. error function), which includes the variables you are interested in and experimental values
   i. Currently, this is a sum of square errors
   ii. Because of the stochasticity of the ABM, currently, each parameter set is executed 3 times and averaged
   d. Output you are interested in tracking after “print”
   e. numpar: number of parameters in your ABM
   f. p1, p2, …, pn: parameter numbers of the parameters you are optimizing
   g. par_s: number of parameters you are optimizing
   h. init: initial values for optimization to begin at
   i. minimize: call to scipy.optimize function, where you can specify the optimization technique you wish to use, as well as any input it requires
3. Edit the submit ABM_optimize_job.sh file as needed
4. Create a folder with the bin and configFiles folders from your ABM, the Sensitivity Analysis.xlsx file with the parameter information, and the code files ABM_optimize.py and ABM_optimize_job.sh
   a. Make sure testRun in bin has execution permissions, “chmod +x testRun”
5. Create a subfolder named output and a subfolder under output named SensitivityAnalysis

## Execution

1. Export your email address for notifications like so:

```bash
export EMAIL="youremail@mail.com"
```

2. Submit the script

```bash
sbatch --mail-user $EMAIL ABM_optimize_job.sh [optional args]
```

3. Results will be outputed in `output`. You should also see a tarball archive of the entire directory.

## Analysis

- Output includes all the variables you included for each ABM execution after the “print” statement, which you can use to monitor how they change with optimization
- Output ends with results of optimization, i.e. optimal parameter values, minimum error, and stopping criteria met
- We are typically most interested in evaluating how much error decreased (absolute and % decrease in error)

## Changing the number of parameters

You can optionally submit the job with a number of parameters and a method for ranking parameter importance. The method must be present as a column title in `Sensitivity Analysis.xlsx`.

The default method is Random Forest with `n=5` parameters.

To change the number of parameters:

```bash
sbatch --mail-user $EMAIL ABM_optimize_job.sh --n [NUMBER HERE]
```

To change the method of selecting parameters (different ranking):

```bash
sbatch --mail-user $EMAIL ABM_optimize_job.sh --method [METHOD HERE]
```

You can customize both the number and method at the same time:

```bash
sbatch --mail-user $EMAIL ABM_optimize_job.sh --n [NUMBER HERE] --method [METHOD HERE]
```

## Testing

TO BE UPDATED

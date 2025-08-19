# Parameter optimization

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

## Saving Snapshots

To enable saving biomarker snapshots, use the `--save-snapshots` flag (do not provide a value):

```bash
sbatch --mail-user $EMAIL ABM_optimize_job.sh --n [NUMBER HERE] --save-snapshots
```

## Testing

It can be time-consuming to run the ABM_optimize_job.sh script as a job.

### Testing dimension fitting

This option can be used to test that the `ABM_optimize.py` script correctly passes inputs from each function to the other.

Inside of ABM_optimize_job.sh, you can uncomment out the lines in the blocks marked as TESTING to run the script in your scratch directory (edit the line to change it to another directory if you want).

Note that if you do this, you should ensure that you run the script with --test True to use the test_ABM method rather than the actual ABM method (which is resource intensive).

```bash
./ABM_optimize_job.sh --test True
```

### Testing `testRun`

`testRun` relies on an available CUDA device, which is only available on the compute clusters.

As a result, this script attempts to run a single iteration of `testRun`, requesting fewer resources than `ABM_optimize.py`. However, it must still be submitted as a job.

```bash
EMAIL=[your email for notifs]
sbatch ABM_light_test.sh
```

You can also run the job in interactive mode with `salloc`, see documentation here: https://slurm.schedmd.com/salloc.html

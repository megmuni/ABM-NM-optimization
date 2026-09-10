#!/bin/bash
#SBATCH --account=rrg-nicoleli
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --gpus=h100:1
#SBATCH --mem=64000M
#SBATCH --mail-user=${EMAIL}
#SBATCH --mail-type=ALL

# =========================================
# ABM Nelder-Mead optimization job.
#
# Usage:
#   export EMAIL="you@mail.com"
#   sbatch --mail-user $EMAIL ABM_optimize_job.sh [args passed to ABM_optimize.py]
#
# Examples:
#   sbatch --mail-user $EMAIL ABM_optimize_job.sh
#   sbatch --mail-user $EMAIL ABM_optimize_job.sh --mode separate
#   sbatch --mail-user $EMAIL ABM_optimize_job.sh --mode joint --iters 3 --maxiter 50
#   sbatch --mail-user $EMAIL ABM_optimize_job.sh --mode separate --conditions low
#
# Every argument is forwarded verbatim to ABM_optimize.py, so
# `python ABM_optimize.py --help` is the authoritative list. This script
# only *reads* --mode and --conditions
#
# NOTE ON WALLTIME: the objective now reads day 21, so each ABM execution
# runs 1009 ticks rather than the 289 the old script assumed -- roughly
# 3.5x longer per run. Total executions are
# (NM evaluations) x (conditions) x (--iters), and --mode separate runs a
# full optimization per condition. Time one ABM execution with
# run_ABM_job.sh and multiply before trusting the 23h above.
# =========================================

# =========================================
# SETUP (virtualenv, modules, etc.)
# =========================================

module load StdEnv/2023 gcc/12.3 cuda/12.2 python/3.11 || { echo "Module load failed"; exit 1; }

# Uncomment the following lines
# to override the default SLURM_TMPDIR
# for testing purposes.

# This will create the virtual environment in the user's
# scratch directory. You can also replace it
# with a different location.

# If you do this, make sure ABM_optimize.py
# is not running the default ABM method, but the test_ABM
# method as the main ABM method is extremely resource-intensive
# and should not be run directly in the file system.

# You may also want to comment out the
# lines at the end of the script that
# remove the environment.

# =========================================
# TESTING ONLY
# =========================================
# export SLURM_TMPDIR="$HOME/scratch"
# echo "Changed SLURM_TMPDIR to $SLURM_TMPDIR for testing purposes."

# if [ -d "$SLURM_TMPDIR/env" ]; then
#       echo "Directory '$SLURM_TMPDIR/env' already exists."
# else
#       echo "Directory '$SLURM_TMPDIR/env' does not exist. Creating it."
#       mkdir "$SLURM_TMPDIR/env" || { echo "Failed to create directory"; exit 1; }
# fi
# =========================================

# ARGUMENTS
#
# All args are forwarded to ABM_optimize.py. --mode and --conditions
# are only used so the tarball name says what the job was for
# =========================================

original_args=("$@")
mode=""
conditions=()

while [[ $# -gt 0 ]]; do
	case $1 in
		--mode=*)
			mode="${1#*=}"
			shift
			;;
		--mode)
			mode="$2"
			shift 2
			;;
		--conditions=*)
			conditions+=("${1#*=}")
			shift
			;;
		--conditions)
			shift
			# --conditions takes one or more values (nargs='+'), so consume
			# every following token until the next flag
			while [[ $# -gt 0 && "$1" != --* ]]; do
				conditions+=("$1")
				shift
			done
			;;
		*)
			shift
			;;
	esac
done

# Defaults here must mirror ABM_optimize.py's argparse defaults
mode=${mode:-joint}
if [ ${#conditions[@]} -eq 0 ]; then
	conditions=(high low)
fi
conditions_safe=$(IFS=-; echo "${conditions[*]}")

echo "Mode: $mode"
echo "Conditions: ${conditions[*]}"

current_date=$(date +"%Y-%m-%d_%H-%M-%S")
tarball_name="../param_opt_${current_date}_${mode}_${conditions_safe}.tar.gz"

# trap to package directory on any exit (success or failure)
trap 'echo "Packaging directory (trap)..."; tar -czf "$tarball_name" . && echo "Packaged directory into: $tarball_name"' EXIT

# =========================================
# Some checks to make the job fail immediately with a clear message:
echo "Checking required input files..."
for required in \
	"ABM.py" \
	"ABM_optimize.py" \
	"simulation_config.template.json" \
	"parameters.xlsx" \
	"experimental_config.csv" \
	"bin/testRun"
do
  if [ ! -e "$required" ]; then
		echo "MISSING REQUIRED INPUT: $required"
		exit 1
	fi
done

if [ ! -x "bin/testRun" ]; then
	echo "bin/testRun is not executable; run 'chmod +x bin/testRun'"
	exit 1
fi

# Directories the Python scripts write into. ABM.run_ABM creates output/
# itself, but output/SensitivityAnalysis (the per-execution ABM
# stdout/stderr logs), configFiles/ (the generated sample config) and
# this script's own output/output.txt must exist up front.
mkdir -p output/SensitivityAnalysis configFiles || { echo "Failed to create directories"; exit 1; }

echo "SLURM_TMPDIR: $SLURM_TMPDIR"
df -h $SLURM_TMPDIR || { echo "Failed to check disk space"; exit 1; }
ls -lah $SLURM_TMPDIR

which python
which virtualenv
python --version
virtualenv --version

echo "Creating virtual environment..."
virtualenv --no-download "$SLURM_TMPDIR/env"
if [ $? -ne 0 ]; then
	echo "Virtualenv creation failed."
	exit 1
fi

echo "Virtual environment created successfully."


if [ -f "$SLURM_TMPDIR/env/bin/activate" ]; then
	source "$SLURM_TMPDIR/env/bin/activate"
	echo "Virtual environment created."
else
	echo "Failed to activate virtual env: $SLURM_TMPDIR/bin/activate not found."
	exit 1
fi

pip list || { echo "pip is not working correctly"; exit 1; }

echo "Installing dependencies..."

pip install --upgrade pip || { echo "Failed to upgrade pip"; exit 1; }
pip install --no-index -r requirements.txt || { echo "Failed to install Python dependencies"; exit 1; }
pip freeze > output/installed_packages.txt

echo "Installed packages:"
cat output/installed_packages.txt

# Confirm the imports the optimizer needs are actually usable, and that
# the parameter spreadsheet and config template agree with each other,
# before committing the walltime to a long run.
python -c "
import ABM
print('Fitted targets:', ABM.FITTED_METRIC_KEYS)
print('Tracked metrics:', ABM.OUTPUT_METRIC_KEYS)
print('Ticks per ABM execution:', ABM.METRICS_NUMTICKS)
paths = ABM.get_template_parameter_paths(ABM.CONFIG_TEMPLATE)
import pandas as pd
df = pd.read_excel(ABM.PARAMETER_FILE)
assert df['Parameter Name'].tolist() == paths, \
    'parameters.xlsx and simulation_config.template.json are out of sync'
print('Parameters in sync:', len(paths))
" || { echo "Preflight import/parameter check failed"; exit 1; }
# =========================================

# ==========================================
# RUN PYTHON SCRIPT
# ==========================================

export OMP_NUM_THREADS=1
export OMP_NESTED=TRUE

if [ -f output/objective_cache.json ]; then
	CACHED=$(python -c "import json;print(len(json.load(open('output/objective_cache.json'))))" 2>/dev/null || echo "?")
	echo "Found an existing objective cache with $CACHED evaluation(s)."
	echo "This run will replay them and continue from where the last one stopped."
	echo "Delete output/objective_cache.json to start fresh."
fi

# The ABM allocates device memory on startup, so a missing GPU means
# every evaluation fails. Fail here instead.
if ! nvidia-smi > /dev/null 2>&1; then
	echo "ERROR: no GPU visible to this job."
	echo "  SLURM_JOB_GPUS      = ${SLURM_JOB_GPUS:-<unset>}"
	echo "  SLURM_GPUS_ON_NODE  = ${SLURM_GPUS_ON_NODE:-<unset>}"
	echo "  hostname            = $(hostname)"
	echo "Check the --gpus request and that the account has a GPU allocation"
	echo "on this cluster."
	exit 1
fi
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo "GPUs on node: ${SLURM_GPUS_ON_NODE:-unknown} (used for --devices)"

echo "Running Python script..."
echo "Arguments: ${original_args[@]}"

# -u keeps stdout unbuffered so output/output.txt stays current and can be
# tailed while the job runs -- useful given how long these take.
python -u ABM_optimize.py "${original_args[@]}" > output/output.txt 2>&1
python_status=$?

if [ $python_status -ne 0 ]; then
	echo "Python script failed with exit status $python_status. Last 40 lines of output/output.txt:"
	tail -n 40 output/output.txt
	echo "Last 20 lines of ABM stderr (output/SensitivityAnalysis/stderr.txt):"
	tail -n 20 output/SensitivityAnalysis/stderr.txt 2>/dev/null
	exit 1
fi

# Echo a summary into the SLURM log so results are visible without
# unpacking the tarball
echo "Optimization finished. Tail of output/output.txt:"
tail -n 25 output/output.txt

if [ -f output/optimization_results.json ]; then
	echo "Results written to output/optimization_results.json"
else
	echo "WARNING: output/optimization_results.json not found."
fi

# =========================================
# TESTING ONLY
# =========================================
# echo "Cleaning up: Deleting the virtual environment..."
# rm -rf "$SLURM_TMPDIR/env" || { echo "Failed to delete virtual environment"; exit 1; }
# echo "Virtual environment deleted successfully."
# =========================================

echo "Script finished successfully."
exit 0

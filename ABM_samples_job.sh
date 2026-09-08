#!/bin/bash
# ABM_samples_job.sh
#
#SBATCH --account=rrg-nicoleli
#SBATCH --time=0-03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --gpus=h100_1g.10gb:1
#SBATCH --mem=8000M
#SBATCH --output=logs/gensamples_%A_%a.out
#SBATCH --error=logs/gensamples_%A_%a.err
#SBATCH --mail-user=${EMAIL}
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-19%20
#
# Sample generation for ABM optimization, structured as a SLURM job array.
# One array task per (sample, condition) pair
# The workflow is 3 steps:
#   1. Write the configs (from login node):
#        python ABM_generate_samples.py write-configs --n-samples 10
#      This draws the parameter sets, writes one JSON config per task
#      into configFiles/samples/, and records the drawn values in
#      sample_parameters.csv. It prints the valid task range
#
#   2. Submit this job array. --array MUST match the task count:
#        NTASKS=$(python ABM_generate_samples.py count --n-samples 10)
#        export EMAIL="you@mail.com"
#        sbatch --array=0-$((NTASKS-1))%20 --mail-user $EMAIL ABM_generate_samples_job.sh
#      The %20 caps concurrency at 20 running tasks
#
#   3. Collect once the array finishes:
#        python ABM_generate_samples.py collect
#      This joins sample_parameters.csv with the per-task metrics into
#      generated_samples_with_outputs.csv, which is what ABM_verify.py
#      reads. Add --strict to fail on any incomplete task
#
# The --array default above (0-19) matches --n-samples 10 x 2 conditions.
# If you change the sample count, override --array at submit time -- a
# mismatch either skips tasks silently or fails them with "no task N".
set -uo pipefail

module load StdEnv/2020 gcc/9.3.0 cuda/11.0 python/3.10 || { echo "Module load failed"; exit 1; }

TASK=${SLURM_ARRAY_TASK_ID:-0}
SUBMIT_DIR=${SLURM_SUBMIT_DIR:-$PWD}
echo "=== Task $TASK (array job ${SLURM_ARRAY_JOB_ID:-local}) ==="
echo "Submit directory: $SUBMIT_DIR"
date

mkdir -p logs

cd "$SUBMIT_DIR" || exit 1

# =========================================
# PREFLIGHT -- step 1 must have run already
# =========================================
for required in \
	"ABM.py" \
	"ABM_generate_samples.py" \
	"sample_parameters.csv" \
	"bin/testRun"
do
	if [ ! -e "$required" ]; then
		echo "MISSING REQUIRED INPUT: $required"
		if [ "$required" = "sample_parameters.csv" ]; then
			echo "Run step 1 first: python ABM_generate_samples.py write-configs --n-samples N"
		fi
		exit 1
	fi
done

if [ ! -x "bin/testRun" ]; then
	echo "bin/testRun is not executable; run 'chmod +x bin/testRun'"
	exit 1
fi

# =========================================
# ENVIRONMENT
#
# Built per task in node-local storage. Redundant across tasks, but a
# shared venv on the project filesystem is a contention point when many
# tasks start at once, and $SLURM_TMPDIR is wiped for us on exit.
# =========================================
echo "Creating virtual environment in $SLURM_TMPDIR/env ..."
virtualenv --no-download "$SLURM_TMPDIR/env" || { echo "Virtualenv creation failed"; exit 1; }
source "$SLURM_TMPDIR/env/bin/activate" || { echo "Failed to activate virtualenv"; exit 1; }

pip install --upgrade pip --quiet || { echo "Failed to upgrade pip"; exit 1; }
pip install --no-index -r requirements.txt || { echo "Failed to install dependencies"; exit 1; }

# =========================================
# PRIVATE WORKING DIRECTORY - so that concurrent
# output files don't conflict. only the finished metrics
# are copied back from this temp dir to the main directory
# =========================================
WORKDIR="$SLURM_TMPDIR/task_${TASK}"
mkdir -p "$WORKDIR/output" || { echo "Failed to create workdir"; exit 1; }

ln -s "$SUBMIT_DIR/bin" "$WORKDIR/bin"
ln -s "$SUBMIT_DIR/configFiles" "$WORKDIR/configFiles"
cp "$SUBMIT_DIR/ABM.py" "$SUBMIT_DIR/ABM_generate_samples.py" "$WORKDIR/"
cp "$SUBMIT_DIR/sample_parameters.csv" "$WORKDIR/"
cp "$SUBMIT_DIR/simulation_config.template.json" "$WORKDIR/" 2>/dev/null

cd "$WORKDIR" || exit 1
echo "Working directory: $WORKDIR"

# =========================================
# RUN THIS TASK'S ABM EXECUTION
# =========================================
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

echo "Running task $TASK ..."
python -u ABM_generate_samples.py run --task "$TASK"
status=$?

if [ $status -ne 0 ]; then
	echo "Task $TASK failed with exit status $status"
	echo "Last 20 lines of ABM stderr:"
	tail -n 20 output/stderr.txt 2>/dev/null
	# Preserve the failed run's raw output for debugging
	mkdir -p "$SUBMIT_DIR/output/failed_task_${TASK}"
	cp output/stdout.txt output/stderr.txt output/Output_Biomarkers.csv \
		"$SUBMIT_DIR/output/failed_task_${TASK}/" 2>/dev/null
	exit 1
fi

# =========================================
# COPY RESULTS BACK
#
# Only the metrics JSON is kept by default -- 20 full biomarker CSVs at
# 1009 rows each is a lot of files for data the metrics already
# summarise. Uncomment below to keep the raw CSVs too.
# =========================================
mkdir -p "$SUBMIT_DIR/output/task_metrics"
cp output/task_metrics/*.json "$SUBMIT_DIR/output/task_metrics/" \
	|| { echo "Failed to copy metrics back to $SUBMIT_DIR"; exit 1; }

# mkdir -p "$SUBMIT_DIR/output/task_biomarkers"
# cp output/Output_Biomarkers.csv "$SUBMIT_DIR/output/task_biomarkers/task_${TASK}.csv"

echo "Task $TASK finished successfully."
date
exit 0

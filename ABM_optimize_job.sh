#!/bin/bash
#SBATCH --account=def-nicoleli
#SBATCH --time=09:00:00
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-node=2
#SBATCH --mem=64000M
#SBATCH --mail-user=${EMAIL}
#SBATCH --mail-type=ALL

# =========================================
# SETUP (virtualenv, modules, etc.)
# =========================================

module load StdEnv/2020 gcc/9.3.0 cuda/11.0 python/3.10 || { echo "Module load failed"; exit 1; }

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

echo "SLURM_TMPDIR: $SLURM_TMPDIR"
df -h $SLURM_TMPDIR || { echo "Failed to check disk space"; exit 1; }
ls -lah $SLURM_TMPDIR

# =============================
# COPY FILES TO SLURM_TMPDIR AND SETUP TRAP
# =============================

mkdir -p "$SLURM_SUBMIT_DIR/finished_runs"
mkdir -p "$SLURM_SUBMIT_DIR/slurm_logs"

# Trap function to package and copy results, capturing exit code
function package_dir() {
	exit_code=$?
	echo "Packaging directory (trap)... (exit code: $exit_code)"
	cd "$SLURM_TMPDIR"
	# Copy the slurm log to the tarball
	if [[ -n "$SLURM_JOB_NAME" && -n "$SLURM_JOB_ID" ]]; then
		cp "$SLURM_SUBMIT_DIR/slurm_logs/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out" . 2>/dev/null
		cp "$SLURM_SUBMIT_DIR/slurm_logs/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err" . 2>/dev/null
	fi
	tar --exclude="env" --exclude="./env" --exclude="$SLURM_SUBMIT_DIR/finished_runs" -czf "$SLURM_SUBMIT_DIR/finished_runs/$tarball_name" .
	echo "Packaged directory into: $SLURM_SUBMIT_DIR/finished_runs/$tarball_name"
	exit $exit_code
}

trap 'package_dir' EXIT

echo "Copying files to SLURM_TMPDIR..."
rsync -av --exclude='slurm_logs' --exclude='finished_runs' "$SLURM_SUBMIT_DIR/" "$SLURM_TMPDIR/" || { echo "Failed to copy files"; exit 1; }
cd "$SLURM_TMPDIR" || { echo "Failed to change to SLURM_TMPDIR"; exit 1; }

# extract arguments for naming output file
for arg in "$@"; do
	case $arg in
		--n=*) n="${arg#*=}" ;;
		--method=*) method="${arg#*=}" ;;
	esac
done

n=""
method=""
while [[ $# -gt 0 ]]; do
	case $1 in
		--n=*)
			n="${1#*=}"
			shift
			;;
		--n)
			n="$2"
			shift 2
			;;
		--method=*)
			method="${1#*=}"
			shift
			;;
		--method)
			method="$2"
			shift 2
			;;
		*)
			shift
			;;
	esac
done

n=${n:-5}
method=${method:-"Random Forest"}
method_safe=$(echo "$method" | tr ' ' '_') # replace spaces with underscores for safe filename

echo "Method: $method_safe"
echo "Parameters: $n"

current_date=$(date +"%Y-%m-%d_%H-%M-%S")
tarball_name="../param_opt_${current_date}_n${n}_${method_safe}.tar.gz"

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
#pip install --no-index numpy scipy pandas xlrd openpyxl || { echo "Failed to install Python dependencies"; exit 1; }
pip install --no-index -r requirements.txt || { echo "Failed to install Python dependencies"; exit 1; }
pip freeze > requirements.txt

echo "Installed packages:"
cat requirements.txt


# ==========================================
# RUN PYTHON SCRIPT
# ==========================================

export OMP_NUM_THREADS=32
export OMP_NESTED=TRUE

# check the CUDA device
nvidia-smi || { echo "Failed to check CUDA device"; exit 1; }

echo "Running Python script..."
python ABM_optimize.py "$@" > output/output.txt 2>&1 || { echo "Python script failed"; exit 1; }

# =========================================
# TESTING ONLY
# =========================================
# echo "Cleaning up: Deleting the virtual environment..."
# rm -rf "$SLURM_TMPDIR/env" || { echo "Failed to delete virtual environment"; exit 1; }
# echo "Virtual environment deleted successfully."
# =========================================

echo "Script finished successfully."
exit 0


#!/bin/bash
#SBATCH --account=def-nicoleli
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-node=2
#SBATCH --mem=16000M
#SBATCH --mail-user=emily.wang10@mail.mcgill.ca
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

echo "Running Python script..."
python test_testRun.py > output/output.txt 2>&1 || { echo "Python script failed"; exit 1; }

# =========================================
# CREATE OUTPUT DIRECTORY
# =========================================

# extract arguments for naming output file
for arg in "$@"; do
    case $arg in
        --n=*) n="${arg#*=}" ;; # extract value for num of parameters
        --method=*) method="${arg#*=}" ;; # extract value for method
    esac
done

# set default values if not provided
n=${n:-5}
method=${method:-"Random_Forest"}

# replace spaces in method with underscores
method_safe=$(echo "$method" | tr ' ' '_')

# get current date (this is used to name the output file)
current_date=$(date +"%Y-%m-%d_%H-%M-%S")


# report the running time of the script
end_time=$(date +%s)
start_time=$(date -d "$SLURM_JOB_START_TIME" +%s)
running_time=$((end_time - start_time))

touch "running_time.txt"
echo "Running time: $((running_time / 3600)) hours $(((running_time % 3600) / 60)) minutes $((running_time % 60)) seconds" > "running_time.txt"

# package the entire directory
echo "Packaging directory..."
tarball_name="../param_opt_${current_date}_n${n}_${method_safe}.tar.gz"

echo "Creating tarball..."
tar -czf "$tarball_name" . || { echo "Failed to create tarball"; exit 1; }

echo "Packaged directory into: $tarball_name"

# =========================================
# TESTING ONLY
# =========================================
# echo "Cleaning up: Deleting the virtual environment..."
# rm -rf "$SLURM_TMPDIR/env" || { echo "Failed to delete virtual environment"; exit 1; }
# echo "Virtual environment deleted successfully."
# =========================================

echo "Script finished successfully."
exit 0


#!/bin/bash
#SBATCH --account=def-nicoleli
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=32
#SBATCH --gpus=h100:2
#SBATCH --mem=60000M
#SBATCH --mail-user=${EMAIL}
#SBATCH --mail-type=ALL

# =========================================
# SETUP (virtualenv, modules, etc.)
# =========================================

module load StdEnv/2020 gcc/9.3.0 cuda/11.0 python/3.10 || { echo "Module load failed"; exit 1; }

current_date=$(date +"%Y-%m-%d_%H-%M-%S")
tarball_name="../param_opt_${current_date}.tar.gz"

trap 'echo "Packaging directory (trap)..."; tar -czf "$tarball_name" . && echo "Packaged directory into: $tarball_name"' EXIT
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

export OMP_NUM_THREADS=32
export OMP_NESTED=TRUE

echo "Starting memory monitoring..."
while true; do
    echo "$(date): $(free -h | grep Mem)" >> memory_usage.log
    sleep 10
done &
MONITOR_PID=$!

# Copy ExampleSample.txt to Sample.txt so it can be fed into the program
cp ExampleSample.txt Sample.txt

./bin/testRun --numticks 289 --inputfile configFiles/config_Scaffold_GH2.txt --wxw 0.6 --wyw 0.6 --wzw 0.6 > output.txt

# Stop memory monitoring
kill $MONITOR_PID 2>/dev/null
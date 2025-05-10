#!/bin/bash
#SBATCH --account=rrg-nicoleli
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-node=2
#SBATCH --mem=64000M
#SBATCH --mail-user=grace.yu@mail.mcgill.ca
#SBATCH --mail-type=ALL

module load cuda
module load python

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r requirements.txt
pip install numpy
pip install scipy
pip install pandas
pip install xlrd
pip install openpyxl

export OMP_NUM_THREADS=32
export OMP_NESTED=TRUE

python ABM_optimize.py > output.txt

# Finish the script
exit 0

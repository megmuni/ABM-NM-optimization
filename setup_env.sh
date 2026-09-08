#!/bin/bash
# setup_env.sh
#
# Creates (or refreshes) a persistent Python environment for this repo,
# for work that runs OUTSIDE a job: the write-configs and collect phases
# of sample generation, ABM_verify smoke tests, and the analysis notebook
#
# Usage (once, from the repo directory):
#   bash setup_env.sh
#
# Then in each new shell:
#   source activate_env.sh
#
# Override the location if you want it elsewhere (e.g. in scratch):
#   ABM_ENV=~/scratch/ABM-env bash setup_env.sh

set -uo pipefail

ABM_ENV="${ABM_ENV:-$HOME/ABM-env}"
MODULES="StdEnv/2023 gcc/12.3 python/3.11"

echo "Loading modules: $MODULES"
module load $MODULES || { echo "Module load failed"; exit 1; }

python --version

if [ -d "$ABM_ENV" ]; then
	echo "Environment already exists at $ABM_ENV -- refreshing packages."
else
	echo "Creating virtual environment at $ABM_ENV ..."
	virtualenv --no-download "$ABM_ENV" || { echo "Virtualenv creation failed"; exit 1; }
fi

source "$ABM_ENV/bin/activate" || { echo "Failed to activate $ABM_ENV"; exit 1; }

# --no-index forces the cluster's local wheelhouse rather than PyPI
pip install --no-index --upgrade pip || { echo "Failed to upgrade pip"; exit 1; }
pip install --no-index -r requirements.txt || { echo "Failed to install dependencies"; exit 1; }

echo
echo "Verifying the imports the scripts need ..."
python -c "
import numpy, pandas, scipy, openpyxl
print('numpy', numpy.__version__)
print('pandas', pandas.__version__)
print('scipy', scipy.__version__)
import ABM
print('ABM imported OK; fitted targets:', ABM.FITTED_METRIC_KEYS)
" || { echo "Import check failed"; exit 1; }

echo
echo "Done. Environment: $ABM_ENV"
echo "In each new shell, run:  source activate_env.sh"

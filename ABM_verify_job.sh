#!/bin/bash
# ABM_verify_job.sh
#
# Usage: bash ABM_verify_job.sh [args for ABM_verify.py]
#
# Examples:
#   bash verify_interactive.sh --limit 2          e.g. for a quick test
#   bash verify_interactive.sh --reps 3
#
# This script requests the allocation and runs verification inside it,
# then exits. You can also control your own interactive job by using the
# salloc line from this script, then running the python command yourself

set -uo pipefail

TIME=${TIME:-0-03:00:00}
ACCOUNT=${ACCOUNT:-rrg-nicoleli}
GPUS=${GPUS:-h100_1g.10gb:1}
MEM=${MEM:-8000M}

echo "Requesting an interactive session (${TIME}, ${GPUS}, ${MEM})..."
echo "This may queue for a few minutes."
echo

# --- what runs inside the allocation ---------------------------------
# Kept in a heredoc so the whole thing is one command. The ABM needs a
# GPU and its own shared libraries, so both are set up here.
read -r -d '' INNER <<'INNER_EOF'
set -uo pipefail

module load StdEnv/2020 gcc/9.3.0 cuda/11.0 python/3.10 || exit 1

if ! nvidia-smi > /dev/null 2>&1; then
	echo "ERROR: no GPU in this allocation. The ABM needs a CUDA device."
	exit 1
fi
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ABM shared libraries (libAgent.so etc.) aren't on the default search path
if [ -n "${ABM_LIB_DIR:-}" ]; then
	LIB_DIRS=("$ABM_LIB_DIR")
else
	LIB_DIRS=("$PWD/bin" "$PWD/lib" "$PWD/build/lib" "$PWD/build/bin")
fi
for candidate in "${LIB_DIRS[@]}"; do
	if compgen -G "$candidate/libAgent.so*" > /dev/null 2>&1; then
		export LD_LIBRARY_PATH="$candidate${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
		echo "ABM libraries: $candidate"
		break
	fi
done
if ldd bin/testRun 2>/dev/null | grep -q "not found"; then
	echo "ERROR: bin/testRun has unresolved libraries:"
	ldd bin/testRun | grep "not found"
	echo "Set ABM_LIB_DIR, or run ./find_abm_libs.sh"
	exit 1
fi

# Python: the persistent env if it exists, else assume a --user install
if [ -f "${ABM_ENV:-$HOME/ABM-env}/bin/activate" ]; then
	source "${ABM_ENV:-$HOME/ABM-env}/bin/activate"
fi
python -c "import numpy, pandas, scipy, openpyxl" || {
	echo "ERROR: Python dependencies missing. See the README setup step:"
	echo "  pip install --no-index --user -r requirements.txt"
	exit 1
}
 
echo
echo "Running verification..."
python -u ABM_verify.py all $VERIFY_ARGS
INNER_EOF
# ---------------------------------------------------------------------

export VERIFY_ARGS="$*"

salloc --account="$ACCOUNT" --time="$TIME" --gpus="$GPUS" \
	--cpus-per-task=1 --mem="$MEM" \
	bash -c "$INNER"

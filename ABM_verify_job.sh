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

TIME=${TIME:-0-01:00:00}
ACCOUNT=${ACCOUNT:-rrg-nicoleli}
GPU_FLAG=${GPU_FLAG:-h100_1g.10gb:1}
MEM=${MEM:-8000M}

echo "Requesting an interactive session (${TIME}, ${GPU_FLAG}, ${MEM})..."
echo "This may queue for a few minutes."
echo

# --- what runs inside the allocation ---------------------------------
# Kept in a heredoc so the whole thing is one command. The ABM needs a
# GPU and its own shared libraries, so both are set up here.
read -r -d '' INNER <<'INNER_EOF'
set -uo pipefail

module load StdEnv/2023 gcc/12.3 cuda/12.2 python/3.11 || exit 1

if ! nvidia-smi > /dev/null 2>&1; then
	echo "ERROR: no GPU visible."
	echo "  SLURM_JOB_GPUS      = ${SLURM_JOB_GPUS:-<unset>}"
	echo "  SLURM_GPUS_ON_NODE  = ${SLURM_GPUS_ON_NODE:-<unset>}"
	echo "  CUDA_VISIBLE_DEVICES= ${CUDA_VISIBLE_DEVICES:-<unset>}"
	echo "  hostname            = $(hostname)"
	echo "  SLURM_NODELIST      = ${SLURM_NODELIST:-<unset>}"
	echo
	echo "If the SLURM_* variables above are unset, the allocation had no GPU"
	echo "attached. Common causes:"
	echo "  - the GPU type in the request doesn't exist on this cluster;"
	echo "    check with 'sinfo -o \"%%N %%G\" | sort -u'"
	echo "  - the account isn't GPU-eligible; GPU allocations are often a"
	echo "    separate account (e.g. def-nicoleli-gpu or an rrg-* RAC)."
	echo "    Check with 'sshare -U' or 'sacctmgr show assoc user=$USER'"
	echo "  - the flag syntax: some clusters want --gres=gpu:h100:1 rather"
	echo "    than --gpus=h100:1"
	exit 1
fi
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader

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
python -u ABM_verify.py $VERIFY_ARGS
INNER_EOF
# ---------------------------------------------------------------------

export VERIFY_ARGS="$*"

# srun (not salloc): the command must execute on the allocated node.
srun --account="$ACCOUNT" --time="$TIME" $GPU_FLAG \
	--cpus-per-task=1 --mem="$MEM" \
	--pty bash -c "$INNER"

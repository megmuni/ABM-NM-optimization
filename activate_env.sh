#!/bin/bash
# activate_env.sh
#
# Loads the modules and activates the persistent environment created by
# setup_env.sh. Must be SOURCED, not executed, so the changes apply to
# your current shell:
#
#   source activate_env.sh
#
# Override the location the same way as setup_env.sh:
#   ABM_ENV=~/scratch/ABM-env source activate_env.sh

ABM_ENV="${ABM_ENV:-$HOME/ABM-env}"

if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
	echo "This script must be sourced, not executed:"
	echo "  source activate_env.sh"
	exit 1
fi

module load StdEnv/2023 gcc/12.3 python/3.11 || return 1

if [ ! -f "$ABM_ENV/bin/activate" ]; then
	echo "No environment at $ABM_ENV -- run 'bash setup_env.sh' first."
	return 1
fi

source "$ABM_ENV/bin/activate" || return 1
echo "Activated $ABM_ENV ($(python --version 2>&1))"

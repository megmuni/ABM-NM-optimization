#!/bin/bash

# This script is used to submit a job to the SLURM scheduler
# For email notifications, set the EMAIL variable.

if [ -z "$1" ]; then
    echo "Usage: $0 your_email@domain.com"
    exit 1
fi

EMAIL="$1"

sbatch --export=EMAIL="$EMAIL" ABM_optimize_job.sh "$@"
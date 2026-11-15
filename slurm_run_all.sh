#!/bin/bash
# Slurm wrapper for the KCL cluster (CREATE) only. On your Mac run:  bash run_all.sh
# Edit the module lines to match your cluster environment.
#SBATCH --job-name=amp_analysis
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=amp_analysis_%j.log
if command -v module >/dev/null 2>&1; then
  module load gromacs      # ASSUMPTION: adjust to the module name you use for the production runs
  module load python       # or: source activate <env with MDAnalysis>
fi
[[ -n "${SLURM_CPUS_PER_TASK:-}" ]] && export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "$SLURM_SUBMIT_DIR"
PARALLEL=1 bash run_all.sh "$@"

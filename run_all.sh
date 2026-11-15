#!/usr/bin/env bash
# Whole pipeline for the systems given (default: every row of systems.tsv).
#   bash run_all.sh                      # all systems
#   bash run_all.sh mel_popc mel_pepg    # selected systems
#   GMX=gmx_mpi DT_ANALYSIS=200 bash run_all.sh
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
bash "$here/00_register_settings.sh" "$@"
bash "$here/01_make_index.sh" "$@"
bash "$here/02_process_trajectories.sh" "$@"
bash "$here/03_run_analysis.sh" "$@"

#!/usr/bin/env bash
# LL-37 2K6O in water, 313 K, 1 us
# Put the run's tpr and xtc (or xtc parts) in:  AMP_MD_processed/ll37/10_ll37_2k6o_313K
# then:  bash pipeline/new_systems/process_ll37_sol313.sh      (TPR=... XTC=... to name the files explicitly)
exec bash "$(dirname "$0")/process_system.sh" ll37_sol313 ll37 solution_313K r1 solution 1 ll37/10_ll37_2k6o_313K "${TPR:-}" "${XTC:-}"

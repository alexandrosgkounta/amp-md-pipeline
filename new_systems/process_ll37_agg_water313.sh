#!/usr/bin/env bash
# Four LL-37 (2K6O) in water, 313 K, 1 us
# Put the run's tpr and xtc (or xtc parts) in:  AMP_MD_processed/ll37/12_ll37x4_water_313K
# then:  bash pipeline/new_systems/process_ll37_agg_water313.sh      (TPR=... XTC=... to name the files explicitly)
exec bash "$(dirname "$0")/process_system.sh" ll37_agg_water313 ll37 agg_water_313K r1 solution 4 ll37/12_ll37x4_water_313K "${TPR:-}" "${XTC:-}"

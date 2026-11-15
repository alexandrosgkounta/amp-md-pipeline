#!/usr/bin/env bash
# Four melittins (2MLT) in water, 300 K, 1 us
# Put the run's tpr and xtc (or xtc parts) in:  AMP_MD_processed/mel/11_mel4_water_300K
# then:  bash pipeline/new_systems/process_mel_agg_water300.sh      (TPR=... XTC=... to name the files explicitly)
exec bash "$(dirname "$0")/process_system.sh" mel_agg_water300 mel agg_water_300K r1 solution 4 mel/11_mel4_water_300K "${TPR:-}" "${XTC:-}"

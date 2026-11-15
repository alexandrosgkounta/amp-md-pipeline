#!/usr/bin/env bash
# LL-37 random coil start A, water, 313 K, 1 us
# Put the run's tpr and xtc (or xtc parts) in:  AMP_MD_processed/ll37/11_ll37_coil_313K/A
# then:  bash pipeline/new_systems/process_ll37_coil313.sh      (TPR=... XTC=... to name the files explicitly)
exec bash "$(dirname "$0")/process_system.sh" ll37_coil313 ll37 coil_313K r1 solution 1 ll37/11_ll37_coil_313K/A "${TPR:-}" "${XTC:-}"

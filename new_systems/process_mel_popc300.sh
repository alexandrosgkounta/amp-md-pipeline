#!/usr/bin/env bash
# Melittin on POPC (r1 pose), 300 K
# Put the run's tpr and xtc (or xtc parts) in:  AMP_MD_processed/mel/12_mel_popc_300K
# then:  bash pipeline/new_systems/process_mel_popc300.sh      (TPR=... XTC=... to name the files explicitly)
exec bash "$(dirname "$0")/process_system.sh" mel_popc300 mel popc_300K r1 membrane 1 mel/12_mel_popc_300K "${TPR:-}" "${XTC:-}"

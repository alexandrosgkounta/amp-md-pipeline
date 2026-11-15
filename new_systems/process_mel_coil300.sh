#!/usr/bin/env bash
# Melittin, random coil start A (mel_coilA_2505ns.pdb), water, 300 K, 1 us
# Put the run's tpr and xtc (or xtc parts) in:  AMP_MD_processed/mel/10_mel_coil_300K/A
# then:  bash pipeline/new_systems/process_mel_coil300.sh      (TPR=... XTC=... to name the files explicitly)
exec bash "$(dirname "$0")/process_system.sh" mel_coil300 mel coil_300K r1 solution 1 mel/10_mel_coil_300K/A "${TPR:-}" "${XTC:-}"

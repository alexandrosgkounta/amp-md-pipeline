#!/usr/bin/env bash
# Melittin, random coil start C (mel_coilC_2847ns.pdb), water, 300 K, 1 us
# Put the run's tpr and xtc (or xtc parts) in:  AMP_MD_processed/mel/10_mel_coil_300K/C
# then:  bash pipeline/new_systems/process_mel_coil300_r3.sh      (TPR=... XTC=... to name the files explicitly)
exec bash "$(dirname "$0")/process_system.sh" mel_coil300_r3 mel coil_300K r3 solution 1 mel/10_mel_coil_300K/C "${TPR:-}" "${XTC:-}"

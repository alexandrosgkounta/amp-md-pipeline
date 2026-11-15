#!/usr/bin/env bash
# Process every new run whose folder already holds a tpr and xtc; the others are listed as skipped.
#   bash pipeline/new_systems/run_all_new.sh          (ANALYSE=0 to process only)
here="$(cd "$(dirname "$0")" && pwd)"
done_=(); skipped=(); failed=()
for f in "$here"/process_*.sh; do
  [[ "$(basename "$f")" == "process_system.sh" ]] && continue
  sid=$(basename "$f" .sh); sid=${sid#process_}
  bash "$f"; rc=$?
  case $rc in 0) done_+=("$sid");; 2) skipped+=("$sid");; *) failed+=("$sid");; esac
done
echo; echo "processed: ${done_[*]:-none}"; echo "skipped (no files yet): ${skipped[*]:-none}"; echo "FAILED: ${failed[*]:-none}"

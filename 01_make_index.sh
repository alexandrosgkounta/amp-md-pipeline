#!/usr/bin/env bash
# Step 1. Index groups for every system, built from the tpr with MDAnalysis
# (peptide-local residue numbers, so 04 with melittin as residues 161-186 works unchanged).
# Writes <dir>/amp_index.ndx (a new name, so your own index.ndx files are never overwritten)
# and and <dir>/amp_index_groups.json (group sizes, check these once).
# Usage: bash 01_make_index.sh [sys_id ...]
source "$(dirname "$0")/lib_common.sh"
while IFS=$'\t' read -r -u 3 sid pep cond rep kind npep dir tpr xtc whole; do
  d="$AMP_ROOT/$dir"
  check_inputs "$sid" "$d" "$tpr" "" || continue
  if [[ "$FORCE" != "1" && -s "$d/amp_index.ndx" ]]; then log "$sid: amp_index.ndx exists, skipped"; continue; fi
  "$PY" "$PIPE_DIR/py/make_index.py" "$sid"
done 3< <(systems_rows "$@")

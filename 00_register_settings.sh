#!/usr/bin/env bash
# Step 0. Register every setting and condition of each run, read from the tpr itself
# (not from the mdp, which may not be the file that was actually used).
#
# Output, per system:  <dir>/settings/tpr_dump.txt      full gmx dump of the run input
#                      <dir>/settings/key_params.txt    the parameters a Methods section needs
#                      <dir>/settings/xtc_check.txt     frame count, spacing and last time
# Combined:            $AMP_ROOT/analysis/system_table.tsv  (composition + key parameters)
#
# Usage: bash 00_register_settings.sh [sys_id ...]
source "$(dirname "$0")/lib_common.sh"
need_gmx
KEYS='integrator|dt|nsteps|tinit|cutoff-scheme|nstlist|coulombtype|coulomb-modifier|rcoulomb|vdw-type|vdwtype|vdw-modifier|rvdw-switch|rvdw|rlist|fourier-spacing|fourierspacing|pme-order|tcoupl|nsttcouple|pcoupl|pcoupltype|tau-p|refcoord-scaling|constraints|constraint-algorithm|lincs-order|nstxout-compressed|compressed-x-precision|comm-mode|nstcomm|ld-seed|gen-vel|gen-temp|gen-seed|continuation|disre|free-energy'

while IFS=$'\t' read -r -u 3 sid pep cond rep kind npep dir tpr xtc whole; do
  d="$AMP_ROOT/$dir"
  check_inputs "$sid" "$d" "$tpr" "$xtc" || continue
  mkdir -p "$d/settings"
  if [[ "$FORCE" != "1" && -s "$d/settings/key_params.txt" && -s "$d/settings/xtc_check.txt" ]]; then
    log "$sid: settings already registered, skipped"; continue
  fi
  log "$sid: dumping $tpr"
  "$GMX" dump -s "$d/$tpr" > "$d/settings/tpr_dump.txt" 2> "$d/settings/tpr_dump.log"
  {
    echo "# $sid  ($dir/$tpr)"
    grep -E "^[[:space:]]+($KEYS)[[:space:]]+=" "$d/settings/tpr_dump.txt" | sed 's/^ *//' | awk '!seen[$1]++'
    # temperature-coupling groups and reference values (grpopts block)
    grep -E "^[[:space:]]+(ref-t|tau-t):" "$d/settings/tpr_dump.txt" | sed 's/^ *//' | head -2
    # pressure: reference and compressibility tensors (3 rows each)
    grep -A3 -E "^[[:space:]]+(ref-p|compressibility) \(3x3\):" "$d/settings/tpr_dump.txt" | sed 's/^ *//'
  } > "$d/settings/key_params.txt"
  log "$sid: checking $xtc (reads the whole file, a few minutes for tens of GB)"
  "$GMX" check -f "$d/$xtc" > "$d/settings/xtc_check.txt" 2>&1 || true
  grep -E "^(Step|Time|Coords|Box|Last frame|Reading frame)" "$d/settings/xtc_check.txt" | tail -4 || true
done 3< <(systems_rows "$@")
log "composition and table"
"$PY" "$PIPE_DIR/py/system_table.py" "$@" || echo "WARNING: system table failed (see the error above); processing continues"

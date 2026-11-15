#!/usr/bin/env bash
# Step 3. Python analysis per system, then the cross-condition figures and summary table.
# Needs: python3 with MDAnalysis >= 2.8 (DSSP module), numpy, pandas, matplotlib.
#   pip install "MDAnalysis>=2.8" pandas matplotlib
#
# Usage: bash 03_run_analysis.sh [sys_id ...]
# Settings: STRIDE=1 (use every frame of dry.xtc), PARALLEL=1 (run systems at the same time),
#           WINDOW_FRAC=0.5 (averages over the last half), START_NS= (or a fixed start, in ns)
source "$(dirname "$0")/lib_common.sh"
STRIDE="${STRIDE:-1}"; PARALLEL="${PARALLEL:-0}"; WINDOW_FRAC="${WINDOW_FRAC:-0.5}"; START_NS="${START_NS:-}"
mkdir -p "$AMP_ROOT/analysis/logs"
peptides=()
while IFS=$'\t' read -r -u 3 sid pep cond rep kind npep dir tpr xtc whole; do
  peptides+=("$pep")
  d="$AMP_ROOT/$dir"
  [[ -s "$d/processed/dry.xtc" ]] || { echo "WARNING $sid: processed/dry.xtc missing, run step 2 first; skipped" >&2; continue; }
  if [[ "$FORCE" != "1" && -s "$AMP_ROOT/analysis/$sid/meta.json" ]]; then log "$sid: analysis exists, skipped"; continue; fi
  logf="$AMP_ROOT/analysis/logs/${sid}.log"
  log "$sid: analysis (log: analysis/logs/${sid}.log)"
  if [[ "$PARALLEL" == "1" ]]; then
    "$PY" "$PIPE_DIR/py/amp_analysis.py" "$sid" --stride "$STRIDE" --end-ns "$(end_ns_for "$sid")" > "$logf" 2>&1 &
  else
    "$PY" "$PIPE_DIR/py/amp_analysis.py" "$sid" --stride "$STRIDE" --end-ns "$(end_ns_for "$sid")" 2>&1 | tee "$logf" | grep -E "done|frames|ERROR|Error" || true
  fi
done 3< <(systems_rows "$@")
wait
for pep in $(printf '%s\n' "${peptides[@]}" | sort -u); do
  log "figures for $pep"
  extra=(); [[ -n "$START_NS" ]] && extra=(--start-ns "$START_NS")
  "$PY" "$PIPE_DIR/py/amp_figures.py" "$pep" --window-frac "$WINDOW_FRAC" "${extra[@]}"
  log "extension check for $pep (is 0 to $END_NS ns enough for every run?)"
  "$PY" "$PIPE_DIR/py/check_extension.py" "$pep" || true
done

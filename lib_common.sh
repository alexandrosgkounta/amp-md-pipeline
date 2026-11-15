#!/usr/bin/env bash
# Shared helpers for the AMP batch scripts. Sourced, not run.
set -euo pipefail

PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export AMP_ROOT="${AMP_ROOT:-$(cd "$PIPE_DIR/.." && pwd)}"      # the amp_simulations folder
export AMP_SYSTEMS="${AMP_SYSTEMS:-$PIPE_DIR/systems.tsv}"
GMX="${GMX:-gmx}"                       # e.g. GMX=gmx_mpi on the cluster
export GMX_MAXBACKUP=-1                 # overwrite instead of keeping #file.1# backups (they cost GBs here)
PY="${PY:-python3}"
DT_ANALYSIS="${DT_ANALYSIS:-200}"       # ps between frames of the dry analysis trajectory (3 us -> 15000 frames)
DT_WET="${DT_WET:-1000}"                # ps between frames of the full system (water) trajectory
DT_VIZ="${DT_VIZ:-1000}"                # ps between frames of the visualisation trajectory
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-0}"

# One pipeline at a time: a second copy started while one is running would write the same files.
LOCK="$AMP_ROOT/.amp_pipeline.lock"
if [[ -z "${AMP_LOCK_HELD:-}" ]]; then
  if mkdir "$LOCK" 2>/dev/null; then
    echo $$ > "$LOCK/pid"; export AMP_LOCK_HELD=1
    trap 'rm -rf "$LOCK"' EXIT
  else
    other=$(cat "$LOCK/pid" 2>/dev/null || echo "?")
    if [[ "$other" != "?" ]] && ! kill -0 "$other" 2>/dev/null; then
      rm -rf "$LOCK"; mkdir "$LOCK"; echo $$ > "$LOCK/pid"; export AMP_LOCK_HELD=1; trap 'rm -rf "$LOCK"' EXIT
    else
      echo "ERROR: the pipeline is already running (pid $other). Stop it first, or wait." >&2; exit 1
    fi
  fi
fi

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

# systems_rows [sys_id ...] -> prints matching TSV rows (all rows if no ids given)
systems_rows() {
  awk -F'\t' -v ids="$*" '
    BEGIN{n=split(ids,a," "); for(i=1;i<=n;i++) want[a[i]]=1}
    /^#/ || NF<10 {next}
    $1=="sys_id" {next}
    (n==0 || ($1 in want)) {print}
  ' "$AMP_SYSTEMS"
}

# Steps are skipped for a system when their output already exists (resume after a stop).
# gsel GROUP... : feed group names to an interactive gmx prompt
gsel() { printf '%s\n' "$@"; }

# check_inputs SID DIR TPR XTC : warn and return 1 if the run input or trajectory is missing
check_inputs() {
  local sid="$1" d="$2" tpr="$3" xtc="$4" ok=0
  [[ -d "$d" ]]       || { echo "WARNING $sid: folder $d not found, skipped" >&2; return 1; }
  [[ -f "$d/$tpr" ]]  || { echo "WARNING $sid: $tpr not found in $d, skipped" >&2; ok=1; }
  if [[ -n "$xtc" && ! -e "$d/$xtc" ]]; then echo "WARNING $sid: $xtc not found in $d, skipped" >&2; ok=1; fi
  return $ok
}
FORCE="${FORCE:-0}"
END_NS="${END_NS:-1000}"                # analysis period 0..END_NS for every system (1 us)
export END_NS
# end_ns_for SYS_ID : END_NS unless analysis_window.tsv extends this system
end_ns_for() {
  local f="$PIPE_DIR/analysis_window.tsv" v=""
  [[ -f "$f" ]] && v=$(awk -F'\t' -v s="$1" '!/^#/ && $1==s {print $2; exit}' "$f")
  echo "${v:-$END_NS}"
}      # FORCE=1 redoes steps whose outputs already exist

need_gmx() {
  command -v "$GMX" >/dev/null 2>&1 || { echo "ERROR: $GMX not found. Load GROMACS (module load gromacs) or set GMX=..." >&2; exit 1; }
}

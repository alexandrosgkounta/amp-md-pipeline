#!/usr/bin/env bash
# =====================================================================================================
# process_system.sh: process ONE new run exactly as the existing runs in AMP_MD_processed were processed,
# using the pipeline in AMP_MD_processed/pipeline (steps 00, 01, 02 and the per-system analysis of 03).
# Normally called by the process_<sys_id>.sh wrappers in this folder, not directly.
#
#   bash process_system.sh SYS_ID PEPTIDE CONDITION REPLICATE KIND N_PEP DIR [TPR] [XTC]
#
# What it does:
#   1. finds the run input and trajectory in AMP_MD_processed/DIR (or uses TPR / XTC if given);
#      several trajectory parts (md.part0002.xtc ..., step7_1.xtc step7_2.xtc ...) are joined first
#   2. registers the run in pipeline/systems.tsv (adds the row, or updates it if it changed)
#   3. 00_register_settings.sh  -> DIR/settings/ (key_params.txt, xtc_check.txt), analysis/system_table_SYS_ID.*
#      (the existing analysis/system_table.* is kept unchanged)
#   4. 01_make_index.sh         -> DIR/amp_index.ndx, amp_index_groups.json
#   5. 02_process_trajectories.sh -> DIR/processed/ (dry.xtc, viz.xtc, viz_start.pdb, dry.tpr; wet.xtc for
#                                  membrane runs), with the same settings as before (DT_ANALYSIS 200 ps,
#                                  DT_WET 1000 ps, DT_VIZ 1000 ps)
#   6. amp_analysis.py SYS_ID   -> analysis/SYS_ID/ (timeseries.csv, perres.npz, ...), 0 to 1000 ns
#      (set ANALYSE=0 to skip; figures are NOT regenerated, see README.md)
# Every step is skipped when its output already exists (FORCE=1 redoes them), so a stopped run resumes.
# Settings: GMX=/path/to/gmx  PY=python3  ANALYSE=1  FORCE=0
# =====================================================================================================
set -euo pipefail
[[ $# -ge 7 ]] || { sed -n 2,25p "$0"; exit 1; }
HERE="$(cd "$(dirname "$0")" && pwd)"
PIPE="$(cd "$HERE/.." && pwd)"
ROOT="$(cd "$PIPE/.." && pwd)"
SID=$1; PEP=$2; COND=$3; REP=$4; KIND=$5; NPEP=$6; DIR=$7; TPR="${8:-}"; XTC="${9:-}"
D="$ROOT/$DIR"
ANALYSE="${ANALYSE:-1}"
export PY="${PY:-python3}"
say() { printf '[%s] %s: %s\n' "$(date +%H:%M:%S)" "$SID" "$*"; }

# ---- GROMACS (Homebrew or a source build) --------------------------------------------------------
if [[ -z "${GMX:-}" ]]; then
  for c in gmx gmx_mpi gmx_d /opt/homebrew/bin/gmx /usr/local/bin/gmx /usr/local/gromacs/bin/gmx; do
    if command -v "$c" >/dev/null 2>&1; then GMX="$(command -v "$c")"; break; fi
  done
fi
[[ -n "${GMX:-}" ]] || { echo "ERROR: gmx not found (brew install gromacs, or GMX=/path/to/gmx bash ...)" >&2; exit 1; }
export GMX
"$PY" -c "import MDAnalysis" 2>/dev/null || { echo "ERROR: $PY has no MDAnalysis (pip install 'MDAnalysis>=2.8' pandas matplotlib)" >&2; exit 1; }

# ---- inputs ---------------------------------------------------------------------------------------
[[ -d "$D" ]] || { echo "SKIP $SID: folder $D not found. Copy the run's tpr and xtc there first." >&2; exit 2; }
cd "$D"
if [[ -z "$TPR" || ! -f "$TPR" ]]; then
  TPR=""
  for c in md.tpr step7.tpr step7_production.tpr step7_1.tpr step5_production.tpr; do [[ -f "$c" ]] && { TPR="$c"; break; }; done
fi
[[ -n "$TPR" ]] || { echo "SKIP $SID: no run input (md.tpr, step7.tpr, step7_1.tpr ...) in $D; set TPR=..." >&2; exit 2; }
if [[ -z "$XTC" || ! -f "$XTC" ]]; then
  XTC=""
  for c in md_joined.xtc md.xtc step7.xtc step7_production.xtc step5_production.xtc; do [[ -f "$c" ]] && { XTC="$c"; break; }; done
fi
# natural order of parts (works with the old bash and sort of macOS): by the number before .xtc
order_parts() { awk '{n=$0; sub(/\.xtc$/,"",n); sub(/.*[^0-9]/,"",n); print (n==""?0:n+0) "\t" $0}' | sort -n | cut -f2; }
PARTS=""
if [[ -n "$XTC" ]]; then
  # continuation parts written with -noappend next to the first file: md.xtc + md.part0002.xtc ...
  PARTS=$(ls -1 "${XTC%.xtc}".part[0-9]*.xtc 2>/dev/null | order_parts || true)
  [[ -n "$PARTS" ]] && PARTS="$XTC"$'\n'"$PARTS"
else
  # only parts: mdrun -noappend (name.partNNNN.xtc) or the CHARMM-GUI loop (step7_1.xtc, step7_2.xtc ...)
  PARTS=$(ls -1 *.part[0-9]*.xtc step7_[0-9]*.xtc 2>/dev/null | order_parts || true)
fi
if [[ -n "$PARTS" && "$XTC" != "md_joined.xtc" ]]; then
  if [[ -s md_joined.xtc && "${FORCE:-0}" != "1" ]]; then
    say "md_joined.xtc exists, using it"
  else
    say "joining $(echo "$PARTS" | wc -l | tr -d ' ') trajectory parts into md_joined.xtc: $(echo $PARTS)"
    # shellcheck disable=SC2086
    "$GMX" trjcat -f $PARTS -o md_joined.xtc > trjcat.log 2>&1 || { echo "ERROR: trjcat failed, see $D/trjcat.log" >&2; exit 1; }
  fi
  XTC=md_joined.xtc
fi
[[ -n "$XTC" ]] || { echo "SKIP $SID: no trajectory (md.xtc, step7.xtc, parts ...) in $D; set XTC=..." >&2; exit 2; }
say "run input $TPR, trajectory $XTC ($(du -h "$XTC" | cut -f1))"

# ---- register in systems.tsv ----------------------------------------------------------------------
TSV="$PIPE/systems.tsv"
ROW="$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' "$SID" "$PEP" "$COND" "$REP" "$KIND" "$NPEP" "$DIR" "$TPR" "$XTC" 0)"
if awk -F'\t' -v s="$SID" '$1==s{f=1} END{exit !f}' "$TSV"; then
  if ! grep -qxF "$ROW" "$TSV"; then
    cp "$TSV" "$TSV.bak_$(date +%Y%m%d_%H%M%S)"
    awk -F'\t' -v s="$SID" -v r="$ROW" '$1==s{print r; next} {print}' "$TSV" > "$TSV.tmp" && mv "$TSV.tmp" "$TSV"
    say "updated its row in systems.tsv (backup kept)"
  fi
else
  grep -q "^# Lower-temperature runs" "$TSV" || printf '# Lower-temperature runs (melittin 300 K, LL-37 313 K), processed locally with pipeline/new_systems\n' >> "$TSV"
  printf '%s\n' "$ROW" >> "$TSV"
  say "added to systems.tsv"
fi

# ---- stale lock from an earlier stopped run -------------------------------------------------------
LOCK="$ROOT/.amp_pipeline.lock"
if [[ -d "$LOCK" ]]; then
  pid=$(cat "$LOCK/pid" 2>/dev/null || echo "")
  if [[ -z "$pid" ]] || ! ps -p "$pid" -o command= 2>/dev/null | grep -qE "pipeline|_register_|_make_index|_process_traj|_run_analysis"; then
    rm -rf "$LOCK"; say "removed a stale pipeline lock (pid ${pid:-?} is not a pipeline job)"
  fi
fi

# ---- pipeline steps, same settings as the existing runs -------------------------------------------
# step 0 rebuilds analysis/system_table.* with this run only; keep the existing table, put this run's beside it
for e in tsv md; do [[ -f "$ROOT/analysis/system_table.$e" ]] && cp "$ROOT/analysis/system_table.$e" "$ROOT/analysis/.system_table.$e.keep"; done
say "step 0: register settings";    bash "$PIPE/00_register_settings.sh"    "$SID"
for e in tsv md; do
  [[ -f "$ROOT/analysis/system_table.$e" ]] && mv "$ROOT/analysis/system_table.$e" "$ROOT/analysis/system_table_$SID.$e"
  [[ -f "$ROOT/analysis/.system_table.$e.keep" ]] && mv "$ROOT/analysis/.system_table.$e.keep" "$ROOT/analysis/system_table.$e"
done
say "step 1: index groups";         bash "$PIPE/01_make_index.sh"           "$SID"
say "step 2: process trajectories"; bash "$PIPE/02_process_trajectories.sh" "$SID"
[[ -s "$D/processed/dry.xtc" ]] || { echo "ERROR: $D/processed/dry.xtc not written, see $D/processed/process.log" >&2; exit 1; }
"$GMX" check -f "$D/processed/dry.xtc" 2>&1 | grep -oE "Last frame +[0-9]+ time +[0-9.]+" | sed "s/^/[$SID] dry.xtc: /" || true

if [[ "$ANALYSE" == "1" ]]; then
  mkdir -p "$ROOT/analysis/logs"
  if [[ "${FORCE:-0}" != "1" && -s "$ROOT/analysis/$SID/meta.json" ]]; then
    say "analysis exists, skipped"
  else
    say "step 3: per-system analysis, 0 to ${END_NS:-1000} ns (log: analysis/logs/$SID.log)"
    AMP_ROOT="$ROOT" "$PY" "$PIPE/py/amp_analysis.py" "$SID" --stride 1 --end-ns "${END_NS:-1000}" \
      > "$ROOT/analysis/logs/$SID.log" 2>&1 || { echo "ERROR: analysis failed, see analysis/logs/$SID.log" >&2; exit 1; }
    tail -2 "$ROOT/analysis/logs/$SID.log"
  fi
fi
say "done: $D/processed  and  analysis/$SID"

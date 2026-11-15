#!/usr/bin/env bash
# Show what the pipeline is doing right now (safe to run any time, in a second terminal tab).
#   bash progress.sh            once
#   watch -n 30 bash progress.sh   (if 'watch' is installed) or: while true; do clear; bash progress.sh; sleep 30; done
here="$(cd "$(dirname "$0")" && pwd)"; root="$(cd "$here/.." && pwd)"
echo "== running processes"
ps -eo pid,etime,command | grep -E "gmx (trjconv|check|dump|convert-tpr)|amp_analysis|amp_figures|make_index" | grep -v grep \
  | sed -E 's#/Volumes/Transcend/MD_LL37_MEL/##g' | cut -c1-160 || true
pgrep -f run_all.sh >/dev/null && echo "run_all.sh: running" || echo "run_all.sh: NOT running"
echo
echo "== processed/ status per system"
awk -F'\t' '!/^#/ && $1!="sys_id"{print $1"\t"$7}' "$here/systems.tsv" | while IFS=$'\t' read -r sid dir; do
  p="$root/$dir/processed"
  if [[ -s "$p/viz.xtc" && -s "$p/dry.xtc" && ! -e "$p/center.xtc" ]]; then st="done"
  elif [[ -d "$p" ]]; then
    newest=$(ls -t "$p" 2>/dev/null | grep -v '^#' | grep -v process.log | head -1)
    last=$(tr '\r' '\n' < "$p/process.log" 2>/dev/null | grep "Reading frame" | tail -1 | awk '{printf "%.2f us", $5/1e6}')
    st="working: writing ${newest:-?}, read up to ${last:-?}"
  else st="waiting"; fi
  [[ -e "$root/analysis/$sid/meta.json" ]] && st="$st, analysis done"
  printf '  %-12s %s\n' "$sid" "$st"
done

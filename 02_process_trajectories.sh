#!/usr/bin/env bash
# Step 2. PBC treatment and trajectory reduction with GROMACS.
# Requires amp_index.ndx from 01_make_index.sh.
#
# All frames are first thinned to one every DT_ANALYSIS ps, so no full-resolution copy is written.
# Every system  : whole.xtc        molecules made whole (skipped if whole=1 in systems.tsv)
# solution      : center.xtc       nojump, then peptide centred, -pbc mol
# membrane/agg  : center.xtc       membrane centred in the box (-pbc res), then peptides made whole again
# then          : dry.xtc          Pep_Memb group (peptide + lipids, no water), every DT_ANALYSIS ps -> Python analysis
#                 wet.xtc          full system every DT_WET ps -> water penetration and water density
#                 viz.xtc/.pdb     Pep_Memb every DT_VIZ ps (solution: fitted on backbone) -> VMD / PyMOL
# whole.xtc and center.xtc are deleted at the end unless KEEP_INTERMEDIATE=1.
# DT_WET and DT_VIZ must be multiples of DT_ANALYSIS.
#
# Usage: bash 02_process_trajectories.sh [sys_id ...]
# Settings: GMX=gmx_mpi DT_ANALYSIS=100 DT_WET=1000 DT_VIZ=1000 KEEP_INTERMEDIATE=0
source "$(dirname "$0")/lib_common.sh"
need_gmx

while IFS=$'\t' read -r -u 3 sid pep cond rep kind npep dir tpr xtc whole; do
  d="$AMP_ROOT/$dir"; out="$d/processed"
  check_inputs "$sid" "$d" "$tpr" "$xtc" || continue
  if [[ "$FORCE" != "1" && -s "$out/dry.xtc" && -s "$out/viz.xtc" ]]; then log "$sid: processed/ exists, skipped"; continue; fi
  mkdir -p "$out"
  ndx="$d/amp_index.ndx"; s="$d/$tpr"
  [[ -f "$ndx" ]] || { echo "WARNING $sid: $ndx missing, run 01_make_index.sh first; skipped" >&2; continue; }
  L="$out/process.log"; : > "$L"
  log "$sid ($kind): start, log in processed/process.log"

  # 1. whole molecules
  if [[ "$whole" == "1" ]]; then
    ln -sf "../$xtc" "$out/whole.xtc"
  else
    gsel System | "$GMX" trjconv -s "$s" -f "$d/$xtc" -n "$ndx" -o "$out/whole.xtc" -pbc whole -dt "$DT_ANALYSIS" >>"$L" 2>&1
  fi

  # 2. centring
  if [[ "$kind" == "solution" ]]; then
    gsel System | "$GMX" trjconv -s "$s" -f "$out/whole.xtc" -n "$ndx" -o "$out/nojump.xtc" -pbc nojump -dt "$DT_ANALYSIS" >>"$L" 2>&1
    gsel Peptide System | "$GMX" trjconv -s "$s" -f "$out/nojump.xtc" -n "$ndx" -o "$out/center.xtc" \
         -center -pbc mol >>"$L" 2>&1
    [[ "$KEEP_INTERMEDIATE" == "1" ]] || rm -f "$out/nojump.xtc"
  else
    # The lipids in the pdb2gmx topologies form ONE moleculetype ("Other"), so -pbc mol would treat
    # the whole bilayer as one molecule. -pbc res wraps lipid by lipid; the second pass re-joins the peptides.
    gsel Membrane System | "$GMX" trjconv -s "$s" -f "$out/whole.xtc" -n "$ndx" -o "$out/center_res.xtc" \
         -center -pbc res -dt "$DT_ANALYSIS" >>"$L" 2>&1
    gsel System | "$GMX" trjconv -s "$s" -f "$out/center_res.xtc" -n "$ndx" -o "$out/center.xtc" -pbc whole >>"$L" 2>&1
    [[ "$KEEP_INTERMEDIATE" == "1" ]] || rm -f "$out/center_res.xtc"
  fi

  # 3. reduced trajectories
  gsel Pep_Memb | "$GMX" trjconv -s "$s" -f "$out/center.xtc" -n "$ndx" -o "$out/dry.xtc" -dt "$DT_ANALYSIS" >>"$L" 2>&1
  if [[ "$kind" != "solution" ]]; then
    gsel System | "$GMX" trjconv -s "$s" -f "$out/center.xtc" -n "$ndx" -o "$out/wet.xtc" -dt "$DT_WET" >>"$L" 2>&1
  fi
  if [[ "$kind" == "solution" ]]; then
    gsel Peptide_BB Pep_Memb | "$GMX" trjconv -s "$s" -f "$out/center.xtc" -n "$ndx" -o "$out/viz.xtc" \
         -fit rot+trans -dt "$DT_VIZ" >>"$L" 2>&1
  else
    gsel Pep_Memb | "$GMX" trjconv -s "$s" -f "$out/center.xtc" -n "$ndx" -o "$out/viz.xtc" -dt "$DT_VIZ" >>"$L" 2>&1
  fi
  gsel Pep_Memb | "$GMX" trjconv -s "$s" -f "$out/viz.xtc" -n "$ndx" -o "$out/viz_start.pdb" -dump 0 >>"$L" 2>&1 || true
  # dry run-input file, so the reduced trajectory can also be used with gmx tools
  gsel Pep_Memb | "$GMX" convert-tpr -s "$s" -n "$ndx" -o "$out/dry.tpr" >>"$L" 2>&1 || true

  "$GMX" check -f "$out/dry.xtc" 2>&1 | grep -E "^(Coords|Last frame)" | tee -a "$L" || true
  if [[ "$KEEP_INTERMEDIATE" != "1" ]]; then
    rm -f "$out/center.xtc"; [[ -L "$out/whole.xtc" ]] || rm -f "$out/whole.xtc"
  fi
  log "$sid: done -> $out"
done 3< <(systems_rows "$@")

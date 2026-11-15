#!/usr/bin/env python3
"""Recompute the RMSD columns (rmsd_bb, rmsd_hn, rmsd_hc) and RMSF of an existing analysis,
with frame 0 of dry.xtc as the reference (the same definitions as amp_analysis.py).

Used for the aggregation runs, which are exactly 1 us long, so their October analysis already
covers 0 to 1 us and only the RMSD reference needed correcting.
Usage: python rmsd_refix.py SYS_ID SRC_ANALYSIS_DIR   (writes to analysis/SYS_ID)
"""
import os, sys, json, shutil, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from MDAnalysis.lib.mdamath import make_whole
from MDAnalysis.analysis.rms import rmsd as _rmsd
from MDAnalysis.analysis.align import rotation_matrix
from amp_config import get_system, PEPTIDES, ANALYSIS_DIR
from amp_analysis import load_dry, copies_of

sid, src = sys.argv[1], sys.argv[2]
s = get_system(sid); pep = PEPTIDES[s["peptide"]]; nres = len(pep["seq"])
out = os.path.join(ANALYSIS_DIR, sid)
if not os.path.exists(out): shutil.copytree(os.path.join(src, sid), out)
full, u, ref, ndx = load_dry(s)
cp = copies_of(u, s["n_pep"], nres); rc = copies_of(ref, s["n_pep"], nres)
for c in rc: make_whole(c)
segs = pep["segments"]
def sel(c):
    R = c.residues
    return (c.select_atoms("backbone"), R[segs[0][1] - 1:segs[0][2]].atoms.select_atoms("backbone"),
            R[segs[-1][1] - 1:segs[-1][2]].atoms.select_atoms("backbone"))
S = [sel(c) for c in cp]; RS = [[g.positions.copy() for g in sel(c)] for c in rc]
ts_df = pd.read_csv(os.path.join(out, "timeseries.csv"))
times = np.sort(ts_df["time_ns"].unique())
idx = np.rint((times * 1000 - u.trajectory[0].time) / u.trajectory.dt).astype(int)
rows = {k: [] for k in range(len(cp))}; bbt = [np.zeros((len(idx), S[k][0].n_atoms, 3), np.float32) for k in range(len(cp))]
for fi, f in enumerate(idx):
    u.trajectory[f]
    for k, c in enumerate(cp):
        make_whole(c)
        rows[k].append([_rmsd(S[k][j].positions, RS[k][j], center=True, superposition=True) / 10 for j in range(3)])
        bbt[k][fi] = S[k][0].positions
for k in range(len(cp)):
    m = ts_df["copy"] == k + 1
    ts_df.loc[m, ["rmsd_bb", "rmsd_hn", "rmsd_hc"]] = np.array(rows[k])
ts_df.to_csv(os.path.join(out, "timeseries.csv"), index=False, float_format="%.4f")
pr = dict(np.load(os.path.join(out, "perres.npz"), allow_pickle=True))
for k in range(len(cp)):
    X = bbt[k].astype(np.float64); ca = S[k][0].names == "CA"; r_ = RS[k][0] - RS[k][0].mean(0)
    for _ in range(2):
        Y = np.empty_like(X)
        for f in range(len(X)):
            xc = X[f] - X[f].mean(0); Rm, _r = rotation_matrix(xc, r_); Y[f] = xc @ Rm.T
        r_ = Y.mean(0)
    pr["rmsf"][k] = np.sqrt(((Y[:, ca] - r_[ca]) ** 2).sum(-1).mean(0)) / 10
np.savez_compressed(os.path.join(out, "perres.npz"), **pr)
meta = json.load(open(os.path.join(out, "meta.json"))); meta["end_ns_applied"] = 1000.0
meta["note"] = "October analysis (0 to 1 us); RMSD and RMSF recomputed with frame 0 of dry.xtc as reference"
json.dump(meta, open(os.path.join(out, "meta.json"), "w"), indent=1, default=str)
print(sid, "done", len(idx), "frames")

"""Amide hydration for the melittin HDX comparison: per residue, the fraction of frames (500 to 1000 ns) in which the
backbone amide H is hydrogen bonded to water (H...O(water) <= 2.4 A, the Best and Vendruscolo H-bond cut-off) and the
mean number of water oxygens within 3.5 A of the amide N. Uses processed/wet.xtc (1 ns). Writes analysis/hdx_compare/hydration_<sid>.csv"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
sys.path.insert(0, os.path.dirname(__file__))
from amp_config import get_system, ANALYSIS_DIR
OUT = os.path.join(ANALYSIS_DIR, "hdx_compare")
for sid in sys.argv[1:]:
    s = get_system(sid); d = s["path"]
    u = mda.Universe(os.path.join(d, s["tpr"]), os.path.join(d, "processed", "wet.xtc"))
    prot = u.select_atoms("protein"); frs = sorted(prot.fragments, key=lambda a: a.ix[0]); ncp = s["n_pep"]
    W = u.select_atoms("resname SOL TIP3 TIP3P HOH WAT and name OW OH2")
    rows = []
    for k, c in enumerate(frs[:ncp]):
        res = c.residues; r0 = res.resids[0]
        H = c.select_atoms("name HN H and not resname PRO"); N = c.select_atoms("name N and not resname PRO")
        hres = H.resids - r0 + 1; nres = N.resids - r0 + 1
        nf = 0; hb = np.zeros(len(H)); nw = np.zeros(len(N))
        for ts in u.trajectory:
            t = ts.time / 1000
            if t < 500 - 1e-6 or t > 1000 + 1e-6: continue
            nf += 1
            p = capped_distance(H.positions, W.positions, 2.4, box=ts.dimensions, return_distances=False)
            if len(p): hb[np.unique(p[:, 0])] += 1
            p = capped_distance(N.positions, W.positions, 3.5, box=ts.dimensions, return_distances=False)
            if len(p): np.add.at(nw, p[:, 0], 1)
        for i, r in enumerate(hres):
            j = np.where(nres == r)[0]
            rows.append(dict(copy=k + 1, resid=int(r), water_hbond_frac=hb[i] / nf, waters_near_N=float(nw[j[0]] / nf) if len(j) else np.nan))
    df = pd.DataFrame(rows); df.to_csv(os.path.join(OUT, f"hydration_{sid}.csv"), index=False)
    print(sid, nf, "frames", df.groupby("resid").water_hbond_frac.mean().round(2).values)

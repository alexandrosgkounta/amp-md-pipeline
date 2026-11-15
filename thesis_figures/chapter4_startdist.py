"""SI check: how far from the bilayer does each single-peptide membrane run start, and when does it first touch lipid?
Minimum heavy-atom distance peptide to lipid (minimum image) and COM height above the proximal phosphate plane,
0 to 50 ns at 0.2 ns. Writes analysis/start_distance.csv and figures/chapter4/SI/Fig_S4_2_start_distance.*"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from amp_config import get_system, ANALYSIS_DIR, FIG_DIR
from amp_analysis import load_dry
from MDAnalysis.lib.distances import distance_array
TMAX = float(os.environ.get("TMAX", 50))
sids = sys.argv[1:]
for sid in sids:
    s = get_system(sid); full, u, ref, ndx = load_dry(s)
    pep = u.select_atoms("protein and not name H*"); lip = u.select_atoms("resname POPC POPE POPG and not name H*")
    out = []
    for ts in u.trajectory:
        if ts.time / 1000 > TMAX + 1e-6: break
        d = distance_array(pep.positions, lip.positions, box=ts.dimensions).min() / 10
        out.append((ts.time / 1000, d))
    df = pd.DataFrame(out, columns=["time_ns", "min_dist_nm"])
    tsr = pd.read_csv(os.path.join(ANALYSIS_DIR, sid, "timeseries.csv"))[["time_ns", "height_com"]]
    df = pd.merge_asof(df.sort_values("time_ns"), tsr.sort_values("time_ns"), on="time_ns", direction="nearest")
    df.to_csv(os.path.join(ANALYSIS_DIR, sid, "start_distance.csv"), index=False)
    fc = df.time_ns[df.min_dist_nm <= 0.4]
    print(f"{sid}: t0 min dist {df.min_dist_nm.iloc[0]:.2f} nm, COM height {df.height_com.iloc[0]:.2f} nm, first contact {fc.iloc[0] if len(fc) else 'none <=%g' % TMAX} ns, frames {len(df)}")

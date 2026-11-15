"""Chapter 4 v2 structure panels: time series in water (aligned on helix C of the representative frame, helix C
along x, helix N towards +y), membrane time series, and the representative bound frame of every POPC run.
Same file format and orientation conventions as make_vmd_pdbs.py. Usage: python3 make_vmd_pdbs_v2.py <job> ..."""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, MDAnalysis as mda
sys.path.insert(0, os.path.dirname(__file__))
import make_vmd_pdbs as M
from amp_config import get_system
from amp_analysis import load_dry
from MDAnalysis.analysis import align
OUT = os.path.join(os.path.dirname(M.OUT), "pdb_v2"); os.makedirs(OUT, exist_ok=True)
HC = {"mel": (13, 20), "ll37": (20, 29)}; HN = {"mel": (2, 10), "ll37": (2, 13)}

def water_series(tag, sid, pep, times):
    s = get_system(sid); full, u, ref, ndx = load_dry(s)
    c = M.peptide_copies(u, 1)[0]; lo, hi = HC[pep]; nlo, nhi = HN[pep]
    frame0 = None
    for t in times:
        M.goto(u, t); X = c.positions.copy(); ca = c.select_atoms("name CA").positions
        hc = ca[lo - 1:hi]
        if frame0 is None:            # axes from the first time given (the representative frame)
            cen = hc.mean(0); v = np.linalg.svd(hc - cen, full_matrices=False)[2][0]
            if np.dot(v, hc[-1] - hc[0]) < 0: v = -v
            w = ca[nlo - 1:nhi].mean(0) - cen; w = w - np.dot(w, v) * v; w /= np.linalg.norm(w)
            R = np.vstack([v, w, np.cross(v, w)]); frame0 = (cen, R, (hc - cen) @ R.T)
            Xn = (X - cen) @ R.T
        else:
            cen, R, ref_hc = frame0
            Rm, _ = align.rotation_matrix(hc - hc.mean(0), ref_hc)
            Xn = (X - hc.mean(0)) @ np.asarray(Rm).T
        keep = np.zeros(u.atoms.n_atoms, bool); keep[c.ix] = True
        F = np.zeros((u.atoms.n_atoms, 3)); F[c.ix] = Xn
        name = f"{tag}_{t:06.1f}ns".replace(".", "p"); M.write_pdb(u, keep, F, [c], pep, os.path.join(OUT, name + ".pdb"))
        print(name)

def mem_frames(tag, sid, pep, times):
    s = get_system(sid); full, u, ref, ndx = load_dry(s)
    for t in times:
        M.goto(u, t); keep, F, copies, flip = M.orient(u, pep, 1)
        name = f"{tag}_{t:06.1f}ns".replace(".", "p"); M.write_pdb(u, keep, F, copies, pep, os.path.join(OUT, name + ".pdb")); print(name, "flipped" if flip else "")

def run_reps(pep):
    for r in ["", "_r2", "_r3", "_r4", "_r5"]:
        sid = f"{pep}_popc{r}"; t, feat = M.pick_frame(sid, pep, need_bound=True)
        mem_frames(f"R_{pep}_popc{r or '_r1'}", sid, pep, [t]); print("  ", sid, t, feat)

JOBS = {
  "W_mel": lambda: water_series("W_mel", "mel_sol", "mel", [844.6, 0.0, 50.0, 200.0, 500.0]),
  "W_ll37": lambda: water_series("W_ll37", "ll37_sol", "ll37", [675.6, 0.0, 50.0, 200.0, 500.0]),
  "P_mel_r2": lambda: mem_frames("P_mel_r2", "mel_popc_r2", "mel", [100.0, 200.0, 300.0, 400.0, 796.0]),
  "P_ll37_r5": lambda: mem_frames("P_ll37_r5", "ll37_popc_r5", "ll37", [100.0, 400.0, 600.0, 800.0, 1000.0]),
  "R_mel": lambda: run_reps("mel"), "R_ll37": lambda: run_reps("ll37"),
}
if __name__ == "__main__":
    for j in sys.argv[1:]: JOBS[j]()

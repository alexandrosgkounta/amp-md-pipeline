"""Representative frames for the Chapter 4 structure figures, written as pre-oriented PDBs for amp_vmd.tcl.

Frame choice (500 to 1000 ns, the averaging window): among frames in the state the text describes (bound, if a
membrane run), the frame that is closest to the window medians of helicity, COM height, tilt, deepest Ca and Rg
and whose per-residue helix pattern (DSSP H) is closest to the window majority pattern (sum of squared robust
z scores). Overlays: frames every 50 ns, aligned on the C-terminal helix Ca.

Orientation (identical for every file, so one camera fits all):  bilayer normal = +y (up), upper phosphate
plane at y = 0, the peptide on the upper leaflet, its helix axis along +x, the peptide COM at x = 0, z = 0.
Lipids are re-wrapped around the peptide and kept within |x| < XW and |z| < ZW (A).
Segment names hold the peptide regions (melittin HN HNG HC; LL-37 HN G14 HC CT), chain = copy (A to D), lipids
chain M, segname MEMB; peptide residues renumbered 1..n per copy.
Usage: python3 make_vmd_pdbs.py <job> [<job> ...]    (jobs listed in JOBS below)"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, MDAnalysis as mda
from MDAnalysis.analysis import align
sys.path.insert(0, os.path.dirname(__file__))
from amp_config import get_system, ANALYSIS_DIR, FIG_DIR, PEPTIDES
from amp_analysis import load_dry
OUT = os.path.join(FIG_DIR, "chapter4", "vmd_pack", "pdb"); os.makedirs(OUT, exist_ok=True)
TMPD = os.path.expanduser("~/vmdtmp"); os.makedirs(TMPD, exist_ok=True)
LOG = os.path.join(FIG_DIR, "chapter4", "vmd_pack", "frames_log_raw.csv")
W0, W1 = 500.0, 1000.0
REG = {"mel": [("HN", 1, 10), ("HNG", 11, 14), ("HC", 15, 26)],
       "ll37": [("HN", 1, 13), ("G14", 14, 16), ("HC", 17, 31), ("CT", 32, 37)]}
HC = {"mel": (15, 25), "ll37": (17, 30)}

def seg_of(pep, r):
    for n, lo, hi in REG[pep]:
        if lo <= r <= hi: return n
    return "PEP"

def pick_frame(sid, pep, need_bound=True):
    ts = pd.read_csv(os.path.join(ANALYSIS_DIR, sid, "timeseries.csv"))
    pr = np.load(os.path.join(ANALYSIS_DIR, sid, "perres.npz"), allow_pickle=True)
    t = pr["time_ns"]; H = (pr["dssp"][:, 0, :] == "H").astype(float)
    ts = ts[ts["copy"] == 1].reset_index(drop=True)
    m = (t >= W0) & (t <= W1)
    memf = os.path.join(ANALYSIS_DIR, sid, "membrane_ts.csv")
    if need_bound and os.path.exists(memf):
        mem = pd.read_csv(memf); m &= mem.lip_contacts_copy1.values > 0
    idx = np.where(m)[0]
    cols = [c for c in ("helicity", "height_com", "tilt", "depth_min", "rg") if c in ts]
    Z = np.zeros(len(idx))
    for c in cols:
        x = ts[c].values[idx]; med = np.nanmedian(x); s = (np.nanpercentile(x, 75) - np.nanpercentile(x, 25)) / 1.35 or 1
        Z += ((x - med) / s) ** 2
    maj = (H[idx].mean(0) >= 0.5).astype(float); ham = np.abs(H[idx] - maj).sum(1)
    Z += (ham / max(ham.std(), 1e-3)) ** 2
    k = idx[np.argmin(Z)]
    return float(t[k]), {c: round(float(ts[c].values[k]), 3) for c in cols}

def peptide_copies(u, n_pep):
    prot = u.select_atoms("protein")
    return [f for f in sorted(prot.fragments, key=lambda a: a.ix[0])][:n_pep]

def goto(u, t):
    dt = u.trajectory.dt; i = int(round((t - u.trajectory[0].time / 1000) * 1000 / dt))
    u.trajectory[i]; assert abs(u.trajectory.time / 1000 - t) < 0.05, (u.trajectory.time, t)

def orient(u, pep, n_pep, center_copies=None, xw=45.0, zw=30.0):
    """Return (atomgroup to write, coordinates in the figure frame)."""
    L = u.dimensions[:3]; X = u.atoms.positions.copy()
    copies = peptide_copies(u, n_pep); lip = u.select_atoms("resname POPC POPE POPG")
    P = u.select_atoms("resname POPC POPE POPG and name P")
    cc = center_copies if center_copies is not None else list(range(n_pep))
    # make the chosen copies one cluster: place each copy at the periodic image (27 candidates) that brings it
    # closest (minimum CA-CA distance) to a copy already placed; then put the cluster COM in the primary cell (z)
    import itertools
    SH = np.array(list(itertools.product((-1, 0, 1), repeat=3))) * L
    placed = [cc[0]]; todo = [i for i in cc[1:]]
    while todo:
        best = None
        for j in todo:
            cj = X[copies[j].select_atoms("name CA").ix]
            for i in placed:
                ci = X[copies[i].select_atoms("name CA").ix]
                for sh in SH:
                    d = np.min(np.linalg.norm(ci[:, None, :] - (cj + sh)[None, :, :], axis=2))
                    if best is None or d < best[0]: best = (d, j, sh)
        d, j, sh = best; X[copies[j].ix] += sh; placed.append(j); todo.remove(j)
    zc0 = X[P.ix, 2].mean(); cz = np.mean([X[copies[i].ix].mean(0) for i in cc], axis=0)
    shz = -L[2] * np.round((cz[2] - zc0) / L[2])
    for i in cc: X[copies[i].ix, 2] += shz
    ref = np.mean([X[copies[i].ix].mean(0) for i in cc], axis=0)
    for k, c in enumerate(copies):          # other copies: nearest image in xy, primary cell in z
        if k in cc: continue
        d = X[c.ix].mean(0) - ref; sh = np.zeros(3); sh[:2] = -L[:2] * np.round(d[:2] / L[:2])
        sh[2] = -L[2] * np.round((X[c.ix, 2].mean() - zc0) / L[2]); X[c.ix] += sh
    cen = ref
    # wrap lipids by residue around the peptide centre (xy)
    for r in lip.residues:
        d = X[r.atoms.ix].mean(0) - cen; sh = np.zeros(3); sh[:2] = -L[:2] * np.round(d[:2] / L[:2]); X[r.atoms.ix] += sh
    zP = X[P.ix, 2]; zc = zP.mean()
    pepcom_z = np.mean([X[copies[i].ix, 2].mean() for i in cc])
    flip = pepcom_z < zc
    if flip:                      # put the peptide on the upper leaflet: rotate 180 deg about x
        X[:, 1] = 2 * cen[1] - X[:, 1]; X[:, 2] = 2 * zc - X[:, 2]
    zP = X[P.ix, 2]; up = zP > zc; zup = zP[up].mean()
    # in-plane direction: helix axis of the (first) chosen copy, or the principal axis of all chosen CA
    ca = np.vstack([X[copies[i].select_atoms("name CA").ix] for i in cc])
    if len(cc) == 1:
        A = PEPTIDES[pep]["axis_all"]; c0 = copies[cc[0]].select_atoms("name CA")
        Y = X[c0.ix][A[0] - 1:A[1]]
    else:
        Y = ca
    Yc = Y[:, :2] - Y[:, :2].mean(0); v = np.linalg.svd(Yc, full_matrices=False)[2][0]
    if len(cc) == 1 and np.dot(v, Y[-1, :2] - Y[0, :2]) < 0: v = -v
    th = np.arctan2(v[1], v[0]); R = np.array([[np.cos(-th), -np.sin(-th)], [np.sin(-th), np.cos(-th)]])
    xy = (X[:, :2] - cen[:2]) @ R.T
    # figure frame: x = along the helix, y = bilayer normal (upper P plane at 0), z = -in-plane perpendicular
    F = np.column_stack([xy[:, 0], X[:, 2] - zup, -xy[:, 1]])
    keep = np.zeros(len(X), bool)
    for c in copies: keep[c.ix] = True
    for r in lip.residues:
        g = F[r.atoms.ix].mean(0)
        if abs(g[0]) < xw and abs(g[2]) < zw: keep[r.atoms.ix] = True
    return keep, F, copies, flip

def write_pdb(u, keep, F, copies, pep, path, chains="ABCD"):
    sel = u.atoms[keep]; ix = sel.ix
    w = mda.Merge(sel); w.atoms.positions = F[ix]
    cmap = {}
    for k, c in enumerate(copies):
        r0 = c.residues.resids[0]
        for r in c.residues: cmap[r.ix] = (chains[k], r.resid - r0 + 1)
    tmp = os.path.join(TMPD, os.path.basename(path))
    w.atoms.write(tmp)
    src = u.atoms[keep]; lines = open(tmp).read().splitlines(); out = []; a = 0; lipres = {}
    for ln in lines:
        if ln.startswith(("ATOM", "HETATM")):
            at = src[a]; a += 1
            if at.residue.ix in cmap:
                ch, rid = cmap[at.residue.ix]; seg = seg_of(pep, rid)
            else:
                ch = "M"; rid = lipres.setdefault(at.residue.ix, 1001 + len(lipres)); seg = "MEMB"
            ln = ln.ljust(80)
            ln = "ATOM  " + ln[6:21] + ch + f"{rid % 10000:4d}" + ln[26:72] + f"{seg:<4s}" + ln[76:]
        elif ln.startswith(("CONECT",)): continue
        out.append(ln.rstrip())
    open(path, "w").write("\n".join(out) + "\n"); os.remove(tmp)

def log(row):
    df = pd.DataFrame([row]); df.to_csv(LOG, mode="a", header=not os.path.exists(LOG), index=False)

def job_single(name, sid, pep, note):
    s = get_system(sid); full, u, ref, ndx = load_dry(s)
    t, feat = pick_frame(sid, pep, need_bound=s["kind"] == "membrane")
    goto(u, t)
    if s["kind"] == "membrane":
        keep, F, copies, flip = orient(u, pep, 1)
    else:
        copies = peptide_copies(u, 1); X = u.atoms.positions
        ca = copies[0].select_atoms("name CA").positions; A = PEPTIDES[pep]["axis_all"]; Y = ca[A[0] - 1:A[1]]
        Yc = Y - Y.mean(0); vh = np.linalg.svd(Yc, full_matrices=False)[2]
        if np.dot(vh[0], Y[-1] - Y[0]) < 0: vh[0] = -vh[0]
        vh[2] = np.cross(vh[0], vh[1]); F = (X - copies[0].center_of_mass()) @ vh.T
        keep = np.zeros(len(X), bool); keep[copies[0].ix] = True; flip = False
    path = os.path.join(OUT, f"{name}.pdb"); write_pdb(u, keep, F, copies, pep, path)
    log(dict(file=f"{name}.pdb", system=sid, time_ns=t, **{k: v for k, v in feat.items()}, note=note))
    print(name, sid, f"{t:.1f} ns", feat, "flipped" if flip else "")

def job_overlay(name, items, pep, note, ref_item=0):
    """items: list of (sid, time_ns or None=representative). One chain per frame, aligned on helix C Ca."""
    ref_ca = None; blocks = {}
    chains = "ABCDEFGHIJKL"
    order = [ref_item] + [i for i in range(len(items)) if i != ref_item]   # axis frame from a helical reference
    for k in order:
        sid, t = items[k]
        s = get_system(sid); full, u, ref, ndx = load_dry(s)
        if t is None: t, _ = pick_frame(sid, pep, need_bound=s["kind"] == "membrane")
        goto(u, t); c = peptide_copies(u, 1)[0]
        lo, hi = HC[pep]; hc = c.select_atoms(f"name CA").positions[lo - 1:hi]
        X = c.positions.copy()
        if ref_ca is None:
            Y = hc - hc.mean(0); vh = np.linalg.svd(Y, full_matrices=False)[2]
            if np.dot(vh[0], Y[-1] - Y[0]) < 0: vh[0] = -vh[0]
            vh[2] = np.cross(vh[0], vh[1]); ref_ca = (hc - hc.mean(0)) @ vh.T
            Xn = (X - hc.mean(0)) @ vh.T
        else:
            Rm, _ = align.rotation_matrix(hc - hc.mean(0), ref_ca)
            Xn = (X - hc.mean(0)) @ np.asarray(Rm).T
        w = mda.Merge(c); w.atoms.positions = Xn
        tmp = os.path.join(TMPD, "_ov.pdb"); w.atoms.write(tmp)
        r0 = c.residues.resids[0]; a = 0; bl = []
        for ln in open(tmp).read().splitlines():
            if ln.startswith(("ATOM", "HETATM")):
                at = c[a]; a += 1; rid = at.residue.resid - r0 + 1; ln = ln.ljust(80)
                bl.append("ATOM  " + ln[6:21] + chains[k] + f"{rid:4d}" + ln[26:72] + f"{seg_of(pep, rid):<4s}" + ln[76:].rstrip())
        blocks[k] = bl + ["TER"]
        log(dict(file=f"{name}.pdb", system=sid, time_ns=t, note=f"{note}; chain {chains[k]}"))
    os.remove(tmp); open(os.path.join(OUT, f"{name}.pdb"), "w").write("\n".join(sum((blocks[k] for k in sorted(blocks)), [])) + "\nEND\n")
    print(name, len(items), "frames")

def job_agg(name, sid, pep, note, xw=75.0, zw=45.0):
    s = get_system(sid); full, u, ref, ndx = load_dry(s)
    cl = pd.read_csv(os.path.join(ANALYSIS_DIR, sid, "agg_clusters.csv")); mem = pd.read_csv(os.path.join(ANALYSIS_DIR, sid, "membrane_ts.csv"))
    ts = pd.read_csv(os.path.join(ANALYSIS_DIR, sid, "timeseries.csv")); hel = ts.groupby("time_ns").helicity.mean()
    nb = sum((mem[f"lip_contacts_copy{k}"].values > 0) for k in range(1, 5))
    t = cl.time_ns.values; m = (t >= W0) & (t <= W1)
    lc_mode = pd.Series(cl.largest_cluster.values[m]).mode()[0]; nb_mode = pd.Series(nb[m]).mode()[0]
    ok = m & (cl.largest_cluster.values == lc_mode) & (nb == nb_mode)
    h = hel.reindex(t, method="nearest").values; target = np.median(h[ok])
    k = np.where(ok)[0][np.argmin(np.abs(h[ok] - target))]; tt = float(t[k])
    goto(u, tt)
    copies = peptide_copies(u, 4)
    # copies in contact with each other form the cluster to centre on
    from MDAnalysis.lib.distances import distance_array
    D = np.array([[distance_array(a.positions, b.positions, box=u.dimensions).min() if i != j else 0 for j, b in enumerate(copies)] for i, a in enumerate(copies)])
    groups = []; seen = set()
    for i in range(4):
        if i in seen: continue
        st = [i]; g = set([i])
        while st:
            a = st.pop()
            for j in range(4):
                if j not in g and D[a, j] < 4.5: g.add(j); st.append(j)
        seen |= g; groups.append(sorted(g))
    big = max(groups, key=len)
    keep, F, cps, flip = orient(u, pep, 4, center_copies=big, xw=xw, zw=zw)
    path = os.path.join(OUT, f"{name}.pdb"); write_pdb(u, keep, F, cps, pep, path)
    bound = [int(mem[f"lip_contacts_copy{j}"].values[k] > 0) for j in range(1, 5)]
    log(dict(file=f"{name}.pdb", system=sid, time_ns=tt, helicity=round(float(h[k]), 3), note=f"{note}; clusters {groups} (copies 0-3 = chains A-D); bound {bound}"))
    print(name, sid, f"{tt:.1f} ns", "clusters", groups, "bound", bound, "largest", lc_mode, "n bound", nb_mode)

JOBS = {
  "A1": lambda: job_single("A1_mel_water", "mel_sol", "mel", "water, representative"),
  "A2": lambda: job_single("A2_ll37_water", "ll37_sol", "ll37", "water, representative"),
  "A3": lambda: job_overlay("A3_mel_water_overlay", [("mel_sol", float(x)) for x in range(500, 1000, 50)], "mel", "water, every 50 ns"),
  "A4": lambda: job_overlay("A4_ll37_water_overlay", [("ll37_sol", float(x)) for x in range(500, 1000, 50)], "ll37", "water, every 50 ns"),
  "B1": lambda: job_single("B1_mel_POPC", "mel_popc", "mel", "POPC reference, bound"),
  "B2": lambda: job_single("B2_mel_PCPEPG", "mel_pcpepg", "mel", "7:2:1, bound"),
  "B3": lambda: job_single("B3_mel_PEPG", "mel_pepg", "mel", "7:3, bound"),
  "B4": lambda: job_single("B4_ll37_POPC", "ll37_popc", "ll37", "POPC reference, bound frames only"),
  "B5": lambda: job_single("B5_ll37_PCPEPG", "ll37_pcpepg", "ll37", "7:2:1, bound"),
  "B6": lambda: job_single("B6_ll37_PEPG", "ll37_pepg", "ll37", "7:3, bound"),
  "C1": lambda: job_single("C1_mel_POPC_opposite_face", "mel_popc_r2", "mel", "POPC opposite face, folded N helix, deep"),
  "C2": lambda: job_single("C2_ll37_POPC_reversed", "ll37_popc_r5", "ll37", "POPC reversed, N-terminal anchor"),
  "C3": lambda: job_overlay("C3_ll37_POPC_five_runs", [("ll37_popc", None), ("ll37_popc_r2", None), ("ll37_popc_r3", None), ("ll37_popc_r4", None), ("ll37_popc_r5", None)], "ll37", "representative frame of each POPC run (A ref, B opposite, C tilt+, D tilt-, E reversed)", ref_item=1),
  "C4": lambda: job_overlay("C4_mel_POPC_five_runs", [("mel_popc", None), ("mel_popc_r2", None), ("mel_popc_r3", None), ("mel_popc_r4", None), ("mel_popc_r5", None)], "mel", "representative frame of each POPC run (A ref, B opposite, C tilt+, D tilt-, E reversed)"),
  "D1": lambda: job_agg("D1_mel_four_peptides", "mel_agg", "mel", "4 peptides, modal cluster state"),
  "D2": lambda: job_agg("D2_ll37_four_peptides", "ll37_agg", "ll37", "4 peptides, modal cluster state"),
}
if __name__ == "__main__":
    for j in sys.argv[1:]: JOBS[j]()

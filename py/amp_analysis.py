#!/usr/bin/env python3
"""Per-system analysis of the AMP simulations (melittin now, LL-37 later).

Reads the reduced trajectories written by 02_process_trajectories.sh:
    <dir>/processed/dry.xtc   peptide + lipids, every DT_ANALYSIS ps
    <dir>/processed/wet.xtc   full system, every DT_WET ps   (membrane systems)
and the ORIGINAL tpr + amp_index.ndx (the dry topology is cut from the tpr with the
Pep_Memb index group, so atom order is guaranteed to match trjconv).

Writes analysis/<sys_id>/
    timeseries.csv        one row per frame and peptide copy (conformation, insertion, orientation)
    membrane_ts.csv       one row per frame (thickness, area per lipid, lipid shell composition)
    water_ts.csv          water in the bilayer core (from wet.xtc)
    agg_pairs.csv         aggregation only: pairwise minimum distance, contact, inter-helix angle
    perres.npz            per-frame per-residue arrays (DSSP, CA depth, lipid contacts, HDX proxy terms)
    blocks.npz            block accumulators (10 equal time blocks) for RMSF-free profiles:
                          order parameters, radial bilayer height, density profiles, contact maps
    meta.json             what was run, on which frames

Usage:  python amp_analysis.py SYS_ID [--stride N] [--max-frames N] [--no-water]
"""
import os, sys, json, argparse, time, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance, distance_array
from MDAnalysis.lib.mdamath import make_whole
from MDAnalysis.analysis.rms import rmsd as _rmsd
from MDAnalysis.analysis.dssp import DSSP
from amp_config import (get_system, PEPTIDES, LIPIDS, WATER, LIPID_HEAD, LIPID_GLYC, SN1, SN2,
                        ANALYSIS_DIR, end_ns_for)

NBLK = 10                     # time blocks for accumulators
CONTACT_CUT = 4.0             # A, heavy atom contact peptide-lipid
SHELL_CUT = 6.0               # A, lipid counted in the peptide's first shell
PEP_PEP_CUT = 4.5             # A, inter-peptide heavy atom contact
BV_NC_CUT, BV_HB_CUT = 6.5, 2.4   # Best and Vendruscolo 2006 (Structure 14:97) phenomenological model
BV_BC, BV_BH = 0.35, 2.0
ZBINS = np.arange(-50.0, 50.01, 1.0)       # A, density profile bins along the normal
RBINS = np.arange(0.0, 50.01, 2.5)         # A, lateral distance from peptide


# ----------------------------------------------------------------------------- helpers
def read_ndx(path):
    groups, cur = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith("["):
            cur = line.strip("[] "); groups[cur] = []
        elif line and cur:
            groups[cur].extend(int(x) - 1 for x in line.split())
    return {k: np.array(v, dtype=int) for k, v in groups.items()}


def axis(X):
    """Principal axis of a set of CA positions, oriented N to C."""
    Xc = X - X.mean(0)
    v = np.linalg.svd(Xc, full_matrices=False)[2][0]
    return v if np.dot(v, X[-1] - X[0]) > 0 else -v


def ang(a, b):
    return np.degrees(np.arccos(np.clip(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b), -1, 1)))


def mic(d, L):
    return d - L * np.round(d / L)


def load_dry(s):
    d = s["path"]
    tpr = os.path.join(d, s["tpr"])
    ndx = read_ndx(os.path.join(d, "amp_index.ndx"))
    full = mda.Universe(tpr)
    sub = full.atoms[ndx["Pep_Memb"]]
    u = mda.Merge(sub)
    u.load_new(os.path.join(d, "processed", "dry.xtc"))
    if u.atoms.n_atoms != u.trajectory.n_atoms:
        raise SystemExit("dry.xtc atom count does not match the Pep_Memb group")
    try:
        has_xyz = full.trajectory is not None      # MDAnalysis >= 2.10 reads tpr coordinates
    except AttributeError:
        has_xyz = False
    if not has_xyz:
        # MDAnalysis <= 2.9 cannot read tpr coordinates: use frame 0 of dry.xtc (t = 0, the same
        # start structure, made whole) as the reference
        ref = mda.Merge(sub); ref.load_new(os.path.join(d, "processed", "dry.xtc")); ref.trajectory[0]
        ref.transfer_to_memory(start=0, stop=1)
        return full, u, ref, ndx
    if full.dimensions is None:          # some MDAnalysis versions return no box for a tpr alone
        full.dimensions = u.trajectory[0].dimensions
        u.trajectory[0]
    return full, u, full, ndx       # reference = tpr coordinates (start of production)


def copies_of(u, n_pep, nres):
    prot = u.select_atoms("protein")
    frs = sorted(prot.fragments, key=lambda a: a.ix[0])
    if len(frs) != n_pep:
        res = prot.residues
        frs = [res[i * nres:(i + 1) * nres].atoms for i in range(n_pep)]
    return frs


# ----------------------------------------------------------------------------- main analysis
def run(sys_id, stride=1, max_frames=None, do_water=True, end_ns=None):
    t0 = time.time()
    s = get_system(sys_id)
    pep = PEPTIDES[s["peptide"]]
    nres = len(pep["seq"])
    out = os.path.join(ANALYSIS_DIR, sys_id); os.makedirs(out, exist_ok=True)
    full, u, ref, ndx = load_dry(s)
    copies = copies_of(u, s["n_pep"], nres)
    rcopies = copies_of(ref, s["n_pep"], nres)
    for c in rcopies:
        make_whole(c)
    ncp = len(copies)
    membrane = s["kind"] != "solution"

    # analysis period: 0 to end_ns (default from amp_config / analysis_window.tsv)
    end_ns = end_ns_for(sys_id) if end_ns is None else end_ns
    t0_ps, dt_ps = u.trajectory[0].time, u.trajectory.dt
    n_keep = min(len(u.trajectory), int(np.floor((end_ns * 1000.0 - t0_ps) / dt_ps + 1e-6)) + 1)
    frames = np.arange(0, n_keep, stride)
    u.trajectory[frames[-1]]
    if u.trajectory.ts.time > end_ns * 1000.0 + 1e-3:
        raise SystemExit(f"{sys_id}: frame spacing is not uniform, cannot cut at {end_ns} ns")
    u.trajectory[0]
    if max_frames: frames = frames[:max_frames]
    nf = len(frames)
    blk_of = np.minimum((np.arange(nf) * NBLK) // max(nf, 1), NBLK - 1)

    # --- peptide selections, per copy
    seg = {n: (lo, hi) for n, lo, hi in pep["segments"]}
    hn_lo, hn_hi = pep["axis_n"]; hc_lo, hc_hi = pep["axis_c"]; ha_lo, ha_hi = pep["axis_all"]
    P = []
    for k, c in enumerate(copies):
        R = c.residues
        d = dict(all=c, heavy=c.select_atoms("not name H*"), bb=c.select_atoms("backbone"),
                 ca=c.select_atoms("name CA"))
        d["bb_n"] = R[seg[pep["segments"][0][0]][0] - 1:seg[pep["segments"][0][0]][1]].atoms.select_atoms("backbone")
        d["bb_c"] = R[pep["segments"][-1][1] - 1:pep["segments"][-1][2]].atoms.select_atoms("backbone")
        # amide N/H of residues 2..n, prolines excluded
        Ns, Hs, rid = [], [], []
        for i, r in enumerate(R):
            if i == 0 or r.resname == "PRO": continue
            n = r.atoms.select_atoms("name N"); h = r.atoms.select_atoms("name HN H")
            if len(n) and len(h):
                Ns.append(n[0].ix); Hs.append(h[0].ix); rid.append(i)
        d["N"] = u.atoms[Ns]; d["H"] = u.atoms[Hs]; d["amide_res"] = np.array(rid)
        d["res_of_atom"] = {a.ix: a.resindex - R[0].resindex for a in c}
        P.append(d)
    rP = [dict(bb=c.select_atoms("backbone"), bb_n=None) for c in rcopies]
    for k, c in enumerate(rcopies):
        R = c.residues
        rP[k]["bb_n"] = R[pep["segments"][0][1] - 1:pep["segments"][0][2]].atoms.select_atoms("backbone")
        rP[k]["bb_c"] = R[pep["segments"][-1][1] - 1:pep["segments"][-1][2]].atoms.select_atoms("backbone")
        rP[k]["bb_pos"] = rP[k]["bb"].positions.copy()
        rP[k]["bbn_pos"] = rP[k]["bb_n"].positions.copy()
        rP[k]["bbc_pos"] = rP[k]["bb_c"].positions.copy()

    allpep = u.select_atoms("protein")
    pep_heavy = allpep.select_atoms("not name H*")
    copy_of_heavy = np.concatenate([np.full(P[k]["heavy"].n_atoms, k) for k in range(ncp)])
    res_of_heavy = np.concatenate([[P[k]["res_of_atom"][a.ix] for a in P[k]["heavy"]] for k in range(ncp)])
    pep_O = allpep.select_atoms("element O or name O* and not name H*")
    copy_of_O = np.array([next(k for k in range(ncp) if a.ix in P[k]["res_of_atom"]) for a in pep_O])
    res_of_O = np.array([P[copy_of_O[j]]["res_of_atom"][a.ix] for j, a in enumerate(pep_O)])

    # --- lipid selections
    if membrane:
        lip = u.select_atoms("resname " + " ".join(LIPIDS))
        lip_heavy = lip.select_atoms("not name H*")
        lip_types = [l for l in LIPIDS if (lip.resnames == l).any()]
        ltype_of_heavy = np.array([lip_types.index(r) for r in lip_heavy.resnames])
        cls = np.full(lip_heavy.n_atoms, 2)      # 0 head+phosphate, 1 glycerol/ester, 2 acyl
        cls[np.isin(lip_heavy.names, LIPID_HEAD)] = 0
        cls[np.isin(lip_heavy.names, LIPID_GLYC)] = 1
        lipres_of_heavy = lip_heavy.resindices
        lip_O = lip.select_atoms("name O*")
        Pat = lip.select_atoms("name P")
        P_res = Pat.resindices
        nlip = Pat.n_atoms
        ltype_of_P = np.array([lip_types.index(r) for r in Pat.resnames])
        # leaflet from the first analysed frame
        u.trajectory[frames[0]]
        zc0 = Pat.positions[:, 2].mean()
        upper = Pat.positions[:, 2] > zc0
        # order parameter C-H pairs
        cpos, chain, Cix, Hix, Lid = [], [], [], [], []
        resi2lip = {r: j for j, r in enumerate(P_res)}
        for ch, names in ((0, SN1), (1, SN2)):
            for pos, nm in enumerate(names):
                for a in lip.select_atoms(f"name {nm}"):
                    for b in a.bonded_atoms:
                        if b.name.startswith("H"):
                            Cix.append(a.ix); Hix.append(b.ix); chain.append(ch); cpos.append(pos)
                            Lid.append(resi2lip[a.resindex])
        Cix, Hix, chain, cpos, Lid = map(np.array, (Cix, Hix, chain, cpos, Lid))
        ncarb = max(len(SN1), len(SN2))
        # density groups (lipids + peptide; water comes from wet.xtc)
        dens_sel = {
            "Phosphate": "name P", "Choline/amine N": "resname POPC POPE and name N",
            "PG glycerol": "resname POPG and name C12 C13 OC2 OC3", "Carbonyl O": "name O22 O32",
            "Terminal methyl": "name C218 C316", "Acyl chains": "resname POPC POPE POPG and name C2* C3* and not name C2 C3 C21 C31",
            "Peptide": "protein and not name H*", "Lys NZ": "protein and resname LYS and name NZ",
            "Arg CZ": "protein and resname ARG and name CZ",
            "Aromatic ring": "protein and resname TRP PHE TYR and name CG CD1 CD2 CE1 CE2 CZ NE1 CE3 CZ2 CZ3 CH2",
        }
        dens_ag = {k: u.select_atoms(v) for k, v in dens_sel.items()}
        dens_ag = {k: v for k, v in dens_ag.items() if v.n_atoms}

    # --- accumulators
    ts_rows, mem_rows, pair_rows = [], [], []
    depth = np.full((nf, ncp, nres), np.nan, np.float32)           # CA height above proximal P plane (A)
    contact = np.zeros((nf, ncp, nres, 4), np.bool_)               # any lipid, head, glycerol, acyl
    contact_type = np.zeros((nf, ncp, nres, 3), np.bool_)          # by lipid type (index in lip_types)
    nc_pep = np.zeros((nf, ncp, nres), np.float32); nh_pep = np.zeros_like(nc_pep)
    nc_env = np.zeros_like(nc_pep); nh_env = np.zeros_like(nc_pep)
    bb_traj = [np.zeros((nf, P[k]["bb"].n_atoms, 3), np.float32) for k in range(ncp)]
    B = {}
    if membrane:
        B["scd_sum"] = np.zeros((NBLK, 4, 2, ncarb)); B["scd_n"] = np.zeros_like(B["scd_sum"])
        B["rad_sum"] = np.zeros((NBLK, 2, len(RBINS) - 1)); B["rad_sq"] = np.zeros_like(B["rad_sum"])
        B["rad_n"] = np.zeros_like(B["rad_sum"])
        B["dens"] = np.zeros((NBLK, len(dens_ag), len(ZBINS) - 1)); B["dens_vol"] = np.zeros(NBLK)
        B["dens_names"] = np.array(list(dens_ag))
    if ncp > 1:
        B["rr_map"] = np.zeros((NBLK, nres, nres)); B["rr_n"] = np.zeros(NBLK)

    print(f"{sys_id}: {nf} frames ({len(u.trajectory)} in dry.xtc, stride {stride}), {ncp} peptide(s), "
          f"{'membrane' if membrane else 'solution'}")
    for fi, frame in enumerate(frames):
        ts = u.trajectory[frame]
        box = ts.dimensions; L = box[:3]
        tns = ts.time / 1000.0
        b = blk_of[fi]
        for k in range(ncp):
            make_whole(P[k]["all"])
        if membrane:
            zP = Pat.positions[:, 2]
            zc = zP.mean()
            zu, zl = zP[upper].mean(), zP[~upper].mean()
        # ---------- per-copy peptide metrics
        dzs, sides = [], []
        for k in range(ncp):
            d = P[k]; ca = d["ca"].positions
            bb_traj[k][fi] = d["bb"].positions
            row = dict(time_ns=tns, copy=k + 1, rg=d["all"].radius_of_gyration() / 10,
                       e2e=np.linalg.norm(ca[-1] - ca[0]) / 10)
            row["rmsd_bb"] = _rmsd(d["bb"].positions, rP[k]["bb_pos"], center=True, superposition=True) / 10
            row["rmsd_hn"] = _rmsd(d["bb_n"].positions, rP[k]["bbn_pos"], center=True, superposition=True) / 10
            row["rmsd_hc"] = _rmsd(d["bb_c"].positions, rP[k]["bbc_pos"], center=True, superposition=True) / 10
            an = axis(ca[hn_lo - 1:hn_hi]); ac = axis(ca[hc_lo - 1:hc_hi]); aa = axis(ca[ha_lo - 1:ha_hi])
            row["kink"] = ang(an, ac)
            d["_ax"] = aa
            if membrane:
                com = d["all"].center_of_mass()
                dz = mic(com[2] - zc, L[2]); side = 1.0 if dz >= 0 else -1.0
                hprox = (zu - zc) if side > 0 else (zc - zl)
                zres = dz + (ca[:, 2] - com[2])
                depth[fi, k] = side * zres - hprox
                nz = np.array([0, 0, side])
                row.update(dz_center=dz / 10, side=int(side), height_com=(side * dz - hprox) / 10,
                           tilt=ang(aa, nz), tilt_n=ang(an, nz), tilt_c=ang(ac, nz),
                           depth_min=np.nanmin(depth[fi, k]) / 10)
                dzs.append(dz); sides.append(side)
            ts_rows.append(row)
        # ---------- peptide-lipid contacts (all copies at once)
        if membrane:
            pr = capped_distance(pep_heavy.positions, lip_heavy.positions, CONTACT_CUT, box=box,
                                 return_distances=False)
            if len(pr):
                pk, pres = copy_of_heavy[pr[:, 0]], res_of_heavy[pr[:, 0]]
                contact[fi, pk, pres, 0] = True
                for c_ in range(3):
                    m = cls[pr[:, 1]] == c_
                    contact[fi, pk[m], pres[m], c_ + 1] = True
                lt = ltype_of_heavy[pr[:, 1]]
                contact_type[fi, pk, pres, lt] = True
            # first-shell lipid composition (all copies pooled)
            sh = capped_distance(pep_heavy.positions, lip_heavy.positions, SHELL_CUT, box=box,
                                 return_distances=False)
            shell = np.unique(lipres_of_heavy[sh[:, 1]]) if len(sh) else np.array([], int)
            mrow = dict(time_ns=tns, zc=zc / 10, thickness_PP=(zu - zl) / 10,
                        apl=L[0] * L[1] / (nlip / 2) / 100, box_x=L[0] / 10, box_z=L[2] / 10,
                        n_shell=len(shell))
            st = np.isin(P_res, shell)
            for j, l in enumerate(lip_types):
                mrow[f"shell_{l}"] = int((st & (ltype_of_P == j)).sum())
            for k in range(ncp):
                mrow[f"lip_contacts_copy{k + 1}"] = int(contact[fi, k, :, 0].sum())
            mem_rows.append(mrow)
            # proximal leaflet = side of the peptides (majority)
            sref = np.sign(np.sum(sides)) or sides[0]
            prox = upper if sref > 0 else ~upper
            # lateral distance of every P to the nearest peptide CA
            caall = np.concatenate([P[k]["ca"].positions for k in range(ncp)])
            dxy = mic(Pat.positions[:, None, :2] - caall[None, :, :2], L[:2])
            rmin = np.sqrt((dxy ** 2).sum(-1)).min(1)
            # radial bilayer height (P height from centre), proximal and distal
            hP = np.abs(mic(zP - zc, L[2]))
            ib = np.digitize(rmin, RBINS) - 1
            for leaf, mask in ((0, prox), (1, ~prox)):
                ok = mask & (ib >= 0) & (ib < len(RBINS) - 1)
                np.add.at(B["rad_sum"][b, leaf], ib[ok], hP[ok])
                np.add.at(B["rad_sq"][b, leaf], ib[ok], hP[ok] ** 2)
                np.add.at(B["rad_n"][b, leaf], ib[ok], 1)
            # order parameters by lipid group: 0 prox <10 A, 1 prox 10-20 A, 2 prox >20 A, 3 distal
            grp = np.where(~prox, 3, np.where(rmin < 10, 0, np.where(rmin < 20, 1, 2)))
            v = u.atoms.positions[Hix] - u.atoms.positions[Cix]
            p2 = 1.5 * (v[:, 2] / np.linalg.norm(v, axis=1)) ** 2 - 0.5
            g = grp[Lid]
            np.add.at(B["scd_sum"][b], (g, chain, cpos), p2)
            np.add.at(B["scd_n"][b], (g, chain, cpos), 1)
            # density profiles, peptide side positive
            for j, (nm, ag) in enumerate(dens_ag.items()):
                zz = sref * mic(ag.positions[:, 2] - zc, L[2])
                B["dens"][b, j] += np.histogram(zz, ZBINS)[0]
            B["dens_vol"][b] += L[0] * L[1] * (ZBINS[1] - ZBINS[0]) / 1000.0     # nm^3 per bin
        # ---------- HDX protection proxy (Best and Vendruscolo), peptide-only and environment-inclusive
        env_heavy = pep_heavy.positions if not membrane else np.concatenate([pep_heavy.positions, lip_heavy.positions])
        nph = pep_heavy.n_atoms
        env_O = pep_O.positions if not membrane else np.concatenate([pep_O.positions, lip_O.positions])
        nO = pep_O.n_atoms
        for k in range(ncp):
            d = P[k]; ar = d["amide_res"]
            pr = capped_distance(d["N"].positions, env_heavy, BV_NC_CUT, box=box, return_distances=False)
            if len(pr):
                i, j = pr[:, 0], pr[:, 1]
                is_pep = j < nph
                jj = np.where(is_pep, j, 0)
                same = is_pep & (copy_of_heavy[jj] == k)
                near_seq = same & (np.abs(res_of_heavy[jj] - ar[i]) <= 2)
                keep = ~near_seq
                np.add.at(nc_env[fi, k], ar[i[keep]], 1)
                np.add.at(nc_pep[fi, k], ar[i[keep & same]], 1)
            hb = capped_distance(d["H"].positions, env_O, BV_HB_CUT, box=box, return_distances=False)
            if len(hb):
                i, j = hb[:, 0], hb[:, 1]
                is_pep = j < nO
                jj = np.where(is_pep, j, 0)
                same = is_pep & (copy_of_O[jj] == k)
                not_self = ~(same & (res_of_O[jj] == ar[i]))
                np.add.at(nh_env[fi, k], ar[i[not_self]], 1)
                np.add.at(nh_pep[fi, k], ar[i[not_self & same]], 1)
        # ---------- aggregation: pairs, clusters, residue contact map
        if ncp > 1:
            pr = capped_distance(pep_heavy.positions, pep_heavy.positions, PEP_PEP_CUT, box=box,
                                 return_distances=False)
            ck = copy_of_heavy[pr[:, 0]], copy_of_heavy[pr[:, 1]]
            inter = ck[0] != ck[1]
            rr = np.zeros((nres, nres), bool)
            rr[res_of_heavy[pr[inter, 0]], res_of_heavy[pr[inter, 1]]] = True
            rr = rr | rr.T                     # residue i of one copy touches residue j of another, either way round
            B["rr_map"][b] += rr; B["rr_n"][b] += 1
            adj = np.zeros((ncp, ncp), bool); adj[ck[0][inter], ck[1][inter]] = True
            for a_ in range(ncp):
                for c_ in range(a_ + 1, ncp):
                    dmin = distance_array(P[a_]["heavy"].positions, P[c_]["heavy"].positions, box=box).min() / 10
                    pair_rows.append(dict(time_ns=tns, pair=f"{a_ + 1}-{c_ + 1}", dmin=dmin,
                                          contact=int(adj[a_, c_]), axis_angle=ang(P[a_]["_ax"], P[c_]["_ax"])))
            # connected components
            lab = list(range(ncp))
            def find(x):
                while lab[x] != x: x = lab[x]
                return x
            for a_ in range(ncp):
                for c_ in range(ncp):
                    if adj[a_, c_]: lab[find(a_)] = find(c_)
            roots = [find(x) for x in range(ncp)]
            sizes = np.bincount(roots)
            pair_rows[-1]["_largest"] = int(sizes.max()); pair_rows[-1]["_nclusters"] = int((sizes > 0).sum())
        if fi % max(1, nf // 20) == 0:
            print(f"  frame {fi + 1}/{nf}  t = {tns:8.1f} ns  ({time.time() - t0:.0f} s)")

    # ---------- DSSP, per copy
    dssp = np.zeros((nf, ncp, nres), "U1")
    for k in range(ncp):
        # DSSP run per frame index list
        D = DSSP(P[k]["all"]).run(frames=frames)
        dssp[:, k] = D.results.dssp
    helix = dssp == "H"

    # ---------- RMSF (CA), two-pass alignment to the average, per copy, all frames
    from MDAnalysis.analysis.align import rotation_matrix
    rmsf = np.zeros((ncp, nres))
    for k in range(ncp):
        X = bb_traj[k].astype(np.float64)
        names = P[k]["bb"].names
        ca_mask = names == "CA"
        ref_ = rP[k]["bb_pos"] - rP[k]["bb_pos"].mean(0)
        for _ in range(2):
            Y = np.empty_like(X)
            for f in range(nf):
                xc = X[f] - X[f].mean(0)
                R, _r = rotation_matrix(xc, ref_)
                Y[f] = xc @ R.T
            ref_ = Y.mean(0)
        rmsf[k] = np.sqrt(((Y[:, ca_mask] - ref_[ca_mask]) ** 2).sum(-1).mean(0)) / 10

    # ---------- write
    ts_df = pd.DataFrame(ts_rows)
    ts_df["helicity"] = helix.reshape(nf * ncp, nres).mean(1) if ncp == 1 else \
        np.concatenate([[helix[f, k].mean() for k in range(ncp)] for f in range(nf)])
    ts_df.to_csv(os.path.join(out, "timeseries.csv"), index=False, float_format="%.4f")
    if membrane:
        pd.DataFrame(mem_rows).to_csv(os.path.join(out, "membrane_ts.csv"), index=False, float_format="%.4f")
    if ncp > 1:
        pdf = pd.DataFrame(pair_rows)
        cl = pdf.dropna(subset=["_largest"])[["time_ns", "_largest", "_nclusters"]]
        cl.columns = ["time_ns", "largest_cluster", "n_clusters"]
        cl.to_csv(os.path.join(out, "agg_clusters.csv"), index=False)
        pdf.drop(columns=["_largest", "_nclusters"]).to_csv(os.path.join(out, "agg_pairs.csv"), index=False,
                                                           float_format="%.4f")
    times = ts_df[ts_df["copy"] == 1]["time_ns"].values
    np.savez_compressed(os.path.join(out, "perres.npz"), time_ns=times, dssp=dssp, depth=depth,
                        contact=contact, contact_type=contact_type,
                        lip_types=np.array(lip_types if membrane else []),
                        nc_pep=nc_pep, nh_pep=nh_pep, nc_env=nc_env, nh_env=nh_env, rmsf=rmsf,
                        seq=np.array(list(pep["seq"])))
    B["block_t0_ns"] = np.array([times[blk_of == i].min() if (blk_of == i).any() else np.nan for i in range(NBLK)])
    B["block_t1_ns"] = np.array([times[blk_of == i].max() if (blk_of == i).any() else np.nan for i in range(NBLK)])
    B["zbins"] = ZBINS; B["rbins"] = RBINS
    np.savez_compressed(os.path.join(out, "blocks.npz"), **B)
    meta = dict(sys_id=sys_id, system=s, frames=int(nf), stride=stride, dry_frames=len(u.trajectory),
                t_first_ns=float(times[0]), t_last_ns=float(times[-1]),
                cutoffs_A=dict(contact=CONTACT_CUT, shell=SHELL_CUT, pep_pep=PEP_PEP_CUT,
                               bv_nc=BV_NC_CUT, bv_hb=BV_HB_CUT), bv_coeffs=[BV_BC, BV_BH],
                lip_types=lip_types if membrane else [], runtime_s=round(time.time() - t0, 1),
                end_ns_applied=float(end_ns))
    json.dump(meta, open(os.path.join(out, "meta.json"), "w"), indent=1, default=str)

    if membrane and do_water:
        water_pass(s, full, ndx, out, sides_ref=None, end_ns=end_ns)
    print(f"{sys_id}: done in {time.time() - t0:.0f} s -> {out}")


def water_pass(s, full, ndx, out, sides_ref=None, end_ns=1e12):
    """Water in the hydrophobic core and water density profile, from wet.xtc (full system)."""
    f = os.path.join(s["path"], "processed", "wet.xtc")
    if not os.path.exists(f):
        print("  wet.xtc not found, water analysis skipped"); return
    u = mda.Universe(os.path.join(s["path"], s["tpr"]), f)
    Pat = u.select_atoms("name P and resname " + " ".join(LIPIDS))
    OW = u.atoms[ndx["Water"]].select_atoms("name OW OH2")
    pepca = u.atoms[ndx["Peptide_CA"]]
    pepall = u.atoms[ndx["Peptide"]]
    rows, prof, vol = [], np.zeros(len(ZBINS) - 1), 0.0
    for ts in u.trajectory:
        if ts.time > end_ns * 1000.0 + 1e-3: break
        L = ts.dimensions[:3]
        zc = Pat.positions[:, 2].mean()
        dzp = mic(pepall.center_of_mass()[2] - zc, L[2]); sref = 1.0 if dzp >= 0 else -1.0
        zw = sref * mic(OW.positions[:, 2] - zc, L[2])
        core = np.abs(zw) < 5.0
        dxy = mic(OW.positions[core, None, :2] - pepca.positions[None, :, :2], L[:2])
        near = (np.sqrt((dxy ** 2).sum(-1)).min(1) < 12.0).sum() if core.any() else 0
        rows.append(dict(time_ns=ts.time / 1000, water_core_10A=int(core.sum()),
                         water_core_near_peptide=int(near), water_within_15A=int((np.abs(zw) < 15).sum())))
        prof += np.histogram(zw, ZBINS)[0]; vol += L[0] * L[1] * (ZBINS[1] - ZBINS[0]) / 1000.0
    pd.DataFrame(rows).to_csv(os.path.join(out, "water_ts.csv"), index=False)
    np.savez_compressed(os.path.join(out, "water_profile.npz"), zbins=ZBINS, density=prof / max(vol, 1e-9))
    print(f"  water: {len(rows)} frames")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("sys_id")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--max-frames", type=int, default=None)
    ap.add_argument("--no-water", action="store_true")
    ap.add_argument("--end-ns", type=float, default=None, help="default: END_NS or analysis_window.tsv")
    a = ap.parse_args()
    run(a.sys_id, a.stride, a.max_frames, not a.no_water, a.end_ns)

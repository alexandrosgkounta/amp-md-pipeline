#!/usr/bin/env python3
"""Is 0 to END_NS (1 us) enough for every run, or does one need a longer analysis period?

Usage:  python check_extension.py mel            (or ll37)

Every system is analysed from 0 to END_NS and averaged over the last half (500 to 1000 ns).
A run is analysed further ONLY if one of three criteria fails, and the report states which:

  C1  bound in time    membrane runs: sustained lipid contact (10 ns) is reached before the
                       averaging window starts (END_NS / 2). Otherwise the window mixes unbound
                       and bound frames.
  C2  stationary       inside the window, the two halves agree: |mean2 - mean1| is smaller than
                       the practical threshold OR within 2 combined block SE. A larger shift is drift.
  C3  nothing new      after END_NS (from processed/full_check.xtc, peptide + lipid P every 1 ns,
                       whole run): the mean over the part after END_NS stays within the practical
                       threshold of the window mean, computed the same way from the same file.
                       Needs at least MIN_AFTER_NS of data after END_NS.

Practical thresholds: helicity 0.10, COM height 0.3 nm, helix tilt 15 deg, Rg 0.15 nm,
largest cluster 0.5 peptides. They are set at roughly the size of a change that would alter a
conclusion in the chapter, not at statistical detectability (block SE from one run is too
small to be a useful test on its own).

Verdict per run:
  ok              all criteria pass, analyse 0 to END_NS
  extend_analysis a criterion fails and the trajectory continues beyond END_NS: analyse to the
                  end of the run (add the line given in the report to analysis_window.tsv)
  extend_run      a criterion fails and the trajectory stops at END_NS: the period is not long
                  enough; extend the simulation or state it as a limitation

Writes analysis/extension_check_<pep>.csv and analysis/extension_check_<pep>.md
"""
import os, sys, json, argparse, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.mdamath import make_whole
from MDAnalysis.lib.distances import capped_distance
from amp_config import read_systems, PEPTIDES, ANALYSIS_DIR, end_ns_for

THR = dict(helicity=0.10, height_com=0.3, tilt=15.0, rg=0.15, largest_cluster=0.5)
UNIT = dict(helicity="", height_com=" nm", tilt=" deg", rg=" nm", largest_cluster=" peptides")
MIN_AFTER_NS = float(os.environ.get("MIN_AFTER_NS", 100))
BIND_SUSTAIN_NS = 10.0


def read_ndx(path):
    g, cur = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith("["):
            cur = line.strip("[] "); g[cur] = []
        elif line and cur:
            g[cur].extend(int(x) - 1 for x in line.split())
    return {k: np.array(v, int) for k, v in g.items()}


def axis(X):
    Xc = X - X.mean(0); v = np.linalg.svd(Xc, full_matrices=False)[2][0]
    return v if np.dot(v, X[-1] - X[0]) > 0 else -v


def ang(a, b):
    return np.degrees(np.arccos(np.clip(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b), -1, 1)))


def mic(d, L):
    return d - L * np.round(d / L)


def block_se(x, nb=5):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 2 * nb: return np.nan
    m = [b.mean() for b in np.array_split(x, nb)]
    return np.std(m, ddof=1) / np.sqrt(nb)


def t_bind(t, ncontacts):
    bound = np.asarray(ncontacts) > 0
    dt = np.median(np.diff(t)) if len(t) > 1 else 1.0
    nw = max(1, int(round(BIND_SUSTAIN_NS / dt)))
    run = np.convolve(bound.astype(int), np.ones(nw, int), "valid") == nw
    return float(t[np.argmax(run)]) if run.any() else np.nan


def full_check_series(s, pep, membrane):
    """Per-frame metrics over the whole run, every ~1 ns.

    Source: processed/full_check.xtc (peptide + lipid P, written by the current step 02 / cluster
    script). If it is missing, the full-length processed/dry.xtc (peptide + lipids, from an earlier
    processing of the whole run) is read instead, every 1 ns. No reprocessing is needed in that case.
    """
    d = s["path"]
    nf_ = os.path.join(d, "amp_index.ndx")
    ndx = read_ndx(nf_ if os.path.exists(nf_) else os.path.join(d, "index.ndx"))
    full = mda.Universe(os.path.join(d, s["tpr"]))
    f = os.path.join(d, "processed", "full_check.xtc")
    if os.path.exists(f) and "Pep_LipidP" in ndx:
        u = mda.Merge(full.atoms[ndx["Pep_LipidP"]]); src = "full_check.xtc"
    else:
        f = os.path.join(d, "processed", "dry.xtc")
        if not os.path.exists(f): return None
        u = mda.Merge(full.atoms[ndx["Pep_Memb"]]); src = "dry.xtc"
    u.load_new(f)
    dt = u.trajectory.dt
    frames = np.arange(0, len(u.trajectory), max(1, int(round(1000.0 / dt))))
    print(f"  {s['sys_id']}: full-length check from {src}, {len(frames)} frames to "
          f"{(u.trajectory[0].time + (len(u.trajectory) - 1) * dt) / 1000:.0f} ns")
    nres = len(pep["seq"]); prot = u.select_atoms("protein")
    frs = sorted(prot.fragments, key=lambda a: a.ix[0])
    if len(frs) != s["n_pep"]:
        frs = [prot.residues[i * nres:(i + 1) * nres].atoms for i in range(s["n_pep"])]
    lo, hi = pep["axis_all"]
    P = u.select_atoms("name P and not protein")
    heavy = [c.select_atoms("not name H*") for c in frs]
    from MDAnalysis.analysis.dssp import DSSP
    hel = [(DSSP(c).run(frames=frames).results.dssp == "H").mean(1) for c in frs]
    rows = []
    for fi, fr in enumerate(frames):
        ts = u.trajectory[fr]
        L = ts.dimensions[:3]; tns = ts.time / 1000.0
        for c in frs: make_whole(c)
        if membrane and P.n_atoms:
            zP = P.positions[:, 2]; zc = zP.mean(); dzP = mic(zP - zc, L[2])
            zu, zl = dzP[dzP > 0].mean(), dzP[dzP < 0].mean()
        for k, c in enumerate(frs):
            ca = c.select_atoms("name CA").positions
            r = dict(time_ns=tns, copy=k + 1, helicity=hel[k][fi], rg=c.radius_of_gyration() / 10)
            if membrane and P.n_atoms:
                dz = mic(c.center_of_mass()[2] - zc, L[2]); side = 1.0 if dz >= 0 else -1.0
                hprox = zu if side > 0 else -zl
                r["height_com"] = (side * dz - hprox) / 10
                r["tilt"] = ang(axis(ca[lo - 1:hi]), np.array([0, 0, side]))
            rows.append(r)
        if len(frs) > 1:   # largest cluster (heavy atoms within 4.5 A)
            par = list(range(len(frs)))
            def find(x):
                while par[x] != x: x = par[x]
                return x
            for a in range(len(frs)):
                for b in range(a + 1, len(frs)):
                    if len(capped_distance(heavy[a].positions, heavy[b].positions, 4.5, box=ts.dimensions,
                                           return_distances=False)):
                        par[find(a)] = find(b)
            roots = [find(x) for x in range(len(frs))]
            big = max(roots.count(x) for x in set(roots))
            for r in rows[-len(frs):]: r["largest_cluster"] = big
    return pd.DataFrame(rows)


def fmt(x, m):
    return "n/a" if not np.isfinite(x) else (f"{x:.2f}" if m in ("helicity", "rg", "height_com", "largest_cluster") else f"{x:.0f}")


def check(s, pep):
    sid = s["sys_id"]; out = os.path.join(ANALYSIS_DIR, sid)
    p = os.path.join(out, "timeseries.csv")
    if not os.path.exists(p): return []
    meta = json.load(open(os.path.join(out, "meta.json")))
    end = float(meta.get("end_ns_applied", end_ns_for(sid)))
    w0 = end / 2.0
    membrane = s["kind"] != "solution"
    ts = pd.read_csv(p)
    mem = pd.read_csv(os.path.join(out, "membrane_ts.csv")) if membrane and os.path.exists(os.path.join(out, "membrane_ts.csv")) else None
    metrics = ["helicity", "height_com", "tilt"] if membrane else ["helicity", "rg"]
    fc = full_check_series(s, pep, membrane)
    if fc is not None and s["n_pep"] > 1: metrics_fc = metrics + ["largest_cluster"]
    else: metrics_fc = metrics
    t_last = float(fc["time_ns"].max()) if fc is not None else float(ts["time_ns"].max())
    rows = []
    for k, g in ts.groupby("copy"):
        t = g["time_ns"].values
        r = dict(sys_id=sid, condition=s["condition"], replicate=s["replicate"], copy=int(k), analysed_to_ns=end,
                 window=f"{w0:.0f}-{end:.0f}", run_length_ns=round(t_last, 1), reasons=[])
        # C1
        if membrane and mem is not None and f"lip_contacts_copy{k}" in mem:
            tb = t_bind(mem["time_ns"].values, mem[f"lip_contacts_copy{k}"].values)
            r["t_bind_ns"] = tb; r["C1_bound_in_time"] = bool(np.isfinite(tb) and tb <= w0)
            if not r["C1_bound_in_time"]:
                r["reasons"].append("no sustained lipid contact before the window starts" if not np.isfinite(tb)
                                    else f"sustained lipid contact only from {tb:.0f} ns, after the window start ({w0:.0f} ns)")
        else:
            r["C1_bound_in_time"] = True
        # C2
        w = g[(t >= w0) & (t <= end)]; h = len(w) // 2
        c2 = True
        for m in metrics:
            if m not in w: continue
            a_, b_ = w[m].values[:h], w[m].values[h:]
            d = b_.mean() - a_.mean(); se = np.sqrt(block_se(a_) ** 2 + block_se(b_) ** 2)
            r[f"C2_{m}_shift"] = d
            drift = abs(d) > THR[m] and not (np.isfinite(se) and abs(d) <= 2 * se)
            if drift:
                c2 = False
                r["reasons"].append(f"{m} drifts inside the window ({fmt(a_.mean(), m)} to {fmt(b_.mean(), m)}{UNIT[m]})")
        r["C2_stationary"] = c2
        # C3
        r["C3_nothing_new"] = None
        if fc is not None:
            gc = fc[fc["copy"] == k]; tc = gc["time_ns"].values
            after = gc[tc > end + 1e-6]; win = gc[(tc >= w0) & (tc <= end)]
            r["ns_after_end"] = float(after["time_ns"].max() - end) if len(after) else 0.0
            if len(after) and r["ns_after_end"] >= MIN_AFTER_NS:
                ok = True
                for m in metrics_fc:
                    d = after[m].mean() - win[m].mean(); r[f"C3_{m}_shift"] = d
                    if abs(d) > THR[m]:
                        ok = False
                        # first time the 50 ns running mean leaves the window mean by the threshold
                        rm = pd.Series(gc[m].values).rolling(max(1, int(50 / max(np.median(np.diff(tc)), 1e-9))),
                                                            center=True, min_periods=1).mean().values
                        idx = np.where((tc > end) & (np.abs(rm - win[m].mean()) > THR[m]))[0]
                        tch = f", from about {tc[idx[0]]:.0f} ns" if len(idx) else ""
                        r["reasons"].append(f"{m} changes after {end:.0f} ns ({fmt(win[m].mean(), m)} in the window, "
                                            f"{fmt(after[m].mean(), m)}{UNIT[m]} after{tch})")
                r["C3_nothing_new"] = ok
        failed = not (r["C1_bound_in_time"] and r["C2_stationary"] and r["C3_nothing_new"] is not False)
        more_data = r.get("ns_after_end", 0) >= MIN_AFTER_NS
        r["verdict"] = "ok" if not failed else ("extend_analysis" if more_data else "extend_run")
        rows.append(r)
    return rows


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("peptide"); a = ap.parse_args()
    pep = PEPTIDES[a.peptide]; rows = []
    for s in read_systems():
        if s["peptide"] != a.peptide: continue
        try:
            rows += check(s, pep)
        except Exception as e:
            print(f"  {s['sys_id']}: check failed ({e})")
    if not rows:
        print("no analysed systems"); return
    df = pd.DataFrame(rows)
    # one verdict per system (worst copy)
    order = {"ok": 0, "extend_analysis": 1, "extend_run": 2}
    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    out = df.copy(); out["reasons"] = out["reasons"].apply("; ".join)
    out.to_csv(os.path.join(ANALYSIS_DIR, f"extension_check_{a.peptide}.csv"), index=False, float_format="%.4g")
    L = [f"# Analysis period check, {pep['name']}", "",
         f"All runs analysed from 0 to their analysis end (default {df['analysed_to_ns'].min():.0f} ns), averages over the last half. "
         f"Criteria and thresholds are in the header of check_extension.py.", "",
         "| system | copy | run length (ns) | t bind (ns) | C1 bound in time | C2 stationary | C3 nothing new after end | verdict |",
         "|---|---|---|---|---|---|---|---|"]
    yn = lambda v: "n/a" if v is None or (isinstance(v, float) and not np.isfinite(v)) else ("yes" if v else "**no**")
    for _, r in df.iterrows():
        tb = r.get("t_bind_ns", np.nan)
        L.append(f"| {r['sys_id']} | {r['copy']} | {r['run_length_ns']:.0f} | {'' if not np.isfinite(tb) else f'{tb:.0f}'} | "
                 f"{yn(r.C1_bound_in_time)} | {yn(r.C2_stationary)} | {yn(r.C3_nothing_new)} | {r.verdict} |")
    L += ["", "## Runs that need more than the default period", ""]
    any_ = False
    for sid, g in df.groupby("sys_id", sort=False):
        v = max(g["verdict"], key=order.get)
        if v == "ok": continue
        any_ = True
        why = "; ".join(sorted({x for rs in g["reasons"] for x in rs}))
        end_full = g["run_length_ns"].max()
        if v == "extend_analysis":
            L.append(f"**{sid}: analyse to {end_full:.0f} ns.** Justification: {why}. "
                     f"Add to pipeline/analysis_window.tsv:  `{sid}\t{end_full:.0f}\t{why}`")
        else:
            L.append(f"**{sid}: the run stops at {end_full:.0f} ns and is not settled.** {why}. "
                     f"Extend the simulation, or report this run's averages as not converged.")
        L.append("")
    if not any_:
        L.append("None. Every run passes C1 to C3, so 0 to 1 us is used for all of them.")
    L += ["", "Notes. C3 needs processed/full_check.xtc and at least "
          f"{MIN_AFTER_NS:.0f} ns after the analysis end; n/a means the run is not longer or the file is missing. "
          "An extension changes that run's averaging window (last half of the longer period), so state it in Methods "
          "and keep the 1 us values of the same run in the supplement for comparison."]
    open(os.path.join(ANALYSIS_DIR, f"extension_check_{a.peptide}.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L[:6 + len(df)]))
    print(f"-> analysis/extension_check_{a.peptide}.md")


if __name__ == "__main__":
    main()

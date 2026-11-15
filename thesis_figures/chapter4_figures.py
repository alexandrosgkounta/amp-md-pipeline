"""Chapter 4 figures: melittin and LL-37 side by side, in four sections.
  Fig 4.A  solution           Fig 4.B  membranes (reference orientation on each bilayer)
  Fig 4.C  POPC starting orientations (five runs per peptide)      Fig 4.D  four-peptide systems
All averages over 500 to 1000 ns of the 0 to 1 us analysis period. Writes figures/chapter4/*.svg|pdf|png
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import amp_style as S
from amp_config import PEPTIDES, ANALYSIS_DIR, FIG_DIR, CONDITIONS
S.apply()
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = os.path.join(FIG_DIR, "chapter4"); os.makedirs(OUT, exist_ok=True)
W0, W1 = 500.0, 1000.0
PEPS = ["mel", "ll37"]
NAME = {"mel": "Melittin", "ll37": "LL-37"}
COND = {"sol": ("Water", "#4d4d4d", ":"), "popc": ("POPC", "#2a78d6", "-"), "pcpepg": ("POPC:POPE:POPG 7:2:1", "#1baf7a", "--"),
        "pepg": ("POPE:POPG 7:3", "#eb6834", "-")}
RUNCOL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
RUNLAB = {"mel": ["reference", "opposite face", "tilted +", "tilted −", "reversed end to end"],
          "ll37": ["reference", "opposite face", "tilted +", "tilted −", "reversed end to end"]}
TEMP = {"mel": "343 K", "ll37": "313 K (water and POPE:POPG 343 K)"}

def load(sid):
    d = os.path.join(ANALYSIS_DIR, sid)
    D = dict(ts=pd.read_csv(os.path.join(d, "timeseries.csv")), pr=dict(np.load(os.path.join(d, "perres.npz"), allow_pickle=True)))
    for n, f in (("mem", "membrane_ts.csv"), ("pairs", "agg_pairs.csv"), ("cl", "agg_clusters.csv")):
        p = os.path.join(d, f); D[n] = pd.read_csv(p) if os.path.exists(p) else None
    p = os.path.join(d, "blocks.npz"); D["blk"] = dict(np.load(p, allow_pickle=True)) if os.path.exists(p) else {}
    return D
def win(t): return (t >= W0) & (t <= W1)
def helix_profile(D, copies=None):
    H = D["pr"]["dssp"] == "H"; m = win(D["pr"]["time_ns"])
    H = H[m] if copies is None else H[m][:, copies]
    return H.mean(axis=(0, 1))
def bound_mask(D, k=1):
    return D["mem"][f"lip_contacts_copy{k}"].values > 0
def rolling(t, x, ns=5.0):
    dt = np.median(np.diff(t)); w = max(1, int(round(ns / dt)))
    return pd.Series(x).rolling(w, center=True, min_periods=1).mean().values
def res_axis(ax, pep, labels=True):
    P = PEPTIDES[pep]; S.residue_axis(ax, P["seq"], P["segments"], label_segments=labels)
def pro_nan(y, pep):
    y = y.astype(float).copy()
    return y

# ======================================================================= Fig 4.A solution
def fig_solution():
    fig = S.figure(height_mm=118)
    gs = fig.add_gridspec(2, 3, hspace=0.75, wspace=0.42, width_ratios=[1.35, 1, 1])
    for i, pep in enumerate(PEPS):
        D = load(f"{pep}_sol"); n = len(PEPTIDES[pep]["seq"])
        a = fig.add_subplot(gs[i, 0]); res_axis(a, pep)
        a.plot(np.arange(1, n + 1), helix_profile(D), "o-", ms=2.2, color="#4d4d4d", lw=1)
        a.set_ylim(0, 1.05); a.set_ylabel("Helix fraction"); a.set_title(NAME[pep], loc="left", fontsize=S.PT_NAME, pad=12)
        b = fig.add_subplot(gs[i, 1]); ts = D["ts"]
        b.plot(ts.time_ns / 1000, rolling(ts.time_ns.values, ts.helicity.values * 100, 5), color="#4d4d4d", lw=0.8)
        b.axvspan(W0 / 1000, W1 / 1000, color="#f2f2f2", zorder=0, lw=0)
        b.set_ylim(0, 100); b.set_xlabel("Time (µs)"); b.set_ylabel("Helicity (%)")
        c = fig.add_subplot(gs[i, 2]); w = ts[win(ts.time_ns)]
        h, x, y = np.histogram2d(w.rg, w.kink, bins=[np.linspace(0.8, 1.8, 26), np.linspace(0, 180, 31)], density=True)
        c.imshow(h.T, origin="lower", aspect="auto", extent=(x[0], x[-1], y[0], y[-1]), cmap=S.SEQ)
        c.set_xlabel("Rg (nm)"); c.set_ylabel("Inter-helix angle (°)")
        for ax, l in zip((a, b, c), "abc" if i == 0 else "def"): S.panel_letter(ax, l)
    return S.save(fig, os.path.join(OUT, "Fig_4A_solution"))

# ======================================================================= Fig 4.B membranes
def fig_membranes():
    fig = S.figure(height_mm=130)
    gs = fig.add_gridspec(2, 3, hspace=0.75, wspace=0.55, width_ratios=[1.35, 1.35, 1])
    for i, pep in enumerate(PEPS):
        n = len(PEPTIDES[pep]["seq"]); x = np.arange(1, n + 1)
        a = fig.add_subplot(gs[i, 0]); b = fig.add_subplot(gs[i, 1]); c = fig.add_subplot(gs[i, 2])
        res_axis(a, pep); res_axis(b, pep, labels=False)
        for key in ["sol", "popc", "pcpepg", "pepg"]:
            lab, col, ls = COND[key]; D = load(f"{pep}_{key}")
            a.plot(x, helix_profile(D), ls=ls, color=col, marker="o", ms=1.8, lw=1, label=lab)
            if key == "sol": continue
            m = win(D["pr"]["time_ns"]) & bound_mask(D)
            dep = np.nanmean(D["pr"]["depth"][m][:, 0], axis=0) / 10
            b.plot(x, dep, ls=ls, color=col, marker="o", ms=1.8, lw=1)
            ts = D["ts"]; w = ts[win(ts.time_ns) & bound_mask(D)]
            hist, e = np.histogram(w.tilt, bins=np.arange(0, 181, 10), density=True)
            c.plot((e[:-1] + e[1:]) / 2, hist, ls=ls, color=col, lw=1)
        a.set_ylim(0, 1.05); a.set_ylabel("Helix fraction"); a.set_title(NAME[pep], loc="left", fontsize=S.PT_NAME, pad=12)
        b.axhline(0, color="#555555", lw=0.5, ls="--"); b.set_ylabel("Cα height above P plane (nm)")
        b.text(n, 0.03, "phosphate plane", ha="right", va="bottom", fontsize=S.PT_FOOT, color="#555555")
        c.axvline(90, color="#555555", lw=0.5, ls="--"); c.set_xlabel("Helix tilt to normal (°)"); c.set_ylabel("Density"); c.set_xlim(0, 180)
        for ax, l in zip((a, b, c), "abc" if i == 0 else "def"): S.panel_letter(ax, l)
    h = [Line2D([], [], color=COND[k][1], ls=COND[k][2], lw=1.2, label=COND[k][0]) for k in COND]
    fig.legend(handles=h, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 0.97), fontsize=S.PT_ANNOT)
    return S.save(fig, os.path.join(OUT, "Fig_4B_membranes"))

# ======================================================================= Fig 4.C POPC starting orientations (v2)
# v2: binding raster instead of overlaid COM traces (a, d); helix per residue at full width so every residue is labelled (b, e)
def fig_replicates():
    from matplotlib.colors import LinearSegmentedColormap
    fig = S.figure(height_mm=178)
    gs = fig.add_gridspec(4, 2, hspace=0.95, wspace=0.62, width_ratios=[2.3, 1], height_ratios=[1, 1.05, 1, 1.05])
    BIND = LinearSegmentedColormap.from_list("bind", ["#0d366b", "#1c5cab", "#3987e5", "#9ec5f4", "#f4f8fd"])
    rows = []
    for i, pep in enumerate(PEPS):
        sids = [f"{pep}_popc"] + [f"{pep}_popc_r{r}" for r in range(2, 6)]
        n = len(PEPTIDES[pep]["seq"]); x = np.arange(1, n + 1)
        a = fig.add_subplot(gs[2 * i, 0]); c = fig.add_subplot(gs[2 * i, 1]); b = fig.add_subplot(gs[2 * i + 1, :]); res_axis(b, pep)
        sol = load(f"{pep}_sol"); b.plot(x, helix_profile(sol), ls=":", color="#4d4d4d", lw=1, label="water")
        M = []
        for j, sid in enumerate(sids):
            D = load(sid); ts = D["ts"]; col = RUNCOL[j]
            t = ts.time_ns.values; edges = np.arange(0, 1001, 2.0)
            idx = np.digitize(t, edges) - 1
            M.append([np.nanmean(ts.height_com.values[idx == k]) if (idx == k).any() else np.nan for k in range(len(edges) - 1)])
            b.plot(x, helix_profile(D), color=col, lw=0.9, marker="o", ms=1.5)
            bm = bound_mask(D); w = win(t)
            bf = bm[w].mean(); hel = ts.helicity.values[w].mean()
            c.scatter(bf, hel, s=20, color=col, edgecolor="black", lw=0.4, zorder=3)
            rows.append(dict(peptide=pep, run=sid, start=RUNLAB[pep][j], bound_fraction=round(bf, 2), helicity=round(hel, 2)))
        im = a.imshow(np.array(M), aspect="auto", cmap=BIND, vmin=-0.5, vmax=4.0, interpolation="nearest", extent=(0, 1, 4.5, -0.5))
        for j in range(5): a.add_patch(plt.Rectangle((-0.035, j - 0.4), 0.022, 0.8, color=RUNCOL[j], clip_on=False, transform=a.transData))
        a.set_yticks(range(5)); a.set_yticklabels(RUNLAB[pep], fontsize=S.PT_FOOT); a.tick_params(axis="y", length=0, pad=9)
        a.axvline(W0 / 1000, color="white", lw=0.8, ls="--")
        a.set_xlabel("Time (µs)"); a.set_title(NAME[pep], loc="left", fontsize=S.PT_NAME, pad=6)
        cb = fig.colorbar(im, ax=a, fraction=0.035, pad=0.015); cb.set_label("COM height (nm)", fontsize=S.PT_FOOT, labelpad=2); cb.ax.tick_params(labelsize=S.PT_FOOT)
        b.set_ylim(0, 1.05); b.set_ylabel("Helix fraction")
        hs = sol["ts"]; c.axhline(hs.helicity[win(hs.time_ns)].mean(), color="#4d4d4d", ls=":", lw=0.8)
        c.set_xlim(-0.05, 1.05); c.set_ylim(0, 1); c.set_xlabel("Fraction of frames bound"); c.set_ylabel("Helicity")
        for ax, l in zip((a, b, c), "abc" if i == 0 else "def"): S.panel_letter(ax, l)
    h = [Line2D([], [], color=RUNCOL[j], lw=1.2, label=RUNLAB["mel"][j]) for j in range(5)] + [Line2D([], [], color="#4d4d4d", ls=":", lw=1, label="water")]
    fig.legend(handles=h, ncol=6, loc="lower center", bbox_to_anchor=(0.5, 0.93), fontsize=S.PT_ANNOT)
    pd.DataFrame(rows).to_csv(os.path.join(ANALYSIS_DIR, "chapter4_popc_runs.csv"), index=False)
    return S.save(fig, os.path.join(OUT, "Fig_4C_popc_starting_orientations"))

# ======================================================================= Fig 4.D four peptides
def fig_aggregation():
    fig = S.figure(height_mm=125)
    gs = fig.add_gridspec(2, 3, hspace=0.75, wspace=0.85, width_ratios=[1, 1.35, 1])
    for i, pep in enumerate(PEPS):
        D = load(f"{pep}_agg"); n = len(PEPTIDES[pep]["seq"]); x = np.arange(1, n + 1)
        a = fig.add_subplot(gs[i, 0]); b = fig.add_subplot(gs[i, 1]); c = fig.add_subplot(gs[i, 2]); res_axis(b, pep)
        cl = D["cl"]; a.step(cl.time_ns / 1000, cl.largest_cluster, where="post", color="#184f95", lw=0.9, label="largest cluster")
        a2 = a.twinx()
        ts = D["ts"]; nb = np.zeros(len(D["mem"]))
        for k in range(1, 5): nb += (D["mem"][f"lip_contacts_copy{k}"].values > 0)
        a2.plot(D["mem"].time_ns / 1000, rolling(D["mem"].time_ns.values, nb, 10), color="#eb6834", lw=0.8)
        a.set_ylim(0.5, 4.5); a.set_yticks([1, 2, 3, 4]); a2.set_ylim(-0.2, 4.2); a2.set_yticks([0, 1, 2, 3, 4])
        a.set_xlabel("Time (µs)"); a.set_ylabel("Largest cluster", color="#184f95"); a2.set_ylabel("Peptides on bilayer", color="#eb6834")
        a2.spines["right"].set_visible(True); a.set_title(NAME[pep], loc="left", fontsize=S.PT_NAME, pad=12)
        b.plot(x, helix_profile(load(f"{pep}_sol")), ls=":", color="#4d4d4d", lw=1, label="1 peptide, water")
        b.plot(x, helix_profile(load(f"{pep}_popc")), color="#2a78d6", lw=1, label="1 peptide, POPC")
        b.plot(x, helix_profile(D), color="#184f95", lw=1.4, ls="-.", marker="o", ms=1.8, label="4 peptides, POPC (mean)")
        b.set_ylim(0, 1.05); b.set_ylabel("Helix fraction")

        blk = D["blk"]; t0 = blk["block_t0_ns"]; sel = np.isfinite(t0) & (t0 >= W0 - 1e-6)
        M = blk["rr_map"][sel].sum(0) / max(blk["rr_n"][sel].sum(), 1)
        im = c.imshow(M, origin="lower", cmap=S.SEQ, vmin=0, vmax=1, extent=(0.5, n + 0.5, 0.5, n + 0.5))
        c.set_xticks(range(5, n + 1, 5)); c.set_yticks(range(5, n + 1, 5))
        c.set_xlabel("Residue, peptide A"); c.set_ylabel("Residue, peptide B")
        cb = fig.colorbar(im, ax=c, fraction=0.046, pad=0.03); cb.set_label("Contact frequency", fontsize=S.PT_FOOT); cb.ax.tick_params(labelsize=S.PT_FOOT)
        for ax, l in zip((a, b, c), "abc" if i == 0 else "def"): S.panel_letter(ax, l)
    h = [Line2D([], [], color="#4d4d4d", ls=":", lw=1, label="1 peptide, water"), Line2D([], [], color="#2a78d6", lw=1, label="1 peptide, POPC"),
         Line2D([], [], color="#184f95", ls="-.", marker="o", ms=2, lw=1.4, label="4 peptides, POPC (mean of copies)")]
    fig.legend(handles=h, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 0.97), fontsize=S.PT_ANNOT)
    return S.save(fig, os.path.join(OUT, "Fig_4D_four_peptides"))

if __name__ == "__main__":
    for f in ([globals()[a] for a in sys.argv[1:]] or (fig_solution, fig_membranes, fig_replicates, fig_aggregation)):
        print("saved", f())

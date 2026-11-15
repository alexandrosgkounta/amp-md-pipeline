#!/usr/bin/env python3
"""Cross-condition figures and summary table for one peptide.

Usage:  python amp_figures.py mel [--window-frac 0.5] [--start-ns T]

Averages (profiles, distributions, summary values) use the production window: by default
the last half of each trajectory (--window-frac 0.5), or everything after --start-ns.
Time series always show the full trajectory with the window start marked.
Errors are block standard errors (5 contiguous blocks in the window) and are descriptive:
each condition is one trajectory until the replicates are added to systems.tsv.

Writes figures/<peptide>/Fig_*.svg|pdf|png and analysis/summary_<peptide>.csv
"""
import os, json, argparse, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import amp_style as S
from amp_config import read_systems, PEPTIDES, CONDITIONS, COPY_COLOURS, ANALYSIS_DIR, FIG_DIR, SN1, SN2

ORDER = ["solution", "popc", "pcpepg", "pepg", "agg_popc", "agg_surf"]


# ----------------------------------------------------------------------------- data
def load(sid):
    d = os.path.join(ANALYSIS_DIR, sid)
    if not os.path.exists(os.path.join(d, "timeseries.csv")):
        return None
    D = dict(sid=sid, ts=pd.read_csv(os.path.join(d, "timeseries.csv")),
             pr=dict(np.load(os.path.join(d, "perres.npz"))), meta=json.load(open(os.path.join(d, "meta.json"))))
    for name, f in (("mem", "membrane_ts.csv"), ("water", "water_ts.csv"), ("pairs", "agg_pairs.csv"),
                    ("clusters", "agg_clusters.csv")):
        p = os.path.join(d, f)
        D[name] = pd.read_csv(p) if os.path.exists(p) else None
    p = os.path.join(d, "blocks.npz")
    D["blk"] = dict(np.load(p, allow_pickle=True)) if os.path.exists(p) else {}
    p = os.path.join(d, "water_profile.npz")
    D["wprof"] = dict(np.load(p)) if os.path.exists(p) else None
    return D


def window_mask(t, frac, start_ns):
    if start_ns is not None:
        return t >= start_ns
    return t >= t.min() + (1 - frac) * (t.max() - t.min())


def block_se(x, nb=5):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < nb * 2:
        return np.nan
    m = [b.mean() for b in np.array_split(x, nb)]
    return np.std(m, ddof=1) / np.sqrt(nb)


def equil_point(t, x, tol=1.0):
    """Reverse cumulative averaging (as in the LARP6 construction guide)."""
    x = np.asarray(x, float)
    if len(x) < 10: return np.nan
    fin = x[len(x) // 2:].mean(); sd = x[len(x) // 2:].std() or 1e-9
    for i in range(len(x)):
        if abs(x[i:].mean() - fin) < tol * sd:
            return float(t[i])
    return float(t[len(x) // 2])


def rolling(t, x, ns=10.0):
    dt = np.median(np.diff(t)) if len(t) > 1 else 1
    w = max(1, int(round(ns / dt)))
    return pd.Series(x).rolling(w, center=True, min_periods=1).mean().values


def blk_window(blk, frac, start_ns):
    """Which accumulator blocks fall inside the window."""
    t0 = blk.get("block_t0_ns")
    if t0 is None: return None
    ok = np.isfinite(t0)
    tmax = np.nanmax(blk["block_t1_ns"]); tmin = np.nanmin(t0)
    thr = start_ns if start_ns is not None else tmin + (1 - frac) * (tmax - tmin)
    sel = ok & (t0 >= thr - 1e-9)
    if not sel.any(): sel = ok
    return sel


def style(cond):
    lab, col, ls = CONDITIONS[cond]
    return dict(color=col, ls=ls, label=lab)


def legend_conditions(fig, conds, y=1.0):
    h = [Line2D([], [], **{k: v for k, v in style(c).items()}, lw=1.4) for c in conds]
    fig.legend(handles=h, loc="lower center", ncol=len(conds), bbox_to_anchor=(0.5, y),
               fontsize=S.PT_ANNOT, handlelength=2.2, columnspacing=1.2)


# ----------------------------------------------------------------------------- figures
def fig_conformation(sets, pep, out, frac, start_ns):
    fig = S.figure(height_mm=118)
    gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.3)
    a = fig.add_subplot(gs[0, :]); b = fig.add_subplot(gs[1, 0]); c = fig.add_subplot(gs[1, 1])
    S.residue_axis(a, pep["seq"], pep["segments"])
    n = len(pep["seq"])
    for cond, Ds in sets.items():
        H = []
        for D in Ds:
            m = window_mask(D["pr"]["time_ns"], frac, start_ns)
            H.append((D["pr"]["dssp"][m] == "H").mean(axis=(0, 1)))
        a.plot(range(1, n + 1), np.mean(H, 0), marker="o", ms=2.2, lw=1.0, **style(cond))
        for D in Ds:
            ts = D["ts"]; t = ts.groupby("time_ns")["helicity"].mean()
            b.plot(t.index / 1000, rolling(t.index.values, t.values) * 100, lw=0.9, **style(cond))
            k = ts[window_mask(ts["time_ns"].values, frac, start_ns)]["kink"].values
            if len(k) > 1:
                hist, edges = np.histogram(k, bins=np.arange(0, 181, 7.5), density=True)
                c.plot((edges[1:] + edges[:-1]) / 2, hist, lw=1.0, **style(cond))
    a.set_ylabel("Helix fraction (DSSP H)"); a.set_ylim(0, 1.05)
    b.set_xlabel("Time (µs)"); b.set_ylabel("Helicity (%)"); b.set_ylim(0, 100)
    c.set_xlabel("Inter-helix bend angle (°)"); c.set_ylabel("Probability density")
    for ax, l in zip((a, b, c), "abc"): S.panel_letter(ax, l, dx=-0.06 if ax is a else -0.15)
    legend_conditions(fig, list(sets), y=0.95)
    return S.save(fig, os.path.join(out, "Fig_M1_conformation"))


def fig_landscape(sets, pep, out, frac, start_ns):
    """Rg against helicity free energy surfaces, one per condition."""
    conds = list(sets); nc = len(conds)
    fig = S.figure(height_mm=42)
    for i, cond in enumerate(conds):
        ax = fig.add_subplot(1, nc, i + 1)
        ts = pd.concat([D["ts"][window_mask(D["ts"]["time_ns"].values, frac, start_ns)] for D in sets[cond]])
        H, xe, ye = np.histogram2d(ts["helicity"] * 100, ts["rg"], bins=[np.linspace(0, 100, 27), np.linspace(0.8, 2.0, 25)])
        with np.errstate(divide="ignore"):
            G = -np.log(H / H.sum()); G -= np.nanmin(G[np.isfinite(G)]) if np.isfinite(G).any() else 0
        G[~np.isfinite(G)] = np.nan
        im = ax.pcolormesh(xe, ye, G.T, cmap=S.SEQ.reversed(), vmin=0, vmax=5, shading="flat")
        ax.set_title(CONDITIONS[cond][0], fontsize=S.PT_FOOT)
        ax.set_xlabel("Helicity (%)")
        if i == 0: ax.set_ylabel("Rg (nm)")
        else: ax.set_yticklabels([])
    cb = fig.colorbar(im, ax=fig.axes, fraction=0.02, pad=0.01); cb.set_label("ΔG (k$_B$T)", fontsize=S.PT_ANNOT)
    return S.save(fig, os.path.join(out, "Fig_M2_conformational_landscape"))


def fig_binding(sets, pep, out, frac, start_ns):
    msets = {c: v for c, v in sets.items() if v[0]["mem"] is not None}
    if not msets: return
    fig = S.figure(height_mm=125)
    gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.32)
    a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[0, 1]); c = fig.add_subplot(gs[1, 0]); d = fig.add_subplot(gs[1, 1])
    n = len(pep["seq"])
    S.residue_axis(b, pep["seq"], pep["segments"], label_segments=False)
    for cond, Ds in msets.items():
        st = style(cond)
        prof = []
        for D in Ds:
            ts = D["ts"]
            for k, g in ts.groupby("copy"):
                lw = 0.9 if g["copy"].iloc[0] == 1 or not cond.startswith("agg") else 0.5
                a.plot(g["time_ns"] / 1000, rolling(g["time_ns"].values, g["height_com"].values, 2), lw=lw,
                       alpha=1 if lw > 0.6 else 0.7, color=st["color"], ls=st["ls"])
            m = window_mask(D["pr"]["time_ns"], frac, start_ns)
            prof.append(np.nanmean(D["pr"]["depth"][m], axis=(0, 1)) / 10)
            w = ts[window_mask(ts["time_ns"].values, frac, start_ns)]
            hist, e = np.histogram(w["tilt"], bins=np.arange(0, 181, 7.5), density=True)
            c.plot((e[1:] + e[:-1]) / 2, hist, **st)
            mem = D["mem"]
            cols = [x for x in mem.columns if x.startswith("lip_contacts")]
            d.plot(mem["time_ns"] / 1000, rolling(mem["time_ns"].values, mem[cols].sum(1).values / len(cols), 5), lw=0.9,
                   color=st["color"], ls=st["ls"])
        p = np.mean(prof, 0)
        b.plot(range(1, n + 1), p, marker="o", ms=2.2, **st)
    for ax in (a, b):
        ax.axhline(0, color="#555555", lw=0.5, ls="--")
    a.text(0.99, 0.0, "phosphate plane", transform=a.get_yaxis_transform(), ha="right", va="bottom", fontsize=S.PT_FOOT, color="#555555")
    a.set_xlabel("Time (µs)"); a.set_ylabel("COM height above P plane (nm)")
    b.set_ylabel("Cα height above P plane (nm)")
    c.set_xlabel("Helix tilt to membrane normal (°)"); c.set_ylabel("Probability density")
    c.axvline(90, color="#555555", lw=0.5, ls="--")
    c.text(90, 1.0, " parallel to surface", transform=c.get_xaxis_transform(), fontsize=S.PT_FOOT, color="#555555", va="top")
    d.set_xlabel("Time (µs)"); d.set_ylabel("Residues touching lipid")
    for ax, l in zip((a, b, c, d), "abcd"): S.panel_letter(ax, l, dx=-0.17)
    legend_conditions(fig, list(msets), y=0.95)
    return S.save(fig, os.path.join(out, "Fig_M3_membrane_binding"))


def fig_lipid_contacts(sets, pep, out, frac, start_ns):
    msets = {c: v for c, v in sets.items() if v[0]["mem"] is not None}
    if not msets: return
    n = len(pep["seq"]); x = np.arange(1, n + 1)
    nrow = len(msets)
    fig = S.figure(height_mm=30 + 30 * nrow)
    gs = fig.add_gridspec(nrow, 2, width_ratios=[3, 1], hspace=0.65, wspace=0.25)
    cls_cols = ["#2a78d6", "#eda100", "#1baf7a"]; cls_lab = ["Head group and phosphate", "Glycerol and ester", "Acyl chains"]
    summ = []
    for i, (cond, Ds) in enumerate(msets.items()):
        ax = fig.add_subplot(gs[i, 0]); S.residue_axis(ax, pep["seq"], pep["segments"], label_segments=(i == 0))
        occ = np.mean([D["pr"]["contact"][window_mask(D["pr"]["time_ns"], frac, start_ns)].mean(axis=(0, 1)) for D in Ds], 0)
        ax.plot(x, occ[:, 0], color="#333333", lw=1.2, label="Any lipid")
        for j in range(3):
            ax.plot(x, occ[:, j + 1], color=cls_cols[j], lw=1.0, marker="o", ms=1.8, label=cls_lab[j])
        ax.set_ylim(0, 1.05); ax.set_ylabel("Contact occupancy")
        ax.text(0.01, 0.97, CONDITIONS[cond][0], transform=ax.transAxes, va="top", fontsize=S.PT_ANNOT, fontweight="bold")
        S.panel_letter(ax, "abcdefghijklmnop"[2 * i], dx=-0.08)
        if i == 0: ax.legend(ncol=4, loc="lower left", bbox_to_anchor=(0, 1.12), fontsize=S.PT_FOOT)
        # first-shell enrichment
        ax2 = fig.add_subplot(gs[i, 1])
        D = Ds[0]; mem = D["mem"]; w = mem[window_mask(mem["time_ns"].values, frac, start_ns)]
        lt = list(D["pr"]["lip_types"])
        nshell = w[[f"shell_{l}" for l in lt]].sum(1).replace(0, np.nan)
        for j, l in enumerate(lt):
            # enrichment = fraction of type l in the first shell / fraction of l in the whole bilayer
            e = (w[f"shell_{l}"] / nshell) / BULK[D["sid"]][l]
            ax2.bar(j, np.nanmean(e), yerr=block_se(e.values), color=["#2a78d6", "#1baf7a", "#eb6834"][LIPS.index(l)],
                    width=0.65, error_kw=dict(lw=0.6, capsize=1.5))
            summ.append(dict(sys_id=D["sid"], lipid=l, bulk_fraction=BULK[D["sid"]][l],
                             enrichment=np.nanmean(e), se=block_se(e.values)))
        ax2.axhline(1, color="#555555", lw=0.5, ls="--")
        ax2.set_xticks(range(len(lt))); ax2.set_xticklabels(lt); ax2.set_ylabel("First-shell enrichment")
        S.panel_letter(ax2, "abcdefghijklmnop"[2 * i + 1], dx=-0.4)
    S.save(fig, os.path.join(out, "Fig_M4_lipid_contacts"))
    return pd.DataFrame(summ)


def fig_perturbation(sets, pep, out, frac, start_ns):
    msets = {c: v for c, v in sets.items() if v[0]["mem"] is not None}
    if not msets: return
    nc = len(msets)
    fig = S.figure(height_mm=125)
    gs = fig.add_gridspec(3, nc, hspace=0.75, wspace=0.3, height_ratios=[1, 1, 1])
    grp_lab = ["Proximal, < 1 nm", "Proximal, 1 to 2 nm", "Proximal, > 2 nm", "Distal leaflet"]
    grp_col = ["#eb6834", "#eda100", "#2a78d6", "#7a7a7a"]; grp_ls = ["-", "--", "-", ":"]
    for i, (cond, Ds) in enumerate(msets.items()):
        blk = Ds[0]["blk"]; sel = blk_window(blk, frac, start_ns)
        ax = fig.add_subplot(gs[0, i])
        ssum = blk["scd_sum"][sel].sum(0); sn = blk["scd_n"][sel].sum(0)
        scd = np.where(sn > 0, -ssum / np.maximum(sn, 1), np.nan)
        for g in range(4):
            ax.plot(np.arange(2, 2 + len(SN1)), scd[g, 0, :len(SN1)], color=grp_col[g], ls=grp_ls[g], lw=1.0,
                    marker="o", ms=1.6, label=grp_lab[g])
        ax.set_title(CONDITIONS[cond][0], fontsize=S.PT_ANNOT)
        ax.set_xlabel("sn-1 carbon"); ax.set_ylim(0, 0.3)
        if i == 0:
            ax.set_ylabel("−S$_{CD}$"); S.panel_letter(ax, "a", dx=-0.35)
            ax.legend(ncol=4, loc="lower left", bbox_to_anchor=(0, 1.2), fontsize=S.PT_FOOT)
        else: ax.set_yticklabels([])
        ax = fig.add_subplot(gs[1, i])
        r = blk["rbins"]; rc = (r[1:] + r[:-1]) / 20
        for leaf, col, ls, lab in ((0, "#eb6834", "-", "Proximal"), (1, "#7a7a7a", ":", "Distal")):
            s_ = blk["rad_sum"][sel, leaf].sum(0); n_ = blk["rad_n"][sel, leaf].sum(0)
            h = np.where(n_ > 20, s_ / np.maximum(n_, 1) / 10, np.nan)
            ax.plot(rc, h, color=col, ls=ls, lw=1.0, marker="o", ms=1.6, label=lab)
        ax.set_xlabel("Lateral distance (nm)")
        if i == 0:
            ax.set_ylabel("P height (nm)"); S.panel_letter(ax, "b", dx=-0.35); ax.legend(fontsize=S.PT_FOOT)
    a = fig.add_subplot(gs[2, :max(1, nc // 2)]); b = fig.add_subplot(gs[2, max(1, nc // 2):]) if nc > 1 else None
    for cond, Ds in msets.items():
        for D in Ds:
            m = D["mem"]; a.plot(m["time_ns"] / 1000, rolling(m["time_ns"].values, m["thickness_PP"].values, 10), lw=0.9, **style(cond))
            if b is not None and D["water"] is not None:
                w = D["water"]; b.plot(w["time_ns"] / 1000, rolling(w["time_ns"].values, w["water_core_10A"].values, 10), lw=0.9, **style(cond))
    a.set_xlabel("Time (µs)"); a.set_ylabel("P–P thickness (nm)"); S.panel_letter(a, "c", dx=-0.15)
    if b is not None:
        b.set_xlabel("Time (µs)"); b.set_ylabel("Core waters (|z| < 0.5 nm)"); S.panel_letter(b, "d", dx=-0.15)
    legend_conditions(fig, list(msets), y=0.0 - 0.07)
    return S.save(fig, os.path.join(out, "Fig_M5_membrane_perturbation"))


def fig_density(sets, pep, out, frac, start_ns):
    msets = {c: v for c, v in sets.items() if v[0]["mem"] is not None}
    if not msets: return
    nc = len(msets)
    fig = S.figure(height_mm=55)
    cols = {"Phosphate": "#eda100", "Carbonyl O": "#1baf7a", "Terminal methyl": "#7a7a7a", "Peptide": "#333333",
            "Lys NZ": "#2a78d6", "Arg CZ": "#184f95", "Aromatic ring": "#eb6834", "PG glycerol": "#e87ba4"}
    for i, (cond, Ds) in enumerate(msets.items()):
        ax = fig.add_subplot(1, nc, i + 1)
        blk = Ds[0]["blk"]; sel = blk_window(blk, frac, start_ns)
        z = blk["zbins"]; zc = (z[1:] + z[:-1]) / 20
        dens = blk["dens"][sel].sum(0) / blk["dens_vol"][sel].sum()
        names = list(blk["dens_names"])
        for j, nm in enumerate(names):
            if nm not in cols: continue
            y = dens[j] / (dens[j].max() or 1)
            ax.plot(zc, y, color=cols[nm], lw=0.9, label=nm, ls="--" if nm in ("Lys NZ", "Arg CZ", "Aromatic ring") else "-")
        if Ds[0]["wprof"] is not None:
            wp = Ds[0]["wprof"]["density"]; ax.plot(zc, wp / (wp.max() or 1), color="#9ec5f4", lw=1.6, label="Water", zorder=0)
        ax.set_title(CONDITIONS[cond][0], fontsize=S.PT_ANNOT); ax.set_xlabel("z from bilayer centre (nm)")
        ax.set_xlim(-4.5, 4.5)
        if i == 0: ax.set_ylabel("Normalised density")
        else: ax.set_yticklabels([])
    hl = {}
    for ax in fig.axes:
        for h, l in zip(*ax.get_legend_handles_labels()): hl.setdefault(l, h)
    fig.legend(list(hl.values()), list(hl), ncol=5, loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=S.PT_FOOT)
    fig.text(0.5, -0.12, "Peptide side is positive z; each curve scaled to its own maximum.", ha="center", fontsize=S.PT_FOOT)
    return S.save(fig, os.path.join(out, "Fig_S2_density_profiles"))


def fig_aggregation(sets, pep, out, frac, start_ns):
    for cond in [c for c in sets if c.startswith("agg")]:
        _fig_aggregation_one(sets[cond][0], cond, pep, out, frac, start_ns)


def _fig_aggregation_one(D, cond, pep, out, frac, start_ns):
    if D["pairs"] is None: return
    fig = S.figure(height_mm=120)
    gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.5)
    a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[0, 1]); c = fig.add_subplot(gs[1, 0]); d = fig.add_subplot(gs[1, 1])
    pal = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#7a7a7a"]
    for j, (pair, g) in enumerate(D["pairs"].groupby("pair")):
        a.plot(g["time_ns"] / 1000, rolling(g["time_ns"].values, g["dmin"].values, 2), color=pal[j % 6], lw=0.9, label=f"Pair {pair}")
    from matplotlib.ticker import FixedLocator, ScalarFormatter, NullLocator
    a.set_yscale("log"); a.axhline(0.45, color="#555555", lw=0.5, ls="--")
    a.yaxis.set_major_locator(FixedLocator([0.2, 0.5, 1, 2, 5, 10])); a.yaxis.set_minor_locator(NullLocator())
    a.yaxis.set_major_formatter(ScalarFormatter())
    a.set_xlabel("Time (µs)"); a.set_ylabel("Minimum distance (nm)")
    a.legend(ncol=3, fontsize=S.PT_FOOT, loc="lower left", bbox_to_anchor=(0, 1.02))
    cl = D["clusters"]
    b.step(cl["time_ns"] / 1000, cl["largest_cluster"], where="post", color="#184f95", lw=1.0)
    b.set_ylim(0.5, 4.5); b.set_yticks([1, 2, 3, 4]); b.set_xlabel("Time (µs)"); b.set_ylabel("Largest peptide cluster")
    blk = D["blk"]; sel = blk_window(blk, frac, start_ns)
    M = blk["rr_map"][sel].sum(0) / max(blk["rr_n"][sel].sum(), 1)
    im = c.imshow(M, origin="lower", cmap=S.SEQ, vmin=0, extent=(0.5, len(pep["seq"]) + 0.5, 0.5, len(pep["seq"]) + 0.5))
    c.set_xticks(range(1, len(pep["seq"]) + 1, 5)); c.set_yticks(range(1, len(pep["seq"]) + 1, 5))
    c.set_xlabel("Residue, peptide A"); c.set_ylabel("Residue, peptide B")
    cb = fig.colorbar(im, ax=c, fraction=0.046, pad=0.03); cb.set_label("Contact frequency", fontsize=S.PT_FOOT)
    for k, g in D["ts"].groupby("copy"):
        d.plot(g["time_ns"] / 1000, rolling(g["time_ns"].values, g["height_com"].values, 2), color=COPY_COLOURS[(k - 1) % 4], lw=0.9, label=f"Peptide {k}")
    d.axhline(0, color="#555555", lw=0.5, ls="--")
    d.set_xlabel("Time (µs)"); d.set_ylabel("COM height above P plane (nm)"); d.legend(fontsize=S.PT_FOOT)
    for ax, l in zip((a, b, c, d), "abcd"): S.panel_letter(ax, l, dx=-0.2)
    name = "Fig_M6_aggregation" if cond == "agg_popc" else f"Fig_M6_aggregation_{cond}"
    return S.save(fig, os.path.join(out, name))


def fig_hdx_proxy(sets, pep, out, frac, start_ns):
    n = len(pep["seq"]); x = np.arange(1, n + 1)
    fig = S.figure(height_mm=110)
    gs = fig.add_gridspec(3, 1, hspace=0.45)
    axs = [fig.add_subplot(gs[i]) for i in range(3)]
    for ax in axs: S.residue_axis(ax, pep["seq"], pep["segments"], label_segments=(ax is axs[0]))
    pro = [i for i, c in enumerate(pep["seq"]) if c == "P"] + [0]
    for cond, Ds in sets.items():
        lp, le, nh = [], [], []
        for D in Ds:
            m = window_mask(D["pr"]["time_ns"], frac, start_ns); pr = D["pr"]
            lp.append((0.35 * pr["nc_pep"][m] + 2 * pr["nh_pep"][m]).mean(axis=(0, 1)))
            le.append((0.35 * pr["nc_env"][m] + 2 * pr["nh_env"][m]).mean(axis=(0, 1)))
            nh.append((pr["nh_pep"][m] > 0).mean(axis=(0, 1)))
        for arr, ax in ((lp, axs[0]), (le, axs[1]), (nh, axs[2])):
            y = np.mean(arr, 0).astype(float); y[pro] = np.nan
            ax.plot(x, y, marker="o", ms=2, **style(cond))
    axs[0].set_ylabel("ln P, peptide only"); axs[1].set_ylabel("ln P, with lipids and\nother peptides")
    axs[2].set_ylabel("Amide H-bond\noccupancy (intra-peptide)"); axs[2].set_ylim(0, 1.05)
    for ax, l in zip(axs, "abc"): S.panel_letter(ax, l, dx=-0.07)
    legend_conditions(fig, list(sets), y=0.92)
    fig.text(0.5, 0.02, "ln P = 0.35 N$_c$ + 2.0 N$_h$ (Best and Vendruscolo 2006). Pro and residue 1 have no amide H.",
             ha="center", fontsize=S.PT_FOOT)
    return S.save(fig, os.path.join(out, "Fig_M7_HDX_protection_proxy"))


def fig_maps(allD, pep, out):
    """SI: DSSP and insertion depth, residue by time, one strip per system and copy."""
    strips = [(D, k) for D in allD for k in range(D["pr"]["dssp"].shape[1])]
    fig = S.figure(height_mm=14 * len(strips) + 15)
    gs = fig.add_gridspec(len(strips), 2, hspace=1.1, wspace=0.12)
    for i, (D, k) in enumerate(strips):
        t = D["pr"]["time_ns"] / 1000; n = len(pep["seq"])
        ax = fig.add_subplot(gs[i, 0])
        ax.pcolormesh(np.r_[t, t[-1] + (t[1] - t[0] if len(t) > 1 else 1e-3)], np.arange(0.5, n + 1),
                      (D["pr"]["dssp"][:, k] == "H").T.astype(float), cmap=S.HELIX_CMAP, vmin=0, vmax=1, shading="flat", rasterized=True)
        lab = CONDITIONS[D["meta"]["system"]["condition"]][0] + (f", peptide {k + 1}" if D["pr"]["dssp"].shape[1] > 1 else "")
        ax.set_title(lab, fontsize=S.PT_FOOT, loc="left", pad=1.5); ax.set_yticks([1, 10, 20]); ax.set_ylabel("Res.")
        if i == len(strips) - 1: ax.set_xlabel("Time (µs)")
        ax2 = fig.add_subplot(gs[i, 1])
        dep = D["pr"]["depth"][:, k] / 10
        if np.isfinite(dep).any():
            im = ax2.pcolormesh(np.r_[t, t[-1] + (t[1] - t[0] if len(t) > 1 else 1e-3)], np.arange(0.5, n + 1), dep.T,
                                cmap=S.DIV, vmin=-1.5, vmax=1.5, shading="flat", rasterized=True)
            ax2.set_yticks([])
            last_im = im
        else:
            ax2.axis("off")
        if i == len(strips) - 1: ax2.set_xlabel("Time (µs)")
    if "last_im" in dir():
        cb = fig.colorbar(last_im, ax=fig.axes[1::2], fraction=0.03, pad=0.02, shrink=0.5)
        cb.set_label("Cα height above P plane (nm)", fontsize=S.PT_FOOT)
    fig.axes[0].text(0, 1.6, "Helix (DSSP H) in blue", transform=fig.axes[0].transAxes, fontsize=S.PT_FOOT)
    fig.axes[1].text(0, 1.6, "Insertion depth: blue in water, orange below the phosphate plane",
                     transform=fig.axes[1].transAxes, fontsize=S.PT_FOOT)
    return S.save(fig, os.path.join(out, "Fig_S1_residue_time_maps"))


def fig_qc(allD, pep, out, frac, start_ns):
    fig = S.figure(height_mm=95)
    gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.3)
    axs = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(4)]
    for D in allD:
        cond = D["meta"]["system"]["condition"]; st = style(cond)
        ts = D["ts"].groupby("time_ns").mean(numeric_only=True)
        t = ts.index.values
        axs[0].plot(t / 1000, rolling(t, ts["rmsd_bb"].values, 5), lw=0.8, **st)
        axs[1].plot(t / 1000, rolling(t, ts["rg"].values, 5), lw=0.8, **st)
        if D["mem"] is not None:
            m = D["mem"]
            axs[2].plot(m["time_ns"] / 1000, m["apl"], lw=0.8, **st)
            axs[3].plot(m["time_ns"] / 1000, m["box_z"], lw=0.8, **st)
        tw = t.min() + (1 - frac) * (t.max() - t.min()) if start_ns is None else start_ns
        for ax in axs[:2]: ax.axvline(tw / 1000, color=st["color"], lw=0.4, ls=":")
    labs = ["Backbone RMSD to start (nm)", "Radius of gyration (nm)", "Box area per lipid (nm$^2$)", "Box z (nm)"]
    for ax, l, yl in zip(axs, "abcd", labs):
        ax.set_ylabel(yl); ax.set_xlabel("Time (µs)"); S.panel_letter(ax, l, dx=-0.17)
    legend_conditions(fig, sorted({D["meta"]["system"]["condition"] for D in allD}, key=ORDER.index), y=0.95)
    return S.save(fig, os.path.join(out, "Fig_S3_QC"))


# ----------------------------------------------------------------------------- summary
def summary(allD, frac, start_ns):
    rows = []
    for D in allD:
        s = D["meta"]["system"]; ts = D["ts"]
        for k, g in ts.groupby("copy"):
            t = g["time_ns"].values; m = window_mask(t, frac, start_ns); w = g[m]
            r = dict(sys_id=D["sid"], condition=s["condition"], replicate=s["replicate"], copy=k,
                     analysed_to_ns=D["meta"].get("end_ns_applied", np.nan),
                     t_total_ns=round(t.max(), 1), window_start_ns=round(t[m].min(), 1), n_frames=int(m.sum()),
                     t_equil_rmsd_ns=equil_point(t, g["rmsd_bb"].values),
                     t_equil_helicity_ns=equil_point(t, g["helicity"].values))
            for col in ["helicity", "rg", "e2e", "kink", "rmsd_bb", "height_com", "tilt", "depth_min"]:
                if col in w:
                    r[col] = w[col].mean(); r[col + "_se"] = block_se(w[col].values)
            if "height_com" in g:
                ct = D["mem"][f"lip_contacts_copy{k}"].values
                bound = ct > 0
                r["bound_fraction"] = bound[m].mean()
                # first time with lipid contact sustained for >= 10 ns
                dt = np.median(np.diff(t)) if len(t) > 1 else 1; nw = max(1, int(10 / dt))
                run = np.convolve(bound.astype(int), np.ones(nw, int), "valid") == nw
                r["t_bind_ns"] = float(t[np.argmax(run)]) if run.any() else np.nan
                r["t_equil_height_ns"] = equil_point(t, g["height_com"].values)
            rows.append(r)
        if D["mem"] is not None:
            mm = D["mem"]; m = window_mask(mm["time_ns"].values, frac, start_ns)
            rows[-1].update(thickness_PP=mm["thickness_PP"][m].mean(), thickness_PP_se=block_se(mm["thickness_PP"][m]),
                            apl=mm["apl"][m].mean(), apl_se=block_se(mm["apl"][m]))
        if D["water"] is not None:
            ww = D["water"]; m = window_mask(ww["time_ns"].values, frac, start_ns)
            rows[-1].update(water_core=ww["water_core_10A"][m].mean(), water_core_near_pep=ww["water_core_near_peptide"][m].mean())
        if D["clusters"] is not None:
            cl = D["clusters"]; m = window_mask(cl["time_ns"].values, frac, start_ns)
            rows[-1].update(largest_cluster_mean=cl["largest_cluster"][m].mean(),
                            frac_any_dimer=(cl["largest_cluster"][m] >= 2).mean())
    return pd.DataFrame(rows)


LIPS = ["POPC", "POPE", "POPG"]
BULK = {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("peptide", default="mel", nargs="?")
    ap.add_argument("--window-frac", type=float, default=0.5)
    ap.add_argument("--start-ns", type=float, default=None)
    a = ap.parse_args()
    S.apply()
    pep = PEPTIDES[a.peptide]
    out = os.path.join(FIG_DIR, a.peptide); os.makedirs(out, exist_ok=True)
    sets, allD = {}, []
    import MDAnalysis as mda
    for s in read_systems():
        if s["peptide"] != a.peptide: continue
        D = load(s["sys_id"])
        if D is None:
            print(f"  {s['sys_id']}: no analysis yet, skipped"); continue
        # bulk lipid fractions from the tpr (for enrichment)
        if D["mem"] is not None:
            u = mda.Universe(os.path.join(s["path"], s["tpr"]))
            rn = u.select_atoms("name P").resnames
            BULK[s["sys_id"]] = {l: (rn == l).mean() for l in LIPS}
        sets.setdefault(s["condition"], []).append(D); allD.append(D)
    sets = {c: sets[c] for c in ORDER if c in sets}
    f = a.window_frac; st = a.start_ns
    fig_conformation(sets, pep, out, f, st)
    fig_landscape(sets, pep, out, f, st)
    fig_binding(sets, pep, out, f, st)
    enr = fig_lipid_contacts(sets, pep, out, f, st)
    fig_perturbation(sets, pep, out, f, st)
    fig_density(sets, pep, out, f, st)
    fig_aggregation(sets, pep, out, f, st)
    fig_hdx_proxy(sets, pep, out, f, st)
    fig_maps(allD, pep, out)
    fig_qc(allD, pep, out, f, st)
    sm = summary(allD, f, st)
    sm.to_csv(os.path.join(ANALYSIS_DIR, f"summary_{a.peptide}.csv"), index=False, float_format="%.4g")
    if enr is not None:
        enr.to_csv(os.path.join(ANALYSIS_DIR, f"lipid_enrichment_{a.peptide}.csv"), index=False, float_format="%.4g")
    print(f"figures -> {out}\nsummary -> {os.path.join(ANALYSIS_DIR, f'summary_{a.peptide}.csv')}")


if __name__ == "__main__":
    main()

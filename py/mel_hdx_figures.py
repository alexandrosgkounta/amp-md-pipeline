import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import amp_style as S; from amp_config import ANALYSIS_DIR, FIG_DIR, PEPTIDES
S.apply(); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
A = os.path.join(ANALYSIS_DIR, "hdx_compare"); O = os.path.join(FIG_DIR, "chapter4", "hdx")
SEQ = PEPTIDES["mel"]["seq"]; n = len(SEQ); x = np.arange(1, n + 1)
SEGS = [("2 to 5", 2, 5), ("6 to 13", 6, 13), ("14 to 16", 14, 16), ("17 to 19", 17, 19), ("20 to 26", 20, 26)]
TS = np.array([3, 10, 30, 60, 120.0])
from hdx_observed import load as _load_obs
OBS, _, _ = _load_obs("melittin", segments=[s[0] for s in SEGS], times=list(TS))
ST = {"water": ("#4d4d4d", ":"), "POPC, five runs": ("#2a78d6", "-"), "POPE:POPG 7:3": ("#eb6834", "--"), "four peptides": ("#184f95", "-.")}
d = pd.read_csv(os.path.join(A, "mel_md_descriptors.csv"))
# ---------------- H1 per residue descriptors with the HDX map
fig = S.figure(height_mm=190); gs = fig.add_gridspec(6, 1, hspace=0.55, height_ratios=[1.1, 1, 1, 1, 1, 1])
ax0 = fig.add_subplot(gs[0]); S.residue_axis(ax0, SEQ, PEPTIDES["mel"]["segments"])
for j, (t, c) in enumerate(zip([10, 30, 60], ["#9ec5f4", "#2a78d6", "#0d366b"])):
    k = list(TS).index(t)
    for si, (nm, lo, hi) in enumerate(SEGS):
        ax0.hlines(OBS[si, k], lo - 0.45, hi + 0.45, color=c, lw=2)
ax0.axhline(-0.199, color="#555555", lw=0.5, ls=":"); ax0.set_ylim(-0.25, 0.02); ax0.set_ylabel("HDX ΔRFU\n(POPC − free)")
ax0.legend(handles=[Line2D([], [], color=c, lw=2, label=f"{t} s") for t, c in zip([10, 30, 60], ["#9ec5f4", "#2a78d6", "#0d366b"])], ncol=3, loc="lower right", fontsize=S.PT_FOOT)
S.panel_letter(ax0, "a")
panels = [("helix", "Helix fraction", (0, 1.05)), ("hb_pep", "Amide H-bonded\nto peptide", (0, 1.05)), ("hb_lipid", "Amide H-bonded\nto lipid", (0, 0.6)),
          ("water_hb", "Amide H-bonded\nto water", (0, 1.05)), ("lnP_bv_env", "BV ln P\n(lipid counted)", (0, 14))]
for i, (col, lab, yl) in enumerate(panels):
    ax = fig.add_subplot(gs[i + 1]); S.residue_axis(ax, SEQ, PEPTIDES["mel"]["segments"], label_segments=False)
    for st, (c, ls) in ST.items():
        if col == "hb_lipid" and st == "water": continue
        y = d[d.state == st].set_index("resid")[col].reindex(x).values
        if col == "water_hb" and st == "water": y = 1 - d[d.state == st].set_index("resid")["hb_pep"].reindex(x).values   # proxy, no wet trajectory
        if col == "hb_lipid" and st == "four peptides": continue
        ax.plot(x, y, color=c, ls=ls, lw=1, marker="o", ms=1.6)
    ax.set_ylim(*yl); ax.set_ylabel(lab); S.panel_letter(ax, "bcdef"[i])
fig.legend(handles=[Line2D([], [], color=c, ls=l, lw=1.2, label=k) for k, (c, l) in ST.items()], ncol=4, loc="lower center", bbox_to_anchor=(0.5, 0.965), fontsize=S.PT_ANNOT)
S.save(fig, os.path.join(O, "Fig_H1_md_descriptors_vs_hdx"))
# ---------------- H2 segment differences, observed against forward models (k_ch as calculated)
sp = pd.read_csv(os.path.join(A, "mel_segment_pred_vs_obs_k1.csv"))
MOD = {"open": ("H-bond open-state model", "#1baf7a", "-"), "bv_pep": ("BV, peptide only", "#eda100", "--"), "bv_env": ("BV, lipid counted", "#e87ba4", "-")}
fig = S.figure(height_mm=58); gs = fig.add_gridspec(1, 5, wspace=0.12)
for si, (nm, _, _) in enumerate(SEGS):
    ax = fig.add_subplot(gs[si]); ax.plot(TS, OBS[si], "o-", color="black", lw=1.3, ms=3)
    for m, (lab, c, ls) in MOD.items():
        v = sp[(sp.state == "POPC, five runs") & (sp.model == m) & (sp.segment == nm)].sort_values("t_s").pred_dRFU.values
        ax.plot(TS, v, color=c, ls=ls, lw=1.1)
    ax.set_xscale("log"); ax.set_xticks(TS); ax.set_xticklabels([3, 10, 30, 60, 120]); ax.set_ylim(-0.33, 0.04); ax.axhline(0, color="#999999", lw=0.5)
    ax.set_title(f"residues {nm}", fontsize=S.PT_ANNOT); ax.set_xlabel("Labelling time (s)")
    if si == 0: ax.set_ylabel("ΔRFU (POPC − free)")
    else: ax.set_yticklabels([])
    S.panel_letter(ax, "abcde"[si])
h = [Line2D([], [], color="black", marker="o", ms=3, lw=1.3, label="HDX-MS")] + [Line2D([], [], color=c, ls=ls, lw=1.1, label=f"MD, {lab}") for lab, c, ls in MOD.values()]
fig.legend(handles=h, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 0.97), fontsize=S.PT_ANNOT)
S.save(fig, os.path.join(O, "Fig_H2_segment_differences"))
# ---------------- H3 free state and protected-amide gain, k_ch as calculated and doubled
fig = S.figure(height_mm=70); gs = fig.add_gridspec(1, 2, wspace=0.35)
a = fig.add_subplot(gs[0]); b = fig.add_subplot(gs[1])
for sc, ls in ((1, "-"), (2, "--")):
    f = pd.read_csv(os.path.join(A, f"mel_free_state_k{sc}.csv")); g = pd.read_csv(os.path.join(A, f"mel_gain_k{sc}.csv"))
    a.plot(TS, f["intrinsic ceiling (P = 1)"], color="#999999", ls=ls, lw=1)
    a.plot(TS, f["MD water, open"], color="#4d4d4d", ls=ls, lw=1)
    a.plot(TS, f["MD four peptides, open"], color="#184f95", ls=ls, lw=1)
    b.plot(TS, g["POPC, five runs, open"], color="#1baf7a", ls=ls, lw=1.1)
    b.plot(TS, g["POPC, five runs, bv_env"], color="#e87ba4", ls=ls, lw=1.1)
    b.plot(TS, g["POPC, five runs vs four peptides, open"], color="#184f95", ls=ls, lw=1.1)
a.plot(TS, f.observed_free, "o-", color="black", lw=1.3, ms=3); b.plot(TS, g.observed_gain, "o-", color="black", lw=1.3, ms=3)
for ax in (a, b): ax.set_xscale("log"); ax.set_xticks(TS); ax.set_xticklabels([3, 10, 30, 60, 120]); ax.set_xlabel("Labelling time (s)")
a.set_ylim(0, 1.02); a.set_ylabel("Mean peptide RFU, free melittin"); b.axhline(0, color="#999999", lw=0.5); b.set_ylabel("Protected amides gained on POPC")
b.set_ylim(-3, 9)
S.panel_letter(a, "a"); S.panel_letter(b, "b")
h = [Line2D([], [], color="black", marker="o", ms=3, lw=1.3, label="HDX-MS"), Line2D([], [], color="#999999", lw=1, label="intrinsic, P = 1"),
     Line2D([], [], color="#4d4d4d", lw=1, label="MD water, open model"), Line2D([], [], color="#184f95", lw=1, label="MD four peptides as free state"),
     Line2D([], [], color="#1baf7a", lw=1, label="MD POPC, open model"), Line2D([], [], color="#e87ba4", lw=1, label="MD POPC, BV lipid counted"),
     Line2D([], [], color="black", ls="-", lw=0.8, label="k_ch as calculated"), Line2D([], [], color="black", ls="--", lw=0.8, label="k_ch × 2")]
fig.legend(handles=h, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 0.97), fontsize=S.PT_FOOT)
S.save(fig, os.path.join(O, "Fig_H3_free_state_and_gain"))
print("ok")

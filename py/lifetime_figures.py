import os, sys, json, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import amp_style as S; from amp_config import ANALYSIS_DIR, FIG_DIR, PEPTIDES
from kint import kch
S.apply(); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
A = os.path.join(ANALYSIS_DIR, "hdx_compare"); O = os.path.join(FIG_DIR, "chapter4", "hdx")
ST = {"water": ("#4d4d4d", ":"), "POPC, five runs": ("#2a78d6", "-"), "POPE:POPG 7:3": ("#eb6834", "--"), "four peptides": ("#184f95", "-.")}
SEGS = {"mel": [("2 to 5", 2, 5), ("6 to 13", 6, 13), ("15 to 19", 15, 19), ("20 to 26", 20, 26)],
        "ll37": [("2 to 13", 2, 13), ("17 to 31", 17, 31), ("32 to 37", 32, 37)]}
SEGC = ["#1baf7a", "#eda100", "#2a78d6", "#e87ba4"]
for pep in sys.argv[1:]:
    seq = PEPTIDES[pep]["seq"]; n = len(seq); x = np.arange(1, n + 1); d = pd.read_csv(os.path.join(A, f"lifetimes_{pep}.csv"))
    dw = np.load(os.path.join(A, f"dwell_{pep}.npz"))
    k = kch(seq); kmax, kmin = np.nanmax(k), np.nanmin(k)
    fig = S.figure(height_mm=150); gs = fig.add_gridspec(3, 2, hspace=0.75, wspace=0.3, height_ratios=[1, 1, 1.1])
    a = fig.add_subplot(gs[0, :]); b = fig.add_subplot(gs[1, :]); c = fig.add_subplot(gs[2, 0]); e = fig.add_subplot(gs[2, 1])
    for ax, col, lab in ((a, "tau_closed_ns", "Mean closed dwell (ns)"), (b, "longest_closed_ns", "Longest closed\nepisode (ns)")):
        S.residue_axis(ax, seq, PEPTIDES[pep]["segments"], label_segments=(ax is a))
        for st, (cc, ls) in ST.items():
            y = d[d.state == st].set_index("resid")[col].reindex(x).values.astype(float); y[y <= 0] = np.nan
            ax.plot(x, y, color=cc, ls=ls, lw=1, marker="o", ms=1.6)
        ax.set_yscale("log"); ax.set_ylabel(lab)
        tk = [0.2, 0.5, 1, 2, 5] if ax is a else [1, 3, 10, 30, 100]
        ax.set_yticks(tk); ax.set_yticklabels([str(t) for t in tk]); ax.set_ylim(tk[0] * 0.8, tk[-1] * 1.25)
    a.set_title(PEPTIDES[pep]["name"], loc="left", fontsize=S.PT_NAME, pad=12)
    # dwell time distributions of closed episodes on POPC, per segment
    for j, (nm, lo, hi) in enumerate(SEGS[pep]):
        v = np.concatenate([dw[f"POPC, five runs|cl|{r}"] for r in range(lo, hi + 1) if f"POPC, five runs|cl|{r}" in dw.files])
        v = np.sort(v); ccdf = 1 - np.arange(len(v)) / len(v)
        c.step(v, ccdf, where="post", color=SEGC[j], lw=1, label=f"residues {nm}")
    c.set_xscale("log"); c.set_yscale("log"); c.set_xlabel("Closed episode length (ns)"); c.set_ylabel("Fraction of episodes ≥ length")
    c.legend(fontsize=S.PT_FOOT, loc="lower left"); c.set_title("POPC, five runs", fontsize=S.PT_ANNOT)
    # timescale separation: open dwell versus 1 / k_ch
    p = d[d.state == "POPC, five runs"].set_index("resid")
    tau_o = p.tau_open_ns.reindex(x).values * 1e-9; inv = 1 / k
    e.scatter(tau_o, inv, s=10, color="#2a78d6", edgecolor="black", lw=0.3, zorder=3)
    e.set_xscale("log"); e.set_yscale("log"); e.set_xlim(1e-10, 1e3); e.set_ylim(1e-10, 1e3)
    e.plot([1e-10, 1e3], [1e-10, 1e3], color="#999999", lw=0.6, ls="--"); e.text(3e-10, 3e-9, "τ_open = 1/k_ch", rotation=37, fontsize=S.PT_FOOT, color="#777777")
    e.set_xlabel("Mean open dwell, MD (s)"); e.set_ylabel("1 / k_ch at 278 K, pD 6 (s)"); e.set_title("EX2 check, POPC", fontsize=S.PT_ANNOT)
    for ax, l in zip((a, b, c, e), "abcd"): S.panel_letter(ax, l)
    fig.legend(handles=[Line2D([], [], color=cc, ls=ls, lw=1.2, label=s) for s, (cc, ls) in ST.items()], ncol=4, loc="lower center", bbox_to_anchor=(0.5, 0.965), fontsize=S.PT_ANNOT)
    S.save(fig, os.path.join(O, f"Fig_H4_hbond_lifetimes_{pep}"))
    print(pep, "1/k_ch range s:", round(1 / kmax, 2), round(1 / kmin, 1), "tau_open POPC range ns:", round(np.nanmin(tau_o) * 1e9, 2), round(np.nanmax(tau_o) * 1e9, 2))

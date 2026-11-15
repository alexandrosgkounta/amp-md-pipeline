"""SI analysis for Chapter 4: is the 500 to 1000 ns averaging window equilibrated?
For every run (0 to 1 us): 50 ns block means of helicity and COM height, the equilibration time t_eq of
Chodera (2016, JCTC 12:1799; start that maximises the number of uncorrelated samples), and the change in the
mean when the window starts at 250, 500 or 750 ns. Writes analysis/equilibration.csv and figures/chapter4/SI/Fig_S4_1_equilibration.*"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import amp_style as S; from amp_config import ANALYSIS_DIR, FIG_DIR
S.apply(); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT = os.path.join(FIG_DIR, "chapter4", "SI"); os.makedirs(OUT, exist_ok=True)
PEPS = {"mel": "Melittin", "ll37": "LL-37"}
RUNS = [("sol", "water", "#4d4d4d", ":"), ("popc", "POPC reference", "#2a78d6", "-"), ("popc_r2", "POPC opposite face", "#eb6834", "-"),
        ("popc_r3", "POPC tilted +", "#1baf7a", "-"), ("popc_r4", "POPC tilted −", "#eda100", "-"), ("popc_r5", "POPC reversed", "#e87ba4", "-"),
        ("pcpepg", "POPC:POPE:POPG 7:2:1", "#1baf7a", "--"), ("pepg", "POPE:POPG 7:3", "#eb6834", "--"), ("agg", "4 peptides, POPC (mean)", "#184f95", "-.")]

def g_ineff(x):
    x = np.asarray(x, float); n = len(x); x = x - x.mean(); v = x.var()
    if v == 0 or n < 4: return 1.0
    f = np.fft.rfft(x, 2 * n); ac = np.fft.irfft(f * np.conj(f))[:n] / (v * np.arange(n, 0, -1))
    g = 1.0
    for t in range(1, n - 1):
        if ac[t] <= 0: break
        g += 2 * ac[t] * (1 - t / n)
    return max(g, 1.0)
def t_equil(t, x, step_ns=10.0):
    best = (0, -1)
    for t0 in np.arange(0, t[-1] * 0.8, step_ns):
        m = t >= t0; neff = m.sum() / g_ineff(x[m])
        if neff > best[1]: best = (t0, neff)
    return best
rows = []
fig = S.figure(height_mm=120); gs = fig.add_gridspec(2, 2, hspace=0.6, wspace=0.3)
for i, (pep, pname) in enumerate(PEPS.items()):
    a = fig.add_subplot(gs[i, 0]); b = fig.add_subplot(gs[i, 1])
    for key, lab, col, ls in RUNS:
        sid = f"{pep}_{key}"; f = os.path.join(ANALYSIS_DIR, sid, "timeseries.csv")
        if not os.path.exists(f): continue
        ts = pd.read_csv(f); ts = ts[ts.time_ns <= 1000.0]
        g = ts.groupby("time_ns").mean(numeric_only=True); t = g.index.values
        # 1 ns sampling for the correlation analysis
        sel = np.isclose(t % 1.0, 0) | (np.diff(np.r_[t, t[-1] + 1]) >= 1)
        obs = {"helicity": g.helicity.values}
        if "height_com" in g and key not in ("agg",): obs["height_com"] = g.height_com.values
        for name, x in obs.items():
            t1, x1 = t[sel], x[sel]
            teq, neff = t_equil(t1, x1)
            w = lambda lo: x1[(t1 >= lo)].mean()
            bl = [x1[(t1 >= lo) & (t1 < lo + 100)].mean() for lo in range(500, 1000, 100)]
            rows.append(dict(system=sid, observable=name, t_eq_ns=round(float(teq)), n_eff_from_t_eq=round(float(neff), 1),
                             mean_250_1000=round(w(250), 3), mean_500_1000=round(w(500), 3), mean_750_1000=round(w(750), 3),
                             sd_100ns_blocks_500_1000=round(float(np.std(bl, ddof=1)), 3),
                             drift_per_100ns_500_1000=round(float(np.polyfit(t1[t1 >= 500], x1[t1 >= 500], 1)[0] * 100), 3),
                             window_after_t_eq=bool(teq <= 500)))
        edges = np.arange(0, 1001, 50); c = (edges[:-1] + edges[1:]) / 2 / 1000
        hb = [g.helicity.values[(t >= lo) & (t < lo + 50)].mean() for lo in edges[:-1]]
        a.plot(c, hb, color=col, ls=ls, lw=0.9, marker="o", ms=1.6)
        if "height_com" in obs:
            zb = [g.height_com.values[(t >= lo) & (t < lo + 50)].mean() for lo in edges[:-1]]
            b.plot(c, zb, color=col, ls=ls, lw=0.9, marker="o", ms=1.6)
    for ax in (a, b):
        ax.axvspan(0.5, 1.0, color="#f2f2f2", zorder=0, lw=0); ax.set_xlim(0, 1); ax.set_xlabel("Time (µs)")
    a.set_ylim(0, 1); a.set_ylabel("Helicity (50 ns blocks)"); a.set_title(pname, loc="left", fontsize=S.PT_NAME)
    b.axhline(0, color="#555555", lw=0.5, ls="--"); b.set_ylabel("COM height above P plane (nm)")
    S.panel_letter(a, "ac"[i]); S.panel_letter(b, "bd"[i])
h = [Line2D([], [], color=c, ls=l, lw=1, label=lab) for _, lab, c, l in RUNS]
fig.legend(handles=h, ncol=5, loc="lower center", bbox_to_anchor=(0.5, 0.965), fontsize=S.PT_FOOT)
S.save(fig, os.path.join(OUT, "Fig_S4_1_equilibration"))
df = pd.DataFrame(rows); df.to_csv(os.path.join(ANALYSIS_DIR, "equilibration.csv"), index=False)
pd.set_option("display.width", 250); print(df.to_string(index=False))

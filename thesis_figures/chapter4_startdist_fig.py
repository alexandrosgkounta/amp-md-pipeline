import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import amp_style as S; from amp_config import ANALYSIS_DIR, FIG_DIR
S.apply(); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT = os.path.join(FIG_DIR, "chapter4", "SI")
RUNS = [("popc", "POPC reference", "#2a78d6", "-"), ("popc_r2", "POPC opposite face", "#eb6834", "-"), ("popc_r3", "POPC tilted +", "#1baf7a", "-"),
        ("popc_r4", "POPC tilted −", "#eda100", "-"), ("popc_r5", "POPC reversed", "#e87ba4", "-"),
        ("pcpepg", "POPC:POPE:POPG 7:2:1", "#1baf7a", "--"), ("pepg", "POPE:POPG 7:3", "#eb6834", "--")]
rows = []
fig = S.figure(height_mm=110); gs = fig.add_gridspec(2, 2, hspace=0.6, wspace=0.3)
for i, (pep, name) in enumerate((("mel", "Melittin"), ("ll37", "LL-37"))):
    a = fig.add_subplot(gs[i, 0]); b = fig.add_subplot(gs[i, 1])
    for key, lab, col, ls in RUNS:
        sid = f"{pep}_{key}"; df = pd.read_csv(os.path.join(ANALYSIS_DIR, sid, "start_distance.csv"))
        an = key in ("pcpepg", "pepg"); lw = 1.3 if an else 0.6; al = 1.0 if an else 0.55; z = 3 if an else 2
        a.plot(df.time_ns, df.min_dist_nm, color=col, ls=ls, lw=lw, alpha=al, zorder=z); b.plot(df.time_ns, df.height_com, color=col, ls=ls, lw=lw, alpha=al, zorder=z)
        fc = df.time_ns[df.min_dist_nm <= 0.4]
        rows.append(dict(system=sid, min_dist_t0_nm=round(df.min_dist_nm.iloc[0], 2), com_height_t0_nm=round(df.height_com.iloc[0], 2),
                         first_contact_ns=float(fc.iloc[0]) if len(fc) else np.nan,
                         frac_contact_0_50ns=round(float((df.min_dist_nm <= 0.4).mean()), 2), com_height_50ns_nm=round(df.height_com.iloc[-1], 2)))
    a.axhline(0.4, color="#555555", lw=0.5, ls="--"); a.text(49, 2.45, "dashed line: contact (0.4 nm)", ha="right", va="top", fontsize=S.PT_FOOT, color="#555555", bbox=dict(fc="white", ec="none", pad=1))
    a.set_ylim(0, 2.5); a.set_xlim(0, 50); a.set_xlabel("Time (ns)"); a.set_ylabel("Minimum distance to lipid (nm)")
    a.set_title(name, loc="left", fontsize=S.PT_NAME)
    b.axhline(0, color="#555555", lw=0.5, ls="--"); b.set_xlim(0, 50); b.set_xlabel("Time (ns)"); b.set_ylabel("COM height above P plane (nm)")
    S.panel_letter(a, "ac"[i]); S.panel_letter(b, "bd"[i])
h = [Line2D([], [], color=c, ls=l, lw=1, label=lab) for _, lab, c, l in RUNS]
fig.legend(handles=h, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 0.965), fontsize=S.PT_FOOT)
S.save(fig, os.path.join(OUT, "Fig_S4_2_start_distance"))
d = pd.DataFrame(rows); d.to_csv(os.path.join(ANALYSIS_DIR, "start_distance_summary.csv"), index=False); print(d.to_string(index=False))

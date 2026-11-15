"""Longest helix and lipid runs include runs cut by the window edges (lower bounds).
Lifetimes of backbone amide hydrogen bonds (and of helix and lipid contact) per residue, 500 to 1000 ns.
Closed = the amide H makes a hydrogen bond (H...O <= 2.4 A) to any peptide or lipid oxygen (nh_env > 0 in perres.npz);
open = no such bond. Frames every 0.2 ns, so events shorter than 0.2 ns are not resolved: dwell times are those of
the 0.2 ns coarse-grained state. Per residue: fraction closed, mean closed and open dwell (run lengths, runs cut by
the window edges excluded), longest closed run, intermittent autocorrelation time of the closed state, and
P_kin = 1 + tau_closed / tau_open. Writes analysis/hdx_compare/lifetimes_<pep>.csv and dwell_<pep>.npz"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from amp_config import ANALYSIS_DIR, PEPTIDES
OUT = os.path.join(ANALYSIS_DIR, "hdx_compare"); W0, W1 = 500.0, 1000.0
def runs(b):
    """lengths of True runs not touching either end"""
    b = np.asarray(b, bool); d = np.diff(np.r_[0, b.astype(int), 0]); s = np.where(d == 1)[0]; e = np.where(d == -1)[0]
    keep = (s > 0) & (e < len(b)); return (e - s)[keep], (e - s)
def acf_time(b, dt, maxlag=500):
    b = b.astype(float); m = b.mean()
    if m in (0, 1): return np.nan
    x = b - m; n = len(x); f = np.fft.rfft(x, 2 * n); c = np.fft.irfft(f * np.conj(f))[:min(maxlag, n)] / (np.arange(n, n - min(maxlag, n), -1) * x.var())
    z = np.where(c <= 0)[0]; k = z[0] if len(z) else len(c)
    return float(np.trapz(c[:k], dx=dt))
STATES = {"mel": {"water": ["mel_sol"], "POPC, five runs": ["mel_popc", "mel_popc_r2", "mel_popc_r3", "mel_popc_r4", "mel_popc_r5"],
                  "POPC:POPE:POPG 7:2:1": ["mel_pcpepg"], "POPE:POPG 7:3": ["mel_pepg"], "four peptides": ["mel_agg"]},
          "ll37": {"water": ["ll37_sol"], "POPC, five runs": ["ll37_popc", "ll37_popc_r2", "ll37_popc_r3", "ll37_popc_r4", "ll37_popc_r5"],
                   "POPC:POPE:POPG 7:2:1": ["ll37_pcpepg"], "POPE:POPG 7:3": ["ll37_pepg"], "four peptides": ["ll37_agg"]}}
for pep in sys.argv[1:]:
    n = len(PEPTIDES[pep]["seq"]); rows = []; dw = {}
    for st, sids in STATES[pep].items():
        acc = {k: [[] for _ in range(n)] for k in ("cl", "op", "hx", "lp", "allcl", "allhx", "alllp")}; frac = np.zeros(n); fh = np.zeros(n); fl = np.zeros(n); tac = [[] for _ in range(n)]; nser = 0
        for sid in sids:
            z = np.load(os.path.join(ANALYSIS_DIR, sid, "perres.npz"), allow_pickle=True); t = z["time_ns"]; m = (t >= W0) & (t <= W1)
            dt = float(np.median(np.diff(t)))
            C = z["nh_env"][m] > 0; H = z["dssp"][m] == "H"; L = z["contact"][m][..., 0]
            for k in range(C.shape[1]):
                nser += 1
                for i in range(n):
                    for key, arr in (("cl", C[:, k, i]), ("op", ~C[:, k, i]), ("hx", H[:, k, i]), ("lp", L[:, k, i])):
                        r, rall = runs(arr); acc[key][i].extend((r * dt).tolist())
                        if key in ("cl", "hx", "lp"): acc["all" + key][i].extend((rall * dt).tolist())
                    frac[i] += C[:, k, i].mean(); fh[i] += H[:, k, i].mean(); fl[i] += L[:, k, i].mean()
                    tac[i].append(acf_time(C[:, k, i], dt))
        for i in range(n):
            mc = np.mean(acc["cl"][i]) if acc["cl"][i] else np.nan; mo = np.mean(acc["op"][i]) if acc["op"][i] else np.nan
            rows.append(dict(state=st, resid=i + 1, aa=PEPTIDES[pep]["seq"][i], frac_closed=frac[i] / nser, helix=fh[i] / nser, lipid_contact=fl[i] / nser,
                             tau_closed_ns=mc, tau_open_ns=mo, n_closed_events=len(acc["cl"][i]),
                             longest_closed_ns=max(acc["allcl"][i]) if acc["allcl"][i] else 0.0,
                             tau_acf_closed_ns=np.nanmean(tac[i]) if np.isfinite(tac[i]).any() else np.nan,
                             P_kin=1 + mc / mo if np.isfinite(mc) and np.isfinite(mo) else np.nan,
                             tau_helix_ns=np.mean(acc["hx"][i]) if acc["hx"][i] else np.nan, longest_helix_ns=max(acc["allhx"][i]) if acc["allhx"][i] else 0.0,
                             tau_lipid_ns=np.mean(acc["lp"][i]) if acc["lp"][i] else np.nan, longest_lipid_ns=max(acc["alllp"][i]) if acc["alllp"][i] else 0.0))
            dw[f"{st}|cl|{i + 1}"] = np.array(acc["allcl"][i], np.float32)
    df = pd.DataFrame(rows); df.to_csv(os.path.join(OUT, f"lifetimes_{pep}.csv"), index=False); np.savez_compressed(os.path.join(OUT, f"dwell_{pep}.npz"), **dw)
    pd.set_option("display.width", 250)
    print(pep); print(df[df.state.isin(["water", "POPC, five runs"])][["state", "resid", "aa", "frac_closed", "tau_closed_ns", "tau_open_ns", "longest_closed_ns", "tau_acf_closed_ns", "P_kin", "tau_helix_ns", "longest_helix_ns", "tau_lipid_ns", "longest_lipid_ns"]].round(2).to_string(index=False))

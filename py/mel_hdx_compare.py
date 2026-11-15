"""Melittin: MD ensembles set against the HDX-MS data (free and with POPC liposomes, 278 K, pD 6.0, 3 to 120 s).
Per-residue descriptors over 500 to 1000 ns and two forward models of deuterium uptake:
  open model : k_obs = k_ch * f_open, f_open = fraction of frames in which the amide H makes no hydrogen bond
               (H...O <= 2.4 A) to peptide or lipid oxygen (EX2, the closed state does not exchange)
  BV model   : ln P = 0.35 Nc + 2.0 Nh (Best and Vendruscolo 2006), averaged over frames; peptide-only and
               environment-inclusive (lipid heavy atoms and oxygens counted) variants
D(t) = 1 - exp(-k_ch t / P). Observed uptake is read from data/hdx/melittin_observed.csv (see data/hdx/README.md for the format);
it is not distributed with this repository. Writes analysis/hdx_compare/*.csv, figures/chapter4/hdx/*"""
import os, sys, json, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import amp_style as S; from amp_config import ANALYSIS_DIR, FIG_DIR, PEPTIDES
from kint import kch
S.apply(); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT = os.path.join(ANALYSIS_DIR, "hdx_compare"); FIGO = os.path.join(FIG_DIR, "chapter4", "hdx"); os.makedirs(FIGO, exist_ok=True)
SEQ = PEPTIDES["mel"]["seq"]; n = len(SEQ); x = np.arange(1, n + 1)
K0 = kch(SEQ); SCALE = float(os.environ.get("KSCALE", 1)); K = K0 * SCALE; TS = np.array([3, 10, 30, 60, 120.0])
SEGS = [("2 to 5", 2, 5), ("6 to 13", 6, 13), ("14 to 16", 14, 16), ("17 to 19", 17, 19), ("20 to 26", 20, 26)]
from hdx_observed import load as _load_obs
OBS, OBS_FREE, OBS_GAIN = _load_obs("melittin", segments=[s[0] for s in SEGS], times=list(TS))
PEPS = [(1,13),(1,16),(4,16),(5,15),(5,16),(6,13),(6,16),(7,13),(7,16),(7,17),(7,18),(7,19),(7,26),(8,16),(9,16),(10,16),(14,26),(17,24),(17,26),(18,25),(19,26),(20,26)]
REP = ~np.isnan(K)                     # 24 reportable amides (residue 1 and Pro14 excluded)
W0, W1 = 500.0, 1000.0

def load(sid, bound_only=False):
    pr = dict(np.load(os.path.join(ANALYSIS_DIR, sid, "perres.npz"), allow_pickle=True)); t = pr["time_ns"]
    m = (t >= W0) & (t <= W1)
    sel = {k: pr[k][m] for k in ("dssp", "depth", "contact", "nc_pep", "nh_pep", "nc_env", "nh_env")}
    if bound_only:
        b = sel["contact"][..., 0].any(axis=2)            # (frames, copies) any residue touching lipid
        sel = {k: v[b] for k, v in sel.items()}            # -> (frames*, residues) flattened over copies
    else:
        sel = {k: v.reshape((-1,) + v.shape[2:]) for k, v in sel.items()}
    return sel
def descriptors(parts):
    """parts: list of dicts from load(); pooled with equal weight per frame"""
    cat = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    hb_pep = cat["nh_pep"] > 0; hb_env = cat["nh_env"] > 0; hb_lip = hb_env & ~hb_pep
    lnP_pep = 0.35 * cat["nc_pep"] + 2.0 * cat["nh_pep"]; lnP_env = 0.35 * cat["nc_env"] + 2.0 * cat["nh_env"]
    d = dict(helix=(cat["dssp"] == "H").mean(0), hb_pep=hb_pep.mean(0), hb_lipid=hb_lip.mean(0), lipid_contact=cat["contact"][..., 0].mean(0),
             depth_nm=np.nanmean(cat["depth"], 0) / 10, f_open=1 - hb_env.mean(0), lnP_bv_pep=lnP_pep.mean(0), lnP_bv_env=lnP_env.mean(0),
             frames=len(cat["dssp"]))
    return d
def uptake(P):
    """residue D(t) for the reportable amides, shape (24..., times)"""
    return 1 - np.exp(-np.outer(K / P, TS))
def seg_mean(D):
    return np.array([[np.nanmean(D[[i for i in range(lo - 1, hi) if REP[i]], j]) for j in range(len(TS))] for _, lo, hi in SEGS])
def pep_mean(D):
    return np.array([[np.nanmean(D[[i for i in range(a, b) if REP[i]], j]) for j in range(len(TS))] for a, b in PEPS])
def models(d):
    Popen = 1 / np.clip(d["f_open"], 1e-3, 1)
    return {"open": uptake(Popen), "bv_pep": uptake(np.exp(d["lnP_bv_pep"])), "bv_env": uptake(np.exp(d["lnP_bv_env"]))}

STATES = {
    "water": [("mel_sol", False)],
    "POPC, five runs": [(f"mel_popc{s}", False) for s in ("", "_r2", "_r3", "_r4", "_r5")],
    "POPC, bound frames": [(f"mel_popc{s}", True) for s in ("", "_r2", "_r3", "_r4", "_r5")],
    "POPE:POPG 7:3": [("mel_pepg", False)],
    "four peptides": [("mel_agg", False)],
}
D = {k: descriptors([load(s, b) for s, b in v]) for k, v in STATES.items()}
# water hydration (wet trajectories); water state has no wet.xtc: proxy 1 - hb_pep (stated in the text)
def hyd(sids):
    v = [pd.read_csv(os.path.join(OUT, f"hydration_{s}.csv")).groupby("resid").water_hbond_frac.mean() for s in sids]
    h = pd.concat(v, axis=1).mean(1); out = np.full(n, np.nan); out[h.index.values - 1] = h.values; return out
D["POPC, five runs"]["water_hb"] = hyd(["mel_popc", "mel_popc_r2", "mel_popc_r3", "mel_popc_r4", "mel_popc_r5"])
D["POPE:POPG 7:3"]["water_hb"] = hyd(["mel_pepg"]); D["four peptides"]["water_hb"] = hyd(["mel_agg"])
rows = []
for st, d in D.items():
    for i in range(n):
        rows.append(dict(state=st, resid=i + 1, aa=SEQ[i], k_ch_s=K[i], **{k: (float(v[i]) if hasattr(v, "__len__") else v) for k, v in d.items()}))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "mel_md_descriptors.csv"), index=False)

M = {st: models(d) for st, d in D.items()}
PRED = {}
for st in ("POPC, five runs", "POPC, bound frames", "POPE:POPG 7:3"):
    for mod in ("open", "bv_pep", "bv_env"):
        PRED[(st, mod)] = seg_mean(M[st][mod]) - seg_mean(M["water"][mod])
# alternative free state: four-peptide assembly
for mod in ("open", "bv_pep", "bv_env"): PRED[("POPC, five runs vs four peptides", mod)] = seg_mean(M["POPC, five runs"][mod]) - seg_mean(M["four peptides"][mod])
srows = []
for (st, mod), v in PRED.items():
    for si, (nm, _, _) in enumerate(SEGS):
        for j, t in enumerate(TS):
            srows.append(dict(state=st, model=mod, segment=nm, t_s=t, pred_dRFU=round(float(v[si, j]), 3), obs_dRFU=OBS[si, j]))
sd = pd.DataFrame(srows); sd.to_csv(os.path.join(OUT, f"mel_segment_pred_vs_obs_k{SCALE:g}.csv"), index=False)
# summary statistics: Pearson r over 25 values, rank of segments at 30 s, ratio of magnitudes
summ = []
for (st, mod), v in PRED.items():
    r = np.corrcoef(v.ravel(), OBS.ravel())[0, 1]
    rk_obs = np.argsort(OBS[:, 2]); rk_pred = np.argsort(v[:, 2])
    from scipy.stats import spearmanr
    rho = spearmanr(v[:, 2], OBS[:, 2]).correlation
    summ.append(dict(state=st, model=mod, pearson_r_25=round(r, 2), spearman_segments_30s=round(rho, 2),
                     most_protected_pred_30s=SEGS[rk_pred[0]][0], least_protected_pred_30s=SEGS[rk_pred[-1]][0],
                     mean_pred_dRFU=round(float(v.mean()), 3), mean_obs_dRFU=round(float(OBS.mean()), 3)))
summ = pd.DataFrame(summ); summ.to_csv(os.path.join(OUT, f"mel_model_summary_k{SCALE:g}.csv"), index=False)
# free state: mean peptide RFU, observed vs ceiling vs water predictions
free = {"intrinsic ceiling (P = 1)": pep_mean(uptake(np.ones(n))).mean(0)}
for mod in ("open", "bv_pep"): free[f"MD water, {mod}"] = pep_mean(M["water"][mod]).mean(0)
free["MD four peptides, open"] = pep_mean(M["four peptides"]["open"]).mean(0)
pd.DataFrame(dict(t_s=TS, observed_free=OBS_FREE, **{k: np.round(v, 3) for k, v in free.items()})).to_csv(os.path.join(OUT, f"mel_free_state_k{SCALE:g}.csv"), index=False)
# protected amide gain (24 amides): POPC minus free
def prot_count(Dm): return np.nansum(1 - Dm[REP], axis=0)
gain = {}
for st in ("POPC, five runs", "POPC, bound frames", "POPE:POPG 7:3"):
    for mod in ("open", "bv_pep", "bv_env"):
        gain[f"{st}, {mod}"] = prot_count(M[st][mod]) - prot_count(M["water"][mod])
gain["POPC, five runs vs four peptides, open"] = prot_count(M["POPC, five runs"]["open"]) - prot_count(M["four peptides"]["open"])
pd.DataFrame(dict(t_s=TS, observed_gain=OBS_GAIN, **{k: np.round(v, 2) for k, v in gain.items()})).to_csv(os.path.join(OUT, f"mel_gain_k{SCALE:g}.csv"), index=False)
pd.set_option("display.width", 250)
print(summ.to_string(index=False)); print(pd.read_csv(os.path.join(OUT, f"mel_free_state_k{SCALE:g}.csv")).to_string(index=False))
print(pd.read_csv(os.path.join(OUT, f"mel_gain_k{SCALE:g}.csv")).to_string(index=False))
for st in ("water", "POPC, five runs", "POPC, bound frames", "four peptides"):
    d = D[st]; print(st, d["frames"], "f_open", np.round(d["f_open"], 2))
json.dump({"k_ch": [None if np.isnan(v) else float(v) for v in K]}, open(os.path.join(OUT, "kch.json"), "w"))

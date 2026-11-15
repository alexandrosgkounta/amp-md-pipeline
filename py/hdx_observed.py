"""Loader for experimental HDX-MS observations.

The comparison scripts need measured deuterium uptake to score the forward models
against. Those measurements are kept out of the repository so that unpublished data
is not distributed with the code.

Expected file: data/hdx/<peptide>_observed.csv, with the columns

    quantity,segment,t_s,value

  quantity : one of  seg_dRFU | free_uptake | protection_gain
  segment  : segment label for seg_dRFU rows (e.g. "6 to 13"); blank otherwise
  t_s      : labelling time in seconds
  value    : the measured number

Override the location with the AMP_HDX_DATA environment variable. See
data/hdx/README.md for the format and data/hdx/melittin_observed.example.csv
for a worked skeleton.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DIR = os.environ.get(
    "AMP_HDX_DATA", os.path.abspath(os.path.join(HERE, "..", "data", "hdx"))
)


def _path(peptide):
    return os.path.join(DEFAULT_DIR, f"{peptide}_observed.csv")


def available(peptide="melittin"):
    """True if measured data for this peptide is present."""
    return os.path.exists(_path(peptide))


def load(peptide="melittin", segments=None, times=None):
    """Return (OBS, OBS_FREE, OBS_GAIN).

    OBS       : (n_segments, n_times) segment-level delta RFU, indexed
                [segment, time] to match the comparison script's OBS[si, j]
    OBS_FREE  : (n_times,) free-state fractional uptake
    OBS_GAIN  : (n_times,) protection gain

    Raises FileNotFoundError with the expected path if the file is absent, so that
    a fresh clone fails loudly rather than silently comparing against nothing.
    """
    p = _path(peptide)
    if not os.path.exists(p):
        raise FileNotFoundError(
            f"No experimental HDX data at {p}.\n"
            "This repository ships the forward models but not the measurements. "
            "Add your own file in the documented format, or copy "
            "data/hdx/melittin_observed.example.csv and replace the placeholder values. "
            "See data/hdx/README.md."
        )
    df = pd.read_csv(p, comment="#")
    df.columns = [c.strip() for c in df.columns]

    seg = df[df["quantity"] == "seg_dRFU"]
    if segments is None:
        segments = list(dict.fromkeys(seg["segment"]))
    if times is None:
        times = sorted(seg["t_s"].unique())

    obs = np.full((len(segments), len(times)), np.nan)
    for si, sname in enumerate(segments):
        for j, t in enumerate(times):
            m = seg[(seg["segment"] == sname) & (seg["t_s"] == t)]
            if len(m):
                obs[si, j] = float(m["value"].iloc[0])
    if np.isnan(obs).any():
        missing = [(segments[a], times[b]) for a, b in zip(*np.where(np.isnan(obs)))]
        raise ValueError(f"Missing seg_dRFU rows for {missing} in {p}")

    def series(q):
        sub = df[df["quantity"] == q].sort_values("t_s")
        return sub["value"].to_numpy(dtype=float)

    return obs, series("free_uptake"), series("protection_gain")

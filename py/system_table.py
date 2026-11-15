#!/usr/bin/env python3
"""Build the system table (composition + run parameters) for the Methods section.

Composition comes from the tpr. Run parameters come from settings/key_params.txt
and settings/xtc_check.txt written by 00_register_settings.sh (if present).
Output: analysis/system_table.tsv and analysis/system_table.md
"""
import os, re, sys, warnings, collections
warnings.filterwarnings("ignore")
import numpy as np, MDAnalysis as mda
from amp_config import read_systems, LIPIDS, WATER, IONS, ANALYSIS_DIR, PEPTIDES


def parse_params(d):
    out = {}
    f = os.path.join(d, "settings", "key_params.txt")
    if os.path.exists(f):
        for line in open(f):
            m = re.match(r"^([\w-]+)\s+=\s+(.*)$", line.strip())
            if m: out[m.group(1)] = m.group(2).strip()
            m = re.match(r"^(ref-t|tau-t):\s+(.*)$", line.strip())
            if m: out[m.group(1)] = " ".join(sorted(set(m.group(2).split())))
    f = os.path.join(d, "settings", "xtc_check.txt")
    if os.path.exists(f):
        txt = open(f).read()
        m = re.findall(r"Last frame\s+(\d+)\s+time\s+([\d.]+)", txt)
        if m: out["xtc_frames"] = int(m[-1][0]) + 1; out["xtc_last_ps"] = float(m[-1][1])
        m = re.search(r"^Coords\s+(\d+)\s+([\d.]+)", txt, re.M)
        if m: out["xtc_dt_ps"] = float(m.group(2))
    return out


def box_from_dump(d):
    f = os.path.join(d, "settings", "tpr_dump.txt")
    if not os.path.exists(f): return None
    rows = []
    with open(f) as fh:
        for line in fh:
            if rows or re.match(r"^\s*box \(3x3\):", line):
                if rows == [] and "box (3x3)" in line: rows.append(None); continue
                m = re.findall(r"[-+]?\d+\.\d+e[-+]\d+", line)
                rows.append(float(m[len(rows) - 1]) * 10.0)      # diagonal, nm -> A
                if len(rows) == 4: return np.array(rows[1:])
    return None


def forcefield(d):
    f = os.path.join(d, "topol.top")
    if not os.path.exists(f): return "check (no topol.top)"
    for line in open(f):
        if line.startswith("#include") and ("forcefield" in line):
            return line.split('"')[1]
    return "check"


def composition(s):
    u = mda.Universe(os.path.join(s["path"], s["tpr"]))
    c = collections.Counter(u.residues.resnames)
    prot = u.select_atoms("protein")
    row = dict(sys_id=s["sys_id"], peptide=PEPTIDES[s["peptide"]]["name"], condition=s["condition"],
               replicate=s["replicate"], n_peptides=s["n_pep"], atoms=u.atoms.n_atoms,
               peptide_charge=round(float(prot.total_charge()) / s["n_pep"], 2),
               net_charge=round(float(u.atoms.total_charge()), 2))
    dims = u.dimensions
    if dims is None:
        # Some MDAnalysis versions return no box for a tpr alone. Take it from the gmx dump written
        # by step 0 instead (opening the raw xtc would index the whole file, minutes per system).
        dims = box_from_dump(s["path"])
    row["box_nm"] = "x".join(f"{b / 10.0:.2f}" for b in dims[:3]) if dims is not None else "n/a"
    nlip = sum(c[l] for l in LIPIDS)
    row["lipids"] = nlip
    P = u.select_atoms("name P and resname " + " ".join(LIPIDS))
    if nlip:
        z = P.positions[:, 2]; zc = np.mean(z)
        up, lo = P[z > zc], P[z <= zc]
        cu, cl = collections.Counter(up.resnames), collections.Counter(lo.resnames)
        row["upper_leaflet"] = " ".join(f"{l}:{cu[l]}" for l in LIPIDS if cu[l])
        row["lower_leaflet"] = " ".join(f"{l}:{cl[l]}" for l in LIPIDS if cl[l])
        tot = {l: c[l] for l in LIPIDS if c[l]}
        row["lipid_mol_percent"] = " ".join(f"{l}:{100*v/nlip:.0f}" for l, v in tot.items())
        row["lipid_charge"] = -c["POPG"]
        row["peptide_to_lipid"] = f"1:{nlip/s['n_pep']:.0f}"
    nw = sum(c[w] for w in WATER)
    row["waters"] = nw
    ions = {i: c[i] for i in IONS if c[i]}
    row["ions"] = " ".join(f"{k}:{v}" for k, v in ions.items())
    anion = sum(v for k, v in ions.items() if k in ("CL", "CLA"))
    cation = sum(v for k, v in ions.items() if k not in ("CL", "CLA"))
    # Salt = ion pairs beyond neutralisation, expressed against water molarity (55.5 M).
    # CHARMM-GUI sets salt from the box volume, so the nominal value is usually lower (about 0.15 M).
    row["salt_pairs"] = min(anion, cation)
    row["salt_M_vs_water"] = round(min(anion, cation) / nw * 55.5, 3) if nw else ""
    row["water_model"] = "/".join(sorted({w for w in WATER if c[w]}))
    row["forcefield_include"] = forcefield(s["path"])
    return row


if __name__ == "__main__":
    ids = sys.argv[1:]
    rows = []
    for s in read_systems():
        if ids and s["sys_id"] not in ids: continue
        if not os.path.exists(os.path.join(s["path"], s["tpr"])):
            print(f"WARNING {s['sys_id']}: {s['tpr']} not found, left out of the table"); continue
        r = composition(s)
        p = parse_params(s["path"])
        r["ref_t_K"] = p.get("ref-t", "run 00_register_settings.sh")
        for k in ["integrator", "dt", "nsteps", "tcoupl", "tau-t", "pcoupl", "pcoupltype", "tau-p",
                  "coulombtype", "rcoulomb", "vdw-modifier", "rvdw-switch", "rvdw", "constraints",
                  "nstxout-compressed", "xtc_dt_ps", "xtc_last_ps", "xtc_frames"]:
            r[k] = p.get(k, "")
        # simulated length from the trajectory itself (nsteps is wrong for extended runs)
        r["length_us"] = round(p["xtc_last_ps"] / 1e6, 3) if "xtc_last_ps" in p else ""
        try:
            r["nsteps_us_in_tpr"] = round(int(p["nsteps"]) * float(p["dt"]) / 1e6, 3)
        except Exception:
            r["nsteps_us_in_tpr"] = ""
        rows.append(r)
    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(ANALYSIS_DIR, "system_table.tsv"), "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows: fh.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")
    with open(os.path.join(ANALYSIS_DIR, "system_table.md"), "w") as fh:
        fh.write("| field | " + " | ".join(r["sys_id"] for r in rows) + " |\n")
        fh.write("|---|" + "---|" * len(rows) + "\n")
        for k in keys[1:]:
            fh.write(f"| {k} | " + " | ".join(str(r.get(k, "")) for r in rows) + " |\n")
    print(open(os.path.join(ANALYSIS_DIR, "system_table.md")).read())

#!/usr/bin/env python3
"""Write a GROMACS index file (amp_index.ndx) for one AMP system, straight from the tpr.

Usage:  python make_index.py SYS_ID            (paths taken from systems.tsv)

Groups are built from peptide-local residue numbers (amp_config.PEPTIDES), so the
same definitions work for every topology. Group names are chosen so that typing
the name at a GROMACS prompt selects it unambiguously.
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import MDAnalysis as mda
from amp_config import get_system, PEPTIDES, LIPIDS, WATER, IONS


def peptide_copies(u, n_pep):
    prot = u.select_atoms("protein")
    copies = prot.fragments if len(prot.fragments) == n_pep else None
    if copies is None:          # fall back on equal splitting by residue count
        res = prot.residues; per = len(res) // n_pep
        copies = [res[i * per:(i + 1) * per].atoms for i in range(n_pep)]
    return sorted(copies, key=lambda ag: ag.ix[0])


def local(copy, lo, hi):
    return copy.residues[lo - 1:hi].atoms


def build_groups(u, sysrow):
    pep = PEPTIDES[sysrow["peptide"]]
    copies = peptide_copies(u, sysrow["n_pep"])
    nres = len(pep["seq"])
    for c in copies:
        if c.n_residues != nres:
            raise SystemExit(f"Copy has {c.n_residues} residues, expected {nres} for {pep['name']}")
    sel = lambda s: u.select_atoms(s)
    allpep = sum(copies[1:], copies[0])
    lipres = " ".join(LIPIDS)
    G = {}
    G["System"] = u.atoms
    G["Peptide"] = allpep
    G["Peptide_heavy"] = allpep.select_atoms("not name H*")
    G["Peptide_BB"] = allpep.select_atoms("backbone")
    G["Peptide_CA"] = allpep.select_atoms("name CA")
    mem = sel(f"resname {lipres}")
    if mem.n_atoms:
        G["Membrane"] = mem
        G["Lipid_P"] = mem.select_atoms("name P")
        for l in LIPIDS:
            ag = sel(f"resname {l}")
            if ag.n_atoms:
                G[f"Lipid_{l}"] = ag
    G["Water"] = sel("resname " + " ".join(WATER))
    G["Ions"] = sel("resname " + " ".join(IONS))
    G["Water_and_ions"] = G["Water"] + G["Ions"]
    G["Pep_Memb"] = allpep + mem if mem.n_atoms else allpep      # the "dry" analysis group
    # peptide segments, pooled over copies
    for name, lo, hi in pep["segments"]:
        key = "Seg_" + name.replace(" ", "")
        G[key] = sum((local(c, lo, hi) for c in copies[1:]), local(copies[0], lo, hi))
        G[key + "_BB"] = G[key].select_atoms("backbone")
    basic = sum((c.residues[i - 1].atoms for c in copies for i in pep["basic"]), u.atoms[[]])
    G["Basic_sidechains"] = basic.select_atoms("not backbone and not name H*")
    arom = sum((c.residues[i - 1].atoms for c in copies for i in pep["aromatic"]), u.atoms[[]])
    G["Aromatic_sidechains"] = arom.select_atoms("not backbone and not name H*")
    if len(copies) > 1:
        for k, c in enumerate(copies, 1):
            G[f"Copy{k}"] = c
            G[f"Copy{k}_BB"] = c.select_atoms("backbone")
    return {k: v for k, v in G.items() if v.n_atoms}


def write_ndx(groups, path):
    with open(path, "w") as fh:
        for name, ag in groups.items():
            fh.write(f"[ {name} ]\n")
            idx = ag.ix + 1
            for i in range(0, len(idx), 15):
                fh.write(" ".join(f"{x:6d}" for x in idx[i:i + 15]) + "\n")
            fh.write("\n")


if __name__ == "__main__":
    s = get_system(sys.argv[1])
    tpr = os.path.join(s["path"], s["tpr"])
    u = mda.Universe(tpr)
    G = build_groups(u, s)
    out = os.path.join(s["path"], "amp_index.ndx")
    write_ndx(G, out)
    summary = {k: int(v.n_atoms) for k, v in G.items()}
    with open(os.path.join(s["path"], "amp_index_groups.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(f"{s['sys_id']}: wrote {out}")
    for k, v in summary.items():
        print(f"   {k:<22s}{v:>8d}")

"""Shared definitions for the AMP (melittin, LL-37) MD pipeline.

Everything system specific lives in systems.tsv. Everything peptide specific
(sequence, segment boundaries, residues of interest) lives in PEPTIDES below.
Residue numbers here are peptide-local (1 = first residue of each copy), so they
work whatever numbering the topology uses (system 04 numbers melittin 161-186).
"""
import csv, os

ROOT = os.environ.get("AMP_ROOT", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
PIPE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SYSTEMS_TSV = os.environ.get("AMP_SYSTEMS", os.path.join(PIPE, "systems.tsv"))
ANALYSIS_DIR = os.environ.get("AMP_ANALYSIS", os.path.join(ROOT, "analysis"))
FIG_DIR = os.environ.get("AMP_FIGURES", os.path.join(ROOT, "figures"))

LIPIDS = ["POPC", "POPE", "POPG"]
WATER = ["SOL", "TIP3", "TIP3P", "HOH", "WAT"]
IONS = ["NA", "CL", "POT", "CLA", "SOD", "K"]

# Lipid atom classes (CHARMM36 names, identical for POPC, POPE and POPG)
LIPID_HEAD = "N C11 C12 C13 C14 C15 OC2 OC3 P O11 O12 O13 O14".split()   # choline/ethanolamine/glycerol + phosphate
LIPID_PHOS = "P O11 O12 O13 O14".split()
LIPID_GLYC = "C1 C2 C3 O21 O22 C21 O31 O32 C31".split()                   # glycerol backbone + ester carbonyls
SN1 = [f"C3{i}" for i in range(2, 17)]    # palmitoyl, C32..C316
SN2 = [f"C2{i}" for i in range(2, 19)]    # oleoyl, C22..C218 (double bond C29=C210)

PEPTIDES = {
    "mel": dict(
        name="Melittin",
        seq="GIGAVLKVLTTGLPALISWIKRKRQQ",
        # Segments used for shading, RMSD and per-segment metrics
        segments=[("Helix N", 1, 10), ("Hinge", 11, 14), ("Helix C", 15, 26)],
        # Residue windows used for helix axes (frayed termini excluded)
        axis_n=(2, 10), axis_c=(15, 25), axis_all=(2, 25),
        basic=[7, 21, 22, 23, 24],          # K7, K21, R22, K23, R24 (+ free N terminus)
        aromatic=[19],                      # W19
        cationic_tail=(21, 24),
    ),
    "ll37": dict(
        # Checked against the start structures (DSSP of the tpr coordinates): one continuous helix
        # 2 to 31 and an unstructured tail 32 to 37, as in 2K6O. Gly14 is the known flexible point.
        name="LL-37",
        seq="LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES",
        segments=[("Helix N", 1, 13), ("G14", 14, 16), ("Helix C", 17, 31), ("C tail", 32, 37)],
        axis_n=(3, 13), axis_c=(17, 30), axis_all=(3, 30),
        basic=[7, 8, 10, 12, 15, 18, 19, 23, 25, 29, 34],   # R7 K8 K10 K12 K15 K18 R19 R23 K25 R29 R34
        aromatic=[5, 6, 17, 27],                            # F5 F6 F17 F27
        cationic_tail=(17, 29),
    ),
}

CONDITIONS = {
    # label, colour (validated categorical slots, light surface), linestyle (secondary encoding)
    "solution": ("Solution", "#4d4d4d", ":"),
    "popc":     ("POPC", "#2a78d6", "-"),
    "pepg":     ("POPE:POPG 7:3", "#eb6834", "-"),
    "pcpepg":   ("POPC:POPE:POPG 7:2:1", "#1baf7a", "--"),
    "agg_popc": ("POPC, 4 peptides", "#184f95", "-."),
    "agg_surf": ("POPC, 4 peptides, surface start", "#eda100", (0, (5, 1, 1, 1))),
}
COPY_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]


def read_systems(path=SYSTEMS_TSV):
    rows = []
    with open(path) as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("#")]
    for r in csv.DictReader(lines, delimiter="\t"):
        r = {k.strip(): v.strip() for k, v in r.items()}
        r["n_pep"] = int(r["n_pep"]); r["whole"] = int(r["whole"])
        r["path"] = os.path.join(ROOT, r["dir"])
        rows.append(r)
    return rows


def get_system(sys_id, path=SYSTEMS_TSV):
    for r in read_systems(path):
        if r["sys_id"] == sys_id:
            return r
    raise KeyError(f"{sys_id} not in {path}")


# Analysis period. Every system is analysed from 0 to END_NS (default 1 us) so that all runs
# average the same period (window = last half of 0..END_NS, i.e. 500 to 1000 ns).
# A system is analysed further only when check_extension.py justifies it and a line is added
# to pipeline/analysis_window.tsv (sys_id, end_ns, reason).
END_NS = float(os.environ.get("END_NS", 1000))
WINDOW_TSV = os.path.join(PIPE, "analysis_window.tsv")


def end_ns_for(sys_id):
    if os.path.exists(WINDOW_TSV):
        for line in open(WINDOW_TSV):
            f = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(f) < 2 or f[0] == "sys_id": continue
            if f[0] == sys_id:
                return float(f[1])
    return END_NS

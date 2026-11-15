# AMP MD pipeline

Reproducible processing and analysis of molecular dynamics simulations of
membrane-active antimicrobial peptides, with a forward-modelling bridge to
hydrogen–deuterium exchange mass spectrometry.

Built for melittin and LL-37 across solution, zwitterionic and anionic bilayers
and multi-peptide aggregation boxes, but the registry design means adding a
peptide or a replicate is one new line, not a new script.

---

## The problem it solves

A set of trajectories recorded under different conditions can only be compared
if every one of them receives the same periodic-boundary treatment, the same
index groups, the same analysis window and the same metrics. Doing that by hand
across a dozen systems produces drift: one trajectory centred differently, one
figure on an older colour scheme, one analysis run over a different time window.

This pipeline makes the system definition the only thing that varies. One
registry file describes every trajectory; four steps read it; every system is
treated identically by construction.

## Design

```
systems.tsv            one row per trajectory: peptide, condition, replicate,
                       kind, number of copies, paths
   |
   v
00_register_settings   run parameters read back from each tpr
01_make_index          index groups built from the topology, by name
02_process_trajectories  PBC treatment chosen by system kind, reduced trajectories
03_run_analysis        per-system analysis, then figures and summary tables
```

Index groups are referenced **by name**, not by number, so group numbering
cannot silently shift between systems.

### PBC treatment by system kind

| kind | treatment | why |
| --- | --- | --- |
| solution | whole, nojump, centre peptide with `-pbc mol`; backbone fit for visualisation only | peptide never split, no jumps across the box |
| membrane, aggregation | whole, centre on membrane with `-pbc res`, then `-pbc whole` again | in these topologies the bilayer is a single moleculetype, so `-pbc mol` would treat every lipid as one molecule; `-pbc res` wraps lipid by lipid and the second pass rejoins the peptides |

That second row is the kind of detail that is invisible until it corrupts an
analysis, which is why it is encoded once rather than remembered each time.

## Installation

```bash
git clone <this repository>
cd amp-md-pipeline
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

GROMACS is a separate dependency and is not installed by pip. The pipeline calls
it through `$GMX` (default `gmx`), so point that at your build:

```bash
export GMX=gmx_mpi
```

## Running it

```bash
cp systems.example.tsv systems.tsv     # then edit paths for your trajectories
bash run_all.sh                        # every system, every step
bash run_all.sh mel_popc mel_pepg      # selected systems
bash 02_process_trajectories.sh mel_pepg   # one step, one system
sbatch slurm_run_all.sh                # on a cluster (edit the module lines first)
```

Settings are environment variables: `GMX`, `DT_ANALYSIS` (default 200 ps),
`DT_WET` (1000 ps), `DT_VIZ` (1000 ps), `STRIDE`, `PARALLEL`, `WINDOW_FRAC`
(default 0.5, the last half of the run) and `START_NS`.

### Cost

At 200 ps per frame a 3 µs membrane run gives 15,000 frames, a reduced
trajectory of roughly 1.5 GB, and five to ten minutes of Python analysis. A
262,000-atom aggregation box takes about 0.2 s per frame. No full-resolution
copy is written: frames are thinned in the first `trjconv` call.

## What it computes

| Question | Metrics |
| --- | --- |
| Does the membrane stabilise the helix and its hinge, relative to water? | per-residue helix fraction (DSSP), helicity against time, bend angle between the N- and C-terminal helices |
| Is the conformational ensemble narrower on the membrane? | free energy surface of helicity against radius of gyration, per condition |
| When, how deep, and in what orientation does the peptide bind? | centre-of-mass height above the proximal phosphate plane, per-residue Cα depth, helix tilt relative to the bilayer normal, residues in lipid contact |
| Which residues make the contacts, and does the peptide recruit anionic lipid? | per-residue contact occupancy split into head group, glycerol/ester and acyl chain; first-shell enrichment per lipid type |
| Does the bound peptide perturb the bilayer locally? | sn-1 order parameters by distance shell and leaflet, phosphate height against lateral distance, P–P thickness, water in the hydrophobic core |
| Does the peptide self-associate at the simulated peptide-to-lipid ratio? | pairwise minimum distances, largest cluster size, inter-peptide residue contact map, per-copy height |

## HDX forward modelling

`py/mel_hdx_compare.py` predicts deuterium uptake from the simulated ensembles
and scores it against measurement. Two independent models are implemented:

- **Open model.** `k_obs = k_ch · f_open`, where `f_open` is the fraction of
  frames in which the amide hydrogen makes no hydrogen bond to peptide or lipid
  oxygen. EX2 regime, closed states do not exchange.
- **Best–Vendruscolo.** `ln P = 0.35 N_c + 2.0 N_h`, averaged over frames, in
  peptide-only and environment-inclusive variants (the latter counting lipid
  heavy atoms and oxygens).

Uptake follows `D(t) = 1 − exp(−k_ch · t / P)`. Intrinsic amide exchange rates
come from `py/kint.py`, implementing Bai et al. 1993 (*Proteins* **17**:75).

**Measurements are not included in this repository.** The models need observed
uptake to score against, and those numbers belong to the experiments that
produced them. Supply your own in `data/hdx/` — the format is documented in
`data/hdx/README.md`, with a skeleton in
`data/hdx/melittin_observed.example.csv`.

## Layout

```
run_all.sh, 0*.sh        pipeline steps
lib_common.sh            shared shell helpers
systems.example.tsv      registry template
analysis_window.tsv      per-system analysis windows beyond the default
py/                      analysis, figures, HDX modelling
  amp_config.py          peptide definitions, conditions, colours, lipid atom classes
  amp_analysis.py        per-system analysis
  amp_figures.py         figure generation
  amp_style.py           house figure style
  kint.py                intrinsic amide exchange rates
  amide_hydration.py     per-residue amide hydration
  mel_hdx_compare.py     forward models and scoring
  hdx_observed.py        loader for experimental measurements
thesis_figures/          chapter-specific figure scripts, not part of the core pipeline
new_systems/             per-system processing helpers
data/hdx/                experimental data location (contents gitignored)
```

## Adding a peptide

Peptides are declared once in `py/amp_config.py`: sequence, segment boundaries
for shading and per-segment metrics, helix axis windows excluding frayed
termini, and the basic and aromatic residues of interest. Residue numbers are
peptide-local, so they are independent of whatever numbering the topology uses.

Adding a system is one row in `systems.tsv`. Adding a replicate is one row with
a different `replicate` value.

## Status and scope

This is research code, published as it was used, with experimental data removed.
It has been exercised on roughly twenty trajectories of two peptides. It is not
a general-purpose package: it assumes CHARMM36 lipid atom naming, GROMACS
output, and the directory conventions described above. The analysis and HDX
modules are the parts most likely to be useful elsewhere.

## Citation and licence

See `LICENSE`. If this code contributes to published work, please cite the
repository and the Bai et al. and Best and Vendruscolo references for the
exchange models.

# Processing the new melittin and LL-37 runs on the Mac

**Use these only if the raw trajectories are on the laptop. For runs kept on the external drive, use the `process.sh` in each run folder of `/Volumes/Transcend/mel_additional` instead.**

One script per run. Each one sends the run through your existing pipeline (`AMP_MD_processed/pipeline`), so the outputs match the runs already in `AMP_MD_processed`: same index groups, same PBC treatment, same frame spacing and the same per-system analysis tables.

## Where each run goes

Copy the production `tpr` and `xtc` of each run into its folder (create the folder). Names are found automatically: `md.tpr`/`md.xtc` (solution builds), `step7.tpr`/`step7.xtc` or `step7_1.tpr` (membrane build). Trajectories in several parts (`md.part0002.xtc` ..., or `step7_1.xtc`, `step7_2.xtc` ...) are joined into `md_joined.xtc` first.

| Script | Run | Folder in AMP_MD_processed | kind, n_pep |
|---|---|---|---|
| `process_mel_coil300.sh` | melittin coil A, water, 300 K | `mel/10_mel_coil_300K/A` | solution, 1 |
| `process_mel_coil300_r2.sh` | melittin coil B | `mel/10_mel_coil_300K/B` | solution, 1 |
| `process_mel_coil300_r3.sh` | melittin coil C | `mel/10_mel_coil_300K/C` | solution, 1 |
| `process_mel_agg_water300.sh` | 4 melittin, water, 300 K | `mel/11_mel4_water_300K` | solution, 4 |
| `process_mel_popc300.sh` | melittin on POPC, 300 K | `mel/12_mel_popc_300K` | membrane, 1 |
| `process_ll37_sol313.sh` | LL-37 2K6O, water, 313 K | `ll37/10_ll37_2k6o_313K` | solution, 1 |
| `process_ll37_coil313.sh` | LL-37 coil A, 313 K | `ll37/11_ll37_coil_313K/A` | solution, 1 |
| `process_ll37_coil313_r2.sh` | LL-37 coil B | `ll37/11_ll37_coil_313K/B` | solution, 1 |
| `process_ll37_coil313_r3.sh` | LL-37 coil C | `ll37/11_ll37_coil_313K/C` | solution, 1 |
| `process_ll37_agg_water313.sh` | 4 LL-37, water, 313 K | `ll37/12_ll37x4_water_313K` | solution, 4 |

## Run

Needs `gmx` (`brew install gromacs`) and Python with `MDAnalysis>=2.8`, `pandas`, `matplotlib` (as for the existing pipeline).

```bash
cd ~/Desktop/AMP_MD_processed
bash pipeline/new_systems/process_mel_popc300.sh          # one run
caffeinate -i bash pipeline/new_systems/run_all_new.sh    # every run whose folder has files; the rest are listed as skipped
```

Options, before `bash`: `GMX=/path/to/gmx`, `TPR=name.tpr XTC=name.xtc` (explicit file names), `ANALYSE=0` (processing only), `FORCE=1` (redo steps that already have output). A stopped run resumes where it stopped. Only one run at a time (the pipeline lock); `run_all_new.sh` runs them one after another.

## What each script does

1. Finds the tpr and xtc, joins parts if needed.
2. Adds the run to `pipeline/systems.tsv` (or updates its row; a dated backup of the file is kept).
3. `00_register_settings.sh`: `settings/key_params.txt`, `settings/xtc_check.txt`, and `analysis/system_table_<sys_id>.*` (your existing `analysis/system_table.*` is left unchanged).
4. `01_make_index.sh`: `amp_index.ndx`, `amp_index_groups.json`.
5. `02_process_trajectories.sh`: `processed/dry.xtc` (every 200 ps), `viz.xtc` and `viz_start.pdb` (every 1 ns), `dry.tpr`, and `wet.xtc` (every 1 ns) for the POPC run.
6. `amp_analysis.py`: `analysis/<sys_id>/` (timeseries, per-residue arrays, blocks; `agg_pairs.csv` and `agg_clusters.csv` for the four-peptide runs), 0 to 1000 ns.

## Decisions to check

* **The water tetramers are `kind = solution`, `n_pep = 4`** (not `aggregation` as suggested in the setup READMEs). In your pipeline every kind other than `solution` is treated as having a membrane, which would fail without lipids. With `solution` and 4 peptides, the analysis still writes the peptide-peptide contact and cluster tables. This path (solution with 4 peptides) has not been run before in your pipeline, so check the first output.
* **Visualisation fit for the tetramers.** For solution runs the pipeline fits `viz.xtc` on the backbone of all peptides. For four separate peptides this fit is not meaningful; use `dry.xtc` (not fitted) in VMD instead.
* **Figures are not regenerated.** `amp_figures.py` needs a label for each condition in `CONDITIONS` in `pipeline/py/amp_config.py`. The new conditions are `coil_300K`, `agg_water_300K`, `popc_300K`, `solution_313K`, `coil_313K`, `agg_water_313K`. Add them when you decide how the 300/313 K runs should appear next to the 343 K runs, then run `bash pipeline/03_run_analysis.sh` as before.
* **`dry.tpr` may be missing** for membrane runs: `gmx convert-tpr` aborts on some CHARMM-GUI tprs, exactly as for 4 of your existing runs. The analysis does not need it.

## Tested here

GROMACS 2023.3 and the pipeline from your folder, on short test trajectories made from the `mel_sol` and `mel_popc_r2` run inputs (solution and membrane branches), with trajectories in two parts each: parts joined without duplicate frames, registration, settings, index, processing and analysis all completed; rerunning skipped every finished step.

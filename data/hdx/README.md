# Experimental HDX-MS observations

The forward models in `py/mel_hdx_compare.py` predict deuterium uptake from the
simulated ensembles. Scoring them needs measured uptake.

`melittin_observed.csv` holds the measurements used in this work: melittin at
278 K, pD 6.0, labelled 3 to 120 s, free peptide and with POPC liposomes.

To analyse a different peptide, add `<peptide>_observed.csv` in the same format,
or point `AMP_HDX_DATA` at a directory elsewhere.

## Format

| column | meaning |
| --- | --- |
| `quantity` | `seg_dRFU`, `free_uptake` or `protection_gain` |
| `segment` | segment label for `seg_dRFU` rows, blank otherwise |
| `t_s` | labelling time in seconds |
| `value` | the measured number |

Segment labels must match the `SEGS` list in the comparison script, and times
must match `TS`. `melittin_observed.example.csv` is a complete skeleton with
placeholder values; copy it, replace the numbers, and drop the `.example`.

`OBS` is indexed `[segment, time]`. Rows in the CSV are grouped by segment and
then by labelling time, so that the file reads in the same order the array is
built. The loader raises if any segment and time combination is missing, rather
than leaving a silent `NaN`.

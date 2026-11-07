# Medication-use and adherence events for oral orphan drugs in FAERS

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23956483.svg)](https://doi.org/10.5281/zenodo.23956483)

Analysis code for a disproportionality study of medication-use, dosing-administration and
adherence-related events reported for oral orphan drugs in the US Food and Drug Administration
Adverse Event Reporting System (FAERS). The rare-disease drugs are compared with common chronic-disease
drugs and with specialty-pharmacy drugs that share the same distribution channel.

Two outcomes are analysed:

- MedDRA preferred terms (PTs) for nonadherence and related categories. These were reviewed manually.
- The Standardised MedDRA Query (SMQ) *Medication errors* (SMQ code 20000224), expanded to its narrow
  and broad scopes.

## Contents

```
scripts/
  expand_smq.py                          expand an SMQ to its PTs and their lowest level terms (LLTs)
  faers_me_feasibility_v7_staged.py      main analysis, run one phase at a time (low memory, resumable)
  faers_me_feasibility_v7_singlepass.py  single-pass reference implementation of the same analysis
  meta_regression.py                     drug-level meta-regression (REML, Hartung–Knapp), Table 4
  make_figures.py                        Figure 1 and Figure S1
  screen_comparator_counts.py            count-only screen of the second specialty set (DEMO and DRUG only)
  summarise_configs.py                   cross-configuration summary of the specialty comparison (version 1.1.0)
configs/                                 comparator configurations for --groups (version 1.1.0)
docs/second_comparator_rules.md          prespecified rules and analysis plan for the second specialty set
me_pt_dictionary_reviewed.csv            reviewed PT dictionary
```

The drug identification rules are in the `CANDIDATES` table near the top of both analysis scripts.
These rules are regular expressions for generic and brand names, exclusions, formulation splits and
comparator groups. Target-indication rules for the indication-restricted analysis are in `INDICATIONS`.

`me_pt_dictionary_reviewed.csv` has the following columns:

- `pt`: PT name.
- `n_reports_all_faers`: database-wide report count at screening.
- `keep`: Y or N.
- `category`: outcome category.
- `note`: reason for exclusion or comment.

Only rows with `keep = Y` define outcomes. Former and current names of renamed PTs are both listed.

## Data version

- FAERS quarterly ASCII files, 2013Q1 to 2025Q4 (52 quarters), downloaded on 18 July 2026.
- MedDRA version 29.1 (English), used for the SMQ expansion.

## Requirements

- Python 3.10 or later, with the packages in `requirements.txt` (`pip install -r requirements.txt`).
- FAERS quarterly ASCII files. These are the `DEMO`, `DRUG`, `REAC` and `INDI` tables, `$`-delimited,
  from the FDA FAERS Quarterly Data Extract Files page. Keep them as the downloaded zip files or unpack
  them. The single-pass script reads unpacked files only.
- Your own licensed MedDRA release. **MedDRA is not distributed here.** `expand_smq.py` reads three
  ASCII files from it: `smq_content.asc`, `pt.asc` and `llt.asc`. The SMQ term lists it generates are
  not part of this repository.

## Running the analysis

Step 1. Expand SMQ 20000224 into narrow and broad term lists:

```bash
python scripts/expand_smq.py --meddra-dir /path/to/MedDRA_29_1_English/MedAscii \
    --smq-code 20000224 --out-dir ./smq
```

Step 2. Run the five analysis phases in order. Each phase writes intermediate arrays under
`out/_work/` and skips quarters that are already done, so an interrupted run can be restarted.

```bash
for phase in demo reac drug indi agg; do
  python scripts/faers_me_feasibility_v7_staged.py \
    --faers-dir /path/to/faers_zips --out ./out --phase $phase \
    --start 2013Q1 --end 2025Q4 \
    --pt-dict me_pt_dictionary_reviewed.csv \
    --smq-file narrow=smq/smq_20000224_narrow.txt \
    --smq-file broad=smq/smq_20000224_broad.txt
done
```

How the script selects and cleans the data:

- Only quarters that contain all of `DEMO`, `DRUG` and `REAC` are used.
- Pre-2012Q4 LAERS files are skipped. Their primary key is `ISR` rather than `primaryid`.
- FDA deleted-case lists are applied before de-duplication. The FDA has published them under four
  file-naming conventions (`DELEnnQn.TXT`, `ADRnnQnDeletedCases.txt`, `nnQnDeletedCases.txt`,
  `AllDeletedCases.txt`); all four are read. Matching only `DELE*` files would miss 133,631 of the
  218,881 listed caseids and change the denominators: 57,984 rather than 112,832 records removed,
  and 16,735,500 rather than 16,687,483 reports after de-duplication.
- De-duplication keeps, for each `caseid`, the record with the latest `fda_dt`. If two records tie,
  it keeps the one with the larger `primaryid`.
- Drugs are attributed only when they are the primary suspect (`role_cod = PS`).

Step 3. Run the meta-regression and draw the figures from the `agg` outputs:

```bash
python scripts/meta_regression.py --drug-level out/7_drug_level.csv --out out/meta_regression_HK.csv \
    --loo-out out/meta_regression_HK_leave_one_out.csv
python scripts/make_figures.py --results ./out --out ./figures
```

Step 4 (version 1.1.0). Re-run the `agg` phase once for each comparator configuration in `configs/`.
Each run writes the same numbered tables to its own folder, plus `config_used.json`, and Step 3 is then
run on that folder. Without `--groups` the outputs are identical to version 1.0.2.

```bash
# the drug phase must be re-run once, because version 1.1.0 adds six candidate drugs;
# the demo, reac and indi caches are reused
rm -rf out/_work/drug5
python scripts/faers_me_feasibility_v7_staged.py --faers-dir /path/to/faers_zips --out ./out \
    --phase drug --start 2013Q1 --end 2025Q4
for cfg in configs/c*.json; do
  name=$(basename "$cfg" .json)
  python scripts/faers_me_feasibility_v7_staged.py \
    --faers-dir /path/to/faers_zips --out ./out_$name --work ./out/_work --phase agg \
    --start 2013Q1 --end 2025Q4 \
    --pt-dict me_pt_dictionary_reviewed.csv \
    --smq-file narrow=smq/smq_20000224_narrow.txt \
    --smq-file broad=smq/smq_20000224_broad.txt \
    --groups "$cfg"
  python scripts/meta_regression.py --drug-level out_$name/7_drug_level.csv \
    --out out_$name/meta_regression_HK.csv --loo-out out_$name/meta_regression_HK_leave_one_out.csv
done
```

The six candidates of the second set (`s2_*`, `s2c_*`) are inactive unless a configuration gives them a
role. A configuration can set the role of each drug or of each original group (`rare`, `specialty`,
`common`, `example` or `off`). It can also drop a drug's reports from a given FAERS quarter onward
(`truncate_from`). The configurations follow section 7 of `docs/second_comparator_rules.md`:

| File | Specialty reference | Purpose |
|---|---|---|
| `c1_specialty2_main.json` | enzalutamide, apalutamide, darolutamide, abemaciclib | Main second set (analyses 1–6) |
| `c2_specialty2_with_cyclic.json` | main second set + palbociclib, ribociclib | Cyclic sensitivity set |
| `c3_specialty_pooled.json` | first set + main second set (nine drugs) | Pooled specialty reference |
| `c4_strict_nonorphan.json` | teriflunomide, tofacitinib + main second set | Sensitivity for C3 (no OOPD designation) |
| `c5_specialty1_truncated.json` | first set, truncated at first US generic entry | Sensitivity for C5 |
| `c6_specialty2_vs_specialty1.json` | first set, with the main second set in the `rare` role | Direct comparison of the two sets |

In `c6`, every column or row labelled `rare` refers to the main second set, and the indication-restricted
rows do not apply.

Step 5 (version 1.1.0). Summarise the configurations and draw Figure 1 with both specialty sets:

```bash
python scripts/summarise_configs.py --main ./out --config ./out_c1_* ./out_c2_* ./out_c3_* ./out_c4_* \
    ./out_c5_* ./out_c6_* --out cross_config_summary.csv
python scripts/make_figures.py --results ./out --out ./figures --second-specialty ./out_c1_specialty2_main
```

`cross_config_summary.csv` gives, for each configuration and outcome:

- the report-level and drug-level (Hartung–Knapp) estimates against the specialty reference;
- the ratio of the drug-level estimate to the main run's estimate for the same outcome, and whether it
  lies in the prespecified range;
- the drug-clustered logistic estimates, with normal-based intervals and intervals using t critical
  values with (number of drug clusters − 1) degrees of freedom;
- the range across report subsets.

Without `--second-specialty`, `make_figures.py` draws Figure 1 exactly as in version 1.0.2.
From version 1.1.1 it also draws the flow of records (`FigureS1_flow.png`, from table 14; the second-set box
is added with `--second-specialty`) and writes the annual figure as `FigureS2_annual_nonadherence.png`.

`faers_me_feasibility_v7_singlepass.py` runs the same analysis in one pass and is kept for reference.
On the full 52-quarter extract it needs far more memory than the staged script.

Both implementations now read all four deleted-case file names. The results reported in the paper
were produced with the staged script on the full 52-quarter extract (2013Q1-2025Q4): 218,881 listed
caseids, 112,832 records removed and 16,687,483 reports after de-duplication.

## Output

The `agg` phase writes numbered CSV tables `0_data_summary.csv` to `14_flow_and_characteristics.csv` (listed in
the script header):

- data summary and per-drug counts;
- PT-level hits with SMQ annotation;
- reporter-type, yearly and sensitivity analyses;
- drug-level Wilson intervals, specialty leave-one-out and DerSimonian–Laird summaries;
- indication coverage;
- the Nefecon worked example;
- the report-level logistic regression with drug-clustered standard errors;
- record counts at each processing step and the characteristics of each group's reports (sex, age,
  country, year of receipt, report type, reporter), for the READUS-PV flow diagram and case description
  (version 1.1.1; this table is new and leaves tables 0–13 unchanged);
- per-drug forest plots.

## Data and licensing notes

- FAERS data are public and de-identified. No data beyond the public quarterly files are used.
- MedDRA® is the international medical terminology developed under the auspices of the International
  Council for Harmonisation of Technical Requirements for Pharmaceuticals for Human Use (ICH). MedDRA
  and the SMQ term lists derived from it are subject to the MedDRA licence and are not redistributed here.
- The code is released under the MIT licence (see `LICENSE`).

## Citation

If you use this code, cite the archived release:

> Xia Yuhang. *Medication-use and adherence events for oral orphan drugs in FAERS* (v1.1.1).
> Zenodo; 2026. doi:10.5281/zenodo.23956483 (concept DOI; the version DOI is listed on the Zenodo record)

The concept DOI 10.5281/zenodo.23956483 always resolves to the most recent version.

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Medication-use and adherence events for oral orphan drugs in FAERS: staged (low-memory) analysis.

Reads FAERS quarterly ASCII files (FAERS format from 2012Q4; DEMO, DRUG, REAC, INDI; '$'-delimited),
either unpacked or as the FDA quarterly zip files. Pre-2012Q4 LAERS files (keyed on ISR) are skipped.
FDA deleted-case lists found under --faers-dir are applied before de-duplication; the file names
recognised are DELE*.TXT, ADR*DeletedCases.txt, nnQnDeletedCases.txt and AllDeletedCases.txt.

The analysis runs in five phases that write intermediate arrays to <out>/_work and can be resumed:
    python faers_me_feasibility_v7_staged.py --faers-dir DIR --out ./out --phase demo|reac|drug|indi|agg \
        --pt-dict me_pt_dictionary_reviewed.csv \
        --smq-file narrow=smq/smq_20000224_narrow.txt --smq-file broad=smq/smq_20000224_broad.txt
Options: --start/--end (e.g. 2013Q1, 2025Q4), --min-drug-n (default 100).
Version 1.1.0: --groups FILE (agg phase) assigns each candidate drug a role (rare, specialty, common,
example, off) and can drop a drug's reports from a given quarter onward; the files used for the paper are
in configs/ (rules: docs/second_comparator_rules.md). Without --groups the six second-set candidates are
off and every output is identical to version 1.0.2. With --groups, config_used.json is also written.

Outputs (CSV, UTF-8 with BOM):
    0_data_summary.csv            quarter range, deleted-case records removed, reports after de-duplication
    14_flow_and_characteristics.csv  record counts at each step and, per group, sex, age, country, year of
                                  receipt, report type and reporter (READUS-PV items 8a and 8b; version 1.1.1)
    1_summary_by_drug.csv         per drug: reports, events per outcome category, crude ROR (95% CI) against
                                  the whole database and both comparator groups, consumer-reported share
    2_pt_by_drug.csv              preferred terms counted for each drug, with SMQ annotation
    3_reporter_by_drug.csv        counts stratified by reporter type (occp_cod) and report type (rept_cod)
    4_yearly_by_drug.csv          counts by year of FDA receipt
    5_keyword_pt_dictionary.csv   every PT matched by the screening keywords, with database-wide counts;
                                  the reviewed version of this table is me_pt_dictionary_reviewed.csv
    6_sensitivity.csv             pooled rare-disease drugs vs each comparator group in report subsets
                                  (periodic reports excluded, health professionals only, peak year excluded,
                                  target indication only, leave-one-drug-out)
    7_drug_level.csv              per-drug event proportions with Wilson 95% CI
    8_specialty_leave_one_out.csv results after dropping each specialty comparator in turn
    9_random_effects.csv          drug-level DerSimonian-Laird pooling within groups
    10_indication_coverage.csv    share of reports whose INDI record matches the target rare disease
    11_worked_example.csv         Nefecon (Tarpeyo) dose-omission signal recomputed against the whole
                                  database and both comparator groups, 2021Q4-2023Q4
    12_multivariable_logistic.csv report-level logistic regression with drug-clustered robust SEs
    13_indication_top.csv         most frequent indication PTs per rare-disease drug
    config_used.json              the --groups file and the roles it resolved to (only with --groups)
    forest_*.png                  per-drug forest plots (requires matplotlib)

Definitions:
    * de-duplication: FDA deleted cases removed first; then, per caseid, the record with the latest
      fda_dt and, on ties, the largest primaryid is kept
    * drug attribution: primary suspect only (role_cod == "PS")
    * outcomes: reviewed PT dictionary (--pt-dict) and, optionally, SMQ term lists (--smq-file)
    * the storage notes in the CANDIDATES table are descriptive only and are not used in any analysis
"""
import argparse
import contextlib
import csv
import os
import re
import sys
import zipfile
from collections import defaultdict

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------
CANDIDATES = [
    ("penicillamine", r"PENICILLAMINE|CUPRIMINE|\bDEPEN\b", None, "rare", "room temperature"),
    ("trientine_all", r"TRIENTINE|SYPRINE|CUVRIOR|CUFENCE", None, "rare", "see the formulation entries"),
    ("trientine_tetrahcl_cuvrior", r"CUVRIOR|TETRAHYDROCHLORIDE", None, "rare", "room temperature (to verify)"),
    ("trientine_dihcl_syprine", r"TRIENTINE|SYPRINE|CUFENCE", r"CUVRIOR|TETRAHYDROCHLORIDE", "rare", "refrigerated (to verify)"),
    ("zinc_acetate_galzin", r"GALZIN", None, "rare", "room temperature; brand name only, to exclude over-the-counter zinc products"),
    ("nitisinone_all", r"NITISINONE|ORFADIN|NITYR|HARLIKU", None, "rare", "see the formulation entries"),
    ("nitisinone_orfadin", r"ORFADIN", None, "rare", "capsule, refrigerated (to verify)"),
    ("nitisinone_nityr", r"NITYR", None, "rare", "tablet, room temperature (to verify)"),
    ("nitisinone_unbranded", r"NITISINONE|HARLIKU", r"ORFADIN|NITYR", "rare", "no brand name, formulation unclear"),
    ("cysteamine_oral", r"CYSTEAMINE|MERCAPTAMINE|CYSTAGON|PROCYSBI", r"CYSTARAN|CYSTADROPS|OPHTH|EYE", "rare", "see the formulation entries"),
    ("cysteamine_cystagon_ir", r"CYSTAGON", None, "rare", "immediate release, every 6 hours"),
    ("cysteamine_procysbi_dr", r"PROCYSBI", None, "rare", "delayed release, every 12 hours; storage to verify"),
    ("cysteamine_unbranded", r"CYSTEAMINE|MERCAPTAMINE", r"CYSTAGON|PROCYSBI|CYSTARAN|CYSTADROPS|OPHTH|EYE", "rare", "no brand name, formulation unclear"),
    ("carglumic_acid", r"CARGLUMIC|CARBAGLU", None, "rare", "refrigerate before opening (to verify)"),
    ("sapropterin", r"SAPROPTERIN|KUVAN|JAVYGTOR|TETRAHYDROBIOPTERIN", None, "rare", "requires dissolution"),
    ("sepiapterin", r"SEPIAPTERIN|SEPHIENCE", None, "rare", "recently approved, likely few reports"),
    ("tiopronin", r"TIOPRONIN|THIOLA", None, "rare", "room temperature"),
    ("miglustat", r"MIGLUSTAT|ZAVESCA|OPFOLDA", None, "rare", "room temperature; OPFOLDA is used with enzyme replacement in Pompe disease, split by indication"),
    ("eliglustat", r"ELIGLUSTAT|CERDELGA", None, "rare", "room temperature; CYP2D6 genotype required"),
    ("sodium_phenylbutyrate", r"SODIUM PHENYLBUTYRATE|BUPHENYL|PHEBURANE|OLPRUVA", r"TAURURSODIOL|RELYVRIO|ALBRIOZA|AMX0035", "rare", "high-dose powder/tablet; the ALS combination products Relyvrio/Albrioza are excluded"),
    ("glycerol_phenylbutyrate", r"GLYCEROL PHENYLBUTYRATE|RAVICTI", None, "rare", "oral liquid"),
    ("betaine_anhydrous", r"BETAINE ANHYDROUS|CYSTADANE", None, "rare", "powder, requires measuring"),
    ("cholic_acid", r"\bCHOLIC ACID\b|CHOLBAM", None, "rare", "capsule"),
    ("spec_teriflunomide", r"TERIFLUNOMIDE|AUBAGIO", None, "specialty", "specialty comparator"),
    ("spec_dimethyl_fumarate", r"DIMETHYL FUMARATE|TECFIDERA", None, "specialty", "specialty comparator"),
    ("spec_fingolimod", r"FINGOLIMOD|GILENYA", None, "specialty", "specialty comparator"),
    ("spec_apremilast", r"APREMILAST|OTEZLA", None, "specialty", "specialty comparator"),
    ("spec_tofacitinib", r"TOFACITINIB|XELJANZ", None, "specialty", "specialty comparator"),
    ("example_nefecon", r"NEFECON|TARPEYO|KINPEYGO", None, "example", "worked example: Wang 2024 Kidney Int Rep"),
    ("ctrl_levothyroxine", r"LEVOTHYROXINE|SYNTHROID|EUTHYROX", None, "common", "comparator"),
    ("ctrl_metformin", r"\bMETFORMIN\b", r"SITAGLIPTIN|SAXAGLIPTIN|LINAGLIPTIN|EMPAGLIFLOZIN|DAPAGLIFLOZIN|CANAGLIFLOZIN|ALOGLIPTIN|ERTUGLIFLOZIN|PIOGLITAZONE|GLIPIZIDE|GLYBURIDE", "common", "comparator (single-ingredient only)"),
    ("ctrl_atorvastatin", r"ATORVASTATIN|LIPITOR", r"AMLODIPINE|EZETIMIBE", "common", "comparator (single-ingredient only)"),
    # Second specialty-pharmacy comparator set (docs/second_comparator_rules.md). Inactive unless a --groups
    # configuration gives them a role, so the default run reproduces version 1.0.2.
    ("s2_enzalutamide", r"ENZALUTAMIDE|XTANDI", None, "inactive", "second set, main: AR inhibitor, continuous"),
    ("s2_apalutamide", r"APALUTAMIDE|ERLEADA", None, "inactive", "second set, main: AR inhibitor, continuous"),
    ("s2_darolutamide", r"DAROLUTAMIDE|NUBEQA", None, "inactive", "second set, main: AR inhibitor, continuous"),
    ("s2_abemaciclib", r"ABEMACICLIB|VERZENIO", None, "inactive", "second set, main: CDK4/6 inhibitor, continuous"),
    ("s2c_palbociclib", r"PALBOCICLIB|IBRANCE", None, "inactive", "second set, cyclic: CDK4/6 inhibitor, 21 days on / 7 off"),
    ("s2c_ribociclib", r"RIBOCICLIB|KISQALI", None, "inactive", "second set, cyclic: CDK4/6 inhibitor, 21 days on / 7 off"),
]
ROLES = ("rare", "specialty", "common", "example", "off")


def load_groups_config(path):
    """Role of every candidate, and optional truncation quarters, from a --groups JSON file.

    Without a file every candidate keeps the group given in CANDIDATES and "inactive" entries are off,
    which reproduces version 1.0.2. The file may contain
        "name":           label written to config_used.json
        "group_roles":    {original group: role}, applied first (e.g. {"rare": "off"})
        "roles":          {candidate key: role}, applied second
        "truncate_from":  {candidate key: "YYYYQn"}; reports from FAERS quarterly files of that quarter
                          onward are dropped for that drug (e.g. from first US generic entry)
    Roles are rare, specialty, common, example and off.
    """
    import json
    orig = {k: g for k, _, _, g, _ in CANDIDATES}
    roles = {k: ("off" if g == "inactive" else g) for k, g in orig.items()}
    cfg = {}
    if path:
        with open(path, encoding="utf-8") as f:
            cfg = json.load(f)
        unknown = set(cfg) - {"name", "description", "group_roles", "roles", "truncate_from"}
        if unknown:
            sys.exit(f"--groups: unknown fields {sorted(unknown)}")
        for g, r in cfg.get("group_roles", {}).items():
            if g not in set(orig.values()) or r not in ROLES:
                sys.exit(f"--groups: bad group_roles entry {g!r}: {r!r}")
            for k in roles:
                if orig[k] == g:
                    roles[k] = r
        for k, r in cfg.get("roles", {}).items():
            if k not in roles or r not in ROLES:
                sys.exit(f"--groups: bad roles entry {k!r}: {r!r}")
            roles[k] = r
    trunc = {}
    for k, qs in cfg.get("truncate_from", {}).items():
        if k not in roles:
            sys.exit(f"--groups: unknown drug in truncate_from: {k!r}")
        yy, q = parse_q(qs)
        trunc[k] = yy * 10 + q
    return roles, trunc, cfg

CATEGORIES = {
    "storage_handling": r"product storage|incorrect (product )?storage|storage of drug|expired (product|drug)|temperature excursion",
    "nonadherence": r"dose omission|noncompliance|non-compliance|intentional underdose",
    "interruption_ambiguous": r"^therapy interrupted$|^therapy cessation$|administration interrupted|^treatment delayed$",
    "dosing_admin_error": r"wrong (technique in (product|drug)|dose|drug|product|patient|route|strength|schedule|rate|dosage)|incorrect (dose|route|drug|product|dosage)|^underdose$|accidental (overdose|underdose)|medication error|dispensing (error|issue)|prescribing error|preparation error|administration error|titration error|monitoring error|extra dose",
    "access_supply": r"availability issue|supply issue|insurance issue|afford",
    "off_label": r"^off label use$",
    "diet_noncompliance": r"^diet noncompliance$",
}
EXCLUDE_PT = r"intentional overdose|self-injur|suicid|intentional self|device"
ANY_ME = ["storage_handling", "nonadherence", "dosing_admin_error", "access_supply"]

# Health-professional occupation codes: the FDA documentation lists MD/PH/OT; HP and RN appear in the data but not in the current documentation, so they are treated as "other/unspecified health professional"
HCP_CODES = {"MD", "PH", "OT", "HP", "RN"}
# HP and RN occur in the data but are not defined in the FAERS documentation; a sensitivity analysis uses only these
HCP_DOCUMENTED = {"MD", "PH", "OT"}
AGE_FACTORS = {"YR": 1.0, "DEC": 10.0, "MON": 1 / 12, "WK": 1 / 52, "DY": 1 / 365, "HR": 1 / 8760}



def characteristic_rows(group, n_grp, sex, age_y, ctry, year, rept, occp):
    """Rows of 14_flow_and_characteristics.csv for one group (READUS-PV item 8b).

    sex, ctry, rept, occp: upper-case strings ("" when missing); age_y: years (NaN when missing);
    year: year of FDA receipt (NaN when missing)."""
    rows = []

    def add(char, level, mask):
        k = int(np.sum(mask))
        rows.append({"section": "characteristics", "group": group, "characteristic": char, "level": level,
                     "n": k, "pct": round(100 * k / n_grp, 1) if n_grp else np.nan})
    add("sex", "female", sex == "F"); add("sex", "male", sex == "M"); add("sex", "missing or unknown", ~np.isin(sex, ["F", "M"]))
    add("age", "<18 years", age_y < 18); add("age", "18-64 years", (age_y >= 18) & (age_y < 65))
    add("age", ">=65 years", age_y >= 65); add("age", "missing", np.isnan(age_y))
    add("country", "United States", ctry == "US"); add("country", "other", (ctry != "US") & (ctry != ""))
    add("country", "missing", ctry == "")
    add("year of FDA receipt", "<=2015", year <= 2015); add("year of FDA receipt", "2016-2018", (year >= 2016) & (year <= 2018))
    add("year of FDA receipt", "2019-2021", (year >= 2019) & (year <= 2021)); add("year of FDA receipt", ">=2022", year >= 2022)
    add("year of FDA receipt", "missing", np.isnan(year))
    add("report type", "expedited (EXP)", rept == "EXP"); add("report type", "periodic (PER)", rept == "PER")
    add("report type", "direct (DIR)", rept == "DIR"); add("report type", "other or missing", ~np.isin(rept, ["EXP", "PER", "DIR"]))
    add("reporter", "consumer (CN)", occp == "CN")
    add("reporter", "health professional, documented code (MD, PH, OT)", np.isin(occp, sorted(HCP_DOCUMENTED)))
    add("reporter", "undocumented code (HP, RN)", np.isin(occp, ["HP", "RN"]))
    add("reporter", "other or missing", ~np.isin(occp, ["CN", "HP", "RN"] + sorted(HCP_DOCUMENTED)))
    return rows

def load_smq_specs(specs):
    """--smq-file may be repeated as label=path; returns {smq_<label>: set(terms)}"""
    out = {}
    for spec in specs or []:
        label, path = spec.split("=", 1) if "=" in spec else ("medication_errors", spec)
        raw = [x.strip() for x in open(path, encoding="utf-8-sig").read().splitlines() if x.strip()]
        if raw and raw[0].lower().startswith("pt"):
            sdf = pd.read_csv(path, encoding="utf-8-sig", dtype=str)
            sdf.columns = [c.strip().lower() for c in sdf.columns]
            terms = {str(x).strip().lower() for x in sdf["pt"].dropna()}
        else:
            terms = {x.lower() for x in raw}
        out[f"smq_{label}"] = terms
    return out

INDICATIONS = {
    "penicillamine": r"hepato-?lenticular|wilson", "trientine": r"hepato-?lenticular|wilson",
    "zinc_acetate_galzin": r"hepato-?lenticular|wilson",
    "nitisinone": r"tyrosin(a)?emia|alkaptonuria|ochronosis",
    "cysteamine": r"cystinosis",
    "carglumic_acid": r"hyperammon|glutamate synthase|propionic|methylmalonic|urea cycle",
    "sapropterin": r"phenylketonuria|hyperphenylalanin|tetrahydrobiopterin|biopterin",
    "sepiapterin": r"phenylketonuria|hyperphenylalanin",
    "tiopronin": r"(?<!homo)cystinuria|cystine",
    "miglustat": r"gaucher|niemann|glycogen storage disease type ii\b|pompe",
    "eliglustat": r"gaucher",
    "sodium_phenylbutyrate": r"urea cycle|hyperammon|ornithine|citrullin|argininosuccin|carbamoyl|arginase|argininaemia",
    "glycerol_phenylbutyrate": r"urea cycle|hyperammon|ornithine|citrullin|argininosuccin|carbamoyl|arginase|argininaemia",
    "betaine_anhydrous": r"homocystin|methylenetetrahydrofolate reductase|cystathionine beta-synthase",
    # MedDRA codes Zellweger syndrome as "cerebrohepatorenal syndrome"; Refsum disease and adrenoleukodystrophy
    # are peroxisomal disorders, covered by the labelled indication
    "cholic_acid": r"bile acid|peroxisom|zellweger|cerebrohepatorenal|refsum|adrenoleukodystrophy",
    "example_nefecon": r"iga nephropathy|iga glomerulo|berger",
}

PARENT_RARE = {"penicillamine", "trientine_all", "zinc_acetate_galzin", "nitisinone_all", "cysteamine_oral",
               "carglumic_acid", "sapropterin", "sepiapterin", "tiopronin", "miglustat", "eliglustat",
               "sodium_phenylbutyrate", "glycerol_phenylbutyrate", "betaine_anhydrous", "cholic_acid"}

FILE_RE = re.compile(r"^(DEMO|DRUG|REAC|INDI)(\d{2})Q([1-4])(?:_NEW)?\.TXT$", re.I)
# the FDA has used four naming conventions for the deleted-case lists
DEL_RE = re.compile(r"^(?:DELE.*|ADR.*DELETEDCASES|ALLDELETEDCASES|\d{2}Q\dDELETEDCASES)\.TXT$", re.I)
UNION_RX = re.compile("|".join(f"(?:{r})" for _, r, _, _, _ in CANDIDATES))
FAR_FUTURE = np.int32(99_999_999)


def indication_pattern(key):
    for k, pat in INDICATIONS.items():
        if key == k or key.startswith(k):
            return pat
    return None


# ----------------------------------------------------------------------------
def find_sources(root):
    found = defaultdict(dict)
    deleted_files = []
    for d, _, files in os.walk(root):
        for f in files:
            full = os.path.join(d, f)
            if f.lower().endswith(".zip"):
                with zipfile.ZipFile(full) as zf:
                    for info in zf.infolist():
                        base = os.path.basename(info.filename)
                        m = FILE_RE.match(base)
                        if m:
                            found[(int(m.group(2)), int(m.group(3)))][m.group(1).upper()] = ("zip", (full, info.filename))
                        elif DEL_RE.match(base):
                            deleted_files.append(("zip", (full, info.filename)))
            else:
                m = FILE_RE.match(f)
                if m:
                    found[(int(m.group(2)), int(m.group(3)))][m.group(1).upper()] = ("file", full)
                elif DEL_RE.match(f):
                    deleted_files.append(("file", full))
    return found, deleted_files


@contextlib.contextmanager
def open_table(src):
    kind, loc = src
    if kind == "file":
        with open(loc, "rb") as f:
            yield f
    else:
        zpath, member = loc
        with zipfile.ZipFile(zpath) as zf, zf.open(member) as f:
            yield f


def read_deleted_caseids(deleted_files):
    ids = set()
    for src in deleted_files:
        with open_table(src) as fh:
            for line in fh:
                m = re.match(rb"\s*(\d+)", line)
                if m:
                    ids.add(int(m.group(1)))
    return ids


def iter_table(src, usecols, chunksize=250_000):
    with open_table(src) as f:
        probe = pd.read_csv(f, sep="$", dtype=str, encoding="latin-1",
                            quoting=csv.QUOTE_NONE, nrows=0)
    have = {c.strip().lower(): c for c in probe.columns}
    if "primaryid" not in have:
        return
    sel = [have[c] for c in usecols if c in have]
    with open_table(src) as f:
        for chunk in pd.read_csv(f, sep="$", dtype=str, encoding="latin-1",
                                 quoting=csv.QUOTE_NONE, on_bad_lines="skip",
                                 low_memory=False, usecols=sel, chunksize=chunksize):
            chunk.columns = [c.strip().lower() for c in chunk.columns]
            for c in usecols:
                if c not in chunk.columns:
                    chunk[c] = np.nan
            yield chunk


def parse_q(s):
    m = re.match(r"^(\d{4})Q([1-4])$", s.upper())
    if not m:
        sys.exit(f"bad quarter format: {s} (expected e.g. 2014Q1)")
    return (int(m.group(1)) % 100, int(m.group(2)))


def ror_ci(a, b, c, d):
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    ror = (a * d) / (b * c)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return ror, np.exp(np.log(ror) - 1.96 * se), np.exp(np.log(ror) + 1.96 * se)


def wilson(a, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = a / n
    den = 1 + z ** 2 / n
    centre = (p + z ** 2 / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def expit(x):
    return 1 / (1 + np.exp(-x))


def dl_pool(events, ns):
    """DerSimonian-Laird random-effects pooling (logit proportion). Returns (mu, se, tau2, I2, k)."""
    ev, n = np.asarray(events, float), np.asarray(ns, float)
    k = len(ev)
    if k == 0:
        return np.nan, np.nan, np.nan, np.nan, 0
    corr = (ev == 0) | (ev == n)
    a = np.where(corr, ev + 0.5, ev)
    b = np.where(corr, n - ev + 0.5, n - ev)
    y = np.log(a / b)
    v = 1 / a + 1 / b
    w = 1 / v
    mu_f = np.sum(w * y) / np.sum(w)
    q = np.sum(w * (y - mu_f) ** 2)
    df = k - 1
    cc = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (q - df) / cc) if (df > 0 and cc > 0) else 0.0
    ws = 1 / (v + tau2)
    mu = np.sum(ws * y) / np.sum(ws)
    se = np.sqrt(1 / np.sum(ws))
    i2 = max(0.0, (q - df) / q) * 100 if q > 0 else 0.0
    return mu, se, tau2, i2, k


def logit_cluster(X, y, groups, names, max_iter=50):
    """Newton-Raphson logistic regression with cluster-robust (CR1) standard errors."""
    beta = np.zeros(X.shape[1])
    for _ in range(max_iter):
        eta = X @ beta
        p = 1 / (1 + np.exp(-eta))
        W = p * (1 - p)
        H = X.T @ (X * W[:, None])
        g = X.T @ (y - p)
        step = np.linalg.solve(H, g)
        beta += step
        if np.max(np.abs(step)) < 1e-8:
            break
    p = 1 / (1 + np.exp(-(X @ beta)))
    Hinv = np.linalg.inv(X.T @ (X * (p * (1 - p))[:, None]))
    score = X * (y - p)[:, None]
    uniq, inv = np.unique(groups, return_inverse=True)
    G, N, K = len(uniq), X.shape[0], X.shape[1]
    S = np.zeros((G, K))
    np.add.at(S, inv, score)
    meat = S.T @ S
    adj = (G / (G - 1)) * ((N - 1) / (N - K))
    cov = adj * Hinv @ meat @ Hinv
    return beta, cov, G


GROUP_STYLE = {
    "rare": ("Rare-disease oral drugs", "#2a78d6", "o"),
    "specialty": ("Specialty-pharmacy comparators", "#eb6834", "s"),
    "common": ("Common chronic-disease comparators", "#1baf7a", "^"),
}
DISPLAY_NAME = {
    "trientine_all": "Trientine (all salts)", "nitisinone_all": "Nitisinone (all products)",
    "cysteamine_oral": "Cysteamine (oral)", "zinc_acetate_galzin": "Zinc acetate (Galzin)",
    "betaine_anhydrous": "Betaine anhydrous", "sodium_phenylbutyrate": "Sodium phenylbutyrate",
    "glycerol_phenylbutyrate": "Glycerol phenylbutyrate", "carglumic_acid": "Carglumic acid",
    "cholic_acid": "Cholic acid", "spec_dimethyl_fumarate": "Dimethyl fumarate",
    "s2_enzalutamide": "Enzalutamide", "s2_apalutamide": "Apalutamide", "s2_darolutamide": "Darolutamide",
    "s2_abemaciclib": "Abemaciclib", "s2c_palbociclib": "Palbociclib", "s2c_ribociclib": "Ribociclib",
}
LABEL_EN = {"nonadherence": "Reports with dose omission / nonadherence PTs (%)",
            "any_me": "Reports with any medication-error / adherence / access PT (%)"}


def forest_plot(df, re_list, cat, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, muted, grid = "#2b2b2a", "#6b6a64", "#e4e3de"
    rows = []
    for g in ("rare", "specialty", "common"):
        sub = df[df["group"] == g].sort_values("pct", ascending=False)
        if sub.empty:
            continue
        rows.append((GROUP_STYLE[g][0], None, None, None, g, "header", ""))
        for _, r in sub.iterrows():
            name = DISPLAY_NAME.get(r["drug"], r["drug"].replace("spec_", "").replace("ctrl_", "")
                                    .replace("_", " ").capitalize())
            rows.append((name, r["pct"], r["pct_CI_low"], r["pct_CI_high"], g, "drug",
                         f"{int(r['events']):,}/{int(r['n_reports']):,}"))
        pooled = [x for x in re_list if x["comparison"] == f"pooled within {g}"]
        if pooled:
            p = pooled[0]
            rows.append(("Random-effects pooled", p["estimate"], p["CI_low"], p["CI_high"], g, "pooled",
                         f"I² {p['I2(%)']:.0f}%"))
    n = len(rows)
    fig, ax = plt.subplots(figsize=(8.2, 0.32 * n + 1.4), dpi=200)
    xmax = max(r[3] for r in rows if r[3] is not None) * 1.08
    for i, (lab, p, lo, hi, g, kind, rt) in enumerate(rows):
        y = n - 1 - i
        col, mk = GROUP_STYLE[g][1], GROUP_STYLE[g][2]
        if kind == "header":
            ax.text(-0.01, y, lab, transform=ax.get_yaxis_transform(), ha="right", va="center",
                    fontsize=8.5, fontweight="bold", color=ink)
            continue
        ax.text(-0.01, y, lab, transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=8, color=ink if kind == "drug" else muted, style="italic" if kind == "pooled" else "normal")
        ax.text(1.01, y, rt, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=7.5, color=muted)
        if kind == "drug":
            ax.plot([lo, hi], [y, y], color=col, lw=1.6, solid_capstyle="round", zorder=2)
            ax.plot(p, y, marker=mk, ms=6, color=col, mec="white", mew=1.2, zorder=3)
        else:
            ax.fill([lo, p, hi, p], [y, y + 0.32, y, y - 0.32], color=col, zorder=3, lw=0)
    ax.set_ylim(-0.8, n - 0.2)
    ax.set_xlim(0, xmax)
    ax.set_yticks([])
    ax.set_xlabel(LABEL_EN.get(cat, cat), fontsize=8.5, color=ink)
    ax.tick_params(axis="x", labelsize=8, colors=muted, length=0)
    ax.grid(axis="x", color=grid, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(grid)
    ax.text(1.01, n - 0.2, "Events/N", transform=ax.get_yaxis_transform(), ha="left", va="bottom",
            fontsize=7.5, color=muted, fontweight="bold")
    h = fig.get_figheight()
    fig.text(0.04, 0.22 / h, "Bars: Wilson 95% CI (comparator CIs are narrower than the markers).\n"
             "Diamonds: DerSimonian–Laird random-effects pooled proportion across drugs.",
             fontsize=6.5, color=muted, linespacing=1.4)
    fig.subplots_adjust(left=0.36, right=0.84, top=1 - 0.3 / h, bottom=0.85 / h)
    fig.savefig(path, facecolor="white")
    plt.close(fig)


def qkey(k):
    return f"20{k[0]:02d}q{k[1]}"


def enc_codes(series, na="NA"):
    s = series.astype("object")
    miss = s.isna()
    s = s.where(~miss, other="").astype(str).str.upper()
    if na != "":
        s = s.where(~miss, other=na)
    v = s.to_numpy().astype(str)
    uniq, inv = np.unique(v, return_inverse=True)
    return inv.astype(np.int32), uniq


def shrink(a):
    if a.dtype != np.int32 or len(a) == 0:
        return a
    mx = int(a.max())
    if mx < 255:
        return a.astype(np.uint8)
    if mx < 32767:
        return a.astype(np.int16)
    return a


# ----------------------------------------------------------------------------
def phase_demo(files, quarters, work, **kw):
    d_dir = os.path.join(work, "demo5")
    os.makedirs(d_dir, exist_ok=True)
    cols = ["primaryid", "caseid", "fda_dt", "occp_cod", "rept_cod",
            "sex", "age", "age_cod", "occr_country"]
    for i, k in enumerate(quarters, 1):
        dst = os.path.join(d_dir, qkey(k) + ".npz")
        if os.path.exists(dst):
            continue
        acc = defaultdict(list)
        for df in iter_table(files[k]["DEMO"], cols, chunksize=500_000):
            pid = pd.to_numeric(df["primaryid"], errors="coerce").to_numpy(dtype="float64")
            cid = pd.to_numeric(df["caseid"], errors="coerce").to_numpy(dtype="float64")
            fda = pd.to_numeric(df["fda_dt"], errors="coerce").to_numpy(dtype="float64")
            m = ~np.isnan(pid) & ~np.isnan(cid)
            acc["pid"].append(pid[m].astype(np.int64))
            acc["cid"].append(cid[m].astype(np.int64))
            acc["fda"].append(np.where(np.isnan(fda[m]), FAR_FUTURE, fda[m]).astype(np.int32))
            acc["age"].append(pd.to_numeric(df["age"][m], errors="coerce").to_numpy(dtype="float32"))
            for c, na in (("occp_cod", "NA"), ("rept_cod", "NA"), ("sex", ""),
                          ("age_cod", ""), ("occr_country", "")):
                acc[c + "_s"].append(df[c][m])
        if not acc["pid"]:
            continue
        out = {"pid": np.concatenate(acc["pid"]), "cid": np.concatenate(acc["cid"]),
               "fda": np.concatenate(acc["fda"]), "age": np.concatenate(acc["age"])}
        for c, na in (("occp_cod", "NA"), ("rept_cod", "NA"), ("sex", ""),
                      ("age_cod", ""), ("occr_country", "")):
            ser = pd.concat(acc[c + "_s"], ignore_index=True) if acc[c + "_s"] else pd.Series(dtype=object)
            codes, uniq = enc_codes(ser, na)
            out[c + "_codes"], out[c + "_uniq"] = codes, uniq
        tmp = dst + ".part.npz"
        np.savez(tmp, **out)
        os.replace(tmp, dst)
        if i % 10 == 0 or i == len(quarters):
            print(f"  DEMO {i}/{len(quarters)}", flush=True)


def phase_reac(files, quarters, work, pt_dict=None, smq_files=None, **kw):
    if pt_dict:
        pt_map = _load_pt_dict(pt_dict)
        cat_names = sorted(set(pt_map.values()))
        print(f"using reviewed PT dictionary: {len(pt_map)} PTs", flush=True)
    else:
        pt_map = None
        cat_names = list(CATEGORIES)
    smq_sets = load_smq_specs(smq_files)
    for lab, terms in smq_sets.items():
        cat_names = cat_names + [lab]
        print(f"SMQ list {lab}: {len(terms)} terms", flush=True)
    cat_idx = {c: i for i, c in enumerate(cat_names)}

    cat_re = {c: re.compile(p) for c, p in CATEGORIES.items()}
    excl_re = re.compile(EXCLUDE_PT)
    cache = {}

    def pt_cats(pt):
        if pt in cache:
            return cache[pt]
        s = str(pt).strip().lower()
        if pt_map is not None:
            cats = (pt_map[s],) if s in pt_map else ()
        else:
            cats = () if excl_re.search(s) else tuple(c for c, r in cat_re.items() if r.search(s))
        for lab, terms in smq_sets.items():
            if s in terms:
                cats = cats + (lab,)
        cache[pt] = cats
        return cats

    r_dir = os.path.join(work, "reac5")
    os.makedirs(r_dir, exist_ok=True)
    for i, k in enumerate(quarters, 1):
        dst = os.path.join(r_dir, qkey(k) + ".npz")
        if os.path.exists(dst):
            continue
        hit_pid_parts, hit_pt_parts, pts = [], [], []
        pt_code = {}
        flag_parts = {c: [] for c in cat_names}
        for df in iter_table(files[k]["REAC"], ["primaryid", "pt"], chunksize=400_000):
            ids = pd.to_numeric(df["primaryid"], errors="coerce").to_numpy(dtype="float64")
            pt = df["pt"].to_numpy()
            ok = ~np.isnan(ids) & pd.notna(pt)
            ids = ids[ok].astype(np.int64)
            pt = pt[ok]
            if len(ids) == 0:
                continue
            pt = pd.Series(pt).str.strip().to_numpy().astype(str)
            uniq_pt, inv = np.unique(pt, return_inverse=True)
            cats_of = [pt_cats(p) for p in uniq_pt]
            keep_uniq = np.array([bool(c) for c in cats_of], dtype=bool)
            if not keep_uniq.any():
                continue
            row_keep = keep_uniq[inv]
            sub_ids, sub_pt = ids[row_keep], inv[row_keep]
            pair = sub_ids * np.int64(len(uniq_pt)) + sub_pt
            u = np.unique(pair)
            sub_ids = u // np.int64(len(uniq_pt))
            sub_pt = (u % np.int64(len(uniq_pt))).astype(np.int64)
            codes = np.empty(len(uniq_pt), dtype=np.int32)
            for j, p in enumerate(uniq_pt):
                c = pt_code.get(p)
                if c is None:
                    c = len(pts)
                    pt_code[p] = c
                    pts.append(p)
                codes[j] = c
            hit_pid_parts.append(sub_ids.astype(np.int64))
            hit_pt_parts.append(codes[sub_pt])
            cmat = np.zeros((len(uniq_pt), len(cat_names)), dtype=bool)
            for j, cs in enumerate(cats_of):
                for c in cs:
                    cmat[j, cat_idx[c]] = True
            for ci, c in enumerate(cat_names):
                sel = sub_ids[cmat[sub_pt, ci]]
                if len(sel):
                    flag_parts[c].append(np.unique(sel))
        tmp = dst + ".part.npz"
        payload = {"hit_pid": np.concatenate(hit_pid_parts) if hit_pid_parts else np.array([], np.int64),
                   "hit_pt": np.concatenate(hit_pt_parts) if hit_pt_parts else np.array([], np.int32),
                   "pts": np.array(pts), "cats": np.array(cat_names)}
        for j, c in enumerate(cat_names):
            payload[f"flagc{j}"] = (np.unique(np.concatenate(flag_parts[c]))
                                    if flag_parts[c] else np.array([], np.int64))
        np.savez(tmp, **payload)
        os.replace(tmp, dst)
        if i % 10 == 0 or i == len(quarters):
            print(f"  REAC {i}/{len(quarters)}", flush=True)


def phase_drug(files, quarters, work, **kw):
    cand = [(k, re.compile(r), re.compile(x) if x else None) for k, r, x, _, _ in CANDIDATES]
    g_dir = os.path.join(work, "drug5")
    os.makedirs(g_dir, exist_ok=True)
    seq_ok = True
    for i, k in enumerate(quarters, 1):
        dst = os.path.join(g_dir, qkey(k) + ".npz")
        if os.path.exists(dst):
            continue
        parts = {j: [] for j in range(len(cand))}
        for df in iter_table(files[k]["DRUG"],
                             ["primaryid", "drug_seq", "role_cod", "drugname", "prod_ai"],
                             chunksize=250_000):
            ps = (df["role_cod"].str.strip().str.upper() == "PS").to_numpy()
            if not ps.any():
                continue
            sub = df[ps]
            ids = pd.to_numeric(sub["primaryid"], errors="coerce").to_numpy(dtype="float64")
            seq = pd.to_numeric(sub["drug_seq"], errors="coerce").to_numpy(dtype="float64")
            ok = ~np.isnan(ids)
            ids = ids[ok].astype(np.int64)
            seqv = np.where(np.isnan(seq[ok]), -1, seq[ok]).astype(np.int32)
            if (seqv >= 4095).any():
                seq_ok = False
            name = sub["drugname"].to_numpy()[ok]
            ai = sub["prod_ai"].to_numpy()[ok]
            if len(ids) == 0:
                continue
            s = pd.Series(
                np.where(pd.isna(name), "", name).astype(str) + " | "
                + np.where(pd.isna(ai), "", ai).astype(str)).str.upper()
            any_hit = s.str.contains(UNION_RX, regex=True).to_numpy()
            if not any_hit.any():
                continue
            s2 = s[any_hit]
            ids2 = ids[any_hit]
            seq2 = seqv[any_hit]
            for j, (key, rx, ex) in enumerate(cand):
                m = s2.str.contains(rx, regex=True).to_numpy()
                if ex is not None:
                    m = m & ~s2.str.contains(ex, regex=True).to_numpy()
                if m.any():
                    parts[j].append(ids2[m] * np.int64(4096) + (seq2[m].astype(np.int64) + 1))
        tmp = dst + ".part.npz"
        payload = {"keys": np.array([k for k, _, _ in cand])}
        for j in range(len(cand)):
            payload[f"c{j}"] = np.unique(np.concatenate(parts[j])) if parts[j] else np.array([], np.int64)
        np.savez(tmp, **payload)
        os.replace(tmp, dst)
        if i % 5 == 0 or i == len(quarters):
            print(f"  DRUG {i}/{len(quarters)}", flush=True)
    if not seq_ok:
        print("  note: records with drug_seq >= 4095 present; the packed key may overflow", flush=True)


def phase_indi(files, quarters, work, **kw):
    """Keep only the INDI rows that belong to candidate drugs (rare/example with an indication pattern)."""
    g_dir = os.path.join(work, "drug5")
    ind_keys = [j for j, (k, _, _, g, _) in enumerate(CANDIDATES)
                if g in ("rare", "example") and indication_pattern(k)]
    keep_pids = []
    for k in quarters:
        p = os.path.join(g_dir, qkey(k) + ".npz")
        if not os.path.exists(p):
            sys.exit("DRUG phase not complete; cannot read the candidate drug list")
        z = np.load(p)
        for j in ind_keys:
            v = z[f"c{j}"]
            if len(v):
                keep_pids.append(np.unique(v // 4096))
    keep = np.unique(np.concatenate(keep_pids)) if keep_pids else np.array([], np.int64)
    print(f"  candidate drug reports: {len(keep):,}", flush=True)

    i_dir = os.path.join(work, "indi")
    os.makedirs(i_dir, exist_ok=True)
    for i, k in enumerate(quarters, 1):
        dst = os.path.join(i_dir, qkey(k) + ".npz")
        if os.path.exists(dst):
            continue
        parts = []
        if "INDI" in files[k]:
            for df in iter_table(files[k]["INDI"], ["primaryid", "indi_drug_seq", "indi_pt"], chunksize=300_000):
                ids = pd.to_numeric(df["primaryid"], errors="coerce").to_numpy(dtype="float64")
                ok = ~np.isnan(ids) & pd.notna(df["indi_pt"])
                ids = ids[ok].astype(np.int64)
                if len(ids) == 0:
                    continue
                pos = np.searchsorted(keep, ids)
                pos_c = np.minimum(pos, max(len(keep) - 1, 0))
                hit = (keep[pos_c] == ids) if len(keep) else np.zeros(len(ids), bool)
                if not hit.any():
                    continue
                seq = pd.to_numeric(df["indi_drug_seq"], errors="coerce").to_numpy(dtype="float64")[ok][hit]
                pt = df["indi_pt"].to_numpy()[ok][hit]
                part = pd.DataFrame({
                    "pid": ids[hit],
                    "seq": np.where(np.isnan(seq), -1, seq).astype(np.int32),
                    "pt": pd.Series(pt).str.strip().str.lower().to_numpy().astype(str),
                }).drop_duplicates()
                parts.append(part)
        if parts:
            allp = pd.concat(parts, ignore_index=True).drop_duplicates()
            payload = {"pid": allp["pid"].to_numpy(np.int64), "seq": allp["seq"].to_numpy(np.int32),
                       "pt": allp["pt"].to_numpy().astype(str)}
        else:
            payload = {"pid": np.array([], np.int64), "seq": np.array([], np.int32),
                       "pt": np.array([], dtype=str)}
        tmp = dst + ".part.npz"
        np.savez(tmp, **payload)
        os.replace(tmp, dst)
        if i % 10 == 0 or i == len(quarters):
            print(f"  INDI {i}/{len(quarters)}", flush=True)


def _load_pt_dict(path):
    pdict = pd.read_csv(path, encoding="utf-8-sig", dtype=str)
    pdict.columns = [c.strip().lower() for c in pdict.columns]
    if not {"keep", "pt", "category"} <= set(pdict.columns):
        sys.exit(f"--pt-dict needs the columns pt / keep / category; found: {list(pdict.columns)}")
    pdict = pdict[pdict["keep"].str.strip().str.upper() == "Y"]
    out = {}
    for p, c in zip(pdict["pt"], pdict["category"]):
        if pd.isna(p) or pd.isna(c):
            continue
        out[p.strip().lower()] = c.strip()
    return out


# ----------------------------------------------------------------------------
def phase_agg(files, quarters, work, out, pt_dict=None, smq_files=None,
              deleted_files=None, min_drug_n=100, groups_cfg=None, **kw):
    qs = [qkey(k) for k in quarters]
    d_dir, r_dir, g_dir, i_dir = (os.path.join(work, x) for x in ("demo5", "reac5", "drug5", "indi"))
    for nm, dd in (("DEMO", d_dir), ("REAC", r_dir), ("DRUG", g_dir)):
        miss = [q for q in qs if not os.path.exists(os.path.join(dd, q + ".npz"))]
        if miss:
            sys.exit(f"{nm} phase not complete; {len(miss)} quarters missing")
    load_groups_config(groups_cfg)  # validate the --groups file before the long reads
    want_keys = [k for k, _, _, _, _ in CANDIDATES]
    for q in qs:
        z = np.load(os.path.join(g_dir, q + ".npz"))
        if "keys" not in z.files or [str(x) for x in z["keys"]] != want_keys:
            sys.exit(f"the DRUG cache for {q} was built for a different candidate list (before version 1.1.0); "
                     f"delete {g_dir} and re-run --phase drug")

    if pt_dict:
        pt_map = _load_pt_dict(pt_dict)
        cat_names = sorted(set(pt_map.values()))
    else:
        pt_map = None
        cat_names = list(CATEGORIES)
    want_cats = list(cat_names)
    smq_sets = load_smq_specs(smq_files)
    for lab in smq_sets:
        want_cats = want_cats + [lab]
        print(f"SMQ category {lab}: {len(smq_sets[lab])} terms", flush=True)

    # ---- DEMO ----
    cols = ["pid", "cid", "fda", "age"]
    keys = ["occp_cod", "rept_cod", "sex", "age_cod", "occr_country"]
    acc = {c: [] for c in cols + keys + ["fq"]}
    maps = {c: ({}, []) for c in keys}
    for q in qs:
        z = np.load(os.path.join(d_dir, q + ".npz"))
        acc["fq"].append(np.full(len(z["pid"]), int(q[2:4]) * 10 + int(q[-1]), dtype=np.int16))
        for c in cols:
            acc[c].append(z[c])
        for c in keys:
            codes, uniq = z[c + "_codes"], z[c + "_uniq"]
            mp, lst = maps[c]
            lut = np.empty(len(uniq), np.int32)
            for j, s in enumerate(uniq):
                key = str(s)
                code = mp.get(key)
                if code is None:
                    code = len(lst)
                    mp[key] = code
                    lst.append(key)
                lut[j] = code
            acc[c].append(lut[codes])
    pid = np.concatenate(acc["pid"]); cid = np.concatenate(acc["cid"])
    fda = np.concatenate(acc["fda"]); age = np.concatenate(acc["age"])
    fq = np.concatenate(acc["fq"])
    occ = shrink(np.concatenate(acc["occp_cod"])); rep = shrink(np.concatenate(acc["rept_cod"]))
    sex = shrink(np.concatenate(acc["sex"])); acod = shrink(np.concatenate(acc["age_cod"]))
    ctry = shrink(np.concatenate(acc["occr_country"]))
    acc = None
    n_read = len(pid)
    print(f"report rows read: {len(pid):,}", flush=True)

    # ---- deleted-case lists ----
    deleted = read_deleted_caseids(deleted_files or [])
    if deleted:
        del_arr = np.fromiter(deleted, dtype=np.int64, count=len(deleted))
        keep_rows = ~np.isin(cid, del_arr)
        n_del = int((~keep_rows).sum())
        pid, cid, fda, age, fq = pid[keep_rows], cid[keep_rows], fda[keep_rows], age[keep_rows], fq[keep_rows]
        occ, rep, sex, acod, ctry = occ[keep_rows], rep[keep_rows], sex[keep_rows], acod[keep_rows], ctry[keep_rows]
        print(f"deleted-case list: {len(deleted):,} caseids, {n_del:,} records removed", flush=True)
    else:
        n_del = 0
        print("no deleted-case list (DELE*.txt) found; nothing removed", flush=True)

    if int(cid.max()) < 2 ** 31:
        cid = cid.astype(np.int32)
    bits = max(32, int(pid.max()).bit_length())
    sort_key = fda.astype(np.int64) * np.int64(1 << bits) + pid
    order = np.lexsort((sort_key, cid))
    del sort_key
    sc = cid[order]
    last = np.empty(len(sc), dtype=bool); last[-1] = True
    last[:-1] = sc[1:] != sc[:-1]
    del sc
    keep_idx = order[last]
    del order, last, cid
    keep_pid = pid[keep_idx].copy(); keep_fda = fda[keep_idx].copy()
    keep_age = age[keep_idx].copy(); keep_fq = fq[keep_idx].copy()
    keep_occ = occ[keep_idx].copy(); keep_rep = rep[keep_idx].copy()
    keep_sex = sex[keep_idx].copy(); keep_acod = acod[keep_idx].copy()
    keep_ctry = ctry[keep_idx].copy()
    del pid, fda, age, fq, occ, rep, sex, acod, ctry, keep_idx
    srt = np.argsort(keep_pid, kind="stable")
    keep_pid, keep_fda, keep_age = keep_pid[srt], keep_fda[srt], keep_age[srt]
    keep_fq = keep_fq[srt]
    keep_occ, keep_rep = keep_occ[srt], keep_rep[srt]
    keep_sex, keep_acod, keep_ctry = keep_sex[srt], keep_acod[srt], keep_ctry[srt]
    del srt
    n_total = len(keep_pid)
    print(f"deduplicated reports: {n_total:,}", flush=True)

    def lookup(ids):
        pos = np.searchsorted(keep_pid, ids)
        pos_c = np.minimum(pos, n_total - 1)
        return pos, keep_pid[pos_c] == ids

    # ---- REAC ----
    hit_pid_parts, hit_pt_parts = [], []
    flag_parts = {c: [] for c in want_cats}
    pts_global, pt_code = [], {}
    for q in qs:
        z = np.load(os.path.join(r_dir, q + ".npz"), allow_pickle=False)
        cats_q = [str(x) for x in z["cats"]]
        missing = [c for c in want_cats if c not in cats_q]
        if missing:
            sys.exit(f"{q} lacks category {missing}; re-run the reac phase with the corresponding --smq-file")
        pts = z["pts"]
        if len(pts):
            lut = np.empty(len(pts), np.int32)
            for j, s in enumerate(pts):
                key = str(s)
                c = pt_code.get(key)
                if c is None:
                    c = len(pts_global); pt_code[key] = c; pts_global.append(key)
                lut[j] = c
            hit_pid_parts.append(z["hit_pid"])
            hit_pt_parts.append(lut[z["hit_pt"]])
        for c in want_cats:
            v = z[f"flagc{cats_q.index(c)}"]
            if len(v):
                flag_parts[c].append(v)
    hit_pid = np.concatenate(hit_pid_parts) if hit_pid_parts else np.array([], np.int64)
    hit_pt = np.concatenate(hit_pt_parts) if hit_pt_parts else np.array([], np.int32)
    hit_pid_parts = hit_pt_parts = None
    pos, keepm = lookup(hit_pid)
    hit_pid, hit_pt = hit_pid[keepm], hit_pt[keepm]
    if len(hit_pid):
        MUL = np.int64(len(pts_global) + 1)
        pair = np.unique(hit_pid * MUL + hit_pt)
        hit_pid = pair // MUL
        hit_pt = (pair % MUL).astype(np.int32)
        del pair
    flag_arr = {}
    for c in want_cats:
        v = np.unique(np.concatenate(flag_parts[c])) if flag_parts[c] else np.array([], np.int64)
        if len(v):
            v = v[lookup(v)[1]]
        flag_arr[c] = v
    flag_parts = None
    present_any = [flag_arr[c] for c in ANY_ME if c in flag_arr]
    flag_arr["any_me"] = np.unique(np.concatenate(present_any)) if present_any else np.array([], np.int64)
    print("event reports in the database: " + ", ".join(f"{c}={len(v):,}" for c, v in flag_arr.items()), flush=True)

    all_cats = want_cats + ["any_me"]

    # ---- DRUG ----
    cand_pairs = {}
    for j in range(len(CANDIDATES)):
        parts = []
        for q in qs:
            v = np.load(os.path.join(g_dir, q + ".npz"))[f"c{j}"]
            if len(v):
                parts.append(v)
        cand_pairs[j] = np.unique(np.concatenate(parts)) if parts else np.array([], np.int64)
    cand_ids = {j: np.unique(v // 4096) for j, v in cand_pairs.items()}
    for j in cand_ids:
        v = cand_ids[j]
        if len(v):
            cand_ids[j] = v[lookup(v)[1]]

    # ---- INDI: indication restriction ----
    ind_keys = [j for j, (k, _, _, g, _) in enumerate(CANDIDATES)
                if g in ("rare", "example") and indication_pattern(k)]
    indi_parts = []
    for q in qs:
        p = os.path.join(i_dir, q + ".npz")
        if not os.path.exists(p):
            sys.exit(f"INDI phase not complete; missing {q}")
        z = np.load(p, allow_pickle=False)
        if len(z["pid"]):
            indi_parts.append(pd.DataFrame({"pid": z["pid"], "seq": z["seq"], "pt": z["pt"]}))
    indi = pd.concat(indi_parts, ignore_index=True) if indi_parts else pd.DataFrame(
        columns=["pid", "seq", "pt"])
    indi = indi[lookup(indi["pid"].to_numpy(np.int64))[1]] if len(indi) else indi
    ind_ok = {}
    indi_top_rows = []
    for j in ind_keys:
        key = CANDIDATES[j][0]
        pairs = cand_pairs[j]
        if len(pairs) == 0:
            ind_ok[key] = np.array([], np.int64)
            continue
        dfp = pd.DataFrame({"pid": pairs // 4096, "seq": (pairs % 4096) - 1})
        dfp = dfp[dfp["seq"] >= 0]
        if len(dfp) == 0 or len(indi) == 0:
            ind_ok[key] = np.array([], np.int64)
            continue
        mg = dfp.merge(indi, on=["pid", "seq"], how="inner")
        pat = re.compile(indication_pattern(key))
        sel = mg.loc[mg["pt"].map(lambda t: bool(pat.search(t))), "pid"]
        ind_ok[key] = np.unique(sel.to_numpy(np.int64))
        top = (mg.drop_duplicates(["pid", "pt"]).groupby("pt")["pid"].nunique()
               .sort_values(ascending=False).head(15))
        for t, cnt in top.items():
            indi_top_rows.append({"drug": key, "indi_pt": t, "n_reports": int(cnt),
                                  "matches_target": bool(pat.search(t))})
    if not indi_parts:
        print("  no INDI file found; skipping the indication-restricted analysis", flush=True)

    # ---- aggregation ----
    roles, trunc, cfg = load_groups_config(groups_cfg)
    orig_group = {k: g for k, _, _, g, _ in CANDIDATES}
    groups = {}
    for j, (key, _, _, g, h) in enumerate(CANDIDATES):
        if roles[key] == "off":
            continue
        ids = cand_ids[j]
        if key in trunc and len(ids):
            ids = ids[keep_fq[lookup(ids)[0]] < trunc[key]]
        groups[key] = (roles[key], h, ids)
    if groups_cfg:
        import json
        with open(os.path.join(out, "config_used.json"), "w", encoding="utf-8") as f:
            json.dump({"file": os.path.abspath(groups_cfg), **cfg,
                       "resolved_roles": {k: r for k, r in roles.items() if r != "off"},
                       "n_reports_after_truncation": {k: int(len(v[2])) for k, v in groups.items()
                                                      if k in trunc}}, f, indent=2)
        print("groups: " + ", ".join(f"{r}={sum(1 for v in groups.values() if v[0] == r)}" for r in ROLES[:4]),
              flush=True)
    rare_ids = np.unique(np.concatenate([v[2] for v in groups.values() if v[0] == "rare"]))
    ctrl_ids = np.unique(np.concatenate([v[2] for v in groups.values() if v[0] == "common"]))
    spec_ids = np.unique(np.concatenate([v[2] for v in groups.values() if v[0] == "specialty"]))
    groups["POOLED_rare_all"] = ("rare", "all rare-disease drugs pooled", rare_ids)
    groups["POOLED_controls"] = ("common", "all comparators pooled", ctrl_ids)
    groups["POOLED_specialty"] = ("specialty", "all specialty comparators pooled", spec_ids)
    n_ctrl = len(ctrl_ids)
    n_spec = len(spec_ids)

    def in_flag(ids, flag):
        if len(flag) == 0:
            return np.zeros(len(ids), dtype=bool)
        p = np.searchsorted(flag, ids)
        p_c = np.minimum(p, len(flag) - 1)
        return flag[p_c] == ids

    n_pt = len(pts_global)
    dict_counts = np.bincount(hit_pt, minlength=n_pt) if len(hit_pid) else np.zeros(n_pt, np.int64)
    code_cn = maps["occp_cod"][0].get("CN", -1)

    rows, pt_rows, rep_rows, yr_rows = [], [], [], []
    for key, (g, h, ids) in groups.items():
        n = int(len(ids))
        rec = {"drug": key, "group": g, "handling(to verify)": h, "n_reports_PS": n}
        if n == 0:
            rows.append(rec)
            continue
        pos, _ = lookup(ids)
        occ_v = keep_occ[pos]
        rep_v = keep_rep[pos]
        rec["pct_consumer_reported"] = round(100 * float(np.mean(occ_v == code_cn)), 1)
        for c in all_cats:
            is_c = in_flag(ids, flag_arr[c])
            a = int(is_c.sum())
            b = n - a
            cc = int(len(flag_arr[c]) - a)
            dd = (n_total - n) - cc
            r, lo, hi = ror_ci(a, b, cc, dd)
            rec[f"n_{c}"] = a
            rec[f"pct_{c}"] = round(100 * a / n, 1)
            rec[f"ROR_vs_allFAERS_{c}"] = f"{r:.2f} ({lo:.2f}-{hi:.2f})"
            if g == "rare" and n_ctrl > 0:
                c1 = int(in_flag(ctrl_ids, flag_arr[c]).sum())
                r2, lo2, hi2 = ror_ci(a, b, c1, n_ctrl - c1)
                rec[f"ROR_vs_controls_{c}"] = f"{r2:.2f} ({lo2:.2f}-{hi2:.2f})"
            if g == "rare" and n_spec > 0:
                c2 = int(in_flag(spec_ids, flag_arr[c]).sum())
                r3, lo3, hi3 = ror_ci(a, b, c2, n_spec - c2)
                rec[f"ROR_vs_specialty_{c}"] = f"{r3:.2f} ({lo3:.2f}-{hi3:.2f})"
        is_me = in_flag(ids, flag_arr["any_me"])
        is_na = in_flag(ids, flag_arr["nonadherence"]) if "nonadherence" in flag_arr else np.zeros(n, dtype=bool)
        if is_me.any():
            rec["pct_consumer_among_me"] = round(100 * float(np.mean(occ_v[is_me] == code_cn)), 1)
        rows.append(rec)

        p = np.searchsorted(ids, hit_pid)
        p_c = np.minimum(p, n - 1)
        sel = ids[p_c] == hit_pid
        if sel.any():
            cnt = np.bincount(hit_pt[sel], minlength=n_pt)
            for code in sorted(np.flatnonzero(cnt), key=lambda z: (-cnt[z], pts_global[z])):
                pt_rows.append({"drug": key, "pt": pts_global[code],
                                "category": "|".join(_cats5(pts_global[code], pt_map, smq_sets)),
                                "n_reports": int(cnt[code])})

        for col, vals, names in (("occp_cod", occ_v, maps["occp_cod"][1]),
                                 ("rept_cod", rep_v, maps["rept_cod"][1])):
            sizes = np.bincount(vals, minlength=len(names))
            sums = np.bincount(vals, weights=is_me, minlength=len(names))
            sums_na = np.bincount(vals, weights=is_na, minlength=len(names))
            for code in sorted(range(len(names)), key=lambda z: names[z]):
                if sizes[code] == 0:
                    continue
                rep_rows.append({"drug": key, "field": col, "level": names[code],
                                 "n_reports": int(sizes[code]), "n_any_me": int(round(sums[code])),
                                 "n_nonadherence": int(round(sums_na[code]))})

        year_v = np.where(keep_fda[pos] == FAR_FUTURE, -1, keep_fda[pos] // 10000).astype(np.int32)
        ok_y = year_v >= 0
        if ok_y.any():
            ys = year_v[ok_y]
            my = int(ys.max())
            sizes = np.bincount(ys, minlength=my + 1)
            masks = {c: in_flag(ids, flag_arr[c]) for c in all_cats}
            sums = {c: np.bincount(ys, weights=masks[c][ok_y], minlength=my + 1) for c in all_cats}
            for y in range(len(sizes)):
                if sizes[y] == 0:
                    continue
                yr_rows.append({"drug": key, "year": y, "n_reports": int(sizes[y]),
                                **{f"n_{c}": int(round(sums[c][y])) for c in all_cats}})

    # ---- sensitivity analyses ----
    parents = [k for k, (g, h, ids) in groups.items()
               if g == "rare" and not k.startswith("POOLED")
               and (k in PARENT_RARE or orig_group.get(k) != "rare")]
    hcp = HCP_CODES
    sens_cats = [c for c in ["nonadherence", "any_me", "dosing_admin_error", "storage_handling",
                             "access_supply"] + sorted(smq_sets) if c in flag_arr]

    def subset_mask(ids, how):
        if how == "all":
            return np.ones(len(ids), dtype=bool)
        pos, _ = lookup(ids)
        rep_v, occ_v = keep_rep[pos], keep_occ[pos]
        if how == "exclude_PER":
            return rep_v != maps["rept_cod"][0].get("PER", -1)
        if how == "HCP_only":
            keep = [maps["occp_cod"][0][c] for c in hcp if c in maps["occp_cod"][0]]
            return np.isin(occ_v, keep) if keep else np.zeros(len(ids), dtype=bool)
        if how == "HCP_strict":
            keep = [maps["occp_cod"][0][c] for c in HCP_DOCUMENTED if c in maps["occp_cod"][0]]
            return np.isin(occ_v, keep) if keep else np.zeros(len(ids), dtype=bool)
        if how == "EXP_and_HCP":
            keep = [maps["occp_cod"][0][c] for c in hcp if c in maps["occp_cod"][0]]
            base = rep_v == maps["rept_cod"][0].get("EXP", -1)
            return base & (np.isin(occ_v, keep) if keep else np.zeros(len(ids), dtype=bool))
        raise ValueError(how)

    def peak_year_ids(keys):
        """Reports from the calendar year with the most nonadherence reports, for each drug in keys."""
        if "nonadherence" not in flag_arr:
            return np.array([], dtype=np.int64)
        parts = []
        for k in keys:
            ids = groups[k][2]
            is_na_f = in_flag(ids, flag_arr["nonadherence"])
            if not is_na_f.any():
                continue
            pos, _ = lookup(ids)
            yv = np.where(keep_fda[pos] == FAR_FUTURE, -1, keep_fda[pos] // 10000).astype(np.int32)
            yn = yv[is_na_f]
            yn = yn[yn >= 0]
            if len(yn) == 0:
                continue
            top = int(np.bincount(yn).argmax())
            parts.append(ids[yv == top])
        return np.unique(np.concatenate(parts)) if parts else np.array([], dtype=np.int64)

    burst_ids = peak_year_ids(parents)
    comparator_keys = [k for k, (g, h, ids) in groups.items()
                       if g in ("specialty", "common") and not k.startswith("POOLED")]
    burst_ref_ids = peak_year_ids(comparator_keys)

    scenarios = [("all", "all reports", rare_ids),
                 ("exclude_PER", "excluding periodic reports (PER)", rare_ids),
                 ("HCP_only", "health-professional reports only", rare_ids),
                 ("EXP_and_HCP", "expedited reports from health professionals only", rare_ids),
                 ("HCP_strict", "health-professional reports only, documented codes (MD, PH, OT)", rare_ids),
                 ("all", "excluding the peak nonadherence year of each drug", np.setdiff1d(rare_ids, burst_ids)),
                 ("all", "excluding the peak nonadherence year of every drug in all groups",
                  np.setdiff1d(rare_ids, burst_ids), burst_ref_ids)]
    restricted_ids = np.unique(np.concatenate([ind_ok[k] for k in parents if len(ind_ok.get(k, []))])) \
        if any(len(ind_ok.get(k, [])) for k in parents) else np.array([], np.int64)
    if len(restricted_ids):
        scenarios.append(("all", "restricted to the rare-disease indication", restricted_ids))
    for k in parents:
        scenarios.append(("all", f"leave-one-out: removing {k}", np.setdiff1d(rare_ids, groups[k][2])))

    sens_rows = []
    for sc in scenarios:
        how, label, r_ids = sc[:3]
        ref_excl = sc[3] if len(sc) > 3 else None
        rs = r_ids[subset_mask(r_ids, how)]
        for ref_name, ref_ids in (("controls", ctrl_ids), ("specialty", spec_ids)):
            if len(ref_ids) == 0:
                continue
            fs = ref_ids[subset_mask(ref_ids, how)]
            if ref_excl is not None:
                fs = np.setdiff1d(fs, ref_excl)
            for c in sens_cats:
                a = int(in_flag(rs, flag_arr[c]).sum()); n1 = int(len(rs))
                c1 = int(in_flag(fs, flag_arr[c]).sum()); n2 = int(len(fs))
                if n1 == 0 or n2 == 0:
                    continue
                r, lo, hi = ror_ci(a, n1 - a, c1, n2 - c1)
                sens_rows.append({"scenario": label, "reference": ref_name, "category": c,
                                  "rare_n": n1, "rare_events": a, "rare_pct": round(100 * a / n1, 2),
                                  "ref_n": n2, "ref_events": c1, "ref_pct": round(100 * c1 / n2, 2),
                                  "ROR": round(r, 2), "CI_low": round(lo, 2), "CI_high": round(hi, 2)})

    out_fn = lambda df, name: df.to_csv(os.path.join(out, name), index=False, encoding="utf-8-sig")
    out_fn(pd.DataFrame(sens_rows), "6_sensitivity.csv")

    # ---- indication coverage ----
    cov_rows = []
    for k in parents + [k for k, (g, h, ids) in groups.items() if g == "example"]:
        n = len(groups[k][2])
        if n and k in ind_ok:
            cov_rows.append({"drug": k, "n_reports_PS": n, "n_target_indication": len(ind_ok[k]),
                             "pct_target_indication": round(100 * len(ind_ok[k]) / n, 1)})
    out_fn(pd.DataFrame(cov_rows), "10_indication_coverage.csv")
    out_fn(pd.DataFrame(indi_top_rows), "13_indication_top.csv")

    # ---- worked example: reproduce on the original period and PT first, then change the comparator ----
    pdoi_code = next((c for p_, c in pt_code.items()
                      if p_.strip().lower() == "product dose omission issue"), None)
    pdoi_ids = hit_pid[hit_pt == pdoi_code] if pdoi_code is not None else np.array([], np.int64)
    # every period restricts the example drug and all three references to the same quarters
    ex_periods = (("2021Q4-2023Q4 (original study period)", keep_pid[(keep_fq >= 214) & (keep_fq <= 234)]),
                  ("2021Q4-2025Q4 (all quarters since approval)", keep_pid[(keep_fq >= 214) & (keep_fq <= 254)]),
                  ("2024Q1-2025Q4 (after the original study period)", keep_pid[(keep_fq >= 241) & (keep_fq <= 254)]))
    ex_outcomes = [("PT: product dose omission issue", pdoi_ids)] + \
                  [(c, flag_arr[c]) for c in ("nonadherence", "any_me") if c in flag_arr]
    ex_rows = []
    for k, (g, h, ids) in groups.items():
        if g != "example" or not len(ids):
            continue
        for period, period_ids in ex_periods:
            for label, eids0 in (("all reports", ids),
                                 ("restricted to target indication", ind_ok.get(k, np.array([], np.int64)))):
                eids = np.intersect1d(eids0, period_ids, assume_unique=True)
                n1 = int(len(eids))
                if n1 == 0:
                    continue
                for oname, oset in ex_outcomes:
                    a = int(in_flag(eids, oset).sum())
                    refs = (("all other reports in the database", np.setdiff1d(period_ids, eids, assume_unique=True)),
                            ("common chronic-disease comparators", np.intersect1d(ctrl_ids, period_ids, assume_unique=True)),
                            ("specialty-pharmacy comparators", np.intersect1d(spec_ids, period_ids, assume_unique=True)))
                    for ref_name, ref in refs:
                        n2 = int(len(ref)); c1 = int(in_flag(ref, oset).sum())
                        r, lo, hi = ror_ci(a, n1 - a, c1, n2 - c1)
                        ex_rows.append({"drug": k, "period": period, "subset": label, "outcome": oname,
                                        "reference": ref_name, "n_reports": n1, "events": a,
                                        "pct": round(100 * a / n1, 2), "ref_n": n2, "ref_events": c1,
                                        "ROR": round(r, 2), "CI_low": round(lo, 2), "CI_high": round(hi, 2)})
    out_fn(pd.DataFrame(ex_rows), "11_worked_example.csv")

    # ---- multivariable logistic regression ----
    unit_map = {}
    for g_name, keys in (("rare", parents),
                         ("specialty", [k for k, (g, h, ids) in groups.items()
                                        if g == "specialty" and not k.startswith("POOLED")]),
                         ("common", [k for k, (g, h, ids) in groups.items()
                                     if g == "common" and not k.startswith("POOLED")])):
        for k in keys:
            for pidv in groups[k][2]:
                unit_map.setdefault(int(pidv), (g_name, k))
    mv_rows = []
    if unit_map:
        mv_ids = np.fromiter(unit_map.keys(), dtype=np.int64, count=len(unit_map))
        pos, _ = lookup(mv_ids)
        grp = np.array([unit_map[int(i)][0] for i in mv_ids])
        drg = np.array([unit_map[int(i)][1] for i in mv_ids])
        repc = keep_rep[pos]
        occc = keep_occ[pos]
        sexc = keep_sex[pos]
        ctryc = keep_ctry[pos]
        acc_c = keep_acod[pos]
        age_v = keep_age[pos].astype(np.float64)
        year_v = np.where(keep_fda[pos] == FAR_FUTURE, np.nan,
                          (keep_fda[pos] // 10000).astype(np.float64))
        def code_of(mp_name, val):
            return maps[mp_name][0].get(val, -1)
        fac_map = {}
        for unit_name, factor in (("YR", 1.0), ("DEC", 10.0), ("MON", 1 / 12), ("WK", 1 / 52),
                                  ("DY", 1 / 365), ("HR", 1 / 8760)):
            fac_map[code_of("age_cod", unit_name)] = factor
        fac = np.full(len(mv_ids), np.nan)
        for code, factor in fac_map.items():
            if code >= 0:
                fac[acc_c == code] = factor
        age_y = age_v * fac
        hcp_codes = {code_of("occp_cod", x) for x in ["CN"] + sorted(HCP_CODES)} - {-1}
        cols_x = {
            "intercept": np.ones(len(mv_ids)),
            "rare (vs specialty)": (grp == "rare").astype(float),
            "common (vs specialty)": (grp == "common").astype(float),
            "periodic report (vs expedited)": (repc == code_of("rept_cod", "PER")).astype(float),
            "direct report (vs expedited)": (repc == code_of("rept_cod", "DIR")).astype(float),
            "consumer (vs health professional)": (occc == code_of("occp_cod", "CN")).astype(float),
            "other/missing reporter": (~np.isin(occc, list(hcp_codes))).astype(float),
            "male (vs female)": (sexc == code_of("sex", "M")).astype(float),
            "sex missing": (~np.isin(sexc, [x for x in (code_of("sex", "M"), code_of("sex", "F")) if x >= 0])).astype(float),
            "age <18 (vs 18-64)": (age_y < 18).astype(float),
            "age >=65": (age_y >= 65).astype(float),
            "age missing": np.isnan(age_y).astype(float),
            "non-US country (vs US)": ((ctryc != code_of("occr_country", "US")) &
                                       (ctryc != code_of("occr_country", ""))).astype(float),
            "country missing": (ctryc == code_of("occr_country", "")).astype(float),
        }
        base_names = list(cols_x)
        X_base = np.column_stack([np.asarray(cols_x[n], dtype=float) for n in base_names])
        YEAR_BINS = [(0, 2015), (2016, 2018), (2019, 2021), (2022, 9999)]

        def year_design(yv):
            """Group years; groups with fewer than 5 events are merged with the adjacent group to avoid complete separation (rare outcomes such as access_supply)."""
            bins = []
            for lo, hi in YEAR_BINS:
                ev = yv[(year_v >= lo) & (year_v <= hi)].sum()
                if bins and bins[-1][2] < 5:
                    bins[-1] = (bins[-1][0], hi, bins[-1][2] + ev)
                else:
                    bins.append((lo, hi, ev))
            if len(bins) > 1 and bins[-1][2] < 5:
                last = bins.pop()
                bins[-1] = (bins[-1][0], last[1], bins[-1][2] + last[2])
            lab = lambda lo, hi: (f"<={hi}" if lo == 0 else (f">={lo}" if hi == 9999 else f"{lo}-{hi}"))
            ref = lab(bins[0][0], bins[0][1])
            cols_y, names_y = [], []
            for lo, hi, _ in bins[1:]:
                cols_y.append(((year_v >= lo) & (year_v <= hi)).astype(float))
                names_y.append(f"year {lab(lo, hi)} (vs {ref})")
            return cols_y, names_y

        for c in [x for x in ("nonadherence", "any_me", "access_supply") if x in flag_arr]:
            yv = in_flag(mv_ids, flag_arr[c]).astype(float)
            if yv.sum() < 10:
                continue
            cols_y, names_y = year_design(yv)
            X = np.column_stack([X_base] + cols_y) if cols_y else X_base
            names = base_names + names_y
            keep_cols = [i for i in range(X.shape[1]) if i < 3 or X[:, i].std() > 0]
            X = X[:, keep_cols]; names = [names[i] for i in keep_cols]
            for adjusted in (False, True):
                idx = list(range(3)) if not adjusted else list(range(len(names)))
                beta, cov, G = logit_cluster(X[:, idx], yv, drg, [names[i] for i in idx])
                ir, ic = 1, 2
                se = np.sqrt(np.diag(cov))
                if np.any(np.abs(beta[1:]) > 12) or np.any(~np.isfinite(se)):
                    print(f"  Warning: the {'adjusted' if adjusted else 'unadjusted'} model for {c} may be subject to complete separation; do not interpret its coefficients", flush=True)
                for lab, est, var in (("rare vs specialty", beta[ir], cov[ir, ir]),
                                      ("rare vs common", beta[ir] - beta[ic],
                                       cov[ir, ir] + cov[ic, ic] - 2 * cov[ir, ic])):
                    sd = np.sqrt(var)
                    mv_rows.append({"category": c, "model": "adjusted" if adjusted else "unadjusted",
                                    "contrast": lab, "OR": round(np.exp(est), 2),
                                    "CI_low": round(np.exp(est - 1.96 * sd), 2),
                                    "CI_high": round(np.exp(est + 1.96 * sd), 2),
                                    "n_reports": len(mv_ids), "n_events": int(yv.sum()), "n_clusters": G})
                if adjusted:
                    for j, nm in enumerate(names):
                        if j > 2:
                            mv_rows.append({"category": c, "model": "adjusted", "contrast": f"covariate: {nm}",
                                            "OR": round(np.exp(beta[j]), 2),
                                            "CI_low": round(np.exp(beta[j] - 1.96 * se[j]), 2),
                                            "CI_high": round(np.exp(beta[j] + 1.96 * se[j]), 2),
                                            "n_reports": len(mv_ids), "n_events": int(yv.sum()), "n_clusters": G})
    out_fn(pd.DataFrame(mv_rows), "12_multivariable_logistic.csv")
    out_fn(pd.DataFrame([{"quarters": f"20{quarters[0][0]:02d}Q{quarters[0][1]}-20{quarters[-1][0]:02d}Q{quarters[-1][1]}",
                          "n_quarters": len(quarters), "deleted_caseids_listed": len(deleted),
                          "records_removed_as_deleted": n_del, "reports_after_dedup": n_total}]),
           "0_data_summary.csv")

    # ---- flow and report characteristics (READUS-PV items 8a and 8b) ----
    def _decode(arr, mp_name):
        lst = np.array([str(x).strip().upper() for x in maps[mp_name][1]] + [""])
        a = np.asarray(arr).astype(np.int64)
        return lst[np.where(a < 0, len(lst) - 1, a)]
    flow = [{"section": "flow", "group": "", "characteristic": c, "level": "", "n": v, "pct": np.nan}
            for c, v in (("DEMO records read", n_read), ("records removed (FDA deleted-case lists)", n_del),
                         ("reports after deduplication", n_total))]
    role_ids = {"rare": groups["POOLED_rare_all"][2], "specialty": groups["POOLED_specialty"][2],
                "common": groups["POOLED_controls"][2]}
    ex_keys = [k for k, (g, h, ids) in groups.items() if g == "example"]
    if ex_keys:
        role_ids["example"] = np.unique(np.concatenate([groups[k][2] for k in ex_keys]))
    char_rows = []
    for g_name, ids in role_ids.items():
        drugs = [k for k, (g, h, x) in groups.items() if g == g_name and not k.startswith("POOLED")
                 and (g_name != "rare" or k in parents)]
        flow.append({"section": "flow", "group": g_name, "characteristic": "primary-suspect reports", "level": "",
                     "n": len(ids), "pct": np.nan})
        flow.append({"section": "flow", "group": g_name, "characteristic": "drugs", "level": "",
                     "n": len(drugs), "pct": np.nan})
        flow.append({"section": "flow", "group": g_name, "characteristic": f"drugs with at least {min_drug_n} reports",
                     "level": "", "n": sum(1 for k in drugs if len(groups[k][2]) >= min_drug_n), "pct": np.nan})
        if len(ids) == 0:
            continue
        pos, _ = lookup(ids)
        unit = _decode(keep_acod[pos], "age_cod")
        fac = np.array([AGE_FACTORS.get(u, np.nan) for u in unit], dtype=float)
        age_y = keep_age[pos].astype(np.float64) * fac
        year = np.where(keep_fda[pos] == FAR_FUTURE, np.nan, (keep_fda[pos] // 10000).astype(np.float64))
        char_rows += characteristic_rows(g_name, len(ids), _decode(keep_sex[pos], "sex"), age_y,
                                         _decode(keep_ctry[pos], "occr_country"), year,
                                         _decode(keep_rep[pos], "rept_cod"), _decode(keep_occ[pos], "occp_cod"))
    out_fn(pd.DataFrame(flow + char_rows), "14_flow_and_characteristics.csv")

    # ---- per-drug analyses ----
    min_n = min_drug_n
    small = [k for k in parents if len(groups[k][2]) < min_n]
    if small:
        print(f"  drugs excluded from the per-drug analysis for having fewer than {min_n} reports: {small} (still shown in table 1)", flush=True)
    unit_keys = {"rare": [k for k in parents if len(groups[k][2]) >= min_n],
                 "specialty": [k for k, (g, h, ids) in groups.items()
                               if g == "specialty" and not k.startswith("POOLED") and len(ids) >= min_n],
                 "common": [k for k, (g, h, ids) in groups.items()
                            if g == "common" and not k.startswith("POOLED") and len(ids) >= min_n]}
    drug_rows, re_rows = [], []
    for c in sens_cats:
        per_group = {}
        for g, keys in unit_keys.items():
            ev, ns = [], []
            for k in keys:
                ids = groups[k][2]
                n = int(len(ids))
                if n == 0:
                    continue
                a = int(in_flag(ids, flag_arr[c]).sum())
                lo, hi = wilson(a, n)
                drug_rows.append({"category": c, "group": g, "drug": k, "n_reports": n, "events": a,
                                  "pct": round(100 * a / n, 2), "pct_CI_low": round(100 * lo, 2),
                                  "pct_CI_high": round(100 * hi, 2)})
                ev.append(a); ns.append(n)
            per_group[g] = dl_pool(ev, ns)
            mu, se, tau2, i2, kk = per_group[g]
            if kk:
                p_, plo, phi = expit(mu), expit(mu - 1.96 * se), expit(mu + 1.96 * se)
                re_rows.append({"category": c, "comparison": f"pooled within {g}", "k_drugs": kk,
                                "estimate": round(100 * p_, 2), "CI_low": round(100 * plo, 2),
                                "CI_high": round(100 * phi, 2), "unit": "pooled event proportion (%)",
                                "tau2": round(tau2, 3), "I2(%)": round(i2, 1)})
        for ref in ("common", "specialty"):
            m1, s1, *_, k1 = per_group["rare"]
            m2, s2, *_, k2 = per_group[ref]
            if k1 and k2:
                d, sd = m1 - m2, np.sqrt(s1 ** 2 + s2 ** 2)
                re_rows.append({"category": c, "comparison": f"rare vs {ref} (drug as unit)",
                                "k_drugs": f"{k1} vs {k2}",
                                "estimate": round(np.exp(d), 2), "CI_low": round(np.exp(d - 1.96 * sd), 2),
                                "CI_high": round(np.exp(d + 1.96 * sd), 2), "unit": "pooled odds ratio (OR)",
                                "tau2": "", "I2(%)": ""})
    drug_df = pd.DataFrame(drug_rows)
    out_fn(drug_df, "7_drug_level.csv")
    out_fn(pd.DataFrame(re_rows), "9_random_effects.csv")

    loo_rows = []
    for k in unit_keys["specialty"]:
        ref_ids = np.setdiff1d(spec_ids, groups[k][2])
        keepk = [x for x in unit_keys["specialty"] if x != k]
        for c in sens_cats:
            a = int(in_flag(rare_ids, flag_arr[c]).sum()); n1 = int(len(rare_ids))
            c1 = int(in_flag(ref_ids, flag_arr[c]).sum()); n2 = int(len(ref_ids))
            r, lo, hi = ror_ci(a, n1 - a, c1, n2 - c1)
            m1, s1, *_ = dl_pool([int(in_flag(groups[x][2], flag_arr[c]).sum())
                                  for x in unit_keys["rare"] if len(groups[x][2])],
                                 [int(len(groups[x][2])) for x in unit_keys["rare"] if len(groups[x][2])])
            m2, s2, *_ = dl_pool([int(in_flag(groups[x][2], flag_arr[c]).sum())
                                  for x in keepk if len(groups[x][2])],
                                 [int(len(groups[x][2])) for x in keepk if len(groups[x][2])])
            d, sd = m1 - m2, np.sqrt(s1 ** 2 + s2 ** 2)
            loo_rows.append({"removed_specialty_drug": k, "category": c,
                             "pooled_ROR": round(r, 2), "pooled_CI": f"{lo:.2f}-{hi:.2f}",
                             "RE_OR_drug_unit": round(np.exp(d), 2),
                             "RE_CI": f"{np.exp(d - 1.96 * sd):.2f}-{np.exp(d + 1.96 * sd):.2f}"})
    out_fn(pd.DataFrame(loo_rows), "8_specialty_leave_one_out.csv")

    for c in ("nonadherence", "any_me"):
        if c in flag_arr and len(drug_df):
            try:
                forest_plot(drug_df[drug_df["category"] == c], [r for r in re_rows if r["category"] == c],
                            c, os.path.join(out, f"forest_{c}.png"))
            except ImportError:
                print("  matplotlib not installed; skipping forest plots (pip install matplotlib)")
                break

    out_fn(pd.DataFrame(rows), "1_summary_by_drug.csv")
    out_fn(pd.DataFrame(pt_rows), "2_pt_by_drug.csv")
    out_fn(pd.DataFrame(rep_rows), "3_reporter_by_drug.csv")
    out_fn(pd.DataFrame(yr_rows), "4_yearly_by_drug.csv")

    hit_only = dict_counts > 0
    dict_df = pd.DataFrame({"pt": [p for p, kk in zip(pts_global, hit_only) if kk],
                            "n_reports_all_faers": dict_counts[hit_only]})
    dict_df["category"] = dict_df["pt"].map(lambda p: "|".join(_cats5(p, pt_map, smq_sets)))
    dict_df = dict_df.sort_values(["n_reports_all_faers", "pt"],
                                  ascending=[False, True]).reset_index(drop=True)
    dict_df["keep_after_review(Y/N)"] = ""
    out_fn(dict_df, "5_keyword_pt_dictionary.csv")

    summary = pd.DataFrame(rows)
    show = [c for c in ["drug", "n_reports_PS", "n_nonadherence", "pct_nonadherence",
                        "ROR_vs_controls_nonadherence", "ROR_vs_specialty_nonadherence",
                        "n_any_me", "pct_consumer_reported"] if c in summary.columns]
    with pd.option_context("display.max_rows", None, "display.width", 220):
        print(summary[show].to_string(index=False))
    print(f"\nResults written to: {os.path.abspath(out)}")


_CAT_RX = {c: re.compile(p) for c, p in CATEGORIES.items()}
_EXCL_RX = re.compile(EXCLUDE_PT)


def pt_name(raw):
    return raw


def _cats(pt, pt_map=None):
    s = str(pt).strip().lower()
    if pt_map is not None:
        return (pt_map[s],) if s in pt_map else ()
    if _EXCL_RX.search(s):
        return ()
    return tuple(c for c, r in _CAT_RX.items() if r.search(s))


def _cats5(pt, pt_map, smq_sets):
    base = _cats(pt, pt_map)
    s = str(pt).strip().lower()
    for lab, terms in (smq_sets or {}).items():
        if s in terms:
            base = base + (lab,)
    return base


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--faers-dir", required=True)
    ap.add_argument("--out", default="./out")
    ap.add_argument("--phase", required=True, choices=["demo", "reac", "drug", "indi", "agg"])
    ap.add_argument("--work", default=None)
    ap.add_argument("--pt-dict", default=None)
    ap.add_argument("--smq-file", action="append", default=[])
    ap.add_argument("--min-drug-n", type=int, default=100)
    ap.add_argument("--start", default="2012Q4")
    ap.add_argument("--end", default="2099Q4")
    ap.add_argument("--groups", default=None,
                    help="JSON file assigning comparator roles (agg phase only); see load_groups_config")
    args = ap.parse_args()
    start, end = parse_q(args.start), parse_q(args.end)
    os.makedirs(args.out, exist_ok=True)
    work = args.work or os.path.join(args.out, "_work")
    os.makedirs(work, exist_ok=True)

    files, deleted_files = find_sources(args.faers_dir)
    quarters = sorted(k for k in files if start <= k <= end
                      and all(t in files[k] for t in ("DEMO", "DRUG", "REAC")))
    if not quarters:
        sys.exit("no complete DEMO/DRUG/REAC quarterly files found; check --faers-dir.")
    print(f"phase {args.phase}: {len(quarters)} quarters, {len(deleted_files)} deleted-case list files", flush=True)
    {"demo": phase_demo, "reac": phase_reac, "drug": phase_drug,
     "indi": phase_indi, "agg": phase_agg}[args.phase](
        files, quarters, work, out=args.out, pt_dict=args.pt_dict, smq_files=args.smq_file,
        deleted_files=deleted_files, min_drug_n=args.min_drug_n, groups_cfg=args.groups)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Medication-use and adherence events for oral orphan drugs in FAERS: single-pass reference implementation.

Same analysis as faers_me_feasibility_v7_staged.py, run in one pass over unpacked quarterly files
(no zip support); it needs considerably more memory on the full extract.
FDA deleted-case lists found under --faers-dir are applied before de-duplication. All four naming
conventions the FDA has used are recognised: DELEnnQn.TXT, ADRnnQnDeletedCases.txt,
nnQnDeletedCases.txt and AllDeletedCases.txt.

Usage:
    python faers_me_feasibility_v7_singlepass.py --faers-dir DIR --out ./out \
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
import csv
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------
# 1. candidate drugs (the regex is applied to the upper-cased concatenation of drugname and prod_ai)
#    group: rare = oral drugs taken lifelong for a rare disease; common = common chronic-disease oral-drug comparators
#    handling: storage / administration burden (to verify against the label)
# ----------------------------------------------------------------------------
# Health-professional occupation codes: the FDA documentation lists MD/PH/OT; HP and RN appear in the data but not in the current documentation, so they are treated as "other/unspecified health professional"
HCP_CODES = {"MD", "PH", "OT", "HP", "RN"}
# HP and RN occur in the data but are not defined in the FAERS documentation; a sensitivity analysis uses only these
HCP_DOCUMENTED = {"MD", "PH", "OT"}
AGE_FACTORS = {"YR": 1.0, "DEC": 10.0, "MON": 1 / 12, "WK": 1 / 52, "DY": 1 / 365, "HR": 1 / 8760}

CANDIDATES = [
    # key, regex, exclude_regex, group, handling
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
    # Non-rare oral drugs distributed through specialty pharmacies with patient support programmes (separates a "rare disease" effect from a "specialty pharmacy / support programme" effect)
    ("spec_teriflunomide", r"TERIFLUNOMIDE|AUBAGIO", None, "specialty", "specialty comparator"),
    ("spec_dimethyl_fumarate", r"DIMETHYL FUMARATE|TECFIDERA", None, "specialty", "specialty comparator"),
    ("spec_fingolimod", r"FINGOLIMOD|GILENYA", None, "specialty", "specialty comparator"),
    ("spec_apremilast", r"APREMILAST|OTEZLA", None, "specialty", "specialty comparator"),
    ("spec_tofacitinib", r"TOFACITINIB|XELJANZ", None, "specialty", "specialty comparator"),
    # Worked example: the oral orphan drug whose dose-omission reports a published study presented as a new adverse event (used only in the worked example; enters no pooled group)
    ("example_nefecon", r"NEFECON|TARPEYO|KINPEYGO", None, "example", "worked example: Wang 2024 Kidney Int Rep"),
    # Common chronic-disease oral-drug comparators
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


# ----------------------------------------------------------------------------
# 2. screening keywords (applied to the lower-cased PT); the analysis itself uses the reviewed dictionary (--pt-dict) and the SMQ lists
# ----------------------------------------------------------------------------
# coarse keywords used when --pt-dict is absent (round 1 showed many false matches; fallback only)
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
# Composite outcome: only these four categories; discontinuation (a physician may stop the drug for an adverse reaction), off-label use and dietary nonadherence are reported separately
ANY_ME = ["storage_handling", "nonadherence", "dosing_admin_error", "access_supply"]

# Indication restriction (indi_pt from the INDI table, lower-cased regex). Sub-entries inherit the parent pattern.
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


def indication_pattern(key):
    for k, pat in INDICATIONS.items():
        if key == k or key.startswith(k):
            return pat
    return None


# The FDA has used four naming conventions for the deleted-case lists.
DEL_RE = re.compile(r"^(?:DELE.*|ADR.*DELETEDCASES|ALLDELETEDCASES|\d{2}Q\dDELETEDCASES)\.TXT$", re.I)


def read_deleted_caseids(root):
    """Collect the caseids listed in the FDA deleted-case files (four naming conventions); take the
    first numeric field of each line."""
    ids = set()
    for d, _, files in os.walk(root):
        for f in files:
            if DEL_RE.match(f):
                with open(os.path.join(d, f), encoding="latin-1") as fh:
                    for line in fh:
                        m = re.match(r"\s*(\d+)", line)
                        if m:
                            ids.add(int(m.group(1)))
    return ids


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


# Leave-one-out and peak-year exclusion use only these parent entries (formulation sub-entries are subsets of them)
PARENT_RARE = {"penicillamine", "trientine_all", "zinc_acetate_galzin", "nitisinone_all", "cysteamine_oral",
               "carglumic_acid", "sapropterin", "sepiapterin", "tiopronin", "miglustat", "eliglustat",
               "sodium_phenylbutyrate", "glycerol_phenylbutyrate", "betaine_anhydrous", "cholic_acid"}

FILE_RE = re.compile(r"^(DEMO|DRUG|REAC|INDI)(\d{2})Q([1-4])(?:_NEW)?\.TXT$", re.I)


def find_files(root):
    found = defaultdict(dict)  # (yy,q) -> {table: path}
    for d, _, files in os.walk(root):
        for f in files:
            m = FILE_RE.match(f)
            if m:
                tbl, yy, q = m.group(1).upper(), int(m.group(2)), int(m.group(3))
                found[(yy, q)][tbl] = os.path.join(d, f)
    return found


def in_range(yyq, start, end):
    return start <= yyq <= end


def parse_q(s):
    m = re.match(r"^(\d{4})Q([1-4])$", s.upper())
    if not m:
        sys.exit(f"bad quarter format: {s} (expected e.g. 2014Q1)")
    return (int(m.group(1)) % 100, int(m.group(2)))


def read_faers(path, usecols):
    df = pd.read_csv(path, sep="$", dtype=str, encoding="latin-1", quoting=csv.QUOTE_NONE,
                     on_bad_lines="skip", low_memory=False)
    df.columns = [c.strip().lower() for c in df.columns]
    if "primaryid" not in df.columns:
        return None  # LAERS format
    keep = [c for c in usecols if c in df.columns]
    df = df[keep].copy()
    for c in usecols:
        if c not in df.columns:
            df[c] = np.nan
    return df


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


GROUP_STYLE = {  # colours checked with a dataviz validator (colour-blind separation passes); shape as the second encoding
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
    rows = []  # (label, pct, lo, hi, group, kind, right_text)
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--faers-dir", required=True)
    ap.add_argument("--out", default="./out")
    ap.add_argument("--start", default="2012Q4")
    ap.add_argument("--end", default="2099Q4")
    ap.add_argument("--min-drug-n", type=int, default=100,
                    help="Minimum number of reports for a drug to enter the per-drug analyses (tables 7/8/9 and the forest plots); default 100")
    ap.add_argument("--smq-file", action="append", default=[],
                    help="PT list for the MedDRA SMQ 'Medication errors' (one PT per line, or a CSV with a pt column). Adds an SMQ outcome as a sensitivity analysis.")
    ap.add_argument("--pt-dict", default=None,
                    help="Reviewed PT dictionary CSV (columns pt, keep, category; only keep=Y rows are used). Coarse keywords are used if omitted.")
    ap.add_argument("--groups", default=None,
                    help="JSON file assigning comparator roles; see load_groups_config")
    args = ap.parse_args()
    start, end = parse_q(args.start), parse_q(args.end)
    os.makedirs(args.out, exist_ok=True)

    files = find_files(args.faers_dir)
    quarters = sorted(k for k in files if in_range(k, start, end)
                      and all(t in files[k] for t in ("DEMO", "DRUG", "REAC")))
    if not quarters:
        sys.exit("no complete DEMO/DRUG/REAC quarterly files found; check --faers-dir.")
    print(f"found {len(quarters)} complete quarters: 20{quarters[0][0]:02d}Q{quarters[0][1]} to 20{quarters[-1][0]:02d}Q{quarters[-1][1]}")

    # ---------- pass 1: DEMO, de-duplication ----------
    demo_parts = []
    for k in quarters:
        df = read_faers(files[k]["DEMO"], ["primaryid", "caseid", "fda_dt", "occp_cod", "rept_cod",
                                           "sex", "age", "age_cod", "occr_country"])
        if df is None:
            print(f"  skipping LAERS-format quarter: {files[k]['DEMO']}")
            continue
        df["fq"] = k[0] * 10 + k[1]
        demo_parts.append(df)
    demo = pd.concat(demo_parts, ignore_index=True)
    demo["primaryid"] = pd.to_numeric(demo["primaryid"], errors="coerce")
    demo["caseid"] = pd.to_numeric(demo["caseid"], errors="coerce")
    demo["fda_dt"] = pd.to_numeric(demo["fda_dt"], errors="coerce")
    demo = demo.dropna(subset=["primaryid", "caseid"])
    n_read = len(demo)
    deleted = read_deleted_caseids(args.faers_dir)
    n_del = int(demo["caseid"].isin(deleted).sum()) if deleted else 0
    if deleted:
        demo = demo[~demo["caseid"].isin(deleted)]
    print(f"deleted-case list: {len(deleted):,} caseids, {n_del:,} records removed" if deleted else "no deleted-case list (DELE*.txt) found; nothing removed")
    demo = demo.sort_values(["caseid", "fda_dt", "primaryid"]).drop_duplicates("caseid", keep="last")
    demo["year"] = (demo["fda_dt"] // 10000).astype("Int64")
    demo = demo.set_index("primaryid")
    keep_ids = demo.index.values
    n_total = len(demo)
    print(f"deduplicated reports: {n_total:,}")

    # ---------- pass 2: REAC, event flags ----------
    cat_re = {c: re.compile(p) for c, p in CATEGORIES.items()}
    excl_re = re.compile(EXCLUDE_PT)
    pt_cat_cache = {}
    pt_map = None
    if args.pt_dict:
        pdict = pd.read_csv(args.pt_dict, encoding="utf-8-sig", dtype=str)
        pdict.columns = [c.strip().lower() for c in pdict.columns]
        pdict = pdict[pdict["keep"].str.strip().str.upper() == "Y"]
        pt_map = {p.strip().lower(): c.strip() for p, c in zip(pdict["pt"], pdict["category"])}
        cat_names = sorted(set(pt_map.values()))
        print(f"using reviewed PT dictionary: {len(pt_map)} PTs, categories {cat_names}")
    else:
        cat_names = list(CATEGORIES)
        print("--pt-dict not supplied; using coarse keywords (magnitude only)")
    smq_sets = {}
    for spec in args.smq_file:
        label, path = spec.split("=", 1) if "=" in spec else ("medication_errors", spec)
        raw = open(path, encoding="utf-8-sig").read().splitlines()
        if raw and raw[0].strip().lower().startswith("pt"):
            smq_df = pd.read_csv(path, encoding="utf-8-sig", dtype=str)
            smq_df.columns = [c.strip().lower() for c in smq_df.columns]
            terms = {x.strip().lower() for x in smq_df["pt"].dropna()}
        else:
            terms = {x.strip().lower() for x in raw if x.strip()}
        smq_sets[f"smq_{label}"] = terms
        cat_names = cat_names + [f"smq_{label}"]
        print(f"SMQ list {label}: {len(terms)} terms")
    missing = [c for c in ANY_ME if c not in cat_names]
    if missing:
        print(f"  note: category {missing} is absent from the dictionary; the composite outcome will lack it")

    def pt_cats(pt):
        if pt in pt_cat_cache:
            return pt_cat_cache[pt]
        s = str(pt).strip().lower()
        if pt_map is not None:
            cats = (pt_map[s],) if s in pt_map else ()
        else:
            cats = () if excl_re.search(s) else tuple(c for c, r in cat_re.items() if r.search(s))
        for lab, terms in smq_sets.items():
            if s in terms:
                cats = cats + (lab,)
        pt_cat_cache[pt] = cats
        return cats

    flag_ids = {c: set() for c in cat_names}
    hit_rows = []          # (primaryid, pt), hits only
    for k in quarters:
        df = read_faers(files[k]["REAC"], ["primaryid", "pt"])
        if df is None:
            continue
        df["primaryid"] = pd.to_numeric(df["primaryid"], errors="coerce")
        df = df[df["primaryid"].isin(keep_ids)].dropna(subset=["pt"])
        df["pt"] = df["pt"].str.strip()
        uniq = df["pt"].unique()
        hit_pts = {p for p in uniq if pt_cats(p)}
        sub = df[df["pt"].isin(hit_pts)].drop_duplicates()
        for pid, pt in zip(sub["primaryid"].values, sub["pt"].values):
            for c in pt_cats(pt):
                flag_ids[c].add(pid)
        hit_rows.append(sub)
    hits = pd.concat(hit_rows, ignore_index=True) if hit_rows else pd.DataFrame(columns=["primaryid", "pt"])
    flag_ids["any_me"] = set().union(*(flag_ids[c] for c in ANY_ME if c in flag_ids))
    print("event reports in the database: " + ", ".join(f"{c}={len(v):,}" for c, v in flag_ids.items()))

    # PT dictionary of keyword hits (for manual review)
    dict_df = (hits.drop_duplicates(["primaryid", "pt"]).groupby("pt").size()
               .rename("n_reports_all_faers").reset_index())
    dict_df["category"] = dict_df["pt"].map(lambda p: "|".join(pt_cats(p)))
    dict_df = dict_df.sort_values("n_reports_all_faers", ascending=False)
    dict_df["keep_after_review(Y/N)"] = ""
    dict_df.to_csv(os.path.join(args.out, "5_keyword_pt_dictionary.csv"), index=False, encoding="utf-8-sig")

    # ---------- pass 3: DRUG, primary suspect drug ----------
    cand = [(k, re.compile(r), re.compile(x) if x else None, g, h) for k, r, x, g, h in CANDIDATES]
    drug_hits = []
    for k in quarters:
        df = read_faers(files[k]["DRUG"], ["primaryid", "drug_seq", "role_cod", "drugname", "prod_ai"])
        if df is None:
            continue
        df = df[df["role_cod"].str.strip().str.upper() == "PS"]
        df["primaryid"] = pd.to_numeric(df["primaryid"], errors="coerce")
        df = df[df["primaryid"].isin(keep_ids)]
        s = (df["drugname"].fillna("") + " | " + df["prod_ai"].fillna("")).str.upper()
        for key, rx, ex, g, h in cand:
            m = s.str.contains(rx, regex=True)
            if ex is not None:
                m = m & ~s.str.contains(ex, regex=True)
            if m.any():
                drug_hits.append(pd.DataFrame({"primaryid": df.loc[m, "primaryid"].values,
                                               "drug_seq": df.loc[m, "drug_seq"].astype(str).str.strip().values,
                                               "drug": key}))
    if not drug_hits:
        sys.exit("no candidate drug matched; check the drugname/prod_ai fields.")
    dh = pd.concat(drug_hits, ignore_index=True).drop_duplicates()

    # ---------- pass 4: INDI, indication restriction ----------
    ind_keys = [k for k, _, _, g, _ in CANDIDATES if g in ("rare", "example") and indication_pattern(k)]
    cand_ids = set(dh.loc[dh["drug"].isin(ind_keys), "primaryid"])
    indi_rows = []
    for k in quarters:
        if "INDI" not in files[k]:
            continue
        df = read_faers(files[k]["INDI"], ["primaryid", "indi_drug_seq", "indi_pt"])
        if df is None:
            continue
        df["primaryid"] = pd.to_numeric(df["primaryid"], errors="coerce")
        df = df[df["primaryid"].isin(cand_ids)].dropna(subset=["indi_pt"])
        df["indi_drug_seq"] = df["indi_drug_seq"].astype(str).str.strip()
        df["indi_pt"] = df["indi_pt"].str.strip().str.lower()
        indi_rows.append(df)
    indi = pd.concat(indi_rows, ignore_index=True) if indi_rows else pd.DataFrame(columns=["primaryid", "indi_drug_seq", "indi_pt"])
    ind_ok = {}
    indi_top_rows = []
    for key in ind_keys:
        sub = dh[dh["drug"] == key].merge(indi, left_on=["primaryid", "drug_seq"],
                                         right_on=["primaryid", "indi_drug_seq"], how="inner")
        pat = re.compile(indication_pattern(key))
        hit = sub["indi_pt"].map(lambda t: bool(pat.search(t)))
        ind_ok[key] = set(sub.loc[hit, "primaryid"])
        top = (sub.drop_duplicates(["primaryid", "indi_pt"]).groupby("indi_pt")["primaryid"].nunique()
               .sort_values(ascending=False).head(15))
        for t, cnt in top.items():
            indi_top_rows.append({"drug": key, "indi_pt": t, "n_reports": int(cnt),
                                  "matches_target": bool(pat.search(t))})
    if not indi_rows:
        print("  no INDI file found; skipping the indication-restricted analysis")

    # ---------- aggregation ----------
    roles, trunc, cfg = load_groups_config(args.groups)
    orig_group = {k: g for k, _, _, g, _ in CANDIDATES}
    groups = {}
    for key, _, _, g, h in CANDIDATES:
        if roles[key] == "off":
            continue
        ids = set(dh.loc[dh["drug"] == key, "primaryid"])
        if key in trunc and ids:
            fqs = demo["fq"].reindex(list(ids))
            ids = set(fqs.index[fqs.notna() & (fqs < trunc[key])])
        groups[key] = (roles[key], h, ids)
    if args.groups:
        import json
        with open(os.path.join(args.out, "config_used.json"), "w", encoding="utf-8") as f:
            json.dump({"file": os.path.abspath(args.groups), **cfg,
                       "resolved_roles": {k: r for k, r in roles.items() if r != "off"},
                       "n_reports_after_truncation": {k: int(len(v[2])) for k, v in groups.items()
                                                      if k in trunc}}, f, indent=2)
        print("groups: " + ", ".join(f"{r}={sum(1 for v in groups.values() if v[0] == r)}" for r in ROLES[:4]))
    rare_ids = set().union(*(v[2] for v in groups.values() if v[0] == "rare"))
    ctrl_ids = set().union(*(v[2] for v in groups.values() if v[0] == "common"))
    spec_ids = set().union(*(v[2] for v in groups.values() if v[0] == "specialty"))
    groups["POOLED_rare_all"] = ("rare", "all rare-disease drugs pooled", rare_ids)
    groups["POOLED_controls"] = ("common", "all comparators pooled", ctrl_ids)
    groups["POOLED_specialty"] = ("specialty", "all specialty comparators pooled", spec_ids)
    n_ctrl = len(ctrl_ids)
    n_spec = len(spec_ids)
    all_cats = cat_names + ["any_me"]

    rows, pt_rows, rep_rows, yr_rows = [], [], [], []
    for key, (g, h, ids) in groups.items():
        n = len(ids)
        rec = {"drug": key, "group": g, "handling(to verify)": h, "n_reports_PS": n}
        if n == 0:
            rows.append(rec)
            continue
        d = demo.loc[list(ids)]
        rec["pct_consumer_reported"] = round(100 * (d["occp_cod"].str.upper() == "CN").mean(), 1)
        for c in all_cats:
            a = len(ids & flag_ids[c])
            b = n - a
            cc = len(flag_ids[c]) - a
            dd = (n_total - n) - cc
            r, lo, hi = ror_ci(a, b, cc, dd)
            rec[f"n_{c}"] = a
            rec[f"pct_{c}"] = round(100 * a / n, 1)
            rec[f"ROR_vs_allFAERS_{c}"] = f"{r:.2f} ({lo:.2f}-{hi:.2f})"
            # Comparison against the common chronic-disease oral-drug comparators (the more appropriate reference)
            if g == "rare" and n_ctrl > 0:
                c1 = len(ctrl_ids & flag_ids[c])
                r2, lo2, hi2 = ror_ci(a, b, c1, n_ctrl - c1)
                rec[f"ROR_vs_controls_{c}"] = f"{r2:.2f} ({lo2:.2f}-{hi2:.2f})"
            if g == "rare" and n_spec > 0:
                c2 = len(spec_ids & flag_ids[c])
                r3, lo3, hi3 = ror_ci(a, b, c2, n_spec - c2)
                rec[f"ROR_vs_specialty_{c}"] = f"{r3:.2f} ({lo3:.2f}-{hi3:.2f})"
        me_ids = ids & flag_ids["any_me"]
        if me_ids:
            dm = demo.loc[list(me_ids)]
            rec["pct_consumer_among_me"] = round(100 * (dm["occp_cod"].str.upper() == "CN").mean(), 1)
        rows.append(rec)

        hp = hits[hits["primaryid"].isin(ids)].drop_duplicates(["primaryid", "pt"])
        for pt, cnt in hp.groupby("pt").size().sort_values(ascending=False).items():
            pt_rows.append({"drug": key, "pt": pt, "category": "|".join(pt_cats(pt)), "n_reports": cnt})

        d2 = d.assign(any_me=d.index.isin(list(flag_ids["any_me"])))
        if "nonadherence" in flag_ids:
            d2["nonadh"] = d2.index.isin(list(flag_ids["nonadherence"]))
        else:
            d2["nonadh"] = False
        for col in ("occp_cod", "rept_cod"):
            t = d2.groupby(d2[col].fillna("NA").str.upper()).agg(
                n_reports=("any_me", "size"), n_any_me=("any_me", "sum"), n_nonadherence=("nonadh", "sum"))
            for lvl, r in t.iterrows():
                rep_rows.append({"drug": key, "field": col, "level": lvl, "n_reports": int(r["n_reports"]),
                                 "n_any_me": int(r["n_any_me"]), "n_nonadherence": int(r["n_nonadherence"])})
        for c in all_cats:
            d2[f"n_{c}"] = d2.index.isin(list(flag_ids[c]))
        yr = d2.groupby("year").agg(n_reports=("any_me", "size"),
                                    **{f"n_{c}": (f"n_{c}", "sum") for c in all_cats}).reset_index()
        yr.insert(0, "drug", key)
        yr_rows.extend(yr.to_dict("records"))

    # ---------- sensitivity analysis: pooled rare-disease drugs vs the two comparator groups ----------
    # Leave-one-out and peak-year exclusion use only the parent entries; formulation sub-entries are subsets of them
    parents = [k for k, (g, h, ids) in groups.items()
               if g == "rare" and not k.startswith("POOLED")
               and (k in PARENT_RARE or orig_group.get(k) != "rare")]
    hcp = HCP_CODES
    sens_cats = [c for c in ["nonadherence", "any_me", "dosing_admin_error", "storage_handling", "access_supply"]
                 + sorted(smq_sets) if c in flag_ids]

    def subset(ids, how):
        if how == "all":
            return ids
        d = demo.loc[list(ids)]
        if how == "exclude_PER":
            return set(d.index[d["rept_cod"].fillna("").str.upper() != "PER"])
        if how == "HCP_only":
            return set(d.index[d["occp_cod"].fillna("").str.upper().isin(hcp)])
        if how == "HCP_strict":
            return set(d.index[d["occp_cod"].fillna("").str.upper().isin(HCP_DOCUMENTED)])
        if how == "EXP_and_HCP":
            m = (d["rept_cod"].fillna("").str.upper() == "EXP") & d["occp_cod"].fillna("").str.upper().isin(hcp)
            return set(d.index[m])
        raise ValueError(how)

    # The year with the most nonadherence reports for each parent drug (used to drop clustering)
    def peak_year_ids(keys):
        out_ids = set()
        if "nonadherence" not in flag_ids:
            return out_ids
        for k in keys:
            ids = groups[k][2]
            dd = demo.loc[list(ids)]
            na = dd[dd.index.isin(list(flag_ids["nonadherence"]))]
            if len(na):
                vc = na["year"].value_counts()
                top_year = vc[vc == vc.max()].index.min()   # ties: earliest year, as in the staged script
                out_ids |= set(dd.index[dd["year"] == top_year])
        return out_ids

    burst_ids = peak_year_ids(parents)
    burst_ref_ids = peak_year_ids([k for k, (g, h, ids) in groups.items()
                                   if g in ("specialty", "common") and not k.startswith("POOLED")])

    scenarios = [("all", "all reports", rare_ids),
                 ("exclude_PER", "excluding periodic reports (PER)", rare_ids),
                 ("HCP_only", "health-professional reports only", rare_ids),
                 ("EXP_and_HCP", "expedited reports from health professionals only", rare_ids),
                 ("HCP_strict", "health-professional reports only, documented codes (MD, PH, OT)", rare_ids),
                 ("all", "excluding the peak nonadherence year of each drug", rare_ids - burst_ids),
                 ("all", "excluding the peak nonadherence year of every drug in all groups",
                  rare_ids - burst_ids, burst_ref_ids)]
    restricted_ids = set().union(*(ind_ok.get(k, set()) for k in parents)) if ind_ok else set()
    if restricted_ids:
        scenarios.append(("all", "restricted to the rare-disease indication", restricted_ids))
    for k in parents:
        scenarios.append(("all", f"leave-one-out: removing {k}", rare_ids - groups[k][2]))

    sens_rows = []
    for sc in scenarios:
        how, label, r_ids = sc[:3]
        ref_excl = sc[3] if len(sc) > 3 else set()
        rs = subset(r_ids, how)
        for ref_name, ref_ids in (("controls", ctrl_ids), ("specialty", spec_ids)):
            if not ref_ids:
                continue
            fs = subset(ref_ids, how) - ref_excl
            for c in sens_cats:
                a = len(rs & flag_ids[c]); n1 = len(rs)
                c1 = len(fs & flag_ids[c]); n2 = len(fs)
                if n1 == 0 or n2 == 0:
                    continue
                r, lo, hi = ror_ci(a, n1 - a, c1, n2 - c1)
                sens_rows.append({"scenario": label, "reference": ref_name, "category": c,
                                  "rare_n": n1, "rare_events": a, "rare_pct": round(100 * a / n1, 2),
                                  "ref_n": n2, "ref_events": c1, "ref_pct": round(100 * c1 / n2, 2),
                                  "ROR": round(r, 2), "CI_low": round(lo, 2), "CI_high": round(hi, 2)})

    out = lambda df, name: df.to_csv(os.path.join(args.out, name), index=False, encoding="utf-8-sig")
    out(pd.DataFrame(sens_rows), "6_sensitivity.csv")

    # ---------- indication coverage ----------
    cov_rows = []
    for k in parents + [k for k, (g, h, ids) in groups.items() if g == "example"]:
        n = len(groups[k][2])
        if n and k in ind_ok:
            cov_rows.append({"drug": k, "n_reports_PS": n, "n_target_indication": len(ind_ok[k]),
                             "pct_target_indication": round(100 * len(ind_ok[k]) / n, 1)})
    out(pd.DataFrame(cov_rows), "10_indication_coverage.csv")
    out(pd.DataFrame(indi_top_rows), "13_indication_top.csv")

    # ---------- worked example: a published dose-omission signal recomputed with different comparators ----------
    # Reproduce first on the original period (2021Q4-2023Q4) and PT (product dose omission issue), then switch comparators; then extend to the full period and the nonadherence composite
    pdoi_ids = set(hits.loc[hits["pt"].str.strip().str.lower() == "product dose omission issue", "primaryid"])
    # every period restricts the example drug and all three references to the same quarters
    ex_periods = (("2021Q4-2023Q4 (original study period)", set(demo.index[(demo["fq"] >= 214) & (demo["fq"] <= 234)])),
                  ("2021Q4-2025Q4 (all quarters since approval)", set(demo.index[(demo["fq"] >= 214) & (demo["fq"] <= 254)])),
                  ("2024Q1-2025Q4 (after the original study period)", set(demo.index[(demo["fq"] >= 241) & (demo["fq"] <= 254)])))
    ex_outcomes = [("PT: product dose omission issue", pdoi_ids)] + \
                  [(c, flag_ids[c]) for c in ("nonadherence", "any_me") if c in flag_ids]
    ex_rows = []
    for k, (g, h, ids) in groups.items():
        if g != "example" or not ids:
            continue
        for period, pids in ex_periods:
            for label, eids0 in (("all reports", ids), ("restricted to target indication", ind_ok.get(k, set()))):
                eids = eids0 & pids
                n1 = len(eids)
                if n1 == 0:
                    continue
                for oname, oset in ex_outcomes:
                    a = len(eids & oset)
                    refs = (("all other reports in the database", pids - eids), ("common chronic-disease comparators", ctrl_ids & pids), ("specialty-pharmacy comparators", spec_ids & pids))
                    for ref_name, ref in refs:
                        n2 = len(ref); c1 = len(ref & oset)
                        r, lo, hi = ror_ci(a, n1 - a, c1, n2 - c1)
                        ex_rows.append({"drug": k, "period": period, "subset": label, "outcome": oname,
                                        "reference": ref_name, "n_reports": n1, "events": a,
                                        "pct": round(100 * a / n1, 2), "ref_n": n2, "ref_events": c1,
                                        "ROR": round(r, 2), "CI_low": round(lo, 2), "CI_high": round(hi, 2)})
    out(pd.DataFrame(ex_rows), "11_worked_example.csv")

    # ---------- multivariable logistic regression (report level, robust SE clustered by drug) ----------
    unit_map = {}
    for g_name, keys in (("rare", parents),
                         ("specialty", [k for k, (g, h, ids) in groups.items() if g == "specialty" and not k.startswith("POOLED")]),
                         ("common", [k for k, (g, h, ids) in groups.items() if g == "common" and not k.startswith("POOLED")])):
        for k in keys:
            for pid in groups[k][2]:
                unit_map.setdefault(pid, (g_name, k))
    mv = demo.loc[list(unit_map.keys())].copy()
    mv["group"] = [unit_map[i][0] for i in mv.index]
    mv["drug"] = [unit_map[i][1] for i in mv.index]
    rc = mv["rept_cod"].fillna("").str.upper()
    oc = mv["occp_cod"].fillna("").str.upper()
    yr = mv["year"].astype("float")
    age = pd.to_numeric(mv["age"], errors="coerce")
    unit = mv["age_cod"].fillna("").str.upper()
    fac = unit.map({"YR": 1, "DEC": 10, "MON": 1 / 12, "WK": 1 / 52, "DY": 1 / 365, "HR": 1 / 8760}).fillna(np.nan)
    age_y = age * fac
    cols = {
        "intercept": np.ones(len(mv)),
        "rare (vs specialty)": (mv["group"] == "rare").astype(float),
        "common (vs specialty)": (mv["group"] == "common").astype(float),
        "periodic report (vs expedited)": (rc == "PER").astype(float),
        "direct report (vs expedited)": (rc == "DIR").astype(float),
        "consumer (vs health professional)": (oc == "CN").astype(float),
        "other/missing reporter": (~oc.isin(["CN"] + sorted(HCP_CODES))).astype(float),
        "male (vs female)": (mv["sex"].fillna("").str.upper() == "M").astype(float),
        "sex missing": (~mv["sex"].fillna("").str.upper().isin(["M", "F"])).astype(float),
        "age <18 (vs 18-64)": (age_y < 18).astype(float),
        "age >=65": (age_y >= 65).astype(float),
        "age missing": age_y.isna().astype(float),
        "non-US country (vs US)": ((mv["occr_country"].fillna("").str.upper() != "US") &
                                   mv["occr_country"].notna()).astype(float),
        "country missing": mv["occr_country"].isna().astype(float),
    }
    base_names = list(cols)
    X_base = np.column_stack([np.asarray(cols[n], dtype=float) for n in base_names])
    yr_arr = yr.to_numpy(dtype=float, na_value=np.nan)
    YEAR_BINS = [(0, 2015), (2016, 2018), (2019, 2021), (2022, 9999)]

    def year_design(yv):
        """Group years; groups with fewer than 5 events are merged with the adjacent group to avoid complete separation (rare outcomes such as access_supply)."""
        bins = []
        for lo, hi in YEAR_BINS:
            ev = yv[(yr_arr >= lo) & (yr_arr <= hi)].sum()
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
            cols_y.append(((yr_arr >= lo) & (yr_arr <= hi)).astype(float))
            names_y.append(f"year {lab(lo, hi)} (vs {ref})")
        return cols_y, names_y

    mv_rows = []
    for c in [x for x in ("nonadherence", "any_me", "access_supply") if x in flag_ids]:
        yv = mv.index.isin(list(flag_ids[c])).astype(float)
        if yv.sum() < 10:
            continue
        cols_y, names_y = year_design(yv)
        X = np.column_stack([X_base] + cols_y) if cols_y else X_base
        names = base_names + names_y
        keep_cols = [i for i in range(X.shape[1]) if i < 3 or X[:, i].std() > 0]
        X = X[:, keep_cols]; names = [names[i] for i in keep_cols]
        for adjusted in (False, True):
            idx = list(range(3)) if not adjusted else list(range(len(names)))
            beta, cov, G = logit_cluster(X[:, idx], yv, mv["drug"].values, [names[i] for i in idx])
            ir, ic = 1, 2
            se = np.sqrt(np.diag(cov))
            if np.any(np.abs(beta[1:]) > 12) or np.any(~np.isfinite(se)):
                print(f"  Warning: the {'adjusted' if adjusted else 'unadjusted'} model for {c} may be subject to complete separation; do not interpret its coefficients")
            for lab, est, var in (("rare vs specialty", beta[ir], cov[ir, ir]),
                                  ("rare vs common", beta[ir] - beta[ic], cov[ir, ir] + cov[ic, ic] - 2 * cov[ir, ic])):
                sd = np.sqrt(var)
                mv_rows.append({"category": c, "model": "adjusted" if adjusted else "unadjusted", "contrast": lab,
                                "OR": round(np.exp(est), 2), "CI_low": round(np.exp(est - 1.96 * sd), 2),
                                "CI_high": round(np.exp(est + 1.96 * sd), 2), "n_reports": len(mv),
                                "n_events": int(yv.sum()), "n_clusters": G})
            if adjusted:
                for j, nm in enumerate(names):
                    if j > 2:
                        mv_rows.append({"category": c, "model": "adjusted", "contrast": f"covariate: {nm}",
                                        "OR": round(np.exp(beta[j]), 2),
                                        "CI_low": round(np.exp(beta[j] - 1.96 * se[j]), 2),
                                        "CI_high": round(np.exp(beta[j] + 1.96 * se[j]), 2),
                                        "n_reports": len(mv), "n_events": int(yv.sum()), "n_clusters": G})
    out(pd.DataFrame(mv_rows), "12_multivariable_logistic.csv")
    out(pd.DataFrame([{"quarters": f"20{quarters[0][0]:02d}Q{quarters[0][1]}-20{quarters[-1][0]:02d}Q{quarters[-1][1]}",
                       "n_quarters": len(quarters), "deleted_caseids_listed": len(deleted),
                       "records_removed_as_deleted": n_del, "reports_after_dedup": n_total}]), "0_data_summary.csv")

    # ---------- flow and report characteristics (READUS-PV items 8a and 8b) ----------
    flow = [{"section": "flow", "group": "", "characteristic": c, "level": "", "n": v, "pct": np.nan}
            for c, v in (("DEMO records read", n_read), ("records removed (FDA deleted-case lists)", n_del),
                         ("reports after deduplication", n_total))]
    role_ids = {"rare": groups["POOLED_rare_all"][2], "specialty": groups["POOLED_specialty"][2],
                "common": groups["POOLED_controls"][2]}
    ex_keys = [k for k, (g, h, ids) in groups.items() if g == "example"]
    if ex_keys:
        role_ids["example"] = set().union(*(groups[k][2] for k in ex_keys))
    char_rows = []
    for g_name, ids in role_ids.items():
        drugs = [k for k, (g, h, x) in groups.items() if g == g_name and not k.startswith("POOLED")
                 and (g_name != "rare" or k in parents)]
        flow.append({"section": "flow", "group": g_name, "characteristic": "primary-suspect reports", "level": "",
                     "n": len(ids), "pct": np.nan})
        flow.append({"section": "flow", "group": g_name, "characteristic": "drugs", "level": "",
                     "n": len(drugs), "pct": np.nan})
        flow.append({"section": "flow", "group": g_name, "characteristic": f"drugs with at least {args.min_drug_n} reports",
                     "level": "", "n": sum(1 for k in drugs if len(groups[k][2]) >= args.min_drug_n), "pct": np.nan})
        if not ids:
            continue
        d = demo.loc[sorted(ids)]
        up = lambda c: d[c].fillna("").astype(str).str.strip().str.upper().to_numpy()
        unit = up("age_cod")
        fac = np.array([AGE_FACTORS.get(u, np.nan) for u in unit], dtype=float)
        age_y = pd.to_numeric(d["age"], errors="coerce").to_numpy(dtype=float) * fac
        year = d["year"].astype("float").to_numpy(dtype=float, na_value=np.nan)
        char_rows += characteristic_rows(g_name, len(ids), up("sex"), age_y, up("occr_country"), year,
                                         up("rept_cod"), up("occp_cod"))
    out(pd.DataFrame(flow + char_rows), "14_flow_and_characteristics.csv")

    # ---------- per-drug analysis: drug-level random-effects pooling ----------
    # The per-drug analysis includes only drugs above the report threshold, so drugs with very few reports do not produce meaningless wide intervals
    min_n = args.min_drug_n
    small = [k for k in parents if len(groups[k][2]) < min_n]
    if small:
        print(f"  drugs excluded from the per-drug analysis for having fewer than {min_n} reports: {small} (still shown in table 1)")
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
                n = len(ids)
                if n == 0:
                    continue
                a = len(ids & flag_ids[c])
                lo, hi = wilson(a, n)
                drug_rows.append({"category": c, "group": g, "drug": k, "n_reports": n, "events": a,
                                  "pct": round(100 * a / n, 2), "pct_CI_low": round(100 * lo, 2),
                                  "pct_CI_high": round(100 * hi, 2)})
                ev.append(a); ns.append(n)
            per_group[g] = dl_pool(ev, ns)
            mu, se, tau2, i2, kk = per_group[g]
            if kk:
                p, plo, phi = expit(mu), expit(mu - 1.96 * se), expit(mu + 1.96 * se)
                re_rows.append({"category": c, "comparison": f"pooled within {g}", "k_drugs": kk,
                                "estimate": round(100 * p, 2), "CI_low": round(100 * plo, 2), "CI_high": round(100 * phi, 2),
                                "unit": "pooled event proportion (%)", "tau2": round(tau2, 3), "I2(%)": round(i2, 1)})
        for ref in ("common", "specialty"):
            m1, s1, *_, k1 = per_group["rare"]
            m2, s2, *_, k2 = per_group[ref]
            if k1 and k2:
                d, sd = m1 - m2, np.sqrt(s1 ** 2 + s2 ** 2)
                re_rows.append({"category": c, "comparison": f"rare vs {ref} (drug as unit)", "k_drugs": f"{k1} vs {k2}",
                                "estimate": round(np.exp(d), 2), "CI_low": round(np.exp(d - 1.96 * sd), 2),
                                "CI_high": round(np.exp(d + 1.96 * sd), 2), "unit": "pooled odds ratio (OR)",
                                "tau2": "", "I2(%)": ""})
    drug_df = pd.DataFrame(drug_rows)
    out(drug_df, "7_drug_level.csv")
    out(pd.DataFrame(re_rows), "9_random_effects.csv")

    # ---------- leave-one-out over the specialty comparators ----------
    loo_rows = []
    for k in unit_keys["specialty"]:
        ref_ids = spec_ids - groups[k][2]
        keep = [x for x in unit_keys["specialty"] if x != k]
        for c in sens_cats:
            a = len(rare_ids & flag_ids[c]); n1 = len(rare_ids)
            c1 = len(ref_ids & flag_ids[c]); n2 = len(ref_ids)
            r, lo, hi = ror_ci(a, n1 - a, c1, n2 - c1)
            m1, s1, *_ = dl_pool([len(groups[x][2] & flag_ids[c]) for x in unit_keys["rare"] if groups[x][2]],
                                 [len(groups[x][2]) for x in unit_keys["rare"] if groups[x][2]])
            m2, s2, *_ = dl_pool([len(groups[x][2] & flag_ids[c]) for x in keep if groups[x][2]],
                                 [len(groups[x][2]) for x in keep if groups[x][2]])
            d, sd = m1 - m2, np.sqrt(s1 ** 2 + s2 ** 2)
            loo_rows.append({"removed_specialty_drug": k, "category": c,
                             "pooled_ROR": round(r, 2), "pooled_CI": f"{lo:.2f}-{hi:.2f}",
                             "RE_OR_drug_unit": round(np.exp(d), 2),
                             "RE_CI": f"{np.exp(d - 1.96 * sd):.2f}-{np.exp(d + 1.96 * sd):.2f}"})
    out(pd.DataFrame(loo_rows), "8_specialty_leave_one_out.csv")

    # ---------- forest plots ----------
    for c in ("nonadherence", "any_me"):
        if c in flag_ids and len(drug_df):
            try:
                forest_plot(drug_df[drug_df["category"] == c], [r for r in re_rows if r["category"] == c],
                            c, os.path.join(args.out, f"forest_{c}.png"))
            except ImportError:
                print("  matplotlib not installed; skipping forest plots (pip install matplotlib)")
                break
    summary = pd.DataFrame(rows)
    out(summary, "1_summary_by_drug.csv")
    out(pd.DataFrame(pt_rows), "2_pt_by_drug.csv")
    out(pd.DataFrame(rep_rows), "3_reporter_by_drug.csv")
    out(pd.DataFrame(yr_rows), "4_yearly_by_drug.csv")

    show = [c for c in ["drug", "n_reports_PS", "n_nonadherence", "pct_nonadherence",
                        "ROR_vs_controls_nonadherence", "ROR_vs_specialty_nonadherence",
                        "n_any_me", "pct_consumer_reported"] if c in summary.columns]
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(summary[show].to_string(index=False))
    print(f"\nResults written to: {os.path.abspath(args.out)}")
    print("Next: review 5_keyword_pt_dictionary.csv manually, mark the PTs that should not count as N, then fix the outcome definition.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Count-only screen for a second specialty-pharmacy comparator set (no outcome data are read).

For each candidate drug, and for the existing comparator and rare-disease groups for reference, this
script counts de-duplicated primary-suspect FAERS reports and describes their report mix. It reads
only DEMO and DRUG files; REAC (reactions) and INDI are never opened, so no outcome is seen before
the comparator set is fixed. The rules are in docs/second_comparator_rules.md.

De-duplication, deleted-case handling, primary-suspect attribution, the quarter range and the year of
FDA receipt are the same as in faers_me_feasibility_v7_staged.py, whose functions are imported.
If <work>/demo5 from an earlier staged run exists for the same quarters it is reused; otherwise the
DEMO phase is run first and its arrays are written to <work>/demo5.

    python screen_comparator_counts.py --faers-dir DIR --out ./screen --work ./out/_work \
        --start 2013Q1 --end 2025Q4

Outputs (CSV, UTF-8 with BOM) in --out:
    S1_screen_counts.csv     per drug: reports, first and last year, largest single-year share,
                             report type, reporter and country mix
    S2_screen_yearly.csv     per drug and year of FDA receipt: reports
    S3_screen_names.csv      the 15 most frequent DRUGNAME | PROD_AI strings matched for each drug,
                             to check the name patterns
    S4_screen_overlap.csv    reports matched to more than one drug in this screen
"""
import argparse
import importlib.util
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("staged", os.path.join(HERE, "faers_me_feasibility_v7_staged.py"))
staged = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(staged)

# key, pattern, exclusion pattern, set, note
SCREEN = [
    ("s2_enzalutamide", r"ENZALUTAMIDE|XTANDI", None, "specialty2", "AR inhibitor; 160 mg once daily, continuous"),
    ("s2_apalutamide", r"APALUTAMIDE|ERLEADA", None, "specialty2", "AR inhibitor; 240 mg once daily, continuous"),
    ("s2_darolutamide", r"DAROLUTAMIDE|NUBEQA", None, "specialty2", "AR inhibitor; 600 mg twice daily, continuous"),
    ("s2_abemaciclib", r"ABEMACICLIB|VERZENIO", None, "specialty2", "CDK4/6 inhibitor; twice daily, continuous"),
    ("s2c_palbociclib", r"PALBOCICLIB|IBRANCE", None, "specialty2_cyclic", "CDK4/6 inhibitor; 21 days on, 7 off"),
    ("s2c_ribociclib", r"RIBOCICLIB|KISQALI", None, "specialty2_cyclic", "CDK4/6 inhibitor; 21 days on, 7 off"),
    ("x_abiraterone", r"ABIRATERONE|ZYTIGA|YONSA", None, "excluded", "excluded a priori: US generics from November 2018"),
]
# the existing groups, for comparison of report mix (parent entries only, as in the main analysis)
for key, rx, ex, grp, _ in staged.CANDIDATES:
    if grp in ("specialty", "common") or key in staged.PARENT_RARE:
        SCREEN.append((key, rx, ex, grp, "existing group"))

UNION = re.compile("|".join(f"(?:{r})" for _, r, _, _, _ in SCREEN))
HCP_DOC = sorted(staged.HCP_DOCUMENTED)


def scan_drug(files, quarters):
    """Primary-suspect report ids per screened drug, and matched name strings (by primaryid)."""
    comp = [(k, re.compile(r), re.compile(x) if x else None) for k, r, x, _, _ in SCREEN]
    ids = {k: [] for k, _, _ in comp}
    names = {k: [] for k, _, _ in comp}
    for i, q in enumerate(quarters, 1):
        for df in staged.iter_table(files[q]["DRUG"], ["primaryid", "role_cod", "drugname", "prod_ai"]):
            ps = (df["role_cod"].str.strip().str.upper() == "PS").to_numpy()
            if not ps.any():
                continue
            sub = df[ps]
            pid = pd.to_numeric(sub["primaryid"], errors="coerce").to_numpy(dtype="float64")
            ok = ~np.isnan(pid)
            pid = pid[ok].astype(np.int64)
            nm = sub["drugname"].to_numpy()[ok]
            ai = sub["prod_ai"].to_numpy()[ok]
            s = pd.Series(np.where(pd.isna(nm), "", nm).astype(str) + " | "
                          + np.where(pd.isna(ai), "", ai).astype(str)).str.upper()
            hit = s.str.contains(UNION, regex=True).to_numpy()
            if not hit.any():
                continue
            s2, pid2 = s[hit].reset_index(drop=True), pid[hit]
            for k, rx, ex in comp:
                m = s2.str.contains(rx, regex=True).to_numpy()
                if ex is not None:
                    m = m & ~s2.str.contains(ex, regex=True).to_numpy()
                if m.any():
                    ids[k].append(pid2[m])
                    names[k].append(pd.DataFrame({"pid": pid2[m], "name": s2[m].to_numpy()}))
        if i % 10 == 0 or i == len(quarters):
            print(f"  DRUG {i}/{len(quarters)}", flush=True)
    ids = {k: (np.unique(np.concatenate(v)) if v else np.array([], np.int64)) for k, v in ids.items()}
    names = {k: (pd.concat(v, ignore_index=True) if v else pd.DataFrame(columns=["pid", "name"]))
             for k, v in names.items()}
    return ids, names


def load_dedup(work, quarters, deleted_files, wanted=None):
    """De-duplicated DEMO records, as in phase_agg of the staged script.

    Only the reports the screen needs are materialised. The string columns are read one
    quarter at a time and kept for the wanted primaryids alone; holding the whole
    de-duplicated database as Python strings needs several GB and does not fit in memory.
    """
    d_dir = os.path.join(work, "demo5")
    qs = [staged.qkey(k) for k in quarters]
    keys = ["occp_cod", "rept_cod", "occr_country"]
    pid, cid, fda, sizes = [], [], [], []
    for q in qs:
        z = np.load(os.path.join(d_dir, q + ".npz"))
        pid.append(z["pid"]); cid.append(z["cid"]); fda.append(z["fda"])
        sizes.append(len(z["pid"]))
    pid = np.concatenate(pid); cid = np.concatenate(cid); fda = np.concatenate(fda)
    offsets = np.concatenate(([0], np.cumsum(sizes)))
    n_read = len(pid)
    deleted = staged.read_deleted_caseids(deleted_files or [])
    n_del = 0
    keep = np.ones(n_read, dtype=bool)
    if deleted:
        keep = ~np.isin(cid, np.fromiter(deleted, dtype=np.int64, count=len(deleted)))
        n_del = int((~keep).sum())
        pid, cid, fda = pid[keep], cid[keep], fda[keep]
    bits = max(32, int(pid.max()).bit_length())
    order = np.lexsort((fda.astype(np.int64) * np.int64(1 << bits) + pid, cid))
    sc = cid[order]
    last = np.empty(len(sc), dtype=bool); last[-1] = True
    last[:-1] = sc[1:] != sc[:-1]
    idx = order[last]
    srt = np.argsort(pid[idx], kind="stable")
    idx = idx[srt]
    n_dedup = len(idx)
    # position of each surviving record in the concatenated per-quarter arrays
    orig = np.flatnonzero(keep)[idx] if deleted else idx
    sel = np.isin(pid[idx], wanted) if wanted is not None else np.ones(n_dedup, dtype=bool)
    pid_s, fda_s, orig = pid[idx][sel], fda[idx][sel], orig[sel]
    q_of = np.searchsorted(offsets, orig, side="right") - 1
    local = orig - offsets[q_of]
    cols = {c: np.empty(len(pid_s), dtype=object) for c in keys}
    for qi in np.unique(q_of):
        m = q_of == qi
        z = np.load(os.path.join(d_dir, qs[qi] + ".npz"))
        for c in keys:
            arr = np.asarray(z[c + "_uniq"]).astype(str)[z[c + "_codes"]]
            cols[c][m] = arr[local[m]]
    out = pd.DataFrame({"pid": pid_s, "fda": fda_s})
    for c in keys:
        out[c] = cols[c]
    out["year"] = np.where(out["fda"] == staged.FAR_FUTURE, -1, out["fda"] // 10000).astype(int)
    print(f"report rows read: {n_read:,}; deleted-case records removed: {n_del:,}; "
          f"deduplicated reports: {n_dedup:,}", flush=True)
    return out


def describe(d):
    n = len(d)
    if n == 0:
        return {"n_reports": 0}
    y = d.loc[d["year"] >= 0, "year"]
    vc = y.value_counts()
    occ, rep, ctry = d["occp_cod"], d["rept_cod"], d["occr_country"]

    def pct(mask):
        return round(100 * float(np.mean(mask)), 1)

    return {
        "n_reports": n,
        "first_year": int(y.min()) if len(y) else None,
        "last_year": int(y.max()) if len(y) else None,
        "n_years_with_reports": int(vc.size),
        "largest_year": int(vc[vc == vc.max()].index.min()) if len(vc) else None,
        "largest_year_share_pct": round(100 * vc.max() / len(y), 1) if len(y) else None,
        "expedited_pct": pct(rep == "EXP"),
        "periodic_pct": pct(rep == "PER"),
        "direct_pct": pct(rep == "DIR"),
        "consumer_pct": pct(occ == "CN"),
        "hcp_documented_pct": pct(occ.isin(HCP_DOC)),
        "hp_rn_pct": pct(occ.isin(["HP", "RN"])),
        "reporter_other_missing_pct": pct(~occ.isin(["CN", "HP", "RN"] + HCP_DOC)),
        "us_pct": pct(ctry == "US"),
        "country_missing_pct": pct(ctry == ""),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--faers-dir", required=True)
    ap.add_argument("--out", default="./screen")
    ap.add_argument("--work", default=None, help="work directory of an earlier staged run, to reuse demo5")
    ap.add_argument("--start", default="2013Q1")
    ap.add_argument("--end", default="2025Q4")
    args = ap.parse_args()
    start, end = staged.parse_q(args.start), staged.parse_q(args.end)
    os.makedirs(args.out, exist_ok=True)
    work = args.work or os.path.join(args.out, "_work")
    os.makedirs(work, exist_ok=True)

    files, deleted_files = staged.find_sources(args.faers_dir)
    quarters = sorted(k for k in files if start <= k <= end and all(t in files[k] for t in ("DEMO", "DRUG")))
    if not quarters:
        sys.exit("no complete DEMO/DRUG quarterly files found; check --faers-dir.")
    print(f"{len(quarters)} quarters, {len(deleted_files)} deleted-case list files; "
          f"REAC and INDI are not read", flush=True)

    staged.phase_demo(files, quarters, work)          # skips quarters already in <work>/demo5
    ids, names = scan_drug(files, quarters)
    hits = [v for v in ids.values() if len(v)]
    wanted = np.unique(np.concatenate(hits)) if hits else None
    demo = load_dedup(work, quarters, deleted_files, wanted=wanted)

    rows, yrows, nrows = [], [], []
    kept = {}
    for key, _, _, grp, note in SCREEN:
        pos = np.searchsorted(demo["pid"].to_numpy(), ids[key])
        pos_c = np.minimum(pos, len(demo) - 1)
        ok = demo["pid"].to_numpy()[pos_c] == ids[key]
        d = demo.iloc[pos_c[ok]]
        kept[key] = d["pid"].to_numpy()
        rows.append({"drug": key, "set": grp, "note": note, **describe(d)})
        for y, n in d.loc[d["year"] >= 0, "year"].value_counts().sort_index().items():
            yrows.append({"drug": key, "set": grp, "year": int(y), "n_reports": int(n)})
        nm = names[key]
        nm = nm[nm["pid"].isin(kept[key])].drop_duplicates()
        for s, n in Counter(nm["name"]).most_common(15):
            nrows.append({"drug": key, "set": grp, "drugname_prod_ai": s, "n_reports": n})

    s2 = [k for k, _, _, g, _ in SCREEN if g.startswith("specialty2") or g == "excluded"]
    owner = defaultdict(list)
    for k in s2:
        for p in kept[k]:
            owner[int(p)].append(k)
    ov = Counter(" + ".join(sorted(v)) for v in owner.values() if len(v) > 1)
    orows = [{"drugs": k, "n_reports": n} for k, n in ov.most_common()]

    def save(df, name):
        df.to_csv(os.path.join(args.out, name), index=False, encoding="utf-8-sig")
        print(f"  wrote {name} ({len(df)} rows)", flush=True)

    counts = pd.DataFrame(rows)
    for c in ("first_year", "last_year", "n_years_with_reports", "largest_year"):
        if c in counts:
            counts[c] = counts[c].astype("Int64")
    save(counts, "S1_screen_counts.csv")
    save(pd.DataFrame(yrows), "S2_screen_yearly.csv")
    save(pd.DataFrame(nrows), "S3_screen_names.csv")
    save(pd.DataFrame(orows, columns=["drugs", "n_reports"]), "S4_screen_overlap.csv")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expand a MedDRA SMQ into term lists that can be matched directly against FAERS REAC.pt (one narrow, one broad).

Why expand to the LLT level: each FAERS quarter was coded with the MedDRA version current at the time,
so PTs that were later renamed (e.g. Drug dose omission, Product dose omission) now sit as LLTs under a newer PT.
Matching on the SMQ's PT names alone would miss every report coded with an old name; after expansion they are caught.

Usage:
    python scripts/expand_smq.py --meddra-dir /path/to/MedAscii --smq-code 20000224 --out-dir ./smq
Outputs:
    smq_<code>_narrow.txt   narrow scope (term_scope=2): the PT name plus all of its LLTs, one per line, lower case
    smq_<code>_broad.txt    broad scope (term_scope=1 or 2, i.e. every scope): as above
    smq_<code>_terms.csv    code, name, scope and status of each PT in the SMQ, with the number of LLTs expanded
Then:
    python scripts/faers_me_feasibility_v7_staged.py ... --smq-file narrow=smq/smq_20000224_narrow.txt --smq-file broad=smq/smq_20000224_broad.txt
"""
import argparse
import csv
import os


def read_asc(path):
    with open(path, encoding="latin-1") as f:
        return [line.rstrip("\r\n").split("$") for line in f if line.strip()]


def find(dirpath, name):
    for f in os.listdir(dirpath):
        if f.lower() == name:
            return os.path.join(dirpath, f)
    raise SystemExit(f"Cannot find {name} in {dirpath}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meddra-dir", required=True, help="MedDRA ASCII directory (containing smq_content.asc, llt.asc and pt.asc)")
    ap.add_argument("--smq-code", default="20000224")
    ap.add_argument("--out-dir", default="./smq")
    ap.add_argument("--include-inactive", action="store_true", help="also include inactive terms (term_status other than A); by default only active terms are used")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    content = read_asc(find(args.meddra_dir, "smq_content.asc"))
    llt = read_asc(find(args.meddra_dir, "llt.asc"))
    pt = read_asc(find(args.meddra_dir, "pt.asc"))
    pt_name = {r[0]: r[1] for r in pt}
    llt_name = {r[0]: r[1] for r in llt}
    llt_by_pt = {}
    for r in llt:
        llt_by_pt.setdefault(r[2], []).append(r[1])

    # smq_content: smq_code, term_code, term_level, term_scope, term_category, term_weight, term_status, ...
    by_smq = {}
    for r in content:
        by_smq.setdefault(r[0], []).append(r)

    def collect(code, seen):
        """Collect recursively (sub-SMQs of a hierarchical SMQ appear with term_level=0)."""
        if code in seen:
            return []
        seen.add(code)
        rows = []
        for r in by_smq.get(code, []):
            level, status = r[2], (r[6] if len(r) > 6 else "A")
            if not args.include_inactive and status.upper() != "A":
                continue
            if level == "0":
                rows += collect(r[1], seen)
            else:
                rows.append(r)
        return rows

    rows = collect(args.smq_code, set())
    if not rows:
        raise SystemExit(f"No valid terms found for SMQ {args.smq_code}")

    lists = {"narrow": set(), "broad": set()}
    term_rows = []
    for r in rows:
        code, level, scope = r[1], r[2], r[3]
        if level == "4":
            names = [pt_name.get(code, "")] + llt_by_pt.get(code, [])
            label = pt_name.get(code, code)
        elif level == "5":
            names = [llt_name.get(code, "")]
            label = llt_name.get(code, code)
        else:
            continue
        names = {n.strip().lower() for n in names if n.strip()}
        lists["broad"] |= names
        if scope == "2":
            lists["narrow"] |= names
        term_rows.append({"term_code": code, "term_name": label, "level": "PT" if level == "4" else "LLT",
                          "scope": {"2": "narrow", "1": "broad"}.get(scope, scope),
                          "status": r[6] if len(r) > 6 else "", "n_names_after_expansion": len(names)})

    for k, v in lists.items():
        with open(os.path.join(args.out_dir, f"smq_{args.smq_code}_{k}.txt"), "w", encoding="utf-8") as f:
            f.write("\n".join(sorted(v)) + "\n")
    with open(os.path.join(args.out_dir, f"smq_{args.smq_code}_terms.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(term_rows[0]))
        w.writeheader(); w.writerows(term_rows)
    n_pt = {s: sum(1 for t in term_rows if t["level"] == "PT" and (s == "broad" or t["scope"] == "narrow"))
            for s in ("narrow", "broad")}
    print(f"SMQ {args.smq_code}: narrow scope {n_pt['narrow']} PTs -> {len(lists['narrow'])} names after expansion; "
          f"broad scope {n_pt['broad']} PTs -> {len(lists['broad'])} names")
    for probe in ("product dose omission issue", "drug dose omission", "product dose omission", "intentional dose omission"):
        print(f"  {probe}: narrow={'yes' if probe in lists['narrow'] else 'no'}, broad={'yes' if probe in lists['broad'] else 'no'}")


if __name__ == "__main__":
    main()

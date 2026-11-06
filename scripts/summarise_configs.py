#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cross-configuration summary of the rare-disease vs specialty comparison (version 1.1.0).

Reads the main run (no --groups) and each --groups run, after meta_regression.py has been run in each
folder, and writes one row per configuration and outcome:

    crude_ROR_rare_vs_specialty   report-level ROR (95% CI) from 1_summary_by_drug.csv, POOLED_rare_all
    HK_OR, HK_CI                  drug-level REML + Hartung-Knapp OR, rare vs specialty (meta_regression_HK.csv)
    OR1_main                      the same estimate from the main run, for the same outcome
    ratio_vs_OR1                  HK_OR / OR1_main, both unrounded; left empty for a run in which the
                                  rare-disease drugs are replaced (configs/c6)
    within_0.67_1.5               whether the ratio lies in the prespecified range (rules, section 8; the
                                  prespecified primary quantity is the nonadherence row of the second set)
    cluster_unadjusted, cluster_adjusted
                                  report-level logistic OR with drug-clustered SEs (12_multivariable_logistic.csv),
                                  with the normal-based 95% CI and the CI using t critical values with
                                  (number of drug clusters - 1) degrees of freedom; the t interval is derived
                                  from the published normal interval, which is rounded to two decimals
    subset_ROR_min, subset_ROR_max
                                  range of the pooled ROR across the report subsets of 6_sensitivity.csv
                                  (leave-one-out rows excluded)
    HK_OR_rare_vs_common          drug-level OR against the common comparators, for reference

Usage:
    python scripts/summarise_configs.py --main ./out --config ./out_c1_specialty2_main ./out_c2_... \
        --out cross_config_summary.csv
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

OUTCOMES = ["nonadherence", "smq_narrow", "smq_broad", "any_me", "dosing_admin_error", "access_supply"]
Z = stats.norm.ppf(0.975)


def read(d, name):
    p = os.path.join(d, name)
    return pd.read_csv(p) if os.path.exists(p) else None


def hk(d):
    m = read(d, "meta_regression_HK.csv")
    if m is None:
        raise SystemExit(f"{d}: meta_regression_HK.csv not found; run meta_regression.py in this folder first")
    return m


def hk_row(m, cat, comparison):
    r = m[(m.category == cat) & (m.comparison == comparison)]
    return None if r.empty else r.iloc[0]


def cluster(m12, cat, model):
    if m12 is None:
        return {}
    r = m12[(m12.category == cat) & (m12.model == model) & (m12.contrast == "rare vs specialty")]
    if r.empty:
        return {}
    r = r.iloc[0]
    se = (np.log(r.CI_high) - np.log(r.CI_low)) / (2 * Z)
    df = int(r.n_clusters) - 1
    t = stats.t.ppf(0.975, df)
    lo, hi = np.exp(np.log(r.OR) - t * se), np.exp(np.log(r.OR) + t * se)
    return {f"cluster_{model}_OR": r.OR, f"cluster_{model}_CI": f"{r.CI_low:.2f}-{r.CI_high:.2f}",
            f"cluster_{model}_CI_t": f"{lo:.2f}-{hi:.2f}", "t_df": df}


def config_of(d):
    p = os.path.join(d, "config_used.json")
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            return json.load(f)
    return {}


def label(d):
    return config_of(d).get("name") or os.path.basename(os.path.normpath(d))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--main", required=True, help="results folder of the run without --groups")
    ap.add_argument("--config", nargs="+", required=True, help="results folders of the --groups runs")
    ap.add_argument("--out", default="cross_config_summary.csv")
    args = ap.parse_args()

    main_hk = hk(args.main)
    rows = []
    for d, name in [(args.main, "main (first specialty set)")] + [(d, label(d)) for d in args.config]:
        m = hk(d)
        summ, sens, m12 = read(d, "1_summary_by_drug.csv"), read(d, "6_sensitivity.csv"), read(d, "12_multivariable_logistic.csv")
        for cat in OUTCOMES:
            r, r1 = hk_row(m, cat, "rare vs specialty"), hk_row(main_hk, cat, "rare vs specialty")
            if r is None:
                continue
            rec = {"config": name, "outcome": cat}
            if summ is not None and f"ROR_vs_specialty_{cat}" in summ:
                rec["crude_ROR_rare_vs_specialty"] = summ.loc[summ.drug == "POOLED_rare_all", f"ROR_vs_specialty_{cat}"].iloc[0]
            rec.update({"HK_OR": round(r.OR, 2), "HK_CI": f"{r.CI_low:.2f}-{r.CI_high:.2f}"})
            # a run in which other drugs take the rare role (c6) is not a rare-disease comparison: no ratio
            rare_replaced = config_of(d).get("group_roles", {}).get("rare") == "off"
            if r1 is not None and not rare_replaced:
                ratio = r.OR / r1.OR
                rec.update({"OR1_main": round(r1.OR, 2), "ratio_vs_OR1": round(ratio, 2),
                            "within_0.67_1.5": bool(0.67 <= ratio <= 1.5)})
            rec.update(cluster(m12, cat, "unadjusted"))
            rec.update(cluster(m12, cat, "adjusted"))
            if sens is not None:
                s = sens[(sens.category == cat) & (sens.reference == "specialty") & ~sens.scenario.str.startswith("leave-one-out")]
                if len(s):
                    rec.update({"subset_ROR_min": s.ROR.min(), "subset_ROR_max": s.ROR.max()})
            rc = hk_row(m, cat, "rare vs common")
            if rc is not None:
                rec["HK_OR_rare_vs_common"] = f"{rc.OR:.2f} ({rc.CI_low:.2f}-{rc.CI_high:.2f})"
            rows.append(rec)
    out = pd.DataFrame(rows)
    out.to_csv(args.out, index=False, encoding="utf-8-sig")
    print(f"wrote {args.out} ({len(out)} rows)")


if __name__ == "__main__":
    main()

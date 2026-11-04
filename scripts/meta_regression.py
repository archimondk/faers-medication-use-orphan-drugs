#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Drug-level meta-regression of event proportions (Table 4 of the article).

For each outcome, the logit of each drug's event proportion is regressed on comparator-group indicators
with a random effect for drug: restricted maximum likelihood (REML) estimate of tau^2 and
Hartung-Knapp adjusted standard errors (scaling factor truncated at 1), t distribution with k - p df.
A 0.5 continuity correction is added to drugs with 0 events or 100% events.

Input:  7_drug_level.csv from faers_me_feasibility_v7_staged.py (agg phase)
Usage:
    python scripts/meta_regression.py --drug-level out/7_drug_level.csv --out out/meta_regression_HK.csv
"""
import argparse

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar

OUTCOMES = ["nonadherence", "any_me", "dosing_admin_error", "storage_handling", "access_supply",
            "smq_narrow", "smq_broad"]
GROUPS = ("rare", "specialty", "common")


def reml_metareg(y, v, X):
    k, p = X.shape

    def negll(log_t2):
        t2 = np.exp(log_t2)
        W = np.diag(1 / (v + t2))
        XtWX = X.T @ W @ X
        b = np.linalg.solve(XtWX, X.T @ W @ y)
        r = y - X @ b
        return 0.5 * (np.sum(np.log(v + t2)) + np.log(np.linalg.det(XtWX)) + r @ W @ r)

    res = minimize_scalar(negll, bounds=(-20, 5), method="bounded")
    t2 = np.exp(res.x)
    if negll(-20) <= res.fun + 1e-9:
        t2 = 0.0
    W = np.diag(1 / (v + t2))
    XtWX = X.T @ W @ X
    b = np.linalg.solve(XtWX, X.T @ W @ y)
    r = y - X @ b
    q = max(1.0, (r @ W @ r) / (k - p))          # Hartung-Knapp scaling, truncated at 1
    se = np.sqrt(np.diag(q * np.linalg.inv(XtWX)))
    tcrit = stats.t.ppf(0.975, k - p)
    pval = 2 * stats.t.sf(np.abs(b / se), k - p)
    return b, se, tcrit, pval, t2, k - p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug-level", required=True, help="7_drug_level.csv")
    ap.add_argument("--out", default="meta_regression_HK.csv")
    ap.add_argument("--min-n", type=int, default=100, help="minimum reports per drug (default 100)")
    args = ap.parse_args()

    d = pd.read_csv(args.drug_level)
    d = d[(d.n_reports >= args.min_n) & d.group.isin(GROUPS)]
    rows = []
    for c in [o for o in OUTCOMES if o in set(d.category)]:
        s = d[d.category == c]
        ev, n = s.events.values.astype(float), s.n_reports.values.astype(float)
        corr = (ev == 0) | (ev == n)
        a = np.where(corr, ev + 0.5, ev)
        bb = np.where(corr, n - ev + 0.5, n - ev)
        y, v = np.log(a / bb), 1 / a + 1 / bb
        for ref in ("specialty", "common"):
            other = "common" if ref == "specialty" else "specialty"
            X = np.column_stack([np.ones(len(s)), (s.group == "rare").astype(float), (s.group == other).astype(float)])
            b, se, tc, pv, t2, df = reml_metareg(y, v, X)
            rows.append({"category": c, "comparison": f"rare vs {ref}", "OR": np.exp(b[1]),
                         "CI_low": np.exp(b[1] - tc * se[1]), "CI_high": np.exp(b[1] + tc * se[1]),
                         "p": pv[1], "tau2": t2, "df": df})
        for g in GROUPS:
            gg = s[s.group == g]
            rows.append({"category": c, "comparison": f"{g} median % (range)", "OR": gg.pct.median(),
                         "CI_low": gg.pct.min(), "CI_high": gg.pct.max(), "p": np.nan, "tau2": np.nan, "df": len(gg)})
    out = pd.DataFrame(rows)
    out.to_csv(args.out, index=False, encoding="utf-8-sig")
    with pd.option_context("display.width", 200):
        print(out.round(3).to_string())


if __name__ == "__main__":
    main()

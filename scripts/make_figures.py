#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 1 (comparator choice across sensitivity analyses), the flow of FAERS records and the annual
nonadherence reports by drug.

Inputs, all from the agg phase of faers_me_feasibility_v7_staged.py:
    1_summary_by_drug.csv   whole-database ROR for the pooled rare-disease group
    6_sensitivity.csv       pooled ROR against each comparator group within report subsets
    7_drug_level.csv        which rare-disease drugs met the minimum report count
    4_yearly_by_drug.csv    yearly counts
    14_flow_and_characteristics.csv  record counts for the flow diagram (version 1.1.1; skipped if absent)
Usage:
    python scripts/make_figures.py --results ./out --out ./figures
    python scripts/make_figures.py --results ./out --out ./figures --second-specialty ./out_c1_specialty2_main
With --second-specialty (version 1.1.0), Figure 1 adds the estimates against the second specialty set,
read from 6_sensitivity.csv of that run; without it the figure is drawn exactly as in version 1.0.2.
Version 1.1.1 adds the flow diagram (FigureS1_flow.png; the second-set box is drawn when
--second-specialty is given) and writes the annual figure as FigureS2_annual_nonadherence.png,
with the figure number removed from its title.
"""
import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Scenario rows of 6_sensitivity.csv shown in Figure 1, in display order. The labels are the ones
# written by the current analysis scripts; a run with older, differently labelled output will not
# match here and those scenarios will be skipped.
SCENARIOS = [
    (("all reports",), "All reports"),
    (("excluding periodic reports (PER)",), "Excluding periodic reports"),
    (("health-professional reports only",), "Health-professional reports"),
    (("health-professional reports only, documented codes (MD, PH, OT)",), "Health-professional, documented codes"),
    (("expedited reports from health professionals only",), "Expedited health-professional reports"),
    (("excluding the peak nonadherence year of each drug",), "Excluding rare drugs' peak year"),
    (("excluding the peak nonadherence year of every drug in all groups",), "Excluding every drug's peak year"),
    (("restricted to the rare-disease indication",), "Target indication only"),
]
NAME = {
    "penicillamine": "Penicillamine", "trientine_all": "Trientine (all salts)", "zinc_acetate_galzin": "Zinc acetate (Galzin)",
    "nitisinone_all": "Nitisinone (all products)", "cysteamine_oral": "Cysteamine, oral (all products)",
    "carglumic_acid": "Carglumic acid", "sapropterin": "Sapropterin", "sepiapterin": "Sepiapterin", "tiopronin": "Tiopronin",
    "miglustat": "Miglustat", "eliglustat": "Eliglustat", "sodium_phenylbutyrate": "Sodium phenylbutyrate",
    "glycerol_phenylbutyrate": "Glycerol phenylbutyrate", "betaine_anhydrous": "Betaine anhydrous", "cholic_acid": "Cholic acid",
}
INK, MUTED, GRID, BLUE = "#2b2b2a", "#6b6a64", "#e4e3de", "#2a78d6"


def parse_ror(text):
    """'2.68 (2.53-2.84)' -> (2.68, 2.53, 2.84)"""
    m = re.match(r"\s*([\d.]+)\s*\(([\d.]+)\s*-\s*([\d.]+)\)", str(text))
    return tuple(float(x) for x in m.groups()) if m else None


def style_axis(ax):
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.set_axisbelow(True)


def figure1(summary, sens, path, sens2=None):
    pooled = summary[summary.drug == "POOLED_rare_all"].iloc[0]
    whole = {c: parse_ror(pooled.get(f"ROR_vs_allFAERS_{c}")) for c in ("nonadherence", "access_supply")}
    style = {"whole": ("Whole database", "#8a8a85", "D"), "controls": ("Common chronic-disease drugs", "#1baf7a", "^"),
             "specialty": ("Specialty-pharmacy drugs", "#eb6834", "s")}
    offs = {"whole": 0.24, "controls": 0.0, "specialty": -0.24}
    sources = {"controls": sens, "specialty": sens}
    if sens2 is not None:
        style["specialty"] = ("Specialty-pharmacy drugs, first set", "#eb6834", "s")
        style["specialty2"] = ("Specialty-pharmacy drugs, second set", "#7b4fc9", "v")
        offs = {"whole": 0.3, "controls": 0.1, "specialty": -0.1, "specialty2": -0.3}
        sources["specialty2"] = sens2
    refs = [r for r in ("controls", "specialty", "specialty2") if r in sources]
    scen = [(labels, en) for labels, en in SCENARIOS if sens.scenario.isin(labels).any()]
    if sens2 is not None:
        scen = [(labels, en) for labels, en in scen if sens2.scenario.isin(labels).any()]
    height = max(3.9, 1.5 + (0.38 if sens2 is None else 0.48) * len(scen))
    fig, axes = plt.subplots(1, 2, figsize=(9.6, height), dpi=200, sharey=True)
    for ax, (cat, title) in zip(axes, [("nonadherence", "A  Nonadherence"), ("access_supply", "B  Access or supply problem")]):
        n = len(scen)
        for i, (labels, _) in enumerate(scen):
            y = n - 1 - i
            if i == 0 and whole[cat]:
                r, lo, hi = whole[cat]
                ax.plot([lo, hi], [y + offs["whole"]] * 2, color=style["whole"][1], lw=1.6)
                ax.plot(r, y + offs["whole"], marker=style["whole"][2], ms=5.5, color=style["whole"][1], mec="white", mew=1)
            for ref in refs:
                src = sources[ref]
                row = src[src.scenario.isin(labels) & (src.reference == ref.replace("specialty2", "specialty"))
                          & (src.category == cat)].iloc[0]
                yy = y + offs[ref]
                ax.plot([row.CI_low, row.CI_high], [yy, yy], color=style[ref][1], lw=1.6)
                ax.plot(row.ROR, yy, marker=style[ref][2], ms=6, color=style[ref][1], mec="white", mew=1)
        ax.axvline(1, color=MUTED, lw=1, ls="--", zorder=0)
        ax.set_xscale("log")
        ax.set_xlim(0.5, 100)
        ax.set_xticks([0.5, 1, 2, 5, 10, 20, 50, 100])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.set_title(title, loc="left", fontsize=9.5, fontweight="bold", color=INK)
        ax.tick_params(axis="x", labelsize=8, colors=MUTED, length=0)
        ax.tick_params(axis="y", length=0)
        ax.minorticks_off()
        ax.grid(axis="x", color=GRID, lw=0.8)
        style_axis(ax)
        ax.set_ylim(-0.6, len(scen) - 0.4)
    axes[0].set_yticks(range(len(scen)))
    axes[0].set_yticklabels([en for _, en in scen][::-1], fontsize=8, color=INK)
    handles = [plt.Line2D([], [], marker=m, color=c, ls="", ms=6, mec="white") for _, c, m in style.values()]
    fig.supxlabel("Reporting odds ratio, rare-disease drugs vs comparator (log scale)", fontsize=8, color=INK, x=0.62, y=0.5 / height + 0.035)
    fig.legend(handles, [l for l, _, _ in style.values()], loc="lower center", ncol=3 if sens2 is None else 2,
               frameon=False, fontsize=8, bbox_to_anchor=(0.6, 0.0))
    fig.text(0.01, 0.005, "Whole-database comparison shown for all reports only. Bars: 95% CI.", fontsize=6.5, color=MUTED)
    fig.subplots_adjust(left=0.27, right=0.98, top=1 - 0.4 / height,
                        bottom=(1.05 if sens2 is None else 1.3) / height, wspace=0.16)
    fig.savefig(path, facecolor="white")
    plt.close(fig)


def figure_s1(yearly, drug_level, path, first_year=2013, last_year=2025):
    keep = sorted(set(drug_level[drug_level.group == "rare"].drug))
    totals = {d: yearly[yearly.drug == d].n_nonadherence.sum() for d in keep}
    keep = sorted(keep, key=lambda d: -totals[d])
    years = list(range(first_year, last_year + 1))
    ncol = 4
    nrow = -(-len(keep) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(10, 2.3 * nrow), dpi=200, sharex=True, squeeze=False)
    for ax, d in zip(axes.flat, keep):
        vals = (yearly[yearly.drug == d][["year", "n_nonadherence"]].set_index("year")
                .reindex(years, fill_value=0).n_nonadherence.values)
        ax.bar(years, vals, width=0.78, color=BLUE, edgecolor="white", linewidth=0.6)
        ax.set_title(f"{NAME.get(d, d)}  (n = {int(totals[d])})", fontsize=8, color=INK, loc="left")
        ax.set_ylim(0, max(4, vals.max() * 1.15))
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3, integer=True))
        ax.tick_params(labelsize=7, colors=MUTED, length=0)
        ax.grid(axis="y", color=GRID, lw=0.6)
        style_axis(ax)
        ax.set_xticks([2013, 2017, 2021, 2025])
    for ax in list(axes.flat)[len(keep):]:
        ax.axis("off")
    for ax in axes.flat:
        ax.tick_params(labelbottom=True)
    for ax in axes[:, 0]:
        ax.set_ylabel("Reports", fontsize=7.5, color=MUTED)
    fig.suptitle("Nonadherence reports per year, rare-disease drugs", x=0.02, ha="left", fontsize=10,
                 color=INK, fontweight="bold")
    fig.text(0.02, 0.005, "Year of FDA receipt. Panels ordered by total nonadherence reports; y-axes differ between panels. "
             "Drugs with at least 100 reports.", fontsize=7, color=MUTED)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    fig.savefig(path, facecolor="white")
    plt.close(fig)


def figure_flow(flow, path, flow2=None):
    """Flow of FAERS records from the quarterly files to the analysed groups (READUS-PV item 8a)."""
    from matplotlib.patches import FancyBboxPatch
    f = flow[flow.section == "flow"]
    val = lambda g, c: int(f[(f.group.fillna("") == g) & (f.characteristic == c)].n.iloc[0])
    n_read, n_del, n_dedup = (val("", "DEMO records read"), val("", "records removed (FDA deleted-case lists)"),
                              val("", "reports after deduplication"))
    minc = [c for c in f.characteristic.unique() if c.startswith("drugs with at least")][0]
    thr = minc.split()[-2]
    groups = [("rare", "Rare-disease drugs", f), ("specialty", "First specialty set" if flow2 is not None else "Specialty-pharmacy drugs", f)]
    if flow2 is not None:
        f2 = flow2[flow2.section == "flow"]
        groups.append(("specialty", "Second specialty set", f2))
    groups += [("common", "Common comparators", f), ("example", "Nefecon (example)", f)]
    INK, MUTED, LINE, FILL = "#222222", "#555555", "#9a9a9a", "#f4f4f2"
    fig, ax = plt.subplots(figsize=(9.6, 7.4), dpi=200)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    def box(x, y, w, h, lines, fill=FILL, size=8.6):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=1.2", fc=fill, ec=LINE, lw=0.9))
        y0 = y + h / 2 + (len(lines) - 1) * 1.275
        for i, ln in enumerate(lines):
            ax.text(x + w / 2, y0 - i * 2.55, ln, ha="center", va="center", fontsize=size,
                    color=INK if i == 0 else MUTED, fontweight="bold" if i == 0 else "normal")

    def arrow(x1, y1, x2, y2):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="-|>", color=LINE, lw=1.0, shrinkA=0, shrinkB=0))
    box(18, 88, 46, 9, ["FAERS quarterly ASCII files", f"{n_read:,} DEMO records"])
    arrow(41, 88, 41, 81.5)
    box(18, 71, 46, 10, ["After removal of deleted cases", f"{n_read - n_del:,} records"])
    box(68, 80.75, 31, 7.5, [f"Removed: {n_del:,} records", "of cases in FDA deleted-case lists"], fill="white", size=7.2)
    arrow(41, 84.5, 67.5, 84.5)
    arrow(41, 71, 41, 64.5)
    box(18, 54, 46, 10, ["After deduplication by CASEID", f"{n_dedup:,} reports"])
    box(68, 63.25, 31, 8.5, [f"Removed: {n_read - n_del - n_dedup:,} records",
                             "earlier versions of the same case", "(latest version kept)"],
        fill="white", size=7.2)
    arrow(41, 67.5, 67.5, 67.5)
    arrow(41, 54, 41, 47.5)
    box(10, 39.5, 62, 8, ["Primary-suspect attribution by generic and brand-name patterns"], size=8.3)
    w = 96 / len(groups) - 2
    for i, (g, lab, src) in enumerate(groups):
        x = 1 + i * (w + 2)
        nd, nr = val_src(src, g, "drugs"), val_src(src, g, "primary-suspect reports")
        nk = val_src(src, g, minc)
        lines = [lab, f"{nd} drug{'s' if nd != 1 else ''}", f"{nr:,} reports", f"{nk} with \u2265{thr} reports"]
        box(x, 12, w, 18, lines, size=7.0)
        arrow(41, 39.5, x + w / 2, 30.5)
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def val_src(f, g, c):
    return int(f[(f.group.fillna("") == g) & (f.characteristic == c)].n.iloc[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, help="directory with the agg-phase CSV files")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--second-specialty", default=None,
                    help="results folder of the run with configs/c1_specialty2_main.json (adds it to Figure 1)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rd = lambda f: pd.read_csv(os.path.join(args.results, f))
    sens2 = pd.read_csv(os.path.join(args.second_specialty, "6_sensitivity.csv")) if args.second_specialty else None
    figure1(rd("1_summary_by_drug.csv"), rd("6_sensitivity.csv"), os.path.join(args.out, "Figure1_comparator_sensitivity.png"),
            sens2=sens2)
    figure_s1(rd("4_yearly_by_drug.csv"), rd("7_drug_level.csv"), os.path.join(args.out, "FigureS2_annual_nonadherence.png"))
    flow_path = os.path.join(args.results, "14_flow_and_characteristics.csv")
    if os.path.exists(flow_path):
        flow2 = None
        if args.second_specialty and os.path.exists(os.path.join(args.second_specialty, "14_flow_and_characteristics.csv")):
            flow2 = pd.read_csv(os.path.join(args.second_specialty, "14_flow_and_characteristics.csv"))
        figure_flow(pd.read_csv(flow_path), os.path.join(args.out, "FigureS1_flow.png"), flow2=flow2)
    else:
        print("14_flow_and_characteristics.csv not found (results from before version 1.1.1); flow diagram skipped")
    print(f"figures written to {args.out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 1 (comparator choice across sensitivity analyses) and Figure S1 (annual nonadherence reports).

Inputs, all from the agg phase of faers_me_feasibility_v7_staged.py:
    1_summary_by_drug.csv   whole-database ROR for the pooled rare-disease group
    6_sensitivity.csv       pooled ROR against each comparator group within report subsets
    7_drug_level.csv        which rare-disease drugs met the minimum report count
    4_yearly_by_drug.csv    yearly counts
Usage:
    python scripts/make_figures.py --results ./out --out ./figures
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
    (("expedited reports from health professionals only",), "Expedited health-professional reports"),
    (("excluding the peak nonadherence year of each drug",), "Excluding each drug's peak year"),
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


def figure1(summary, sens, path):
    pooled = summary[summary.drug == "POOLED_rare_all"].iloc[0]
    whole = {c: parse_ror(pooled.get(f"ROR_vs_allFAERS_{c}")) for c in ("nonadherence", "access_supply")}
    style = {"whole": ("Whole database", "#8a8a85", "D"), "controls": ("Common chronic-disease drugs", "#1baf7a", "^"),
             "specialty": ("Specialty-pharmacy drugs", "#eb6834", "s")}
    offs = {"whole": 0.24, "controls": 0.0, "specialty": -0.24}
    scen = [(labels, en) for labels, en in SCENARIOS if sens.scenario.isin(labels).any()]
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9), dpi=200, sharey=True)
    for ax, (cat, title) in zip(axes, [("nonadherence", "A  Nonadherence"), ("access_supply", "B  Access or supply problem")]):
        n = len(scen)
        for i, (labels, _) in enumerate(scen):
            y = n - 1 - i
            if i == 0 and whole[cat]:
                r, lo, hi = whole[cat]
                ax.plot([lo, hi], [y + offs["whole"]] * 2, color=style["whole"][1], lw=1.6)
                ax.plot(r, y + offs["whole"], marker=style["whole"][2], ms=5.5, color=style["whole"][1], mec="white", mew=1)
            for ref in ("controls", "specialty"):
                row = sens[sens.scenario.isin(labels) & (sens.reference == ref) & (sens.category == cat)].iloc[0]
                yy = y + offs[ref]
                ax.plot([row.CI_low, row.CI_high], [yy, yy], color=style[ref][1], lw=1.6)
                ax.plot(row.ROR, yy, marker=style[ref][2], ms=6, color=style[ref][1], mec="white", mew=1)
        ax.axvline(1, color=MUTED, lw=1, ls="--", zorder=0)
        ax.set_xscale("log")
        ax.set_xlim(0.5, 50)
        ax.set_xticks([0.5, 1, 2, 5, 10, 20, 50])
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
    fig.supxlabel("Reporting odds ratio, rare-disease drugs vs comparator (log scale)", fontsize=8, color=INK, x=0.62, y=0.13)
    fig.legend(handles, [l for l, _, _ in style.values()], loc="lower center", ncol=3, frameon=False, fontsize=8,
               bbox_to_anchor=(0.6, 0.0))
    fig.text(0.01, 0.005, "Whole-database comparison shown for all reports only. Bars: 95% CI.", fontsize=6.5, color=MUTED)
    fig.subplots_adjust(left=0.27, right=0.98, top=0.9, bottom=0.27, wspace=0.16)
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
    fig.suptitle("Figure S1. Nonadherence reports per year, rare-disease drugs", x=0.02, ha="left", fontsize=10,
                 color=INK, fontweight="bold")
    fig.text(0.02, 0.005, "Year of FDA receipt. Panels ordered by total nonadherence reports; y-axes differ between panels. "
             "Drugs with at least 100 reports.", fontsize=7, color=MUTED)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    fig.savefig(path, facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, help="directory with the agg-phase CSV files")
    ap.add_argument("--out", default="figures")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rd = lambda f: pd.read_csv(os.path.join(args.results, f))
    figure1(rd("1_summary_by_drug.csv"), rd("6_sensitivity.csv"), os.path.join(args.out, "Figure1_comparator_sensitivity.png"))
    figure_s1(rd("4_yearly_by_drug.csv"), rd("7_drug_level.csv"), os.path.join(args.out, "FigureS1_annual_nonadherence.png"))
    print(f"figures written to {args.out}")


if __name__ == "__main__":
    main()

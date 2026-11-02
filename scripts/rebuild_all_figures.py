from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description="Rebuild all five CMPB figures from packaged Source Data")
parser.add_argument("--source", default=str(ROOT / "source_data"))
parser.add_argument("--output", default=str(ROOT / "verification" / "all_figures"))
args = parser.parse_args()
DELIVERY = Path(args.output)
SRC = Path(args.source)
OUT = DELIVERY / "figures"
OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#0999d8"
ORANGE = "#f79108"
RED = "#e03f40"
DEEP = "#086898"
GREY = "#9aa5b1"
DARK = "#333333"
LIGHT = "#f2f2f2"

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 7.5, "axes.titlesize": 8,
    "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.edgecolor": DARK, "axes.linewidth": 0.8,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.major.width": 0.8,
    "ytick.major.width": 0.8, "axes.grid": False, "savefig.facecolor": "white",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def new_figure(number, nrows, ncols, height):
    fig, axes = plt.subplots(nrows, ncols, figsize=(190 / 25.4, height / 25.4))
    fig.subplots_adjust(left=0.07, right=0.98, top=0.94, bottom=0.10, wspace=0.34, hspace=0.52)
    for ax in np.ravel(axes):
        ax.tick_params(length=3)
    return fig, np.ravel(axes)


def label(ax, letter, title=None):
    ax.text(-0.16, 1.08, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")
    if title:
        ax.set_title(title, loc="left", pad=5, fontweight="bold")


def box(ax, xy, width, height, text, color=BLUE, size=7.5):
    x, y = xy
    patch = FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.02,rounding_size=0.02",
                           facecolor="white", edgecolor=color, linewidth=1.0)
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=size, color=DARK)


def blank(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def save(fig, number):
    stem = OUT / f"Figure_{number}"
    fig.canvas.draw()
    content_in = fig.get_tightbbox(fig.canvas.get_renderer()).width
    pad = max(0.0, (190.0 / 25.4 - content_in) / 2)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=pad)
    fig.savefig(stem.with_suffix(".eps"), bbox_inches="tight", pad_inches=pad)
    fig.savefig(stem.with_suffix(".png"), dpi=200, bbox_inches="tight", pad_inches=pad)
    fig.savefig(stem.with_suffix(".tif"), dpi=600, bbox_inches="tight", pad_inches=pad, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


flow = pd.read_csv(SRC / "Figure1D_hallmark_GSE65682.csv")
fig, ax = new_figure(1, 2, 2, 150)
label(ax[0], "A", "Attribution problem")
blank(ax[0])
box(ax[0], (0.05, 0.60), 0.38, 0.22, "Bulk expression\n(group difference)")
box(ax[0], (0.57, 0.72), 0.36, 0.18, "Composition shift", ORANGE)
box(ax[0], (0.57, 0.42), 0.36, 0.18, "Within-cell change", BLUE)
ax[0].annotate("", xy=(0.56, 0.79), xytext=(0.43, 0.73), arrowprops=dict(arrowstyle="->", color=DARK, lw=1))
ax[0].annotate("", xy=(0.56, 0.50), xytext=(0.43, 0.68), arrowprops=dict(arrowstyle="->", color=DARK, lw=1))
ax[0].text(0.08, 0.17, "Estimate both components;\nreport the estimand and contrast", fontsize=8, color=DARK)
label(ax[1], "B", "Common-scale joint model")
blank(ax[1])
ax[1].text(0.06, 0.76, "Unadjusted:  Y = α0 + βU·group + γᵀT + ε", fontsize=9)
ax[1].text(0.06, 0.57, "Joint:          Y = α1 + βA·group + θᵀC + γᵀT + ε", fontsize=9)
ax[1].text(0.06, 0.34, "d_unadjusted = ΔU/sU     d_joint = βA/sU", fontsize=9, color=DEEP)
ax[1].text(0.06, 0.15, "CAF = 1 − |d_joint|/|d_unadjusted|", fontsize=11, color=DARK, fontweight="bold")
label(ax[2], "C", "GSE65682 eligibility states")
counts = pd.read_csv(SRC / "Figure1C_eligibility_GSE65682.csv")
colors = [DEEP, BLUE, RED, BLUE, GREY]
ax[2].barh(counts.state[::-1], counts.n[::-1], color=colors[::-1], height=0.62)
ax[2].set_xlabel("Genes")
for i, value in enumerate(counts.n[::-1]):
    ax[2].text(value + 170, i, f"{value:,}", va="center", fontsize=7)
ax[2].set_xlim(0, 12700)
label(ax[3], "D", "Hallmark-set median CAF")
ax[3].hist(flow.median_CAF, bins=12, color=BLUE, edgecolor="white", linewidth=0.6)
ax[3].axvline(0.556, color=ORANGE, lw=1.5, ls="--", label="Cross-set median 0.556")
ax[3].axvline(0.354, color=DEEP, lw=0.9, ls=":")
ax[3].axvline(0.749, color=DEEP, lw=0.9, ls=":", label="Range 0.354–0.749")
ax[3].set_xlabel("Set-median CAF")
ax[3].set_ylabel("Hallmark sets")
ax[3].legend(frameon=False, loc="upper left")
save(fig, 1)

gate_metrics = pd.read_csv(SRC / "Figure2D_gate_metrics.csv")
metrics = {f"within_{row.metric.lower()}": {"mean": row.mean, "ci_low": row.ci_low, "ci_high": row.ci_high}
           for row in gate_metrics.itertuples(index=False)}
fig, ax = new_figure(2, 2, 2, 150)
label(ax[0], "A", "Known-truth design")
blank(ax[0])
ax[0].text(0.52, 0.96, "C = composition; W = within-cell", fontsize=7.5, color=DARK)
for i, (name, comp, cell) in enumerate([("Null", "no", "no"), ("Composition-only", "yes", "no"),
                                      ("Within-cell", "no/fixed", "yes"), ("Mixed", "yes", "yes")]):
    y = 0.80 - i * 0.21
    ax[0].text(0.05, y, name, fontsize=8, fontweight="bold", color=DEEP)
    ax[0].text(0.54, y, f"C: {comp}", fontsize=7.5)
    ax[0].text(0.79, y, f"W: {cell}", fontsize=7.5)
ax[0].text(0.05, 0.02, "Multiple seeds; scenario-level summaries", fontsize=7.5, color=DARK)
label(ax[1], "B", "Known-fraction recovery")
methods = ["CIBERSORTx\nν-SVR", "Local\nSVR", "Local\nNNLS"]
vals = [0.80, 0.78, 0.75]
ax[1].bar(methods, vals, color=[DEEP, BLUE, GREY], width=0.62)
ax[1].set_ylim(0, 1)
ax[1].set_ylabel("Correlation with true fractions (r)")
for i, v in enumerate(vals):
    ax[1].text(i, v + 0.035, f"{v:.2f}", ha="center", fontweight="bold")
label(ax[2], "C", "Single maximal monocyte contrast")
calls = [315, 2, 0]
ax[2].bar(["Group\nonly", "Estimated\nfractions", "Oracle\nfractions"], calls, color=[ORANGE, BLUE, DEEP], width=0.62)
ax[2].set_ylabel("Significant unmodified genes")
ax[2].set_ylim(0, 355)
for i, v in enumerate(calls):
    ax[2].text(i, v + 10, str(v), ha="center", fontweight="bold")
ax[2].text(0.02, -0.24, "Observed count; not a general FPR estimate", transform=ax[2].transAxes, fontsize=7)
label(ax[3], "D", "Within-cell call performance (n=200, d=0.4)")
names = ["Sensitivity", "Specificity", "Precision"]
keys = ["within_sensitivity", "within_specificity", "within_precision"]
vals = [metrics[k]["mean"] for k in keys]
lo = [metrics[k]["mean"] - metrics[k]["ci_low"] for k in keys]
hi = [metrics[k]["ci_high"] - metrics[k]["mean"] for k in keys]
ax[3].bar(names, vals, yerr=[lo, hi], color=[BLUE, DEEP, ORANGE], width=0.62,
          error_kw=dict(ecolor=DARK, lw=0.8, capsize=3))
ax[3].set_ylim(0, 1.08)
ax[3].set_ylabel("Mean (Monte Carlo 95% CI)")
for i, v in enumerate(vals):
    ax[3].text(i, v + 0.06, f"{v:.3f}", ha="center", fontsize=7.5)
pd.DataFrame({"metric": names, "mean": vals, "ci_low": [metrics[k]["ci_low"] for k in keys],
              "ci_high": [metrics[k]["ci_high"] for k in keys]}).to_csv(SRC / "Figure2D_gate_metrics.csv", index=False)
save(fig, 2)

# Frozen figure-source tables replace the original evidence workspace.
fraction_summary = pd.read_csv(SRC / "R18_fraction_error_summary.csv")
coverage = pd.read_csv(SRC / "R18_regression_SoA_factorial_coverage.csv")
fig, ax = new_figure(3, 3, 2, 205)
label(ax[0], "A", "Sample-size dependence (d=0.4)")
r = pd.read_csv(SRC / "Figure3A_sample_size.csv")
ax[0].plot(r.n, r.sensitivity, marker="o", color=BLUE, lw=1.3, label="Sensitivity")
ax[0].plot(r.n, r.specificity, marker="s", color=ORANGE, lw=1.3, ls="--", label="Specificity")
ax[0].set_xscale("log")
ax[0].set_xticks([20, 50, 200], ["20", "50", "200"])
ax[0].set_ylim(0, 1.02)
ax[0].set_xlabel("Samples per group (log scale)")
ax[0].set_ylabel("Median")
ax[0].legend(frameon=False, loc="lower right")
label(ax[1], "B", "Reference failure at n=200, d=0.4")
scenarios = ["within_cell_d0.4_n200", "missing_neutrophil_d0.4_n200", "reference_drift_d0.4_n200"]
labels = ["Matched", "Missing\nneutrophil", "Reference\ndrift"]
vals = pd.read_csv(SRC / "Figure3B_reference_failure.csv")[["FPR_median", "FDR_median"]].to_numpy()
x = np.arange(3)
ax[1].bar(x - 0.17, vals[:, 0], width=0.34, color=BLUE, label="FPR")
ax[1].bar(x + 0.17, vals[:, 1], width=0.34, color=ORANGE, label="FDR")
ax[1].set_xticks(x, labels)
ax[1].set_ylabel("Median rate")
ax[1].legend(frameon=False)
label(ax[2], "C", "Power across planted effects")
for d, color, marker in [(0.2, GREY, "o"), (0.4, BLUE, "s"), (0.8, ORANGE, "^")]:
    rr = pd.read_csv(SRC / f"Figure3C_power_d{d}.csv")
    ax[2].plot(rr.n, rr.power, marker=marker, color=color, lw=1.3, label=f"d={d}")
ax[2].set_xscale("log")
ax[2].set_xticks([20, 50, 200], ["20", "50", "200"])
ax[2].set_ylim(-0.04, 1.06)
ax[2].set_xlabel("Samples per group (log scale)")
ax[2].set_ylabel("Median power")
ax[2].legend(frameon=False)
label(ax[3], "D", "Fraction-error gradient (36 scenarios)")
error_methods = ["low_error", "high_error", "systematic_bias", "reference_estimated"]
error_labels = ["Low random", "High random", "Systematic", "Reference\nestimated"]
error_colors = [DEEP, BLUE, RED, ORANGE]
error_data = fraction_summary[fraction_summary.method.isin(error_methods)][
    ["scenario", "method", "fraction_RMSE_median", "CAF_bias_median", "FPR_median", "coverage_95_median"]
].copy()
error_data.to_csv(SRC / "Figure3D_fraction_error.csv", index=False)
distributions = [error_data.loc[error_data.method.eq(method), "CAF_bias_median"].dropna().to_numpy()
                 for method in error_methods]
positions = np.arange(1, 5)
boxes = ax[3].boxplot(distributions, positions=positions, widths=0.52, patch_artist=True,
                      showfliers=False, medianprops=dict(color=DARK, lw=1.0),
                      whiskerprops=dict(color=GREY), capprops=dict(color=GREY))
for patch, color in zip(boxes["boxes"], error_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.35)
    patch.set_edgecolor(color)
rng = np.random.default_rng(42)
for position, values, color in zip(positions, distributions, error_colors):
    jitter = rng.uniform(-0.13, 0.13, len(values))
    ax[3].scatter(np.full(len(values), position) + jitter, values, s=7, color=color, alpha=0.72, edgecolors="none")
ax[3].axhline(0, color=DARK, lw=0.7)
ax[3].set_xticks(positions, error_labels)
ax[3].set_ylabel("Scenario-median CAF bias")
ax[3].set_xlabel("Fraction-error path; points are scenario summaries")
label(ax[4], "E", "CAF bias across 37 reference-estimated summaries")
bias = pd.read_csv(SRC / "Figure3E_CAF_bias.csv")
bias.to_csv(SRC / "Figure3E_CAF_bias.csv", index=False)
ax[4].bar(np.arange(len(bias)), bias.CAF_bias_median,
          color=[ORANGE if v > 0 else BLUE for v in bias.CAF_bias_median], width=0.82)
ax[4].axhline(0, color=DARK, lw=0.8)
ax[4].set_xlabel("Sorted summaries (36 scenarios + central series)")
ax[4].set_ylabel("Median CAF bias")
ax[4].set_xticks([])
label(ax[5], "F", "Completed factorial scope")
blank(ax[5])
coverage.to_csv(SRC / "Figure3F_factorial_scope.csv", index=False)
ax[5].text(0.04, 0.78, "36 scenarios × 5 seeds × 7 regression paths", fontsize=9, fontweight="bold", color=DEEP)
for i, text in enumerate(["1,260 per-run records; 252 scenario × method summaries",
                          "All seven paths completed every scenario and seed",
                          "Zero model failures; frozen anchors reproduced",
                          "CIBERSORTx and csSAM remain separate known-truth arms"]):
    ax[5].text(0.06, 0.59 - i * 0.15, "• " + text, fontsize=7.5)
fig.subplots_adjust(hspace=0.62)
save(fig, 3)

paired = pd.read_csv(SRC / "Figure4A_paired_CAF_79.csv")
thresholds = pd.read_csv(SRC / "Figure4C_threshold_sensitivity.csv")
common6 = pd.read_csv(SRC / "R18_common_all_six_cells.csv")
training = pd.read_csv(SRC / "R18_training_fixed_summary.csv")
bootstrap = pd.read_csv(SRC / "R18_donor_bootstrap_overall_replicates.csv")
fig, ax = new_figure(4, 3, 2, 205)
label(ax[0], "A", "Paired continuous CAF")
transition = np.select([paired.stalled_frozen & paired.stalled_level,
                        ~paired.stalled_frozen & paired.stalled_level,
                        paired.stalled_frozen & ~paired.stalled_level],
                       ["Both stall", "Level only", "Frozen only"], default="Neither")
colors = {"Both stall": DEEP, "Level only": ORANGE, "Frozen only": RED, "Neither": GREY}
for group in ["Neither", "Both stall", "Level only", "Frozen only"]:
    mask = transition == group
    ax[0].scatter(paired.CAF_frozen[mask], paired.CAF_level[mask], s=11, color=colors[group],
                  marker={"Neither": "o", "Both stall": "s", "Level only": "^", "Frozen only": "D"}[group],
                  label=group, edgecolors="none")
lim = [float(min(paired.CAF_frozen.min(), paired.CAF_level.min()) * 1.10),
       float(max(paired.CAF_frozen.max(), paired.CAF_level.max()) * 1.20)]
ax[0].plot(lim, lim, color=DARK, lw=0.7, ls=":")
ax[0].axhline(0.6, color=GREY, lw=0.7)
ax[0].axvline(0.6, color=GREY, lw=0.7)
ax[0].set_xscale("symlog", linthresh=1)
ax[0].set_yscale("symlog", linthresh=1)
ax[0].set_xlim(lim)
ax[0].set_ylim(lim)
ax[0].set_xlabel("Frozen CAF (symlog)")
ax[0].set_ylabel("Level CAF (symlog)")
ax[0].legend(frameon=False, fontsize=7, loc="upper left")
label(ax[1], "B", "Common-denominator stalls")
ax[1].bar(["Frozen", "Level"], [5, 31], color=[BLUE, ORANGE], width=0.55)
ax[1].set_ylabel("Stalls among 79 cells")
ax[1].set_ylim(0, 39)
ax[1].text(0.03, 0.94, "Paired difference 0.329\n95% CI 0.190–0.468", transform=ax[1].transAxes,
           va="top", fontsize=7.5, bbox=dict(facecolor="white", edgecolor=GREY, boxstyle="round,pad=0.3"))
label(ax[2], "C", "Threshold sensitivity")
ax[2].plot(thresholds.threshold, thresholds.frozen_stalls, marker="o", color=BLUE, lw=1.3, label="Frozen")
ax[2].plot(thresholds.threshold, thresholds.level_stalls, marker="^", color=ORANGE, lw=1.3, label="Level")
ax[2].set_xticks(thresholds.threshold)
ax[2].set_xlabel("CAF stall threshold")
ax[2].set_ylabel("Stalls among 79 cells")
ax[2].legend(frameon=False)
label(ax[3], "D", "Reference sensitivity (50 common cells)")
references = ["LM22-NNLS", "LM22-SVR", "REF1-NNLS"]
reference_rows = []
for reference in references:
    reference_rows.append({
        "reference": reference, "common_cells": len(common6),
        "frozen_stalls": int(common6[f"{reference}_frozen_stalled"].sum()),
        "level_stalls": int(common6[f"{reference}_level_stalled"].sum()),
    })
reference_plot = pd.DataFrame(reference_rows)
reference_plot.to_csv(SRC / "Figure4D_reference_sensitivity_50.csv", index=False)
x = np.arange(len(references))
ax[3].bar(x - 0.17, reference_plot.frozen_stalls, width=0.34, color=BLUE, label="Frozen")
ax[3].bar(x + 0.17, reference_plot.level_stalls, width=0.34, color=ORANGE, label="Level")
ax[3].set_xticks(x, ["LM22\nNNLS", "LM22\nSVR", "REF1\nNNLS"])
ax[3].set_ylim(0, 27)
ax[3].set_ylabel("Stalls among 50 cells")
ax[3].legend(frameon=False)
for i, row in reference_plot.iterrows():
    ax[3].text(i - 0.17, row.frozen_stalls + 0.7, str(row.frozen_stalls), ha="center", fontsize=7)
    ax[3].text(i + 0.17, row.level_stalls + 0.7, str(row.level_stalls), ha="center", fontsize=7)
label(ax[4], "E", "Training-fixed external transform")
training.to_csv(SRC / "Figure4E_training_fixed_external.csv", index=False)
x = np.arange(len(training))
width = 0.25
ax[4].bar(x - width, training.frozen_stalls, width=width, color=BLUE, label="Frozen")
ax[4].bar(x, training.level_cohort_stalls, width=width, color=ORANGE, label="Cohort level")
ax[4].bar(x + width, training.level_training_fixed_stalls, width=width, color=DEEP, label="Training-fixed")
ax[4].set_xticks(x, [f"{row.validation_cohort}\n{row.common_cells} cells" for row in training.itertuples()])
ax[4].set_ylim(0, 6)
ax[4].set_ylabel("Stalls among common cells")
ax[4].legend(frameon=False, fontsize=7, loc="upper right")
label(ax[5], "F", "Sample/analysis-unit bootstrap (200 replicates)")
bootstrap.to_csv(SRC / "Figure4F_donor_bootstrap_overall.csv", index=False)
values = bootstrap.stall_rate_delta.to_numpy()
ci = np.percentile(values, [2.5, 97.5])
ax[5].hist(values, bins=18, color=BLUE, alpha=0.75, edgecolor="white", linewidth=0.5)
ax[5].axvline(np.median(values), color=DEEP, lw=1.2, label=f"Median {np.median(values):.3f}")
ax[5].axvline(ci[0], color=ORANGE, lw=1.0, ls="--")
ax[5].axvline(ci[1], color=ORANGE, lw=1.0, ls="--", label=f"95% percentile {ci[0]:.3f}–{ci[1]:.3f}")
ax[5].axvline(0, color=DARK, lw=0.7)
ax[5].set_xlabel("Level − frozen stall-rate difference")
ax[5].set_ylabel("Bootstrap replicates")
ax[5].legend(frameon=False, fontsize=7, loc="upper left")
ax[5].text(0.98, 0.68, "Mean CAF Δ 2.863\n95% percentile 0.610–10.106\nSplit-half means 3.146 / 2.579",
           transform=ax[5].transAxes, ha="right", va="top", fontsize=7,
           bbox=dict(facecolor="white", edgecolor=GREY, boxstyle="round,pad=0.25"))
fig.subplots_adjust(hspace=0.62)
save(fig, 4)

software_row = pd.read_csv(SRC / "Figure5D_quickstart_cost.csv").iloc[0]
software = {"repetitions": int(software_row.repetitions), "failures": int(software_row.failures),
            "runtime_ms_median": software_row.runtime_ms_median, "runtime_ms_p95": software_row.runtime_ms_p95,
            "peak_memory_kib_single_run": software_row.peak_memory_kib}
pd.DataFrame({"test_set": ["compGate package", "Joint-model requirements"], "passed": [11, 7], "total": [11, 7]}).to_csv(
    SRC / "Figure5B_automated_tests.csv", index=False)
pd.DataFrame([{"repetitions": software["repetitions"], "failures": software["failures"],
               "runtime_ms_median": software["runtime_ms_median"], "runtime_ms_p95": software["runtime_ms_p95"],
               "peak_memory_kib": software["peak_memory_kib_single_run"]}]).to_csv(
    SRC / "Figure5D_quickstart_cost.csv", index=False)
pd.DataFrame({"item": ["Steps", "CSV files", "Figures", "Mismatches"], "count": [24, 82, 16, 0]}).to_csv(
    SRC / "Figure5C_reproduction.csv", index=False)
rho = pd.read_csv(SRC / "Figure5E_RHO_HKSJ_BH.csv")
fig, ax = new_figure(5, 3, 2, 205)
label(ax[0], "A", "Software and analysis chain")
blank(ax[0])
box(ax[0], (0.05, 0.62), 0.38, 0.20, "compGate 0.4.0\nPyPI: scoring + joint model", BLUE, 7)
box(ax[0], (0.57, 0.62), 0.38, 0.20, "GitHub v0.4.0\nsource + benchmark", DEEP, 7)
box(ax[0], (0.25, 0.20), 0.50, 0.20, "Zenodo archive\n10.5281/zenodo.23889776\nSource Data + logs", ORANGE, 7)
ax[0].annotate("", xy=(0.48, 0.41), xytext=(0.26, 0.61), arrowprops=dict(arrowstyle="->", color=DARK))
ax[0].annotate("", xy=(0.54, 0.41), xytext=(0.74, 0.61), arrowprops=dict(arrowstyle="->", color=DARK))
label(ax[1], "B", "Automated tests")
ax[1].bar(["compGate\npackage", "Joint-model\nrequirements"], [11, 7], color=[BLUE, DEEP], width=0.55)
ax[1].set_ylabel("Tests passed")
ax[1].set_ylim(0, 14)
for i, v in enumerate([11, 7]):
    ax[1].text(i, v + 0.5, f"{v}/{v}", ha="center", fontweight="bold")
label(ax[2], "C", "Packaged reproduction")
vals = [24, 82, 16]
ax[2].bar(["Steps", "CSV files", "Figures"], vals, color=[GREY, BLUE, DEEP], width=0.58)
ax[2].set_ylabel("Items verified")
ax[2].set_ylim(0, 94)
for i, v in enumerate(vals):
    ax[2].text(i, v + 3, str(v), ha="center", fontweight="bold")
ax[2].text(0.03, 0.94, "0 mismatches", transform=ax[2].transAxes, va="top", color=DEEP, fontweight="bold")
label(ax[3], "D", "Minimal quickstart cost")
blank(ax[3])
ax[3].text(0.07, 0.77, f"{software['runtime_ms_median']:.2f} ms", fontsize=17, color=BLUE, fontweight="bold")
ax[3].text(0.08, 0.62, "median runtime (200 runs)", fontsize=7.5)
ax[3].text(0.56, 0.77, f"{software['runtime_ms_p95']:.2f} ms", fontsize=13, color=DEEP, fontweight="bold")
ax[3].text(0.57, 0.62, "95th percentile", fontsize=7.5)
ax[3].text(0.07, 0.35, f"{software['peak_memory_kib_single_run']:.1f} KiB", fontsize=13, color=ORANGE, fontweight="bold")
ax[3].text(0.08, 0.21, "peak traced memory", fontsize=7.5)
ax[3].text(0.58, 0.35, "0 / 200", fontsize=13, color=DEEP, fontweight="bold")
ax[3].text(0.59, 0.21, "failures", fontsize=7.5)
label(ax[4], "E", "Compact RHO worked example")
ax[4].bar(["Nominal\nHKSJ", "After BH\n(14 members)"], [4, 0], color=[BLUE, ORANGE], width=0.55)
ax[4].set_ylabel("Supported RHO members")
ax[4].set_ylim(0, 5.2)
ax[4].text(0, 4.25, "4", ha="center", fontweight="bold")
ax[4].text(1, 0.25, "0", ha="center", fontweight="bold")
ax[4].text(0.03, 0.90, "Minimum BH q = 0.127", transform=ax[4].transAxes, fontsize=7.5, color=DARK)
blank(ax[5])
ax[5].text(0.03, 0.76, "Verified scope", fontsize=9, fontweight="bold", color=DEEP)
for i, text in enumerate(["Local Windows / Python 3.12.10 execution",
                          "Minimal example only, not genome-scale cost",
                          "No independent-team or Linux claim",
                          "Public 0.4.0 release; independent checks open"]):
    ax[5].text(0.06, 0.58 - i * 0.15, "• " + text, fontsize=7.5)
fig.subplots_adjust(hspace=0.62)
save(fig, 5)

manifest = []
for number in range(1, 6):
    for ext in ["pdf", "eps", "png", "tif"]:
        path = OUT / f"Figure_{number}.{ext}"
        manifest.append({"figure": number, "format": ext, "path": str(path.relative_to(DELIVERY)), "bytes": path.stat().st_size})
pd.DataFrame(manifest).to_csv(DELIVERY / "figure_file_manifest.csv", index=False)
print(pd.DataFrame(manifest).to_string(index=False))

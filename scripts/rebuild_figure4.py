from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description="Rebuild CMPB Figure 4 from packaged frozen source data")
parser.add_argument("--source", default=str(ROOT / "examples" / "figure4_source"))
parser.add_argument("--output", default=str(ROOT / "verification" / "figure4"))
args = parser.parse_args()
SRC = Path(args.source)
OUT = Path(args.output)
OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#0999d8"
ORANGE = "#f79108"
RED = "#e03f40"
DEEP = "#086898"
GREY = "#9aa5b1"
DARK = "#333333"
mpl.rcParams.update({
    "font.family": "Arial", "font.size": 7.5, "axes.titlesize": 8,
    "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.edgecolor": DARK, "axes.linewidth": 0.8,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.major.width": 0.8,
    "ytick.major.width": 0.8, "axes.grid": False, "savefig.facecolor": "white",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def new_figure():
    fig, axes = plt.subplots(3, 2, figsize=(190 / 25.4, 205 / 25.4))
    fig.subplots_adjust(left=0.07, right=0.98, top=0.94, bottom=0.10, wspace=0.34, hspace=0.52)
    for ax in np.ravel(axes):
        ax.tick_params(length=3)
    return fig, np.ravel(axes)


def label(ax, letter, title=None):
    ax.text(-0.16, 1.08, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")
    if title:
        ax.set_title(title, loc="left", pad=5, fontweight="bold")


paired = pd.read_csv(SRC / "Figure4A_paired_CAF_79.csv")
thresholds = pd.read_csv(SRC / "Figure4C_threshold_sensitivity.csv")
common6 = pd.read_csv(SRC / "R18_common_all_six_cells.csv")
training = pd.read_csv(SRC / "R18_training_fixed_summary.csv")
bootstrap = pd.read_csv(SRC / "R18_donor_bootstrap_overall_replicates.csv")
assert len(paired) == 79 and len(common6) == 50 and len(training) == 2 and len(bootstrap) == 200

fig, ax = new_figure()
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
        "reference": reference,
        "common_cells": len(common6),
        "frozen_stalls": int(common6[f"{reference}_frozen_stalled"].sum()),
        "level_stalls": int(common6[f"{reference}_level_stalled"].sum()),
    })
reference_plot = pd.DataFrame(reference_rows)
reference_plot.to_csv(OUT / "Figure4D_reference_sensitivity_50.csv", index=False)
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
training.to_csv(OUT / "Figure4E_training_fixed_external.csv", index=False)
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
bootstrap.to_csv(OUT / "Figure4F_donor_bootstrap_overall.csv", index=False)
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
fig.canvas.draw()
content_in = fig.get_tightbbox(fig.canvas.get_renderer()).width
pad = max(0.0, (190.0 / 25.4 - content_in) / 2)
fig.savefig(OUT / "Figure_4.pdf", bbox_inches="tight", pad_inches=pad)
fig.savefig(OUT / "Figure_4.png", dpi=200, bbox_inches="tight", pad_inches=pad)
plt.close(fig)
anchors = {
    "paired_rows": len(paired),
    "common_six_arm_rows": len(common6),
    "reference_stalls": reference_rows,
    "training_fixed_rows": len(training),
    "bootstrap_replicates": len(bootstrap),
    "bootstrap_stall_rate_median": float(np.median(values)),
    "bootstrap_stall_rate_percentile_95": ci.tolist(),
    "outputs": ["Figure_4.png", "Figure_4.pdf"],
}
(OUT / "figure4_anchors.json").write_text(json.dumps(anchors, indent=2), encoding="utf-8")
print(json.dumps(anchors, indent=2))

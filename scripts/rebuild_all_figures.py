# -*- coding: utf-8 -*-
"""Rebuild all five CMPB figures from the frozen packaged Source Data.

This is the single authoritative figure-generation script for compGate 0.4.0.
It was merged on 2026-11-06 from the editorial-revision build: Figures 1, 2, 3
and 5 now use the higher-density panels requested by the advisor, while Figure 4
is unchanged.  Output is deliberately limited to two formats:

    vector master = PDF (single page, no raster objects, Arial embedded)
    preview       = PNG

EPS and TIFF were discontinued as redundant (2026-11-06).  The script is fully
self-contained: it only reads ``source_data/*.csv`` and writes inside the
``--output`` directory.  It never writes to Source Data, never depends on a
workspace path, and never contacts the network.

Usage:
    python scripts/rebuild_all_figures.py \
        --source source_data --output verification/all_figures
"""
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Patch
from matplotlib.colors import LogNorm, TwoSlopeNorm
from matplotlib.lines import Line2D
import fitz

ROOT = Path(__file__).resolve().parents[1]
PARSER = argparse.ArgumentParser(description="Rebuild all five CMPB figures from packaged Source Data")
PARSER.add_argument("--source", default=str(ROOT / "source_data"))
PARSER.add_argument("--output", default=str(ROOT / "verification" / "all_figures"))
ARGS = PARSER.parse_args()
SRC = Path(ARGS.source)
DELIVERY = Path(ARGS.output)
OUT = DELIVERY / "figures"
OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#0999d8"
ORANGE = "#f79108"
RED = "#e03f40"
DEEP = "#086898"
GREY = "#9aa5b1"
DARK = "#333333"
LIGHT = "#eef2f7"
WHITE = "#ffffff"

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 7.5, "axes.titlesize": 8,
    "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.edgecolor": DARK, "axes.linewidth": 0.8,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.major.width": 0.8,
    "ytick.major.width": 0.8, "axes.grid": False, "savefig.facecolor": "white",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

LEDGER = {"revision": "editorial_revision_merged_20261106",
          "source_dir": "source_data", "output_dir": "verification/all_figures",
          "formats": ["pdf", "png"], "figures": {}}


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def snapshot():
    return {p.name: sha(p) for p in sorted(SRC.glob("*.csv"))}


SNAP_BEFORE = snapshot()
assert len(SNAP_BEFORE) == 50, f"expected 50 Source Data CSVs, got {len(SNAP_BEFORE)}"


def new_figure(nrows, ncols, height, wspace=0.34, hspace=0.52):
    fig, axes = plt.subplots(nrows, ncols, figsize=(190 / 25.4, height / 25.4))
    fig.subplots_adjust(left=0.07, right=0.98, top=0.94, bottom=0.10,
                        wspace=wspace, hspace=hspace)
    for ax in np.ravel(axes):
        ax.tick_params(length=3)
    return fig, np.ravel(axes)


def label(ax, letter, title=None):
    ax.text(-0.16, 1.08, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top")
    if title:
        ax.set_title(title, loc="left", pad=5, fontweight="bold")


def blank(ax, xlim=(0, 1), ylim=(0, 1)):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.axis("off")


def box(ax, xy, width, height, text, color=BLUE, size=7.0, face=WHITE):
    x, y = xy
    patch = FancyBboxPatch((x, y), width, height,
                           boxstyle="round,pad=0.0,rounding_size=0.015",
                           facecolor=face, edgecolor=color, linewidth=1.0)
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height / 2, text, ha="center", va="center",
            fontsize=size, color=DARK)


def arrow(ax, x0, y0, x1, y1, color=DARK, lw=1.0):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="->", color=color, lw=lw,
                                shrinkA=0, shrinkB=0))


def vector_colorbar(fig, ax, cmap, norm, label_text, ticks):
    """Colour scale drawn as Rectangle patches (matplotlib's default Colorbar
    emits a small raster gradient into the PDF)."""
    fig.canvas.draw()
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + 0.010, pos.y0 + pos.height * 0.10,
                        0.011, pos.height * 0.80])
    lo, hi = float(norm.vmin), float(norm.vmax)
    is_log = isinstance(norm, LogNorm)
    values = (np.logspace(np.log10(lo), np.log10(hi), 120) if is_log
              else np.linspace(lo, hi, 120))
    cmap_obj = mpl.colormaps.get_cmap(cmap)
    n = len(values)
    for i, v in enumerate(values):
        cax.add_patch(Rectangle((0, i / n), 1, 1 / n + 0.002,
                                facecolor=cmap_obj(norm(v)), edgecolor="none"))
    cax.set_xlim(0, 1)
    cax.set_ylim(0, 1)
    cax.axis("off")
    for t in ticks:
        frac = ((np.log10(t) - np.log10(lo)) / (np.log10(hi) - np.log10(lo))
                if is_log else (t - lo) / (hi - lo))
        cax.plot([0.85, 1.25], [frac, frac], color=DARK, lw=0.6,
                 transform=cax.transAxes, clip_on=False)
        cax.text(1.35, frac, f"{t:g}", fontsize=6, va="center", ha="left",
                 transform=cax.transAxes, clip_on=False)
    cax.text(4.6, 0.5, label_text, fontsize=6, ha="center", va="center",
             rotation=90, transform=cax.transAxes, clip_on=False)


def vector_heatmap(ax, matrix, cmap, norm, xlabels=None, ylabels=None,
                   xgroup=None, show_y_every=None, cbar_label=None,
                   cb_ticks=None, grid_lines=True):
    m = np.asarray(matrix, dtype=float)
    nrow, ncol = m.shape
    cmap_obj = mpl.colormaps.get_cmap(cmap)
    for i in range(nrow):
        for j in range(ncol):
            v = m[i, j]
            color = cmap_obj(norm(v)) if np.isfinite(v) else "#dddddd"
            ax.add_patch(Rectangle((j, i), 1, 1, facecolor=color, edgecolor="none"))
    if grid_lines and ncol <= 40:
        for j in range(ncol + 1):
            ax.axvline(j, color="white", lw=0.25)
    ax.set_xlim(0, ncol)
    ax.set_ylim(nrow, 0)
    if xlabels is not None:
        ax.set_xticks(np.arange(ncol) + 0.5, xlabels, rotation=45, ha="right", fontsize=6)
    elif xgroup is not None:
        centres, names = xgroup
        ax.set_xticks(centres, names, rotation=30, ha="right", fontsize=6.5)
    else:
        ax.set_xticks([])
    if ylabels is not None:
        if show_y_every:
            idx = list(range(0, nrow, show_y_every))
            ax.set_yticks([i + 0.5 for i in idx], [ylabels[i] for i in idx], fontsize=6)
        else:
            ax.set_yticks(np.arange(nrow) + 0.5, ylabels, fontsize=6)
    else:
        ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    if cbar_label:
        vector_colorbar(ax.figure, ax, cmap, norm, cbar_label, cb_ticks or [])


def save(fig, number, panel_letters, forbidden, creation_date=None):
    stem = OUT / f"Figure_{number}"
    fig.canvas.draw()
    content_in = fig.get_tightbbox(fig.canvas.get_renderer()).width
    pad = max(0.0, (190.0 / 25.4 - content_in) / 2)
    # A fixed date makes repeated builds byte-reproducible.  Figure 4 is the
    # unchanged frozen master, so it keeps the 2026-11-03 date verbatim and the
    # rebuilt PDF is byte-identical to the delivered vector master.
    meta = {"CreationDate": creation_date} if creation_date else {"CreationDate": None}
    pdf = stem.with_suffix(".pdf")
    fig.savefig(pdf, bbox_inches="tight", pad_inches=pad, metadata=meta)
    fig.savefig(stem.with_suffix(".png"), dpi=200, bbox_inches="tight", pad_inches=pad)
    plt.close(fig)

    d = fitz.open(str(pdf))
    assert d.page_count == 1, f"Figure {number}: {d.page_count} pages"
    text = d[0].get_text()
    imgs = d[0].get_images(full=True)
    d.close()
    assert len(imgs) == 0, f"Figure {number}: raster image embedded ({len(imgs)})"
    for letter in panel_letters:
        assert letter in text, f"Figure {number}: panel letter {letter} missing"
    for bad in forbidden:
        assert bad not in text, f"Figure {number}: forbidden text present: {bad}"
    LEDGER["figures"][str(number)] = {
        "pdf": pdf.name, "pdf_bytes": pdf.stat().st_size, "pdf_sha256": sha(pdf),
        "png": stem.with_suffix(".png").name,
        "png_bytes": stem.with_suffix(".png").stat().st_size,
        "png_sha256": sha(stem.with_suffix(".png")),
        "vector_pdf_no_raster_images": True, "panel_letters": list(panel_letters),
    }
    print(f"Figure {number}: PDF {pdf.stat().st_size} B, PNG "
          f"{stem.with_suffix('.png').stat().st_size} B, panels {panel_letters}")


# =====================================================================
# Figure 1
# =====================================================================
refconv = pd.read_csv(SRC / "R18_reference_convention_long.csv")
counts = pd.read_csv(SRC / "Figure1C_eligibility_GSE65682.csv")
flow = pd.read_csv(SRC / "Figure1D_hallmark_GSE65682.csv")

fig, ax = new_figure(2, 2, 150)
label(ax[0], "A", "Attribution problem")
blank(ax[0])
box(ax[0], (0.02, 0.62), 0.30, 0.20, "Bulk group\neffect", DEEP, 7)
box(ax[0], (0.66, 0.78), 0.32, 0.16, "Composition shift", ORANGE, 7)
box(ax[0], (0.66, 0.48), 0.32, 0.16, "Within-cell change", BLUE, 7)
arrow(ax[0], 0.33, 0.74, 0.65, 0.86)
arrow(ax[0], 0.33, 0.70, 0.65, 0.56)
ax[0].text(0.35, 0.94, "two components", fontsize=6.5, color=DARK)
med = float(flow.median_CAF.median())
ax[0].text(0.02, 0.35, r"$\mathrm{CAF}_g$ = attenuated share of $|d_{\mathrm{unadj},g}|$",
           fontsize=6.8, color=DARK)
ax[0].add_patch(Rectangle((0.02, 0.16), 0.96 * med, 0.11, facecolor=ORANGE, edgecolor="none"))
ax[0].add_patch(Rectangle((0.02 + 0.96 * med, 0.16), 0.96 * (1 - med), 0.11,
                          facecolor=GREY, edgecolor="none"))
ax[0].text(0.02, 0.06, f"Cross-set median CAF {med:.3f} (50 Hallmark sets)", fontsize=6.5, color=DARK)
ax[0].text(0.02 + 0.96 * med / 2, 0.215, f"{med:.3f}", ha="center", va="center",
           fontsize=7, fontweight="bold", color=WHITE)
ax[0].text(0.02 + 0.96 * med + 0.96 * (1 - med) / 2, 0.215, f"{1 - med:.3f}",
           ha="center", va="center", fontsize=7, fontweight="bold", color=WHITE)

label(ax[1], "B", "Common-scale estimands")
a = np.abs(refconv.d_unadj.to_numpy())
b = np.abs(refconv.d_joint.to_numpy())
ok = refconv.interpretable.to_numpy().astype(bool)
ax[1].scatter(a[ok], b[ok], s=4.5, color=BLUE, alpha=0.55, edgecolors="none",
              label=f"Interpretable (n={ok.sum()})")
ax[1].scatter(a[~ok], b[~ok], s=4.5, color=GREY, alpha=0.55, edgecolors="none",
              label=f"Not interpretable (n={(~ok).sum()})")
ax[1].set_xscale("symlog", linthresh=0.5)
ax[1].set_yscale("symlog", linthresh=0.5)
ax[1].set_xlim(0, 6)
ax[1].set_ylim(0, 8)
xx = np.linspace(0, 6, 100)
ax[1].plot(xx, xx, color=DARK, lw=0.8, ls="--",
           label=r"$|d_{\mathrm{joint},g}|=|d_{\mathrm{unadj},g}|$")
ax[1].plot(xx, 0.4 * xx, color=ORANGE, lw=0.9, ls=":",
           label=r"$\mathrm{CAF}_g=0.6$ boundary")
ax[1].set_xlabel(r"$|d_{\mathrm{unadj},g}|$ (symlog)")
ax[1].set_ylabel(r"$|d_{\mathrm{joint},g}|$ (symlog)")
ax[1].legend(frameon=False, fontsize=6, loc="upper left")
ax[1].text(0.97, 0.03,
           r"$d_{\mathrm{unadj},g}=\beta_{U,g}/s_{\mathrm{ref},g}$" + "\n" +
           r"$d_{\mathrm{joint},g}=\beta_{A,g}/s_{\mathrm{ref},g}$" + "\n" +
           r"$\mathrm{CAF}_g=1-|d_{\mathrm{joint},g}|/|d_{\mathrm{unadj},g}|$",
           transform=ax[1].transAxes, ha="right", va="bottom", fontsize=6.5,
           color=DEEP, bbox=dict(facecolor="white", edgecolor=GREY,
                                 boxstyle="round,pad=0.25", alpha=0.9))

label(ax[2], "C", "GSE65682 eligibility states")
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
save(fig, 1, "ABCD", ["Estimate both components", "report the estimand"])
LEDGER["figures"]["1"]["provenance"] = {
    "A_partition": {"file": "Figure1D_hallmark_GSE65682.csv", "column": "median_CAF",
                    "statistic": "median of 50 set medians", "value": med},
    "B_scatter": {"file": "R18_reference_convention_long.csv",
                  "columns": ["d_unadj", "d_joint", "interpretable"], "rows": len(refconv)},
}

# =====================================================================
# Figure 2
# =====================================================================
gate_metrics = pd.read_csv(SRC / "Figure2D_gate_metrics.csv")
metrics = {f"within_{row.metric.lower()}": {"mean": row.mean, "ci_low": row.ci_low,
                                            "ci_high": row.ci_high}
           for row in gate_metrics.itertuples(index=False)}

fig, ax = new_figure(2, 2, 150)
label(ax[0], "A", "Known-truth design matrix")
blank(ax[0], (0, 2), (0, 4))
design = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 1.0], [1.0, 1.0]])
cmap_design = {0.0: LIGHT, 0.5: ORANGE, 1.0: BLUE}
for i in range(4):
    for j in range(2):
        v = design[i, j]
        ax[0].add_patch(Rectangle((j + 0.08, i + 0.10), 0.84, 0.80,
                                  facecolor=cmap_design[v], edgecolor=WHITE, linewidth=1.2))
        ax[0].text(j + 0.5, i + 0.5, {0.0: "absent", 0.5: "fixed", 1.0: "present"}[v],
                   ha="center", va="center", fontsize=7,
                   color=DARK if v != 1.0 else WHITE)
for i, name in enumerate(["Null", "Composition-only", "Within-cell", "Mixed"]):
    ax[0].text(-0.06, i + 0.5, name, ha="right", va="center", fontsize=7, clip_on=False)
for j, name in enumerate(["C  composition", "W  within-cell"]):
    ax[0].text(j + 0.5, 4.16, name, ha="center", va="top", fontsize=7, clip_on=False)
ax[0].invert_yaxis()

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
ax[2].bar(["Group\nonly", "Estimated\nfractions", "Oracle\nfractions"], calls,
          color=[ORANGE, BLUE, DEEP], width=0.62)
ax[2].set_ylabel("Significant unmodified genes")
ax[2].set_ylim(0, 355)
for i, v in enumerate(calls):
    ax[2].text(i, v + 10, str(v), ha="center", fontweight="bold")
ax[2].text(0.02, -0.24, "Observed count; not a general FPR estimate",
           transform=ax[2].transAxes, fontsize=7)

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
save(fig, 2, "ABCD", ["Multiple seeds; scenario-level summaries"])

# =====================================================================
# Figure 3
# =====================================================================
fraction_summary = pd.read_csv(SRC / "R18_fraction_error_summary.csv")
error_methods = ["low_error", "high_error", "systematic_bias", "reference_estimated"]
error_data = fraction_summary[fraction_summary.method.isin(error_methods)][
    ["scenario", "method", "fraction_RMSE_median", "CAF_bias_median",
     "FPR_median", "coverage_95_median"]].copy()
frozen_3d = pd.read_csv(SRC / "Figure3D_fraction_error.csv")
recomputed = error_data.sort_values(["scenario", "method"]).reset_index(drop=True)
frozen_sorted = frozen_3d.sort_values(["scenario", "method"]).reset_index(drop=True)
assert len(recomputed) == len(frozen_sorted) == 144
assert np.allclose(recomputed.CAF_bias_median, frozen_sorted.CAF_bias_median)
assert np.allclose(recomputed.fraction_RMSE_median, frozen_sorted.fraction_RMSE_median)

per_run = pd.read_csv(SRC / "R18_fraction_error_per_run.csv")
assert len(per_run) == 1260 and (per_run.model_state == "OK").all()
method_order = ["group_only", "expression_residual", "oracle", "low_error",
                "high_error", "systematic_bias", "reference_estimated"]
scenario_order = sorted(per_run.scenario.unique())
seed_order = sorted(per_run.replicate.unique())
cols = [(m, s) for m in method_order for s in seed_order]
lookup = {(r.scenario, r.method, r.replicate): r.seconds for r in per_run.itertuples(index=False)}
wall = np.array([[lookup[(sc, m, s)] for (m, s) in cols] for sc in scenario_order])
assert wall.shape == (36, 35) and np.isfinite(wall).all()

fig, ax = new_figure(3, 2, 205)
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
vals = pd.read_csv(SRC / "Figure3B_reference_failure.csv")[["FPR_median", "FDR_median"]].to_numpy()
x = np.arange(3)
ax[1].bar(x - 0.17, vals[:, 0], width=0.34, color=BLUE, label="FPR")
ax[1].bar(x + 0.17, vals[:, 1], width=0.34, color=ORANGE, label="FDR")
ax[1].set_xticks(x, ["Matched", "Missing\nneutrophil", "Reference\ndrift"])
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
error_labels = ["Low random", "High random", "Systematic", "Reference\nestimated"]
error_colors = [DEEP, BLUE, RED, ORANGE]
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
    ax[3].scatter(np.full(len(values), position) + jitter, values, s=7, color=color,
                  alpha=0.72, edgecolors="none")
ax[3].axhline(0, color=DARK, lw=0.7)
ax[3].set_xticks(positions, error_labels)
ax[3].set_ylabel("Scenario-median CAF bias")
ax[3].set_xlabel("Fraction-error path; points are scenario summaries")

label(ax[4], "E", "CAF bias across 37 reference-estimated summaries")
bias = pd.read_csv(SRC / "Figure3E_CAF_bias.csv")
ax[4].bar(np.arange(len(bias)), bias.CAF_bias_median,
          color=[ORANGE if v > 0 else BLUE for v in bias.CAF_bias_median], width=0.82)
ax[4].axhline(0, color=DARK, lw=0.8)
ax[4].set_xlabel("Sorted summaries (36 scenarios + central series)")
ax[4].set_ylabel("Median CAF bias")
ax[4].set_xticks([])

label(ax[5], "F", "Completed factorial: 36 x 5 x 7 = 1,260 fits (all OK)")
centres = [2.5 + 5 * i for i in range(7)]
short_names = ["Group\nonly", "Expr.\nresidual", "Oracle", "Low\nerror", "High\nerror",
               "Systematic", "Ref.\nestimated"]
vector_heatmap(ax[5], wall, "viridis",
               LogNorm(vmin=float(wall.min()), vmax=float(wall.max())),
               xgroup=(centres, short_names), cbar_label="Wall time (s)",
               cb_ticks=[0.002, 0.01, 0.05, 0.2], grid_lines=False)
for j in [5, 10, 15, 20, 25, 30]:
    ax[5].axvline(j, color="white", lw=0.8)
ax[5].set_ylabel("36 scenarios (sorted)")
fig.subplots_adjust(hspace=0.68)
save(fig, 3, "ABCDEF", ["Zero model failures", "CIBERSORTx and csSAM remain separate",
                        "Completed factorial scope", "All seven paths completed"])

# =====================================================================
# Figure 4 - design unchanged (advisor did not flag this figure)
# =====================================================================
paired = pd.read_csv(SRC / "Figure4A_paired_CAF_79.csv")
thresholds = pd.read_csv(SRC / "Figure4C_threshold_sensitivity.csv")
common6 = pd.read_csv(SRC / "R18_common_all_six_cells.csv")
training = pd.read_csv(SRC / "R18_training_fixed_summary.csv")
bootstrap = pd.read_csv(SRC / "R18_donor_bootstrap_overall_replicates.csv")
fig, ax = new_figure(3, 2, 205)
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
values = bootstrap.stall_rate_delta.to_numpy()
ci = np.percentile(values, [2.5, 97.5])
ax[5].hist(values, bins=18, color=BLUE, alpha=0.75, edgecolor="white", linewidth=0.5)
ax[5].axvline(values.mean(), color=DEEP, lw=1.2, label=f"Mean {values.mean():.3f}")
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
save(fig, 4, "ABCDEF", [], creation_date="D:20261103131258+08'00'")

# =====================================================================
# Figure 5 - 4 panels A-D (editorial revision)
# =====================================================================
software_row = pd.read_csv(SRC / "Figure5D_quickstart_cost.csv").iloc[0]
tests = pd.read_csv(SRC / "Figure5B_automated_tests.csv")
repro = pd.read_csv(SRC / "Figure5C_reproduction.csv")
rho = pd.read_csv(SRC / "Figure5E_RHO_HKSJ_BH.csv").sort_values("HKSJ_BH_q").reset_index(drop=True)
n_tests = int(tests["passed"].sum())
assert n_tests == 18 and int(tests["total"].sum()) == 18
assert repro.set_index("item").loc["Mismatches", "count"] == 0
assert int(software_row.repetitions) == 200 and int(software_row.failures) == 0
assert int(rho.nominal_HKSJ_supported.sum()) == 4 and int(rho.BH_FDR_supported.sum()) == 0
assert np.isclose(rho.HKSJ_BH_q.min(), 0.127218)

cb = frozen_3d.pivot(index="scenario", columns="method", values="CAF_bias_median")
cb = cb.reindex(index=sorted(cb.index), columns=error_methods)
assert cb.shape == (36, 4) and np.isfinite(cb.to_numpy()).all()

fig = plt.figure(figsize=(190 / 25.4, 200 / 25.4))
gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 2.2],
                      left=0.07, right=0.98, top=0.94, bottom=0.09,
                      wspace=0.28, hspace=0.62)
A = fig.add_subplot(gs[0, :])
B = fig.add_subplot(gs[1, :])
C = fig.add_subplot(gs[2, 0])
D = fig.add_subplot(gs[2, 1])

label(A, "A", "Compact compGate 0.4.0 workflow")
blank(A)
stage_x = [0.005, 0.205, 0.405, 0.605, 0.805]
stage_w = 0.17
stages = [("Expression matrix\n+ cell fractions", DEEP),
          ("Scoring\nfrozen / level", BLUE),
          ("Joint model\nFWL diagnostics", BLUE),
          ("Attribution\neligibility / CAF / stall", ORANGE),
          ("Reports\nSource Data + figures", DEEP)]
for x, (text, color) in zip(stage_x, stages):
    box(A, (x, 0.30), stage_w, 0.44, text, color, 6.6)
for i in range(4):
    arrow(A, stage_x[i] + stage_w + 0.005, 0.52, stage_x[i + 1] - 0.005, 0.52)
A.text(0.005, 0.88, "compGate 0.4.0  |  PyPI / GitHub v0.4.0  |  Zenodo 10.5281/zenodo.23889776",
       fontsize=6.8, color=DARK)
A.text(0.005, 0.08, "single scoring-to-report chain; every stage emits a frozen Source Data record",
       fontsize=6.5, color=DARK)

label(B, "B", "Verification chain and quickstart cost over 200 runs")
blank(B)
chain = [("18 / 18", "automated tests", BLUE),
         ("24", "reproduction steps", DEEP),
         ("82 + 16", "CSV files + figures", BLUE),
         ("0", "mismatches", DEEP)]
for i, (num, text, color) in enumerate(chain):
    x = 0.005 + i * 0.252
    box(B, (x, 0.52), 0.22, 0.42, num + "\n" + text, color, 6.8)
    if i < 3:
        arrow(B, x + 0.225, 0.73, x + 0.247, 0.73)
med_ms = float(software_row.runtime_ms_median)
p95_ms = float(software_row.runtime_ms_p95)
mem = float(software_row.peak_memory_kib)
B.plot([0.005, p95_ms], [0.22, 0.22], color=GREY, lw=3.2, solid_capstyle="round")
B.plot([0.005, med_ms], [0.22, 0.22], color=BLUE, lw=3.2, solid_capstyle="round")
B.plot([med_ms, med_ms], [0.17, 0.27], color=DARK, lw=1.2)
B.plot([p95_ms, p95_ms], [0.17, 0.27], color=ORANGE, lw=1.2)
B.text(med_ms, 0.30, f"median {med_ms:.2f} ms", ha="center", fontsize=6.5, color=DEEP)
B.text(p95_ms, 0.30, f"p95 {p95_ms:.2f} ms", ha="center", fontsize=6.5, color=ORANGE)
B.set_xlim(0, 1.65)
B.set_ylim(0, 1)
B.text(0.005, 0.06, f"0 / {int(software_row.repetitions)} failures   |   {mem:.1f} KiB peak traced memory   |   runtime (ms)",
       fontsize=6.5, color=DARK)

label(C, "C", f"BH-adjusted q values (14 RHO members; min q = {rho.HKSJ_BH_q.min():.3f})")
y = np.arange(len(rho))
qcol = [ORANGE if s else DEEP for s in rho.nominal_HKSJ_supported]
C.barh(y, rho.HKSJ_BH_q, color=qcol, height=0.62)
C.axvline(0.05, color=RED, lw=1.0, ls="--")
C.set_yticks(y, rho.gene, fontsize=6.4)
C.invert_yaxis()
C.set_xlabel("HKSJ BH-adjusted q")
C.set_xlim(0, 0.9)
C.legend(handles=[Patch(facecolor=ORANGE, label="Nominal HKSJ support (n=4)"),
                  Patch(facecolor=DEEP, label="No nominal support (n=10)"),
                  Line2D([0], [0], color=RED, lw=1.0, ls="--", label="q = 0.05")],
         frameon=False, fontsize=6, loc="upper right")

label(D, "D", "CAF bias across the 36-scenario factorial")
norm = TwoSlopeNorm(vmin=float(cb.to_numpy().min()), vcenter=0.0,
                    vmax=float(cb.to_numpy().max()))
n_neg = int((cb.to_numpy() < 0).sum())
vector_heatmap(D, cb.to_numpy(), "RdBu_r", norm,
               xlabels=["Low\nrandom", "High\nrandom", "Systematic", "Ref.\nestimated"],
               cbar_label="CAF bias", cb_ticks=[-1.5, -0.5, 0.0], grid_lines=True)
D.set_ylabel("36 scenarios (sorted)")
D.text(0.99, -0.34, f"{cb.size} scenario x path summaries; {n_neg} negative",
       transform=D.transAxes, ha="right", fontsize=6.5, color=DEEP)
save(fig, 5, "ABCD", ["Local Windows / Python 3.12.10 execution",
                      "Minimal example only", "No independent-team or Linux claim",
                      "Verified scope", "Automated tests"])

# =====================================================================
# Manifest + invariants
# =====================================================================
assert snapshot() == SNAP_BEFORE, "Source Data was modified by the figure build"
rows = []
for number in range(1, 6):
    for ext in ("pdf", "png"):
        p = OUT / f"Figure_{number}.{ext}"
        rows.append({"figure": number, "format": ext,
                     "path": str(p.relative_to(DELIVERY)).replace("\\", "/"),
                     "bytes": p.stat().st_size, "sha256": sha(p)})
pd.DataFrame(rows).to_csv(DELIVERY / "figure_file_manifest.csv", index=False)
assert len(rows) == 10
LEDGER["source_data_unchanged"] = True
LEDGER["figure_file_manifest"] = rows
(DELIVERY / "figure_build_ledger.json").write_text(
    json.dumps(LEDGER, ensure_ascii=False, indent=2), encoding="utf-8")
print(pd.DataFrame(rows)[["figure", "format", "bytes"]].to_string(index=False))
print("Source Data unchanged: 50/50 CSVs")
print("Output formats: PDF (vector master) + PNG (preview); EPS/TIFF discontinued")

# -*- coding: utf-8 -*-
"""RE-CAPT six-gate audit and the CAR statistic (reference implementation).

The numeric conventions replicate the frozen analyses of the R10.4 package:
Cohen's d (pooled SD), Welch two-sided P, residualization by OLS (lstsq), CAF
= 1 - |d_adj| / |d_unadj|. Gate T (temporal kinetics) was added in R10.6
(2026-08-18) and implements the rise-and-return characterization used for the
MAI axis in manuscript section 6.2.
"""
import numpy as np
import pandas as pd
from scipy import stats


def cohen_d(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) /
                 (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / sp if sp > 0 else np.nan


def recap_gate(matrix, meta, case="Sepsis", ctrl="Healthy", fractions=None,
               marker_scores=None):
    """Run gates 2-3 (residualization) for every gene.

    matrix : DataFrame genes x samples
    meta   : DataFrame samples x {group}
    fractions : optional DataFrame samples x cell types (gate 3)
    marker_scores : optional DataFrame samples x signatures (gate 2)

    Returns DataFrame [gene, d_unadj, d_marker_adj, CAF_marker, d_frac_adj, CAF_frac].
    """
    y = (meta["group"] == case).astype(float).reindex(matrix.columns)
    rows = []
    for g in matrix.index:
        v = matrix.loc[g].values
        du = cohen_d(v[y.values == 1], v[y.values == 0])
        dm = ca_m = df_ = ca_f = np.nan
        if marker_scores is not None:
            X = marker_scores.reindex(matrix.columns).values
            r = v - X @ np.linalg.lstsq(X, v, rcond=None)[0]
            dm = cohen_d(r[y.values == 1], r[y.values == 0])
            ca_m = 1 - abs(dm) / abs(du) if abs(du) > 1e-9 else np.nan
        if fractions is not None:
            F = fractions.reindex(matrix.columns).values
            r = v - F @ np.linalg.lstsq(F, v, rcond=None)[0]
            df_ = cohen_d(r[y.values == 1], r[y.values == 0])
            ca_f = 1 - abs(df_) / abs(du) if abs(du) > 1e-9 else np.nan
        rows.append(dict(gene=g, d_unadj=du, d_marker_adj=dm, CAF_marker=ca_m,
                         d_frac_adj=df_, CAF_frac=ca_f))
    return pd.DataFrame(rows)


def gate_T(times, per_subject_values, return_fraction=0.5):
    """Gate T (temporal kinetics): characterize rise-and-return dynamics.

    times : 1-D array of measurement times (h), ascending.
    per_subject_values : list of 1-D per-subject arrays aligned to times.
    return_fraction : fraction of the peak deviation below which the final
        timepoint counts as a "return" (default 0.5).

    Returns dict(peak_h, peak_delta, return_delta, temporal_support); the
    caller states the biological interpretation (composition drift versus
    state-level or within-cell dynamics).
    """
    t = np.asarray(times, float)
    base = np.nanmean([v[0] for v in per_subject_values])
    means = np.array([np.nanmean([v[i] for v in per_subject_values])
                      for i in range(len(t))])
    ipeak = int(np.argmax(np.abs(means - base)))
    peak_delta = means[ipeak] - base
    return_delta = means[-1] - base
    if abs(peak_delta) < 1e-12:
        support = "flat"
    elif abs(return_delta) < return_fraction * abs(peak_delta):
        support = "rise-and-return"
    else:
        support = "sustained-shift"
    return dict(peak_h=t[ipeak], peak_delta=peak_delta,
                return_delta=return_delta, temporal_support=support)


def gate_report(per_gene, thresholds=None):
    """Six-gate checklist report from a per-gene recap_gate output.

    thresholds: dict with keys caf (default 0.5).
    """
    thresholds = thresholds or {}
    caf = thresholds.get("caf", 0.5)
    n = len(per_gene)
    return {
        "G1_replication": "pass input = cross-cohort direction-consistent gene list",
        "G2_marker_residualization": f"{100 * (per_gene.CAF_marker >= caf).mean():.1f}% of genes with CAF_marker >= {caf} ({n} genes)",
        "G3_reference_attribution": f"{100 * (per_gene.CAF_frac >= caf).mean():.1f}% of genes with CAF_frac >= {caf} ({n} genes)",
        "G4_powered_within_cell": "supply donor-level within-cell tests separately (MDE must be stated)",
        "G5_specificity_margin": "supply a non-septic inflammatory contrast separately",
        "G6_temporal_kinetics": "supply a time-course trajectory with stated kinetics (gate_T) separately",
    }


ATTRIBUTION_FIELDS = ("contrast", "reference", "convention", "estimator")


def attribution_report(claim, per_gene=None, thresholds=None):
    """One-click attribution record: the four parameters every
    composition-attribution claim must state, plus the CAF/stall summary.

    The quadruple follows the four reporting practices of the CompositionWall
    paper (contrast, reference, scoring convention, estimator form); a CAF is
    only interpretable inside all four.

    claim : dict with keys
        contrast    : str, e.g. "sepsis vs healthy" or "ATB vs LTBI"
        reference   : str, e.g. "LM22-NNLS" / "REF1-NNLS" / "LM22-SVR"
        convention  : "frozen" (within-sample-centred mean-Z) or
                      "level"  (per-gene standardization across samples)
        estimator   : str, e.g. "unsigned equal-weight mean-Z" or
                      "signed weighted score as deployed"
        optional    : platform, cohort, caf_floor (|d_unadj| floor, default 0.1)
    per_gene     : optional DataFrame from recap_gate
    thresholds   : dict with keys caf (stall cut-off, default 0.6) and
                   floor (|d_unadj| interpretability floor, default 0.1)

    Returns an ordered dict: the quadruple verbatim, completeness of the four
    declared parameters, the stall rule actually applied, and (when per_gene is
    given) the resulting counts. No value is invented: an absent parameter is
    reported as missing, never imputed.
    """
    thresholds = thresholds or {}
    caf_cut = float(thresholds.get("caf", 0.6))
    floor = float(claim.get("caf_floor", thresholds.get("floor", 0.1)))
    conv = claim.get("convention")
    if conv is not None and str(conv).lower() not in ("frozen", "level"):
        raise ValueError(f"convention must be 'frozen' or 'level', got {conv!r}")

    quad = {k: (claim.get(k) or None) for k in ATTRIBUTION_FIELDS}
    missing = [k for k, v in quad.items() if not v]
    out = {
        "attribution_quadruple": quad,
        "platform": claim.get("platform"),
        "cohort": claim.get("cohort"),
        "scoring_convention": conv,
        "complete": len(missing) == 0,
        "missing_parameters": missing,
        "stall_rule": {
            "caf_cut_off": caf_cut,
            "interpretability_floor_abs_d_unadj": floor,
            "note": ("the cut-off is a reporting convention, not a calibrated "
                     "threshold; it is exposed as a parameter"),
        },
    }
    if per_gene is not None:
        g = per_gene
        col = "CAF_frac" if "CAF_frac" in g.columns else (
            "CAF_marker" if "CAF_marker" in g.columns else None)
        if col is None:
            raise ValueError("per_gene must carry CAF_frac or CAF_marker")
        du = g["d_unadj"].abs() if "d_unadj" in g.columns else None
        interpretable = (du >= floor) if du is not None else g[col].notna()
        caf = g[col].where(interpretable)
        n_int = int(caf.notna().sum())
        n_stall = int((caf >= caf_cut).sum())
        out["result"] = {
            "n_genes": int(len(g)),
            "n_caf_interpretable": n_int,
            "n_stalled": n_stall,
            "stall_fraction": (n_stall / n_int) if n_int else float("nan"),
            "mean_CAF": float(caf.mean()) if n_int else float("nan"),
            "caf_column_used": col,
        }
    return out


def CAR(replicating_genes, per_gene, caf=0.5, col="CAF_frac"):
    """Composition-Adjusted Replication statistic.

    replicating_genes : iterable of gene symbols that replicate across cohorts
    per_gene          : DataFrame from recap_gate (gene, CAF_frac, ...)
    Returns (CAR, n_replicating, n_surviving).
    """
    sub = per_gene[per_gene.gene.isin(set(replicating_genes))]
    n_rep = len(sub)
    surv = sub[sub[col] < caf] if n_rep else sub
    n_surv = len(surv)
    return (n_surv / n_rep if n_rep else np.nan), n_rep, n_surv

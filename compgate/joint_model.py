from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg


def design_diagnostics(frac, y, fraction_names=None):
    y = np.asarray(y, dtype=float)
    values = np.asarray(frac, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    if fraction_names is None:
        fraction_names = [f"fraction_{i}" for i in range(values.shape[1])]
    finite_rows = np.isfinite(values).all(axis=1) & np.isfinite(y)
    active = np.nanvar(values, axis=0) > 1e-12
    zero_names = [name for name, keep in zip(fraction_names, active) if not keep]
    active_names = [name for name, keep in zip(fraction_names, active) if keep]
    reference = active_names[-1] if active_names else None
    keep = np.flatnonzero(active)
    reduced = keep[:-1] if len(keep) else keep
    full = np.column_stack([np.ones(len(y)), y, values[:, keep]])
    reduced_design = np.column_stack([np.ones(len(y)), y, values[:, reduced]])
    candidate_names = ["_intercept", "_group"] + active_names[:-1]
    raw_rank = int(np.linalg.matrix_rank(full))
    fit_rank = int(np.linalg.matrix_rank(reduced_design))
    selected = np.arange(reduced_design.shape[1])
    if fit_rank < reduced_design.shape[1]:
        _, r, piv = linalg.qr(reduced_design, mode="economic", pivoting=True)
        tolerance = np.finfo(float).eps * max(reduced_design.shape) * abs(r[0, 0])
        fit_rank = int(np.sum(np.abs(np.diag(r)) > tolerance))
        selected = np.sort(piv[:fit_rank])
        reduced_design = reduced_design[:, selected]
    additional_dropped = [name for i, name in enumerate(candidate_names) if i not in set(selected)]
    group_position = int(np.flatnonzero(selected == 1)[0]) if np.any(selected == 1) else -1
    fit_rank = int(np.linalg.matrix_rank(reduced_design))
    return {
        "values": values, "finite_rows": finite_rows, "active": active,
        "active_names": active_names, "zero_names": zero_names,
        "reference_dropped": reference, "full": full, "reduced": reduced_design,
        "raw_rank": raw_rank, "fit_rank": fit_rank,
        "raw_dof": int(finite_rows.sum() - raw_rank),
        "fit_dof": int(finite_rows.sum() - fit_rank),
        "n_columns": int(reduced_design.shape[1]),
        "constraint_applied": bool(reference is not None),
        "additional_collinearity_dropped": "|".join(additional_dropped),
        "group_position": group_position,
    }


def pooled_effect(values, y):
    values = np.asarray(values, dtype=float)
    y = np.asarray(y, dtype=float)
    a, b = values[y == 1], values[y == 0]
    n1, n0 = len(a), len(b)
    if n1 < 2 or n0 < 2:
        return np.nan, np.nan
    su = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n0 - 1) * b.var(ddof=1)) / (n1 + n0 - 2))
    du = (a.mean() - b.mean()) / su if su > 0 else np.nan
    return du, su


def joint_effect_diagnostics(score, frac, y, fraction_names=None):
    score = np.asarray(score, dtype=float)
    y = np.asarray(y, dtype=float)
    d = design_diagnostics(frac, y, fraction_names)
    du, su = pooled_effect(score, y)
    state = "OK"
    if len(score) != len(y) or not d["finite_rows"].all():
        state = "MISSING_INPUT"
    elif d["group_position"] < 0:
        state = "GROUP_NONIDENTIFIABLE"
    elif d["fit_rank"] < d["n_columns"]:
        state = "RANK_DEFICIENT_UNRESOLVED"
    elif d["fit_dof"] <= 0:
        state = "INSUFFICIENT_DOF"
    elif not np.isfinite(su) or su <= 0:
        state = "ZERO_SD"
    if state == "OK":
        beta = np.linalg.lstsq(d["reduced"], score, rcond=None)[0][d["group_position"]]
        da = beta / su
    else:
        da = np.nan
    floor = bool(np.isfinite(du) and abs(du) >= 0.1)
    no_reversal = bool(floor and np.isfinite(da) and ((du * da > 0) or da == 0))
    caf = (1 - abs(da) / abs(du)) if no_reversal else np.nan
    interpretable = bool(state == "OK" and no_reversal and np.isfinite(caf))
    return {
        "d_unadj": du, "d_joint": da, "CAF": caf, "floor_pass": floor,
        "no_sign_reversal": no_reversal, "interpretable": interpretable,
        "stalled": bool(interpretable and caf >= 0.6), "model_state": state,
        "convergence": state == "OK", "raw_rank": d["raw_rank"],
        "fit_rank": d["fit_rank"], "raw_dof": d["raw_dof"], "fit_dof": d["fit_dof"],
        "n_design_columns": d["n_columns"], "reference_dropped": d["reference_dropped"],
        "zero_variance_dropped": "|".join(d["zero_names"]),
        "additional_collinearity_dropped": d["additional_collinearity_dropped"],
        "constraint_applied": d["constraint_applied"], "sd_unadjusted": su,
    }


def joint_matrix_diagnostics(mix, frac, y, fraction_names=None):
    y = np.asarray(y, dtype=float)
    values = np.asarray(mix, dtype=float)
    d = design_diagnostics(frac, y, fraction_names)
    a, b = values[:, y == 1], values[:, y == 0]
    n1, n0 = a.shape[1], b.shape[1]
    su = np.sqrt(((n1 - 1) * a.var(axis=1, ddof=1) + (n0 - 1) * b.var(axis=1, ddof=1)) / (n1 + n0 - 2))
    with np.errstate(divide="ignore", invalid="ignore"):
        du = (a.mean(axis=1) - b.mean(axis=1)) / su
    state = "OK"
    if not d["finite_rows"].all():
        state = "MISSING_INPUT"
    elif d["group_position"] < 0:
        state = "GROUP_NONIDENTIFIABLE"
    elif d["fit_rank"] < d["n_columns"]:
        state = "RANK_DEFICIENT_UNRESOLVED"
    elif d["fit_dof"] <= 0:
        state = "INSUFFICIENT_DOF"
    if state == "OK":
        beta = np.linalg.lstsq(d["reduced"], values.T, rcond=None)[0][d["group_position"]]
        with np.errstate(divide="ignore", invalid="ignore"):
            da = beta / su
    else:
        da = np.full(values.shape[0], np.nan)
    floor = np.isfinite(du) & (np.abs(du) >= 0.1)
    no_reversal = floor & np.isfinite(da) & ((du * da > 0) | (da == 0))
    with np.errstate(divide="ignore", invalid="ignore"):
        caf = 1 - np.abs(da) / np.abs(du)
    interpretable = (state == "OK") & no_reversal & np.isfinite(caf)
    index = mix.index if isinstance(mix, pd.DataFrame) else [f"gene_{i}" for i in range(values.shape[0])]
    return pd.DataFrame({
        "gene": index, "d_unadj": du, "d_joint": da, "CAF": caf,
        "floor_pass": floor, "no_sign_reversal": no_reversal,
        "interpretable": interpretable, "stalled": interpretable & (caf >= 0.6),
        "model_state": state, "convergence": state == "OK",
        "raw_rank": d["raw_rank"], "fit_rank": d["fit_rank"],
        "raw_dof": d["raw_dof"], "fit_dof": d["fit_dof"],
        "n_design_columns": d["n_columns"], "reference_dropped": d["reference_dropped"],
        "zero_variance_dropped": "|".join(d["zero_names"]),
        "additional_collinearity_dropped": d["additional_collinearity_dropped"],
        "constraint_applied": d["constraint_applied"], "sd_unadjusted": su,
    })


def fwl_group_check(mix, frac, y, fraction_names=None):
    values = np.asarray(mix, dtype=float)
    y = np.asarray(y, dtype=float)
    d = design_diagnostics(frac, y, fraction_names)
    reduced_beta = np.linalg.lstsq(d["reduced"], values.T, rcond=None)[0][d["group_position"]]
    naive = np.column_stack([np.ones(len(y)), y, np.asarray(frac, dtype=float)])
    naive_beta = np.linalg.lstsq(naive, values.T, rcond=None)[0][1]
    fractions = np.column_stack([np.ones(len(y)), d["reduced"][:, 2:]])
    residual_y = y - fractions @ np.linalg.lstsq(fractions, y, rcond=None)[0]
    residual_x = values.T - fractions @ np.linalg.lstsq(fractions, values.T, rcond=None)[0]
    fwl_beta = residual_y @ residual_x / (residual_y @ residual_y)
    return {
        "max_reduced_vs_full": float(np.max(np.abs(reduced_beta - naive_beta))),
        "max_reduced_vs_fwl": float(np.max(np.abs(reduced_beta - fwl_beta))),
    }


def score_frozen(data, genes):
    measured = [gene for gene in genes if gene in data.index]
    sub = data.loc[measured]
    z = (sub - sub.mean(axis=0)) / (sub.std(axis=0) + 1e-12)
    return z.mean(axis=0, skipna=True).to_numpy(float), measured


def score_level_cohort(data, genes):
    measured = [gene for gene in genes if gene in data.index]
    sub = data.loc[measured]
    mu = sub.mean(axis=1).to_numpy()[:, None]
    sd = sub.std(axis=1).to_numpy()[:, None]
    z = (sub - mu) / (sd + 1e-12)
    return z.mean(axis=0, skipna=True).to_numpy(float), measured


def score_level_training_fixed(training, application, genes):
    measured = [gene for gene in genes if gene in training.index and gene in application.index]
    train = training.loc[measured]
    apply = application.loc[measured]
    mu = train.mean(axis=1).to_numpy()[:, None]
    sd = train.std(axis=1).to_numpy()[:, None]
    z = (apply - mu) / (sd + 1e-12)
    score = z.mean(axis=0, skipna=True).to_numpy(float)
    parameters = pd.DataFrame({"gene": measured, "training_mean": mu[:, 0], "training_sd": sd[:, 0]})
    return score, measured, parameters


def bh(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    q = ranked * len(pvalues) / np.arange(1, len(pvalues) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(len(pvalues))
    out[order] = np.minimum(q, 1)
    return out

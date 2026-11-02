# -*- coding: utf-8 -*-
"""Deterministic tests for compGate 0.3.0 (run: python tests/test_recap.py OR pytest)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, pandas as pd
from compgate import (recap_gate, gate_T, gate_report, attribution_report, CAR,
                      __version__)

def test_version():
    assert __version__ == "0.4.0"

def test_cohen_d_known():
    from compgate.recap import cohen_d
    a = np.array([1.0, 2.0, 3.0, 4.0]); b = np.array([3.0, 4.0, 5.0, 6.0])
    assert abs(cohen_d(a, b) + 1.54919296) < 1e-6  # pooled sample-SD convention (ddof=1)

def test_recap_gate_composition_only_collapses():
    rng = np.random.default_rng(42)
    n = 60
    # composition differs by group (sepsis neutrophil-expanded) and the gene rides on it
    neut = np.concatenate([rng.uniform(0.60, 0.90, 30), rng.uniform(0.20, 0.50, 30)])
    frac = pd.DataFrame({"Neut": neut}, index=[f"S{i}" for i in range(n)])
    meta = pd.DataFrame({"group": ["Sepsis"] * 30 + ["Healthy"] * 30}, index=frac.index)
    g = neut * 2 + rng.normal(0, 0.05, n)
    matrix = pd.DataFrame([g], index=["COMP_ONLY"], columns=frac.index)
    res = recap_gate(matrix, meta, fractions=frac)
    row = res.iloc[0]
    assert abs(row.d_unadj) > 1.0 and abs(row.d_frac_adj) < 0.5 and row.CAF_frac > 0.8

def test_recap_gate_within_cell_survives():
    rng = np.random.default_rng(7)
    n = 60
    idx = [f"S{i}" for i in range(n)]
    frac = pd.DataFrame({"Neut": rng.uniform(0.2, 0.8, n)}, index=idx)
    meta = pd.DataFrame({"group": ["Sepsis"] * 30 + ["Healthy"] * 30}, index=idx)
    shift = np.concatenate([np.full(30, 1.0), np.full(30, 0.0)])          # group shift, fraction-independent
    g = shift + rng.normal(0, 0.2, n)
    matrix = pd.DataFrame([g], index=["WITHIN"], columns=idx)
    res = recap_gate(matrix, meta, fractions=frac)
    row = res.iloc[0]
    assert row.CAF_frac < 0.3 and abs(row.d_frac_adj) > 0.8

def test_CAR():
    res = pd.DataFrame({"gene": ["a", "b", "c", "d"],
                        "d_unadj": [1.5, -1.2, 0.8, 0.05],
                        "CAF_frac": [0.9, 0.2, 0.95, np.nan]})
    car, n_rep, n_surv = CAR(("a", "b", "c"), res, caf=0.5, col="CAF_frac")
    assert abs(car - 1 / 3) < 1e-9 and (n_rep, n_surv) == (3, 1)

def test_gate_T_rise_and_return_vs_sustained():
    t = np.array([0.0, 4.0, 24.0])
    rise_ret = [np.array([1.0, 2.0, 1.05]), np.array([1.1, 2.1, 1.0])]
    sustained = [np.array([1.0, 2.0, 1.9]), np.array([1.1, 2.2, 2.0])]
    assert gate_T(t, rise_ret)["temporal_support"] == "rise-and-return"
    assert gate_T(t, sustained)["temporal_support"] == "sustained-shift"

def test_gate_report_smoke():
    res = pd.DataFrame({"gene": ["a"], "d_unadj": [1.0], "d_marker_adj": [0.1],
                        "CAF_marker": [0.9], "d_frac_adj": [0.2], "CAF_frac": [0.8]})
    rep = gate_report(res)
    assert isinstance(rep, dict) and "G6_temporal_kinetics" in rep

def test_attribution_report_quadruple_and_completeness():
    rec = attribution_report(dict(contrast="sepsis vs healthy", reference="LM22-NNLS",
                                  convention="level", estimator="unsigned equal-weight mean-Z",
                                  platform="GPL570", cohort="GSE65682"))
    q = rec["attribution_quadruple"]
    assert tuple(q) == ("contrast", "reference", "convention", "estimator")
    assert q == {"contrast": "sepsis vs healthy", "reference": "LM22-NNLS",
                 "convention": "level", "estimator": "unsigned equal-weight mean-Z"}
    assert rec["complete"] is True and rec["missing_parameters"] == []
    assert rec["stall_rule"]["caf_cut_off"] == 0.6

def test_attribution_report_flags_missing_parameters():
    rec = attribution_report(dict(contrast="ATB vs LTBI", reference="LM22-NNLS"))
    assert rec["complete"] is False
    assert rec["missing_parameters"] == ["convention", "estimator"]

def test_attribution_report_stall_counts():
    res = pd.DataFrame({"gene": list("abcd"),
                        "d_unadj": [1.5, -1.2, 0.8, 0.05],
                        "CAF_frac": [0.90, 0.20, 0.95, np.nan]})
    rec = attribution_report(dict(contrast="sepsis vs healthy", reference="LM22-NNLS",
                                  convention="level", estimator="mean-Z"),
                             per_gene=res, thresholds={"caf": 0.6})
    r = rec["result"]
    # the |d_unadj| = 0.05 row is below the 0.1 interpretability floor -> excluded
    assert r["n_genes"] == 4 and r["n_caf_interpretable"] == 3 and r["n_stalled"] == 2
    assert abs(r["stall_fraction"] - 2 / 3) < 1e-12

def test_attribution_report_rejects_unknown_convention():
    try:
        attribution_report(dict(contrast="x", reference="y", convention="z-score",
                                estimator="mean-Z"))
    except ValueError:
        return
    raise AssertionError("unknown convention must raise")

if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"all {len(fns)} tests passed")

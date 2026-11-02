# -*- coding: utf-8 -*-
"""compgate audit — CLI wrapper around recap_gate (gates 2-3) for CSV inputs."""
import argparse
import pandas as pd
from .recap import recap_gate, gate_report


def main(argv=None):
    ap = argparse.ArgumentParser(prog="compgate", description="RE-CAPT six-gate audit reference implementation")
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("audit", help="run gates 2-3 (marker / fraction residualization) on a matrix")
    a.add_argument("--matrix", required=True, help="CSV genes x samples (first column = gene)")
    a.add_argument("--meta", required=True, help="CSV samples x group (first column = sample)")
    a.add_argument("--fractions", help="CSV samples x cell types (gate 3)")
    a.add_argument("--markers", help="CSV samples x marker scores (gate 2)")
    a.add_argument("--case", default="Sepsis")
    a.add_argument("--ctrl", default="Healthy")
    a.add_argument("--out", default="compgate_audit_report.csv")
    r = sub.add_parser("report", help="emit the one-click attribution record (contrast / reference / convention / estimator)")
    r.add_argument("--contrast", required=True, help='e.g. "sepsis vs healthy"')
    r.add_argument("--reference", required=True, help='e.g. "LM22-NNLS"')
    r.add_argument("--convention", required=True, choices=["frozen", "level"])
    r.add_argument("--estimator", required=True, help='e.g. "unsigned equal-weight mean-Z"')
    r.add_argument("--platform", help="array platform, e.g. GPL570")
    r.add_argument("--cohort", help="cohort accession, e.g. GSE65682")
    r.add_argument("--genes", help="CSV from compgate audit (per-gene table) to summarise")
    r.add_argument("--caf", type=float, default=0.6, help="stall cut-off (default 0.6)")
    args = ap.parse_args(argv)

    if args.cmd == "report":
        from .recap import attribution_report
        genes = pd.read_csv(args.genes) if args.genes else None
        rec = attribution_report(
            dict(contrast=args.contrast, reference=args.reference,
                 convention=args.convention, estimator=args.estimator,
                 platform=args.platform, cohort=args.cohort),
            per_gene=genes, thresholds={"caf": args.caf})
        for k, v in rec.items():
            print(f"{k}: {v}")
        return 0


    matrix = pd.read_csv(args.matrix, index_col=0)
    meta = pd.read_csv(args.meta, index_col=0)
    fractions = pd.read_csv(args.fractions, index_col=0) if args.fractions else None
    markers = pd.read_csv(args.markers, index_col=0) if args.markers else None
    res = recap_gate(matrix, meta, case=args.case, ctrl=args.ctrl,
                     fractions=fractions, marker_scores=markers)
    res.to_csv(args.out, index=False)
    print(gate_report(res))
    print(f"[compgate] per-gene table written to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

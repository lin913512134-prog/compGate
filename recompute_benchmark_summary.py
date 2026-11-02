from pathlib import Path
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--source-data", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()
source = Path(args.source_data)
out = Path(args.output)
out.mkdir(parents=True, exist_ok=True)
frame = pd.read_csv(source / "R18_fraction_error_per_run.csv")
frame["truth_class"] = frame.truth_class.fillna("null")
assert len(frame) == 1260 and not frame.duplicated(["scenario", "replicate", "method"]).any()
assert frame.model_state.eq("OK").all()
coverage = frame.groupby("scenario").agg(n_replicates=("replicate", "nunique"), n_methods=("method", "nunique"), n_rows=("method", "size"), all_models_ok=("model_state", lambda x: bool(x.eq("OK").all()))).reset_index()
coverage["factorial_status"] = "COMPLETE_7_METHODS_X_5_SEEDS"
coverage.to_csv(out / "R18_regression_SoA_factorial_coverage.csv", index=False)
metrics = ["FPR", "FDR", "power", "gate_sensitivity", "gate_specificity", "gate_precision", "fraction_RMSE", "fraction_bias", "CAF_bias", "coverage_95", "seconds"]
summary = frame.groupby(["scenario", "truth_class", "n_per_group", "d_true_target", "method"])[metrics].agg(["median", "mean", "min", "max"]).reset_index()
summary.columns = ["_".join(str(y) for y in col if y).rstrip("_") for col in summary.columns]
assert len(summary) == 252
summary.to_csv(out / "R18_fraction_error_summary.csv", index=False)
print({"per_run_rows": len(frame), "summary_rows": len(summary), "coverage_rows": len(coverage), "model_failures": int((~frame.model_state.eq("OK")).sum())})

from pathlib import Path
import argparse
import json
import time
import tracemalloc
import numpy as np
import pandas as pd
import compgate
from compgate import joint_model as J

parser = argparse.ArgumentParser()
parser.add_argument("--output", default="verification/quickstart.json")
args = parser.parse_args()
rng = np.random.default_rng(42)
n = 60
idx = [f"S{i}" for i in range(n)]
frac = pd.DataFrame({"Neut": rng.uniform(0.2, 0.8, n), "Other": rng.uniform(0.2, 0.8, n)}, index=idx)
frac = frac.div(frac.sum(axis=1), axis=0)
meta = pd.DataFrame({"group": ["Sepsis"] * 30 + ["Healthy"] * 30}, index=idx)
y = np.r_[np.ones(30), np.zeros(30)]
shift = np.r_[np.ones(30), np.zeros(30)]
matrix = pd.DataFrame([shift + rng.normal(0, 0.2, n), rng.normal(0, 0.2, n)], index=["WITHIN", "NOISE"], columns=idx)
claim = dict(contrast="sepsis vs healthy", reference="LM22-NNLS", convention="level", estimator="unsigned equal-weight mean-Z", platform="synthetic", cohort="minimal")


def run_once():
    result = compgate.recap_gate(matrix, meta, fractions=frac)
    report = compgate.attribution_report(claim, per_gene=result)
    joint = J.joint_matrix_diagnostics(matrix, frac, y, frac.columns)
    return report, joint


tracemalloc.start()
report, joint = run_once()
peak = tracemalloc.get_traced_memory()[1]
tracemalloc.stop()
times = []
failures = 0
for _ in range(200):
    started = time.perf_counter()
    try:
        run_once()
    except Exception:
        failures += 1
    times.append((time.perf_counter() - started) * 1000)
payload = {
    "package_version": compgate.__version__,
    "scope": "two genes x 60 samples",
    "attribution_complete": report["complete"],
    "joint_states": joint.model_state.value_counts().to_dict(),
    "joint_rows": len(joint),
    "repetitions": len(times),
    "failures": failures,
    "runtime_ms_median": float(np.median(times)),
    "runtime_ms_p95": float(np.percentile(times, 95)),
    "peak_memory_kib_single_run": peak / 1024,
}
out = Path(args.output)
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
print(json.dumps(payload, indent=2))
assert compgate.__version__ == "0.4.0" and report["complete"] and len(joint) == 2 and failures == 0

from pathlib import Path
import argparse
import json
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--source-data", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()
source = Path(args.source_data)
out = Path(args.output)
out.mkdir(parents=True, exist_ok=True)
long = pd.read_csv(source / "R18_reference_convention_long.csv")
references = ["LM22-NNLS", "LM22-SVR", "REF1-NNLS"]
assert len(long) == 1308
reference_rows = []
for reference in references:
    subset = long[long.reference == reference]
    wide = subset.pivot_table(index=["cohort", "signature", "counts_toward_30"], columns="convention", values=["d_unadj", "d_joint", "CAF", "interpretable", "stalled"], aggfunc="first")
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    wide["common_interpretable"] = wide.interpretable_frozen & wide.interpretable_level
    wide.to_csv(out / f"R18_paired_cells_{reference}.csv", index=False)
    common = wide[wide.common_interpretable & (wide.counts_toward_30 == "yes")].copy()
    per_signature = common.groupby("signature").agg(common_cells=("common_interpretable", "size"), frozen_stalls=("stalled_frozen", "sum"), level_stalls=("stalled_level", "sum")).reset_index()
    per_signature["frozen_majority"] = per_signature.frozen_stalls > per_signature.common_cells / 2
    per_signature["level_majority"] = per_signature.level_stalls > per_signature.common_cells / 2
    per_signature.to_csv(out / f"R18_per_signature_{reference}.csv", index=False)
    reference_rows.append({"reference": reference, "common_cells": len(common), "frozen_stalls": int(common.stalled_frozen.sum()), "level_stalls": int(common.stalled_level.sum()), "signature_majority_frozen": int(per_signature.frozen_majority.sum()), "signature_majority_level": int(per_signature.level_majority.sum()), "n_signatures": len(per_signature)})
pd.DataFrame(reference_rows).to_csv(out / "R18_reference_paired_summary.csv", index=False)
all_wide = long.pivot_table(index=["cohort", "signature", "counts_toward_30"], columns=["reference", "convention"], values=["CAF", "interpretable", "stalled"], aggfunc="first")
all_wide.columns = [f"{reference}_{convention}_{metric}" for metric, reference, convention in all_wide.columns]
all_wide = all_wide.reset_index()
interpretable = [f"{ref}_{conv}_interpretable" for ref in references for conv in ["frozen", "level"]]
all_wide["common_all_six"] = all_wide[interpretable].all(axis=1)
common6 = all_wide[all_wide.common_all_six & (all_wide.counts_toward_30 == "yes")].copy()
common6.to_csv(out / "R18_common_all_six_cells.csv", index=False)
anchors = {"long_rows": len(long), "common_all_six": len(common6), "reference_rows": reference_rows}
(out / "paired_anchors.json").write_text(json.dumps(anchors, indent=2), encoding="utf-8")
print(json.dumps(anchors, indent=2))
assert len(common6) == 50
assert (reference_rows[0]["common_cells"], reference_rows[0]["frozen_stalls"], reference_rows[0]["level_stalls"]) == (79, 5, 31)

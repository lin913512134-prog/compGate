"""compGate — RE-CAPT six-gate composition-aware interpretation framework.

Gates: Replication -> Effect survival after composition residualization -> Composition attribution by
reference-based deconvolution (>=2 cohorts) -> Powered within-cell test ->
composition-informed Specificity margin -> Temporal kinetics (gate_T).
Also implements the CAR
(Composition-Adjusted Replication) statistic: among genes that replicate
across cohorts, the proportion whose effect survives composition adjustment.
"""
from .recap import (recap_gate, gate_T, gate_report, attribution_report,
                    ATTRIBUTION_FIELDS, CAR)
from .scoring import mean_z
from . import joint_model

__all__ = ["recap_gate", "gate_T", "gate_report", "attribution_report",
           "ATTRIBUTION_FIELDS", "CAR", "mean_z", "joint_model"]
__version__ = "0.4.0"

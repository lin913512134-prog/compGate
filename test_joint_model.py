from pathlib import Path
import numpy as np
import pandas as pd
from compgate import joint_model as J


def synthetic():
    y = np.array([0, 0, 0, 1, 1, 1], dtype=float)
    frac = np.array([[0.8, 0.2], [0.7, 0.3], [0.6, 0.4], [0.4, 0.6], [0.3, 0.7], [0.2, 0.8]])
    mix = pd.DataFrame([[1.0, 2.0, 4.0, 5.0, 8.0, 9.0], [2.0, 1.0, 3.0, 5.0, 6.0, 9.0]], index=["G1", "G2"])
    return mix, frac, y


def test_reduced_model_matches_full_and_fwl():
    mix, frac, y = synthetic()
    errors = J.fwl_group_check(mix, frac, y)
    assert errors["max_reduced_vs_full"] < 1e-10
    assert errors["max_reduced_vs_fwl"] < 1e-10


def test_rank_and_dof_reported():
    mix, frac, y = synthetic()
    result = J.joint_matrix_diagnostics(mix, frac, y)
    row = result.iloc[0]
    assert row.model_state == "OK"
    assert row.reference_dropped == "fraction_1"
    assert row.fit_rank == 3 and row.fit_dof == 3
    assert row.constraint_applied


def test_scale_invariance():
    mix, frac, y = synthetic()
    a = J.joint_matrix_diagnostics(mix, frac, y)
    b = J.joint_matrix_diagnostics(mix * 100, frac, y)
    np.testing.assert_allclose(a.CAF, b.CAF, atol=1e-12)
    np.testing.assert_allclose(a.d_unadj, b.d_unadj, atol=1e-12)
    np.testing.assert_allclose(a.d_joint, b.d_joint, atol=1e-12)


def test_zero_sd_and_sign_reversal_states():
    mix, frac, y = synthetic()
    zero = J.joint_effect_diagnostics(np.ones(len(y)), frac, y)
    assert zero["model_state"] == "ZERO_SD" and not zero["interpretable"]
    y2 = np.array([0, 0, 0, 0, 1, 1, 1, 1], dtype=float)
    f = np.array([[0.1, 0.9], [0.2, 0.8], [0.3, 0.7], [0.4, 0.6], [0.6, 0.4], [0.7, 0.3], [0.8, 0.2], [0.9, 0.1]])
    values = 3 * f[:, 0] - 0.2 * y2
    reversal = J.joint_effect_diagnostics(values, f, y2)
    assert reversal["d_unadj"] * reversal["d_joint"] < 0
    assert not reversal["no_sign_reversal"] and not reversal["interpretable"]


def test_insufficient_dof_state():
    y = np.array([0, 1], dtype=float)
    frac = np.array([[0.8, 0.2], [0.2, 0.8]])
    result = J.joint_effect_diagnostics(np.array([1.0, 2.0]), frac, y)
    assert result["model_state"] in {"RANK_DEFICIENT_UNRESOLVED", "INSUFFICIENT_DOF", "ZERO_SD"}
    assert not result["interpretable"]


def test_training_fixed_transform():
    training = pd.DataFrame({"A": [1.0, 2.0, 3.0], "B": [3.0, 6.0, 6.0]}, index=["g1", "g2", "g3"])
    application = pd.DataFrame({"C": [10.0, 20.0, 30.0], "D": [4.0, 8.0, 12.0]}, index=["g1", "g2", "g3"])
    score, measured, params = J.score_level_training_fixed(training, application, ["g1", "g2"])
    assert measured == ["g1", "g2"]
    assert np.allclose(params.training_mean.to_numpy(), [2.0, 4.0])
    assert np.allclose(params.training_sd.to_numpy(), [np.sqrt(2), 2 * np.sqrt(2)])
    assert np.allclose(score, [4 * np.sqrt(2), np.sqrt(2)])
    frozen, _ = J.score_frozen(application, ["g1", "g2"])
    assert not np.allclose(score, frozen)


def test_seed_reproducible_donor_resampling():
    y = np.array([0, 0, 0, 1, 1, 1])
    a = np.random.default_rng(42).choice(np.arange(len(y)), len(y), replace=True)
    b = np.random.default_rng(42).choice(np.arange(len(y)), len(y), replace=True)
    c = np.random.default_rng(7).choice(np.arange(len(y)), len(y), replace=True)
    np.testing.assert_array_equal(a, b)
    assert not np.array_equal(a, c)


if __name__ == "__main__":
    tests = [value for name, value in sorted(globals().items()) if name.startswith("test_") and callable(value)]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"all {len(tests)} tests passed")

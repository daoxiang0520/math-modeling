"""Validate the Question 3 result and figure contracts."""
from pathlib import Path
import os

import numpy as np
import pandas as pd


ROOT = Path(os.getenv("MODELING_OUTPUT_DIR", Path(__file__).resolve().parent))
RESULT = ROOT / "results" / "output.csv"
FIGURES = ROOT / "figures"
COLS = ["group", "bmi_low", "bmi_high", "n", "median_bmi", "t_g", "ci_low", "ci_high",
        "pi_g_at_tg", "median_uncensored", "n_unsolved", "r_cens", "t_star", "delta_t_sigma_tech"]
RANGES = {
    "bmi_low": (15, 60), "bmi_high": (15, 60), "n": (1, 300), "median_bmi": (15, 60),
    "t_g": (10, 25), "ci_low": (10, 25), "ci_high": (10, 25), "pi_g_at_tg": (0, 1),
    "median_uncensored": (10, 25), "n_unsolved": (0, 300), "r_cens": (0, 1),
    "t_star": (10, 25), "delta_t_sigma_tech": (-5, 5),
}
FIG_IDS = [
    "fig_q3_data_profile", "fig_q3_covariate_selection", "fig_q3_prob_curves",
    "fig_q3_individual_tp_hist", "fig_q3_group_t_tradeoff", "fig_q3_equivalence_check",
    "fig_q3_monotone_diagnostic", "fig_q3_error_shift_sigma", "fig_q3_bmi_boundary_bootstrap",
    "fig_q3_p_sensitivity", "fig_q3_risk_sensitivity", "fig_q3_model_structure_sensitivity",
    "fig_q3_model_coef_dualchannel",
]


def main() -> None:
    d = pd.read_csv(RESULT)
    assert 2 <= len(d) <= 6, "row count outside [2, 6]"
    assert all(c in d.columns for c in COLS), "missing contract column"
    assert not d[COLS].isna().any().any(), "NaN in contract columns"
    assert d.group.nunique() == len(d), "group labels must be unique"
    assert d.t_g.nunique() > 1, "t_g must differ or unified-time evidence must be supplied"
    assert np.all(d.bmi_low < d.bmi_high), "invalid BMI interval"
    assert np.all(d.ci_low <= d.t_g) and np.all(d.t_g <= d.ci_high), "CI does not cover estimate"
    assert np.array_equal(d.n.astype(int).to_numpy(), d.n.to_numpy()), "n must be integer"
    assert np.array_equal(d.n_unsolved.astype(int).to_numpy(), d.n_unsolved.to_numpy()), "n_unsolved must be integer"
    for c, (lo, hi) in RANGES.items():
        assert d[c].between(lo, hi, inclusive="both").all(), f"{c} outside [{lo}, {hi}]"
    assert np.allclose(d.r_cens, d.n_unsolved / d.n, atol=1e-12), "censoring rate mismatch"
    for stem in FIG_IDS:
        for ext in ("png", "pdf", "svg"):
            p = FIGURES / f"{stem}.{ext}"
            assert p.exists() and p.stat().st_size > 0, f"missing figure: {p.name}"
    print(f"PASS: {len(d)} output rows, 14 contract columns, {len(FIG_IDS)} figures × 3 formats")


if __name__ == "__main__":
    main()

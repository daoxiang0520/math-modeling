"""Validate the Question 4 result, resampling, runtime, and figure contracts."""
from pathlib import Path
import json
import os

import numpy as np
import pandas as pd

ROOT = Path(os.getenv("MODELING_OUTPUT_DIR", Path(__file__).resolve().parent))
RES = ROOT / "results"; FIG = ROOT / "figures"
COLS = ["level", "chrom", "rule", "tau", "coverage", "sens", "spec", "ppv", "npv", "f1", "auc", "ci_low", "ci_high"]
FIG_IDS = ["fig_q4_data_profile", "fig_q4_gate_curve", "fig_q4_z_roc", "fig_q4_cost_sens",
           "fig_q4_merge_compare", "fig_q4_calibration", "fig_q4_gc_sens", "fig_q4_w_trunc",
           "fig_q4_cv_bootstrap", "fig_male_female_auc", "fig_q4_time_sens"]


def main() -> None:
    q = pd.read_csv(RES / "output.csv", keep_default_na=False)
    assert 12 <= len(q) <= 200 and all(c in q for c in COLS)
    assert not q[COLS].isna().any().any() and np.isfinite(q.select_dtypes("number")).all().all()
    assert set(q.level) == {"record", "pregnant"}; assert set(q.chrom.astype(str)) == {"13", "18", "21"}
    assert set(q[q.level == "pregnant"].rule) == {"conservative", "majority", "max_risk"}
    assert len(q[q.level == "record"]) == 3 and len(q[q.level == "pregnant"]) == 9
    assert q.tau.between(-5, 10).all()
    for c in ["coverage", "sens", "spec", "ppv", "npv", "f1", "auc", "ci_low", "ci_high"]:
        assert q[c].between(0, 1).all(), c
    assert (q.ci_low <= q.ci_high).all()
    assert (RES / "q4.csv").read_bytes() == (RES / "output.csv").read_bytes()
    boot = pd.read_csv(RES / "q4_bootstrap_metrics.csv")
    assert boot.replicate.nunique() == 200 and len(boot) == 2400
    perm = pd.read_csv(RES / "q4_permutation_distribution.csv")
    assert perm.replicate.nunique() == 100 and len(perm) == 300
    summary = json.loads((RES / "q4_summary.json").read_text(encoding="utf-8"))
    assert summary["runtime_contract_passed"] and summary["runtime_seconds"] <= 90
    for stem in FIG_IDS:
        for ext in ("png", "pdf", "svg"):
            p = FIG / f"{stem}.{ext}"; assert p.exists() and p.stat().st_size > 0, p
    print(f"PASS: {len(q)} contract rows, 200 bootstraps, 100 permutations, {len(FIG_IDS)} figures × 3 formats")


if __name__ == "__main__":
    main()

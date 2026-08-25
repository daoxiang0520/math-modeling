"""Q2 Plan 1 joint selection of K, BMI cuts, and group NIPT times."""
from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd

from joint_grouping import assign_groups, bootstrap_selection, select_one_se, solve_joint
from solution_q2_plan1 import (
    BOOTSTRAP_ACTUAL, BOOTSTRAP_B, BOOTSTRAP_MC_DRAWS, BOOTSTRAP_STEP,
    GA_MAX, GA_MIN, GA_STEP, MC_DRAWS, P_LEVELS, P_MAIN, SEED, SIGMA_TECH,
    decision_eta, fit_model, grouped_cv, inverse_time, load_data, misclassification,
    probability_lookup, probability_matrix, refine_minimum, resample_mothers, resolve_paths,
)

K_MAX = 5
N_MIN = 30
GROUP_COVERAGE = 0.80
EARLY_WEIGHT = 4.0
LATE_WEIGHT = 1.0
MIN_TIME_GAP = 0.5
SELECTION_BOOTSTRAP = 300


def group_rows(tp: np.ndarray, censored: np.ndarray, bmi: np.ndarray,
               boundaries: tuple[float, ...], coverage: float = GROUP_COVERAGE,
               early_weight: float = EARLY_WEIGHT) -> list[dict[str, float]]:
    gid = assign_groups(bmi, boundaries)
    lows = [float(np.min(bmi))] + list(boundaries)
    highs = list(boundaries) + [float(np.max(bmi))]
    rows = []
    for g in range(len(boundaries) + 1):
        mask = gid == g
        values = np.asarray(tp[mask], float)
        group_censored = np.asarray(censored[mask], bool)
        decision_values = np.where(group_censored, GA_MAX + GA_STEP, values)
        try:
            raw = float(np.quantile(decision_values, coverage, method="higher"))
        except TypeError:
            raw = float(np.quantile(decision_values, coverage, interpolation="higher"))
        week = float(np.clip(np.ceil((raw - 1e-12) / GA_STEP) * GA_STEP, GA_MIN, GA_MAX))
        loss = (early_weight * np.maximum(decision_values - week, 0) +
                LATE_WEIGHT * np.maximum(week - decision_values, 0))
        uncensored_values = values[~group_censored]
        rows.append({
            "group": f"G{g + 1}", "bmi_low": lows[g], "bmi_high": highs[g],
            "n": int(mask.sum()), "median_bmi": float(np.median(bmi[mask])),
            "recommendation": week,
            "median_uncensored": (float(np.median(uncensored_values))
                                  if uncensored_values.size else float("nan")),
            "n_uncensored": int(uncensored_values.size),
            "coverage": float(np.mean((values <= week + 1e-12) & ~group_censored)),
            "expected_asymmetric_loss": float(np.mean(loss)),
            "n_unsolved": int(np.sum(censored[mask])),
            "r_cens": float(np.mean(censored[mask])),
            "n_at_lower": int(np.sum(values <= GA_MIN + 1e-12)),
        })
    return rows


def main() -> None:
    started = time.time()
    data_path, out = resolve_paths()
    results = out / "results"
    df = load_data(data_path).reset_index(drop=True)
    units = df.groupby("mother_id").agg(
        b_i=("b_i", "first"), age=("age", "median"), ivf=("ivf", lambda s: int(s.mode().iloc[0]))
    ).reset_index()
    bmi = units["b_i"].to_numpy(float)
    t = np.round(np.arange(GA_MIN, GA_MAX + 0.001, GA_STEP), 1)

    linear_cv = grouped_cv(df, "linear")
    piecewise_cv = grouped_cv(df, "piecewise")
    model, names = fit_model(df, "linear", full_random=True)
    if model.params[1] <= 0:
        raise RuntimeError(f"单调主模型失败：beta_ga={model.params[1]:.6g} <= 0")
    p, lookup_error = probability_matrix(model, t, bmi)
    min_delta = float(np.min(np.diff(p, axis=0)))
    if min_delta < -1e-10:
        raise RuntimeError(f"单调主模型失败：min delta P={min_delta:.6g}")
    tp80, cens80 = inverse_time(p, t, P_MAIN)

    tp_decision = np.where(cens80, GA_MAX + GA_STEP, tp80)
    full_solutions = solve_joint(
        tp_decision, bmi, k_max=K_MAX, n_min=N_MIN, coverage=GROUP_COVERAGE,
        early_weight=EARLY_WEIGHT, late_weight=LATE_WEIGHT, min_time_gap=MIN_TIME_GAP,
    )
    selection_boot = bootstrap_selection(
        tp_decision, bmi, repeats=SELECTION_BOOTSTRAP, seed=SEED + 17, k_max=K_MAX,
        n_min=N_MIN, coverage=GROUP_COVERAGE, early_weight=EARLY_WEIGHT,
        late_weight=LATE_WEIGHT, min_time_gap=MIN_TIME_GAP,
    )
    selected_k, k_table = select_one_se(full_solutions, selection_boot)
    selected = full_solutions[selected_k]
    boundaries = selected.boundaries
    groups = group_rows(tp80, cens80, bmi, boundaries)
    k_table["full_data_loss"] = k_table["K"].map({k: v.mean_loss for k, v in full_solutions.items()})
    k_table["boundaries"] = k_table["K"].map(
        {k: "|".join(f"{x:.1f}" for x in v.boundaries) for k, v in full_solutions.items()}
    )
    k_table["times"] = k_table["K"].map(
        {k: "|".join(f"{x:.1f}" for x in v.times) for k, v in full_solutions.items()}
    )
    k_table.to_csv(results / "tab_q2_joint_selection.csv", index=False)
    k_table.to_csv(results / "tab_q2_k_selection.csv", index=False)
    selection_boot.to_csv(results / "q2_plan1_joint_selection_bootstrap.csv", index=False)

    # Convert the selected-K bootstrap partitions to a long table for stability plots.
    boundary_rows = []
    zsel = selection_boot[selection_boot["K"] == selected_k]
    for _, row in zsel.iterrows():
        values = [] if not row["boundaries"] else [float(x) for x in str(row["boundaries"]).split("|")]
        for slot, value in enumerate(values, 1):
            boundary_rows.append({"replicate": int(row["replicate"]), "K": selected_k,
                                  "boundary_slot": slot, "selected_boundary": value})
    boundary_df = pd.DataFrame(boundary_rows)
    boundary_df.to_csv(results / "q2_plan1_boundary_bootstrap.csv", index=False)

    # Guarantee-level sensitivity using the selected BMI partition.
    sensitivity_p = []
    for guarantee in P_LEVELS:
        tp, cens = inverse_time(p, t, guarantee)
        for row in group_rows(tp, cens, bmi, boundaries):
            sensitivity_p.append({"p_guarantee": guarantee, **row})
    pd.DataFrame(sensitivity_p).to_csv(results / "tab_q2_sensitivity_p.csv", index=False)

    # Joint asymmetric-loss curves and early-risk-weight sensitivity.
    gid = assign_groups(bmi, boundaries)
    loss_rows, rho_rows = [], []
    for g, row in enumerate(groups):
        values = np.where(cens80[gid == g], GA_MAX + GA_STEP, tp80[gid == g])
        for factor in (0.5, 1.0, 2.0):
            early = EARLY_WEIGHT * factor
            try:
                q = early / (early + LATE_WEIGHT)
                week = float(np.quantile(values, q, method="higher"))
            except TypeError:
                week = float(np.quantile(values, q, interpolation="higher"))
            week = float(np.clip(np.ceil((week - 1e-12) / GA_STEP) * GA_STEP, GA_MIN, GA_MAX))
            rho_rows.append({"rho": factor, "group": row["group"], "t_star": week})
        for week in t:
            loss = EARLY_WEIGHT * np.maximum(values - week, 0) + LATE_WEIGHT * np.maximum(week - values, 0)
            loss_rows.append({"ga": week, "group": row["group"], "loss": float(np.mean(loss)),
                              "t_star": row["recommendation"]})
    pd.DataFrame(loss_rows).to_csv(results / "q2_plan1_loss_curves.csv", index=False)
    pd.DataFrame(rho_rows).to_csv(results / "q2_plan1_rho_sensitivity.csv", index=False)
    pd.DataFrame(rho_rows).to_csv(results / "tab_q2_sensitivity_rho.csv", index=False)

    # Technical-error channel, keeping the selected BMI partition fixed.
    base_rec = {r["group"]: r["recommendation"] for r in groups}
    sigma_rows, error_rows = [], []
    for factor in (0.0, 0.5, 1.0, 2.0):
        sigma = factor * SIGMA_TECH
        ps, _ = probability_matrix(model, t, bmi, sigma=sigma, seed=SEED)
        ts, cs = inverse_time(ps, t, P_MAIN)
        for row in group_rows(ts, cs, bmi, boundaries):
            delta = row["recommendation"] - base_rec[row["group"]]
            sigma_rows.append({"sigma_factor": factor, "sigma": sigma, "group": row["group"],
                               "t_p0.80_p80": row["recommendation"], "delta_t": delta,
                               "n_unsolved": row["n_unsolved"]})
            fnr, fpr = misclassification(
                model, base_rec[row["group"]], row["median_bmi"], sigma,
                seed=SEED + int(factor * 100) + int(row["group"][1:]),
            )
            error_rows.append({"sigma_factor": factor, "sigma": sigma, "group": row["group"],
                               "FNR": fnr, "FPR": fpr})
    sigma_df = pd.DataFrame(sigma_rows)
    error_df = pd.DataFrame(error_rows)
    sigma_df.to_csv(results / "q2_plan1_sigma_sensitivity.csv", index=False)
    sigma_df.merge(error_df, on=["sigma_factor", "sigma", "group"]).to_csv(
        results / "tab_q2_sensitivity_sigma.csv", index=False
    )
    error_df.to_csv(results / "q2_plan1_error_classification.csv", index=False)

    # Cluster bootstrap: refit the model and estimate fixed-partition group-time CIs.
    rng = np.random.default_rng(SEED + 100)
    boot_t = np.round(np.arange(GA_MIN, GA_MAX + 0.001, BOOTSTRAP_STEP), 1)
    group_bmi = {r["group"]: r["median_bmi"] for r in groups}
    display_bmi = np.array([group_bmi[f"G{i+1}"] for i in range(selected_k)])
    boot_rows, boot_curves = [], []
    completed = 0
    for rep in range(BOOTSTRAP_ACTUAL):
        sub = resample_mothers(df, rng)
        try:
            bm, _ = fit_model(sub, "linear", full_random=False)
            if bm.params[1] <= 0:
                continue
            bu = sub.groupby("mother_id")["b_i"].first().reset_index()
            bb = bu["b_i"].to_numpy(float)
            eta_ind = decision_eta(bm, boot_t, bb)
            eta_curve = decision_eta(bm, t, display_bmi)
            combined, _ = probability_lookup(
                bm, np.concatenate([eta_ind.ravel(), eta_curve.ravel()]),
                draws=BOOTSTRAP_MC_DRAWS, seed=SEED + 1000 + rep, nodes=512,
            )
            split = eta_ind.size
            bp = combined[:split].reshape(eta_ind.shape)
            pc = combined[split:].reshape(eta_curve.shape)
            if np.min(np.diff(bp, axis=0)) < -1e-10:
                continue
            btp, bcens = inverse_time(bp, boot_t, P_MAIN)
            recs = group_rows(btp, bcens, bb, boundaries)
            if any(r["n"] == 0 for r in recs):
                continue
            for row in recs:
                boot_rows.append({"replicate": rep, "group": row["group"],
                                  "recommendation": row["recommendation"],
                                  "coverage": row["coverage"], "n_unsolved": row["n_unsolved"]})
            for j in range(selected_k):
                boot_curves.extend({"replicate": rep, "ga": week, "group": f"G{j+1}", "p_marg": value}
                                   for week, value in zip(t, pc[:, j]))
            completed += 1
        except Exception:
            continue
    if completed < 60:
        raise RuntimeError(f"有效bootstrap仅{completed}次，低于最低60次")
    boot_df = pd.DataFrame(boot_rows)
    boot_curve_df = pd.DataFrame(boot_curves)
    boot_df.to_csv(results / "q2_plan1_bootstrap_times.csv", index=False)

    # Probability curves at each selected group's median BMI.
    main_curves, _ = probability_matrix(model, t, display_bmi)
    prob_rows = []
    for j in range(selected_k):
        group = f"G{j+1}"
        pivot = boot_curve_df[boot_curve_df["group"] == group].pivot(
            index="replicate", columns="ga", values="p_marg"
        )
        lo = pivot.quantile(0.025, axis=0).reindex(t).to_numpy()
        hi = pivot.quantile(0.975, axis=0).reindex(t).to_numpy()
        prob_rows.extend({"ga": week, "group": group, "median_bmi": display_bmi[j],
                          "p_marg": main_curves[k, j], "ci_low": lo[k], "ci_high": hi[k]}
                         for k, week in enumerate(t))
    pd.DataFrame(prob_rows).to_csv(results / "q2_plan1_prob_curves.csv", index=False)

    # Main result table.  t_p0.80_p80 is the jointly optimized group time
    # (80th percentile of within-group t_p0.80 under the coverage-0.80 constraint);
    # median_uncensored reports the uncensored-median metric required by the LTM
    # for groups with >20% right censoring (G5: 22.2 weeks, n=23).
    delta_sigma = sigma_df[np.isclose(sigma_df["sigma_factor"], 1.0)].set_index("group")["delta_t"].to_dict()
    distinct = selected_k > 1 and bool(np.all(np.diff([r["recommendation"] for r in groups]) >= MIN_TIME_GAP))
    main_rows = []
    for row in groups:
        vals = boot_df.loc[boot_df["group"] == row["group"], "recommendation"].to_numpy()
        main_rows.append({
            "group": row["group"], "bmi_low": row["bmi_low"], "bmi_high": row["bmi_high"],
            "n": row["n"], "median_bmi": row["median_bmi"],
            "t_p0.80_p80": row["recommendation"],
            "median_uncensored": row["median_uncensored"],
            "ci_low": float(np.quantile(vals, 0.025)), "ci_high": float(np.quantile(vals, 0.975)),
            "t_star": row["recommendation"], "distinct_required": distinct,
            "delta_t_sigma_tech": float(delta_sigma[row["group"]]),
            "n_unsolved": row["n_unsolved"],
        })
    main_df = pd.DataFrame(main_rows)
    main_df.to_csv(results / "q2.csv", index=False)
    main_df.to_csv(results / "tab_q2_main_results.csv", index=False)

    individual = units[["mother_id", "b_i"]].rename(columns={"b_i": "bmi"}).copy()
    individual["t_p0.80"] = tp80
    individual["censored"] = cens80
    individual["group"] = [f"G{x+1}" for x in assign_groups(bmi, boundaries)]
    individual.to_csv(results / "q2_plan1_individual_t80.csv", index=False)
    pd.DataFrame(groups).to_csv(results / "tab_q2_sensitivity_bins.csv", index=False)

    # Internal calibration.
    first = df.sort_values("ga").groupby("mother_id", as_index=False).first()
    xfirst = __import__("solution_q2_plan1").design(first, "linear")[0]
    eta_first = np.sum(xfirst * model.params[None, :], axis=1)
    pred_first, _ = probability_lookup(model, eta_first, seed=SEED + 77)
    first["pred"] = pred_first
    first["observed"] = (first["y"] >= 0.04).astype(float)
    first["group"] = [f"G{x+1}" for x in assign_groups(first["b_i"].to_numpy(), boundaries)]
    first["ga_band"] = pd.cut(first["ga"], [11, 12, 13, 15, 20], right=False,
                              labels=["11-12", "12-13", "13-15", "15-20"])
    cal_group = first.groupby("group", observed=True).agg(
        n=("observed", "size"), p_observed=("observed", "mean"), p_model=("pred", "mean")
    ).reset_index()
    cal_group["calibration_type"] = "BMI组"
    cal_band = first.dropna(subset=["ga_band"]).groupby("ga_band", observed=True).agg(
        n=("observed", "size"), p_observed=("observed", "mean"), p_model=("pred", "mean")
    ).reset_index().rename(columns={"ga_band": "group"})
    cal_band["calibration_type"] = "孕周带"
    calibration = pd.concat([cal_group, cal_band], ignore_index=True)
    calibration["error"] = calibration["p_observed"] - calibration["p_model"]
    calibration.to_csv(results / "q2_plan1_calibration.csv", index=False)
    calibration.to_csv(results / "tab_q2_calibration.csv", index=False)

    # Main-model coefficient and monotonicity diagnostics.
    piecewise, _ = fit_model(df, "piecewise", full_random=False)
    diag_rows = []
    for form_model, form_name in ((model, "linear_monotone"), (piecewise, "piecewise_diagnostic")):
        curves, _ = probability_matrix(form_model, t, display_bmi, seed=SEED + 91)
        for j in range(selected_k):
            diag_rows.extend({"ga": week, "group": f"G{j+1}", "model_form": form_name, "p_marg": value}
                             for week, value in zip(t, curves[:, j]))
    pd.DataFrame(diag_rows).to_csv(results / "q2_plan1_monotone_diagnostic.csv", index=False)
    pd.DataFrame([
        {"model_form": "linear", "cv_rmse": linear_cv[0], "cv_mae": linear_cv[1], "role": "主决策模型"},
        {"model_form": "piecewise", "cv_rmse": piecewise_cv[0], "cv_mae": piecewise_cv[1], "role": "单调性诊断"},
    ]).to_csv(results / "tab_q2_model_selection.csv", index=False)
    pd.DataFrame({"term": names, "estimate": model.params, "std_error": model.std_errors,
                  "p_value": model.pvalues}).to_csv(results / "tab_q2_model_coefficients_plan1.csv", index=False)

    exact_boundary_stability = 1.0 if selected_k == 1 else float(
        zsel["boundaries"].eq("|".join(f"{x:.1f}" for x in boundaries)).mean()
    )
    summary = {
        "records": int(len(df)), "mothers": int(len(units)),
        "model": "linear Beta + REML random intercept",
        "selection": "joint K, half-BMI boundaries, and asymmetric-loss group times",
        "selected_K": selected_k, "boundaries": list(boundaries),
        "group_coverage_target": GROUP_COVERAGE, "early_weight": EARLY_WEIGHT,
        "late_weight": LATE_WEIGHT, "n_min": N_MIN, "min_time_gap": MIN_TIME_GAP,
        "selection_bootstrap": SELECTION_BOOTSTRAP,
        "rho_semantics": "multiplier on the early-risk weight (4*rho : 1) in the "
                         "asymmetric group-time loss; replaces the former rho=1 t* risk-view loss",
        "exact_boundary_vector_stability": exact_boundary_stability,
        "beta_ga": float(model.params[1]), "beta_ga_p": float(model.pvalues[1]),
        "phi": model.phi, "random_intercept_variance": model.var_b0,
        "min_delta_p_marg": min_delta, "p_guarantee_main": P_MAIN,
        "bootstrap_target": BOOTSTRAP_B, "bootstrap_actual": BOOTSTRAP_ACTUAL,
        "bootstrap_completed": completed, "bootstrap_grid_step": BOOTSTRAP_STEP,
        "bootstrap_mc_draws": BOOTSTRAP_MC_DRAWS, "mc_draws": MC_DRAWS,
        "lookup_nodes_main": 1024, "lookup_nodes_bootstrap": 512,
        "max_lookup_error": lookup_error, "distinct_required": distinct,
        "runtime_seconds": float(time.time() - started),
        "method_limit": "two-stage Beta+REML; selection bootstrap is conditional on the fitted probability model",
    }
    (results / "q2_plan1_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

"""Q2 Plan 1 revised: monotone linear Beta random-intercept decision model.

This module provides the shared data/model/bootstrap primitives used by the
joint optimizer (solution_q2_joint.py), which is the authoritative Plan 1
pipeline and produces the outputs in results/.  The clinical concentration
threshold is 4%; p=0.80 is a decision guarantee level and is explicitly
sensitivity-tested.  The former fixed K=2 / median-recommendation pipeline in
main_legacy_deprecated() is superseded and kept only as historical reference.
"""
from __future__ import annotations

import json
import os
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.stats import beta as beta_dist
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GroupKFold
import statsmodels.api as sm
from statsmodels.othermod.betareg import BetaModel

from solution_q2 import load_data

Y_THR = 0.04
P_MAIN = 0.80
P_LEVELS = (0.75, 0.80, 0.85, 0.90)
GA_MIN, GA_MAX, GA_STEP = 10.0, 25.0, 0.1
SIGMA_TECH = 0.133
MC_DRAWS = 1000
BOOTSTRAP_B = 100
BOOTSTRAP_ACTUAL = 60
BOOTSTRAP_MC_DRAWS = 500
BOOTSTRAP_STEP = 0.2
SEED = 2025
BOUNDARY = 30.0
BOUNDARY_CANDIDATES = (24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0)


@dataclass
class Model:
    params: np.ndarray
    pvalues: np.ndarray
    std_errors: np.ndarray
    phi: float
    var_b0: float
    form: str
    age_ref: float
    ivf_ref: int
    random_method: str


def resolve_paths() -> tuple[Path, Path]:
    raw = os.environ.get("MODELING_DATA_PATH")
    if not raw:
        paths = json.loads(os.environ.get("MODELING_DATA_PATHS", "[]"))
        raw = paths[0] if paths else None
    if not raw:
        local = Path.cwd() / "附件.xlsx"
        raw = str(local) if local.exists() else None
    if not raw:
        raise FileNotFoundError("请通过 MODELING_DATA_PATH 或 MODELING_DATA_PATHS 传入附件.xlsx")
    out = Path(os.environ.get("MODELING_OUTPUT_DIR", Path.cwd()))
    (out / "results").mkdir(parents=True, exist_ok=True)
    return Path(raw), out


def design(df: pd.DataFrame, form: str = "linear", gc: bool = False) -> tuple[np.ndarray, list[str]]:
    ga = df["ga"].to_numpy(float)
    cols = [ga]
    names = ["ga"]
    if form == "piecewise":
        cols.extend([np.maximum(ga - 12.5, 0), np.maximum(ga - 20.0, 0)])
        names.extend(["hinge_12_5", "hinge_20"])
    cols.extend([df["b_i"].to_numpy(float), df["age"].to_numpy(float), df["ivf"].to_numpy(float)])
    names.extend(["bmi", "age", "ivf"])
    if gc:
        cols.append(df["gc"].to_numpy(float))
        names.append("gc")
    return sm.add_constant(np.column_stack(cols), has_constant="add"), ["const"] + names


def decision_eta(model: Model, t: np.ndarray, b: np.ndarray) -> np.ndarray:
    tt, bb = np.meshgrid(np.asarray(t, float), np.asarray(b, float), indexing="ij")
    cols = [tt.ravel()]
    if model.form == "piecewise":
        cols.extend([np.maximum(tt.ravel() - 12.5, 0), np.maximum(tt.ravel() - 20.0, 0)])
    cols.extend([bb.ravel(), np.full(tt.size, model.age_ref), np.full(tt.size, model.ivf_ref)])
    x = sm.add_constant(np.column_stack(cols), has_constant="add")
    return (x @ model.params).reshape(len(t), len(b))


def fit_beta(y: np.ndarray, x: np.ndarray):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for method in ("bfgs", "newton", "nm"):
            try:
                result = BetaModel(y, x).fit(method=method, maxiter=350, disp=False)
                if np.all(np.isfinite(result.params)):
                    return result
            except Exception:
                pass
    raise RuntimeError("Beta regression failed")


def quick_random_intercept_var(df: pd.DataFrame, x: np.ndarray, beta: np.ndarray) -> float:
    resid = df["logit_y"].to_numpy() - x @ beta
    key = "boot_id" if "boot_id" in df else "mother_id"
    tmp = pd.DataFrame({"key": df[key].to_numpy(), "r": resid})
    stats = tmp.groupby("key")["r"].agg(["mean", "count"])
    within = tmp.assign(centered=tmp["r"] - tmp.groupby("key")["r"].transform("mean"))["centered"].var(ddof=1)
    corrected = stats["mean"].var(ddof=1) - within * np.mean(1.0 / stats["count"].to_numpy())
    return float(max(corrected, 1e-8))


def fit_model(df: pd.DataFrame, form: str = "linear", full_random: bool = True) -> tuple[Model, list[str]]:
    x, names = design(df, form)
    fit = fit_beta(df["y"].to_numpy(), x)
    p = x.shape[1]
    params = np.asarray(fit.params)[:p]
    var_b0 = None
    method = "empirical residual random-intercept variance"
    if full_random:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                mixed = sm.MixedLM(df["logit_y"], x, groups=df["mother_id"]).fit(
                    reml=True, method="lbfgs", maxiter=220, disp=False
                )
                var_b0 = float(np.asarray(mixed.cov_re)[0, 0])
                method = "REML random intercept"
            except Exception:
                var_b0 = None
    if var_b0 is None or not np.isfinite(var_b0) or var_b0 <= 0:
        var_b0 = quick_random_intercept_var(df.reset_index(drop=True), x, params)
    units = df.groupby("mother_id").agg(age=("age", "median"), ivf=("ivf", lambda s: int(s.mode().iloc[0])))
    model = Model(
        params=params,
        pvalues=np.asarray(fit.pvalues)[:p],
        std_errors=np.asarray(fit.bse)[:p],
        phi=float(np.exp(np.clip(np.asarray(fit.params)[p], -15, 15))),
        var_b0=float(var_b0),
        form=form,
        age_ref=float(units["age"].median()),
        ivf_ref=int(units["ivf"].mode().iloc[0]),
        random_method=method,
    )
    return model, names


def grouped_cv(df: pd.DataFrame, form: str) -> tuple[float, float]:
    y = df["y"].to_numpy()
    pred = np.empty(len(df))
    groups = df["mother_id"].to_numpy()
    for train, test in GroupKFold(5).split(df, y, groups):
        xtr, _ = design(df.iloc[train], form)
        xte, _ = design(df.iloc[test], form)
        fit = fit_beta(y[train], xtr)
        pred[test] = expit(xte @ np.asarray(fit.params)[:xtr.shape[1]])
    return float(np.sqrt(mean_squared_error(y, pred))), float(mean_absolute_error(y, pred))


def probability_lookup(model: Model, eta: np.ndarray, sigma: float = 0.0,
                       draws: int = MC_DRAWS, seed: int = SEED,
                       nodes: int = 1024) -> tuple[np.ndarray, float]:
    """MC-integrate random intercept and logit measurement error on a 1-D eta lookup."""
    rng = np.random.default_rng(seed)
    b0 = rng.normal(0, np.sqrt(model.var_b0), draws)
    eps = rng.normal(0, sigma, draws) if sigma > 0 else np.zeros(draws)
    lo = float(np.min(eta) + np.min(b0)) - 0.05
    hi = float(np.max(eta) + np.max(b0)) + 0.05
    grid = np.linspace(lo, hi, nodes)
    vals = np.empty(nodes)
    threshold = expit(logit(Y_THR) - eps)
    for start in range(0, nodes, 64):
        z = grid[start:start + 64, None] + b0[None, :]
        mu = np.clip(expit(z), 1e-8, 1 - 1e-8)
        vals[start:start + 64] = beta_dist.sf(
            threshold[None, :], mu * model.phi, (1 - mu) * model.phi
        ).mean(axis=1)
    out = np.interp(np.asarray(eta).ravel(), grid, vals).reshape(np.asarray(eta).shape)
    # Direct audit at evenly spaced eta values using the same MC samples.
    audit_eta = np.linspace(float(np.min(eta)), float(np.max(eta)), 15)
    z = audit_eta[:, None] + b0[None, :]
    mu = np.clip(expit(z), 1e-8, 1 - 1e-8)
    direct = beta_dist.sf(threshold[None, :], mu * model.phi, (1 - mu) * model.phi).mean(axis=1)
    approx = np.interp(audit_eta, grid, vals)
    return np.clip(out, 0, 1), float(np.max(np.abs(direct - approx)))


def probability_matrix(model: Model, t: np.ndarray, b: np.ndarray, sigma: float = 0.0,
                       draws: int = MC_DRAWS, seed: int = SEED) -> tuple[np.ndarray, float]:
    return probability_lookup(model, decision_eta(model, t, b), sigma, draws, seed)


def inverse_time(p: np.ndarray, t: np.ndarray, guarantee: float) -> tuple[np.ndarray, np.ndarray]:
    p = np.maximum.accumulate(np.asarray(p), axis=0)
    solved = np.any(p >= guarantee, axis=0)
    idx = np.argmax(p >= guarantee, axis=0)
    tp = np.full(p.shape[1], GA_MAX, float)
    for j in np.flatnonzero(solved):
        k = int(idx[j])
        if k == 0:
            tp[j] = t[0]
        else:
            den = p[k, j] - p[k - 1, j]
            frac = 0.0 if den <= 0 else (guarantee - p[k - 1, j]) / den
            tp[j] = t[k - 1] + np.clip(frac, 0, 1) * (t[k] - t[k - 1])
    return tp, ~solved


def group_recommendations(tp: np.ndarray, censored: np.ndarray, b: np.ndarray,
                          boundary: float = BOUNDARY) -> list[dict[str, float]]:
    result = []
    for label, mask in (("G1", b < boundary), ("G2", b >= boundary)):
        n = int(mask.sum())
        n_cens = int(censored[mask].sum())
        rate = n_cens / n
        eligible = mask & ~censored if rate > 0.20 and np.any(mask & ~censored) else mask
        result.append({
            "group": label,
            "n": n,
            "median_bmi": float(np.median(b[mask])),
            "recommendation": float(np.median(tp[eligible])),
            "n_unsolved": n_cens,
            "r_cens": rate,
            "n_at_lower": int(np.sum(tp[mask] <= GA_MIN + 1e-9)),
        })
    return result


def refine_minimum(t: np.ndarray, y: np.ndarray) -> float:
    k = int(np.argmin(y))
    if k == 0 or k == len(t) - 1:
        return float(t[k])
    coef = np.polyfit(t[k - 1:k + 2], y[k - 1:k + 2], 2)
    if coef[0] <= 0:
        return float(t[k])
    return float(np.clip(-coef[1] / (2 * coef[0]), t[k - 1], t[k + 1]))


def boundary_choice(tp: np.ndarray, b: np.ndarray) -> float:
    best = None
    for boundary in BOUNDARY_CANDIDATES:
        left, right = b < boundary, b >= boundary
        if left.sum() < 20 or right.sum() < 20:
            continue
        score = np.abs(tp[left] - np.median(tp[left])).sum() + np.abs(tp[right] - np.median(tp[right])).sum()
        item = (float(score), abs(boundary - BOUNDARY), boundary)
        if best is None or item < best:
            best = item
    return float(best[2]) if best else np.nan


def resample_mothers(df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    ids = df["mother_id"].unique()
    selected = rng.choice(ids, size=len(ids), replace=True)
    frames = []
    for k, mid in enumerate(selected):
        part = df[df["mother_id"] == mid].copy()
        part["boot_id"] = f"B{k:03d}"
        part["mother_id"] = part["boot_id"]
        frames.append(part)
    return pd.concat(frames, ignore_index=True)


def misclassification(model: Model, week: float, bmi: float, sigma: float,
                      draws: int = 30000, seed: int = SEED) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    eta = float(decision_eta(model, np.array([week]), np.array([bmi]))[0, 0])
    b0 = rng.normal(0, np.sqrt(model.var_b0), draws)
    mu = np.clip(expit(eta + b0), 1e-8, 1 - 1e-8)
    y_true = rng.beta(mu * model.phi, (1 - mu) * model.phi)
    y_obs = expit(logit(np.clip(y_true, 1e-9, 1 - 1e-9)) + rng.normal(0, sigma, draws))
    fnr = np.mean((y_true >= Y_THR) & (y_obs < Y_THR))
    fpr = np.mean((y_true < Y_THR) & (y_obs >= Y_THR))
    return float(fnr), float(fpr)


def main_legacy_deprecated() -> None:
    """DEPRECATED: superseded fixed K=2 / median-recommendation pipeline.

    Kept only as historical reference.  The active Plan 1 contract jointly
    selects K, half-BMI boundaries and 80%-coverage group times; run
    solution_q2_joint.py (or this module as __main__) instead.
    """
    start = time.time()
    data_path, out = resolve_paths()
    results = out / "results"
    df = load_data(data_path).reset_index(drop=True)
    units = df.groupby("mother_id").agg(
        b_i=("b_i", "first"), age=("age", "median"), ivf=("ivf", lambda s: int(s.mode().iloc[0]))
    ).reset_index()
    b = units["b_i"].to_numpy(float)
    t = np.round(np.arange(GA_MIN, GA_MAX + 0.001, GA_STEP), 1)

    linear_cv = grouped_cv(df, "linear")
    piecewise_cv = grouped_cv(df, "piecewise")
    model, names = fit_model(df, "linear", full_random=True)
    if model.params[1] <= 0:
        raise RuntimeError(f"单调主模型失败：beta_ga={model.params[1]:.6g} <= 0")
    p, lookup_error = probability_matrix(model, t, b)
    min_delta = float(np.min(np.diff(p, axis=0)))
    if min_delta < -1e-10:
        raise RuntimeError(f"单调主模型失败：min delta P={min_delta:.6g}")
    tp80, cens80 = inverse_time(p, t, P_MAIN)
    groups = group_recommendations(tp80, cens80, b)

    # Guarantee-level sensitivity.
    sens_p = []
    for guarantee in P_LEVELS:
        tp, cens = inverse_time(p, t, guarantee)
        for row in group_recommendations(tp, cens, b):
            sens_p.append({"p_guarantee": guarantee, **row})
    pd.DataFrame(sens_p).to_csv(results / "tab_q2_sensitivity_p.csv", index=False)

    # Risk-view t* and rho sensitivity at group median BMI.
    loss_rows, rho_rows = [], []
    tstars = {}
    for row in groups:
        mask = b < BOUNDARY if row["group"] == "G1" else b >= BOUNDARY
        curve, _ = probability_matrix(model, t, np.array([np.median(b[mask])]), seed=SEED + 11)
        curve = curve[:, 0]
        for rho in (0.5, 1.0, 2.0):
            loss = (t - 10.0) / 17.0 + rho * (1 - curve)
            star = refine_minimum(t, loss)
            rho_rows.append({"rho": rho, "group": row["group"], "t_star": star})
            if rho == 1.0:
                tstars[row["group"]] = star
                loss_rows.extend({"ga": w, "group": row["group"], "p_marg": pp, "delay_risk": rr,
                                  "loss": ll, "t_star": star}
                                 for w, pp, rr, ll in zip(t, curve, (t - 10) / 17, loss))
    pd.DataFrame(loss_rows).to_csv(results / "q2_plan1_loss_curves.csv", index=False)
    pd.DataFrame(rho_rows).to_csv(results / "q2_plan1_rho_sensitivity.csv", index=False)
    pd.DataFrame(rho_rows).to_csv(results / "tab_q2_sensitivity_rho.csv", index=False)

    # Measurement-error sensitivity and classification.
    sigma_rows, error_rows = [], []
    base_rec = {r["group"]: r["recommendation"] for r in groups}
    for factor in (0.0, 0.5, 1.0, 2.0):
        sigma = factor * SIGMA_TECH
        # Common random numbers: sigma=0 exactly reproduces the main probability baseline.
        ps, _ = probability_matrix(model, t, b, sigma=sigma, seed=SEED)
        ts, cs = inverse_time(ps, t, P_MAIN)
        for row in group_recommendations(ts, cs, b):
            delta = row["recommendation"] - base_rec[row["group"]]
            sigma_rows.append({"sigma_factor": factor, "sigma": sigma, "group": row["group"],
                               "t_p0.80_median": row["recommendation"], "delta_t": delta,
                               "n_unsolved": row["n_unsolved"]})
            fnr, fpr = misclassification(model, base_rec[row["group"]], row["median_bmi"], sigma,
                                         seed=SEED + int(factor * 100) + (1 if row["group"] == "G1" else 2))
            error_rows.append({"sigma_factor": factor, "sigma": sigma, "group": row["group"],
                               "FNR": fnr, "FPR": fpr})
    sigma_df = pd.DataFrame(sigma_rows)
    error_df = pd.DataFrame(error_rows)
    sigma_df.to_csv(results / "q2_plan1_sigma_sensitivity.csv", index=False)
    sigma_df.merge(error_df, on=["sigma_factor", "sigma", "group"]).to_csv(
        results / "tab_q2_sensitivity_sigma.csv", index=False
    )
    error_df.to_csv(results / "q2_plan1_error_classification.csv", index=False)

    # Bootstrap parameter/sampling uncertainty; fixed seeds make the run reproducible.
    rng = np.random.default_rng(SEED + 100)
    boot_rows, boot_boundaries, boot_curves = [], [], []
    boot_t = np.round(np.arange(GA_MIN, GA_MAX + 0.001, BOOTSTRAP_STEP), 1)
    group_bmi = {r["group"]: r["median_bmi"] for r in groups}
    completed = 0
    for rep in range(BOOTSTRAP_ACTUAL):
        sub = resample_mothers(df, rng)
        try:
            bm, _ = fit_model(sub, "linear", full_random=False)
            if bm.params[1] <= 0:
                continue
            bu = sub.groupby("mother_id")["b_i"].first().reset_index()
            bb = bu["b_i"].to_numpy(float)
            # One shared lookup covers the bootstrap individuals and the two display curves.
            eta_ind = decision_eta(bm, boot_t, bb)
            eta_curve = decision_eta(bm, t, np.array([group_bmi["G1"], group_bmi["G2"]]))
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
            recs = group_recommendations(btp, bcens, bb)
            if any(r["n"] == 0 for r in recs):
                continue
            for row in recs:
                boot_rows.append({"replicate": rep, "group": row["group"],
                                  "recommendation": row["recommendation"], "n_unsolved": row["n_unsolved"]})
            boot_boundaries.append({"replicate": rep, "selected_boundary": boundary_choice(btp, bb)})
            for j, group in enumerate(("G1", "G2")):
                boot_curves.extend({"replicate": rep, "ga": week, "group": group, "p_marg": value}
                                   for week, value in zip(t, pc[:, j]))
            completed += 1
        except Exception:
            continue
    if completed < 60:
        raise RuntimeError(f"有效bootstrap仅{completed}次，低于最低60次")
    boot_df = pd.DataFrame(boot_rows)
    boundary_df = pd.DataFrame(boot_boundaries)
    boot_curve_df = pd.DataFrame(boot_curves)
    boot_df.to_csv(results / "q2_plan1_bootstrap_times.csv", index=False)
    boundary_df.to_csv(results / "q2_plan1_boundary_bootstrap.csv", index=False)

    # Main probability curves with bootstrap 95% bands.
    main_curves, _ = probability_matrix(model, t, np.array([group_bmi["G1"], group_bmi["G2"]]))
    prob_rows = []
    for j, group in enumerate(("G1", "G2")):
        subset = boot_curve_df[boot_curve_df["group"] == group].pivot(index="replicate", columns="ga", values="p_marg")
        lo = subset.quantile(0.025, axis=0).reindex(t).to_numpy()
        hi = subset.quantile(0.975, axis=0).reindex(t).to_numpy()
        prob_rows.extend({"ga": week, "group": group, "median_bmi": group_bmi[group],
                          "p_marg": main_curves[k, j], "ci_low": lo[k], "ci_high": hi[k]}
                         for k, week in enumerate(t))
    pd.DataFrame(prob_rows).to_csv(results / "q2_plan1_prob_curves.csv", index=False)

    # Bootstrap CIs and main contract.
    delta_sigma = sigma_df[np.isclose(sigma_df["sigma_factor"], 1.0)].set_index("group")["delta_t"].to_dict()
    main_rows = []
    for row in groups:
        vals = boot_df.loc[boot_df["group"] == row["group"], "recommendation"].to_numpy()
        mask = b < BOUNDARY if row["group"] == "G1" else b >= BOUNDARY
        main_rows.append({
            "group": row["group"],
            "bmi_low": float(np.min(b[mask])),
            "bmi_high": BOUNDARY if row["group"] == "G1" else float(np.max(b[mask])),
            "n": row["n"],
            "median_bmi": row["median_bmi"],
            "t_p0.80_median": row["recommendation"],
            "ci_low": float(np.quantile(vals, 0.025)),
            "ci_high": float(np.quantile(vals, 0.975)),
            "t_star": tstars[row["group"]],
            "distinct_required": bool(abs(groups[1]["recommendation"] - groups[0]["recommendation"]) >= 0.5),
            "delta_t_sigma_tech": float(delta_sigma[row["group"]]),
            "n_unsolved": row["n_unsolved"],
        })
    main_df = pd.DataFrame(main_rows)
    main_df.to_csv(results / "q2.csv", index=False)
    main_df.to_csv(results / "tab_q2_main_results.csv", index=False)

    # Individual source and boundary sensitivity.
    individual = units[["mother_id", "b_i"]].rename(columns={"b_i": "bmi"}).copy()
    individual["t_p0.80"] = tp80
    individual["censored"] = cens80
    individual["group"] = np.where(b < BOUNDARY, "G1", "G2")
    individual.to_csv(results / "q2_plan1_individual_t80.csv", index=False)
    boundary_sens = []
    for boundary in (28.0, 30.0, 32.0):
        for row in group_recommendations(tp80, cens80, b, boundary):
            boundary_sens.append({"boundary": boundary, **row})
    pd.DataFrame(boundary_sens).to_csv(results / "tab_q2_sensitivity_bins.csv", index=False)

    # Internal calibration using each mother's first observation.
    first = df.sort_values("ga").groupby("mother_id", as_index=False).first()
    eta_first = np.sum(design(first, "linear")[0] * model.params[None, :], axis=1)
    pred_first, _ = probability_lookup(model, eta_first, seed=SEED + 77)
    first["pred"] = pred_first
    first["observed"] = (first["y"] >= Y_THR).astype(float)
    first["group"] = np.where(first["b_i"] < BOUNDARY, "G1", "G2")
    first["ga_band"] = pd.cut(first["ga"], [11, 12, 13, 15, 20], right=False,
                              labels=["11-12", "12-13", "13-15", "15-20"])
    cal_group = first.groupby("group", observed=True).agg(n=("observed", "size"),
        p_observed=("observed", "mean"), p_model=("pred", "mean")).reset_index()
    cal_group["calibration_type"] = "BMI组"
    cal_band = first.dropna(subset=["ga_band"]).groupby("ga_band", observed=True).agg(
        n=("observed", "size"), p_observed=("observed", "mean"), p_model=("pred", "mean")).reset_index()
    cal_band = cal_band.rename(columns={"ga_band": "group"})
    cal_band["calibration_type"] = "孕周带"
    calibration = pd.concat([cal_group, cal_band], ignore_index=True)
    calibration["error"] = calibration["p_observed"] - calibration["p_model"]
    calibration.to_csv(results / "q2_plan1_calibration.csv", index=False)
    calibration.to_csv(results / "tab_q2_calibration.csv", index=False)

    # Linear versus piecewise monotonicity diagnostic.
    piecewise, piece_names = fit_model(df, "piecewise", full_random=False)
    diag_rows = []
    for form_model, form_name in ((model, "linear_monotone"), (piecewise, "piecewise_diagnostic")):
        curves, _ = probability_matrix(form_model, t, np.array([group_bmi["G1"], group_bmi["G2"]]), seed=SEED + 91)
        for j, group in enumerate(("G1", "G2")):
            diag_rows.extend({"ga": week, "group": group, "model_form": form_name, "p_marg": value}
                             for week, value in zip(t, curves[:, j]))
    pd.DataFrame(diag_rows).to_csv(results / "q2_plan1_monotone_diagnostic.csv", index=False)

    pd.DataFrame([
        {"model_form": "linear", "cv_rmse": linear_cv[0], "cv_mae": linear_cv[1], "role": "主决策模型"},
        {"model_form": "piecewise", "cv_rmse": piecewise_cv[0], "cv_mae": piecewise_cv[1], "role": "单调性诊断"},
    ]).to_csv(results / "tab_q2_model_selection.csv", index=False)
    pd.DataFrame({"term": names, "estimate": model.params, "std_error": model.std_errors,
                  "p_value": model.pvalues}).to_csv(results / "tab_q2_model_coefficients_plan1.csv", index=False)

    delta_boot = (boot_df.pivot(index="replicate", columns="group", values="recommendation").dropna())
    delta_vals = (delta_boot["G2"] - delta_boot["G1"]).to_numpy()
    all_boot = []
    for rep, sub in boot_df.groupby("replicate"):
        # A size-weighted pooled median cannot be reconstructed from group medians; use midpoint only as summary absent raw tp.
        all_boot.append(float(np.mean(sub["recommendation"])))
    summary = {
        "records": int(len(df)), "mothers": int(len(units)), "model": "linear Beta + REML random intercept",
        "beta_ga": float(model.params[1]), "beta_ga_p": float(model.pvalues[1]), "phi": model.phi,
        "random_intercept_variance": model.var_b0, "age_ref": model.age_ref, "ivf_ref": model.ivf_ref,
        "min_delta_p_marg": min_delta, "p_guarantee_main": P_MAIN,
        "bootstrap_target": BOOTSTRAP_B, "bootstrap_actual": BOOTSTRAP_ACTUAL,
        "bootstrap_completed": completed, "bootstrap_grid_step": BOOTSTRAP_STEP,
        "bootstrap_mc_draws": BOOTSTRAP_MC_DRAWS, "mc_draws": MC_DRAWS,
        "lookup_nodes_main": 1024, "lookup_nodes_bootstrap": 512, "max_lookup_error": lookup_error,
        "fixed_boundary": BOUNDARY, "boundary_selection_rule": "minimum within-group absolute deviation; n>=20 each side",
        "delta_t_point": float(groups[1]["recommendation"] - groups[0]["recommendation"]),
        "delta_t_ci": [float(np.quantile(delta_vals, 0.025)), float(np.quantile(delta_vals, 0.975))],
        "distinct_required": bool(main_df["distinct_required"].iloc[0]),
        "unified_t_point": float(np.median(tp80)),
        "runtime_seconds": float(time.time() - start),
        "method_limit": "two-stage Beta fixed effects plus REML random intercept; bootstrap uses empirical residual variance",
    }
    (results / "q2_plan1_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    # Active Plan 1 contract: jointly selects K, BMI cuts and group times.
    from solution_q2_joint import main as joint_main
    joint_main()

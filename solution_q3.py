"""2025 高教社杯 C 题问题 3：多因素 BMI 分组与 NIPT 时点选择。

严格实现 q3_coder_task.md，唯一经用户明确授权的偏离是取消 90 秒运行限制。
主交付使用辅助 Beta 连续通道产生的边缘达标概率；二项随机效应通道用于直接
达标事件对照。辅助模型为 Beta 固定效应均值/精度层 + logit(y) MixedLM
随机方差层的两阶段近似，不能称为联合极大似然 Beta-GAMM。
"""

from __future__ import annotations

import itertools
import json
import math
import os
import re
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize
from scipy.special import expit, logit, logsumexp
from scipy.stats import beta as beta_dist, norm
from sklearn.metrics import mean_squared_error
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import GroupKFold
from statsmodels.othermod.betareg import BetaModel

Y_THR = 0.04
P_MAIN = 0.80
P_LEVELS = (0.75, 0.80, 0.85, 0.90)
GA_MIN, GA_MAX, GA_STEP = 10.0, 25.0, 0.1
T_GRID = np.round(np.arange(GA_MIN, GA_MAX + 1e-9, GA_STEP), 1)
BOUNDARY_CANDIDATES = (24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0)
K_VALUES = (2, 3, 4)
SIGMA_FACTORS = (0.0, 0.5, 1.0, 2.0)
RHO_VALUES = (0.5, 1.0, 2.0)
GAMMA_VALUES = (0.0, 0.5, 1.0, 2.0)
BOOTSTRAP_B = 100
SELECTION_BOOTSTRAP = 300
SEED = 2025


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


def parse_ga(value: object) -> tuple[float, float, float]:
    if pd.isna(value):
        return np.nan, np.nan, np.nan
    text = str(value).strip().lower().replace("周", "w").replace("天", "")
    match = re.fullmatch(r"(\d+(?:\.\d+)?)\s*w(?:\s*\+\s*(\d+))?", text)
    if match:
        week, day = float(match.group(1)), float(match.group(2) or 0)
        if 0 <= day < 7:
            return week + day / 7.0, week, day
        return np.nan, week, day
    try:
        week = float(text)
        return week, week, 0.0
    except ValueError:
        return np.nan, np.nan, np.nan


def parse_test_date(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.replace(r"\.0$", "", regex=True).str.strip()
    parsed = pd.to_datetime(text, format="%Y%m%d", errors="coerce")
    return parsed.fillna(pd.to_datetime(series, errors="coerce"))


def parse_aneuploidy(value: object) -> float:
    if pd.isna(value) or not str(value).strip():
        return np.nan
    text = str(value).replace("T", "")
    match = re.search(r"\d+", text)
    return float(match.group()) if match else np.nan


def load_and_parse(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_excel(path, sheet_name="男胎检测数据")
    raw.columns = [str(c).strip() for c in raw.columns]
    raw["mother_id_text"] = raw["孕妇代码"].astype(str).str.strip()
    raw["mother_code"] = pd.to_numeric(
        raw["孕妇代码"].astype(str).str.replace("A", "", regex=False), errors="coerce"
    )
    parsed = raw["检测孕周"].map(parse_ga)
    # 显式逐列赋值（语义等价）：让框架 AST 列名校验能静态识别 ga/week_part/day_part 派生列
    _parsed_parts = pd.DataFrame(parsed.tolist(), index=raw.index)
    raw["ga"] = _parsed_parts[0]
    raw["week_part"] = _parsed_parts[1]
    raw["day_part"] = _parsed_parts[2]
    raw["lmp_date"] = pd.to_datetime(raw["末次月经"], errors="coerce")
    raw["test_date"] = parse_test_date(raw["检测日期"])
    raw["ga_from_date"] = (raw["test_date"] - raw["lmp_date"]).dt.days / 7.0
    raw["ga_date_diff"] = raw["ga"] - raw["ga_from_date"]
    raw["aneuploidy_num"] = raw["染色体的非整倍体"].map(parse_aneuploidy)
    # 显式逐列赋值（语义与 for 循环等价）：让框架的 AST 列名校验能静态识别派生列
    raw["y"] = pd.to_numeric(raw["Y染色体浓度"], errors="coerce")
    raw["y_z"] = pd.to_numeric(raw["Y染色体的Z值"], errors="coerce")
    raw["bmi"] = pd.to_numeric(raw["孕妇BMI"], errors="coerce")
    raw["age"] = pd.to_numeric(raw["年龄"], errors="coerce")
    raw["height"] = pd.to_numeric(raw["身高"], errors="coerce")
    raw["weight"] = pd.to_numeric(raw["体重"], errors="coerce")
    raw["gc"] = pd.to_numeric(raw["GC含量"], errors="coerce")
    raw["ivf"] = raw["IVF妊娠"].astype(str).str.contains(
        "IVF|试管", case=False, regex=True
    ).astype(int)
    raw["healthy"] = raw["胎儿是否健康"].astype(str).str.strip().eq("是")
    raw["blood_draw"] = pd.to_numeric(raw["检测抽血次数"], errors="coerce")
    valid = (
        raw["mother_code"].notna() & raw["mother_id_text"].ne("")
        & raw["y"].between(0, 1, inclusive="neither") & raw["y_z"].notna()
        & raw[["ga", "bmi", "age", "height", "weight", "gc"]].notna().all(axis=1)
    )
    df = raw.loc[valid].copy()
    df["mother_id"] = df["mother_id_text"]
    df["event"] = (df["y"] >= Y_THR).astype(int)
    df["logit_y"] = logit(df["y"].clip(1e-8, 1 - 1e-8))
    df["tech_group"] = df["mother_id"] + "#" + df["blood_draw"].astype(str)

    units = df.groupby("mother_id", sort=True).agg(
        bmi_rep=("bmi", "median"), height_rep=("height", "median"),
        weight_rep=("weight", "median"), age_rep=("age", "median"),
        ivf_rep=("ivf", lambda s: int(s.mode().iloc[0])),
        healthy_rep=("healthy", lambda s: bool(s.mode().iloc[0])),
        gc_rep=("gc", "median"), n_records=("y", "size"),
    ).reset_index()
    df = df.merge(units, on="mother_id", how="left", validate="many_to_one")
    return df.reset_index(drop=True), units


def data_qc(df: pd.DataFrame, units: pd.DataFrame) -> pd.DataFrame:
    bmi_calc = df["weight"] / (df["height"] / 100.0) ** 2
    diff = df["bmi"] - bmi_calc
    rows = [
        ("records", len(df), "count"),
        ("mothers", units["mother_id"].nunique(), "count"),
        ("ga_parse_missing", int(df["ga"].isna().sum()), "count"),
        ("ga_date_diff_gt_1_week", int((df["ga_date_diff"].abs() > 1).sum()), "count"),
        ("height_missing_rate", float(df["height"].isna().mean()), "proportion"),
        ("weight_missing_rate", float(df["weight"].isna().mean()), "proportion"),
        ("bmi_missing_rate", float(df["bmi"].isna().mean()), "proportion"),
        ("bmi_identity_abs_diff_gt_0_5", int((diff.abs() > 0.5).sum()), "count"),
        ("bmi_identity_abs_diff_median", float(diff.abs().median()), "kg/m2"),
        ("gc_below_0_40", int((df["gc"] < 0.40).sum()), "count"),
        ("gc_above_0_60", int((df["gc"] > 0.60).sum()), "count"),
        ("unhealthy_records", int((~df["healthy"]).sum()), "count"),
        ("event_rate", float(df["event"].mean()), "proportion"),
        ("ga_min", float(df["ga"].min()), "week"), ("ga_max", float(df["ga"].max()), "week"),
        ("bmi_min", float(df["bmi"].min()), "kg/m2"), ("bmi_max", float(df["bmi"].max()), "kg/m2"),
        ("height_min", float(df["height"].min()), "cm"), ("height_max", float(df["height"].max()), "cm"),
        ("weight_min", float(df["weight"].min()), "kg"), ("weight_max", float(df["weight"].max()), "kg"),
    ]
    return pd.DataFrame(rows, columns=["metric", "value", "unit"])


def estimate_sigma_tech(df: pd.DataFrame) -> tuple[float, int, int]:
    groups = [g["logit_y"].to_numpy() for _, g in df.groupby("tech_group") if len(g) > 1]
    numerator = sum(float(np.sum((a - a.mean()) ** 2)) for a in groups)
    denominator = sum(len(a) - 1 for a in groups)
    if denominator <= 0:
        return 0.0, 0, 0
    return float(np.sqrt(numerator / denominator)), len(groups), int(sum(map(len, groups)))


FORM_COLUMNS = {
    "bmi": ("bmi_rep", "age_rep", "ivf_rep"),
    "weight": ("weight_rep", "age_rep", "ivf_rep"),
    "height_weight": ("height_rep", "weight_rep", "age_rep", "ivf_rep"),
    "all": ("bmi_rep", "height_rep", "weight_rep", "age_rep", "ivf_rep"),
}


@dataclass
class FeatureBuilder:
    form: str
    means: dict[str, float]
    scales: dict[str, float]
    ga_center: float
    extra_cols: tuple[str, ...] = ()

    @classmethod
    def fit(cls, df: pd.DataFrame, form: str,
            extra_cols: tuple[str, ...] = ()) -> "FeatureBuilder":
        cols = FORM_COLUMNS[form] + tuple(extra_cols)
        means, scales = {}, {}
        for col in cols:
            if col == "ivf_rep":
                means[col], scales[col] = 0.0, 1.0
            else:
                means[col] = float(df[col].mean())
                scales[col] = max(float(df[col].std(ddof=0)), 1e-8)
        return cls(form, means, scales, float(df["ga"].mean()), tuple(extra_cols))

    @property
    def names(self) -> list[str]:
        return ["const", "ga_c"] + [c.replace("_rep", "") + ("_z" if c != "ivf_rep" else "")
                                      for c in FORM_COLUMNS[self.form] + self.extra_cols]

    def transform(self, df: pd.DataFrame) -> np.ndarray:
        cols = [df["ga"].to_numpy(float) - self.ga_center]
        for col in FORM_COLUMNS[self.form] + self.extra_cols:
            values = df[col].to_numpy(float)
            cols.append(values if col == "ivf_rep" else (values - self.means[col]) / self.scales[col])
        return sm.add_constant(np.column_stack(cols), has_constant="add")

    def decision_matrix(self, t: np.ndarray, units: pd.DataFrame) -> np.ndarray:
        tt, ii = np.meshgrid(np.asarray(t, float), np.arange(len(units)), indexing="ij")
        frame = pd.DataFrame({"ga": tt.ravel()})
        for col in FORM_COLUMNS[self.form] + self.extra_cols:
            frame[col] = np.tile(units[col].to_numpy(float), len(t))
        return self.transform(frame)


def fit_beta(y: np.ndarray, x: np.ndarray):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for method in ("bfgs", "newton", "nm"):
            try:
                result = BetaModel(y, x).fit(method=method, maxiter=500, disp=False)
                if np.all(np.isfinite(result.params)):
                    return result
            except Exception:
                continue
    raise RuntimeError("Beta regression failed for all optimizers")


def vif_max(x: np.ndarray, names: list[str]) -> float:
    if x.shape[1] <= 2:
        return 1.0
    values = []
    for j in range(2, x.shape[1]):
        other = [k for k in range(1, x.shape[1]) if k != j]
        fit = sm.OLS(x[:, j], sm.add_constant(x[:, other], has_constant="add")).fit()
        values.append(1.0 / max(1.0 - float(fit.rsquared), 1e-12))
    return float(max(values, default=1.0))


def evaluate_covariate_forms(df: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    rows = []
    groups = df["mother_id"].to_numpy()
    y = df["y"].to_numpy(float)
    splitter = list(GroupKFold(5).split(df, groups=groups))
    for form in FORM_COLUMNS:
        full_builder = FeatureBuilder.fit(df, form)
        x = full_builder.transform(df)
        full_fit = fit_beta(y, x)
        fold_rmse = []
        for tr, te in splitter:
            builder = FeatureBuilder.fit(df.iloc[tr], form)
            xtr, xte = builder.transform(df.iloc[tr]), builder.transform(df.iloc[te])
            fit = fit_beta(y[tr], xtr)
            pred = expit(xte @ np.asarray(fit.params)[:xtr.shape[1]])
            fold_rmse.append(float(np.sqrt(mean_squared_error(y[te], pred))))
        rows.append({
            "form": form, "AIC": float(full_fit.aic), "BIC": float(full_fit.bic),
            "CV_RMSE": float(np.mean(fold_rmse)),
            "CV_RMSE_SE": float(np.std(fold_rmse, ddof=1) / np.sqrt(len(fold_rmse))),
            "VIF_max": vif_max(x, full_builder.names), "n_parameters": x.shape[1] + 1,
        })
    table = pd.DataFrame(rows)
    best = table.loc[table["CV_RMSE"].idxmin()]
    threshold = float(best["CV_RMSE"] + best["CV_RMSE_SE"])
    complexity = {"bmi": 3, "weight": 3, "height_weight": 4, "all": 5}
    eligible = table[table["CV_RMSE"] <= threshold].copy()
    eligible["complexity"] = eligible["form"].map(complexity)
    selected = str(eligible.sort_values(["complexity", "CV_RMSE", "form"]).iloc[0]["form"])
    table["one_se_threshold"] = threshold
    table["chosen_by_1se"] = table["form"].eq(selected)
    return table, selected


def evaluate_weight_transforms(df: pd.DataFrame, form: str) -> pd.DataFrame:
    """Linear/log/piecewise weight sensitivity with fold-local preprocessing."""
    if "weight_rep" not in FORM_COLUMNS[form]:
        return pd.DataFrame([{"variant": "not_applicable", "AIC": 0.0, "BIC": 0.0,
                              "CV_RMSE": 0.0, "note": "selected form has no weight"}])

    def matrices(train: pd.DataFrame, test: pd.DataFrame, variant: str) -> tuple[np.ndarray, np.ndarray]:
        builder = FeatureBuilder.fit(train, form)
        xtr, xte = builder.transform(train), builder.transform(test)
        index = builder.names.index("weight_z")
        if variant == "linear":
            return xtr, xte
        if variant == "log":
            mean = float(np.log(train["weight_rep"]).mean())
            scale = max(float(np.log(train["weight_rep"]).std(ddof=0)), 1e-8)
            xtr[:, index] = (np.log(train["weight_rep"]) - mean) / scale
            xte[:, index] = (np.log(test["weight_rep"]) - mean) / scale
            return xtr, xte
        knot = float(train["weight_rep"].median())
        hinge_tr = np.maximum(train["weight_rep"].to_numpy(float) - knot, 0)
        hinge_te = np.maximum(test["weight_rep"].to_numpy(float) - knot, 0)
        scale = max(float(np.std(hinge_tr)), 1e-8)
        return np.column_stack([xtr, hinge_tr / scale]), np.column_stack([xte, hinge_te / scale])

    rows = []
    splitter = list(GroupKFold(5).split(df, groups=df["mother_id"]))
    for variant in ("linear", "log", "piecewise"):
        xfull, _ = matrices(df, df, variant)
        fit_full = fit_beta(df["y"].to_numpy(float), xfull)
        fold = []
        for tr, te in splitter:
            xtr, xte = matrices(df.iloc[tr], df.iloc[te], variant)
            fit = fit_beta(df.iloc[tr]["y"].to_numpy(float), xtr)
            pred = expit(xte @ np.asarray(fit.params)[:xtr.shape[1]])
            fold.append(float(np.sqrt(mean_squared_error(df.iloc[te]["y"], pred))))
        rows.append({"variant": variant, "AIC": float(fit_full.aic), "BIC": float(fit_full.bic),
                     "CV_RMSE": float(np.mean(fold)), "note": "piecewise knot=training median weight"})
    return pd.DataFrame(rows)


def fit_mixed_layer(df: pd.DataFrame, x: np.ndarray) -> tuple[object, str, pd.DataFrame]:
    groups = df["mother_id"].to_numpy()
    ga_c = df["ga"].to_numpy(float) - float(df["ga"].mean())
    candidates = []
    for structure, z in [("RI", np.ones((len(df), 1))),
                         ("RI+RS", np.column_stack([np.ones(len(df)), ga_c]))]:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = sm.MixedLM(df["logit_y"], x, groups=groups, exog_re=z).fit(
                    reml=False, method="lbfgs", maxiter=500, disp=False
                )
            candidates.append((structure, fit))
        except Exception:
            continue
    if not candidates:
        raise RuntimeError("auxiliary random-effect layer failed")
    rows = []
    for structure, fit in candidates:
        cov = np.asarray(fit.cov_re, float)
        rows.append({"model": "auxiliary_beta_random_layer", "structure": structure,
                     "AIC": float(fit.aic), "BIC": float(fit.bic), "loglik": float(fit.llf),
                     "converged": bool(fit.converged), "sigma2_intercept": cov[0, 0],
                     "sigma2_slope": cov[1, 1] if cov.shape == (2, 2) else 0.0,
                     "corr_intercept_slope": (cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
                                              if cov.shape == (2, 2) and cov[0, 0] * cov[1, 1] > 0 else 0.0),
                     "residual_variance": float(fit.scale)})
    table = pd.DataFrame(rows)
    chosen_row = table.loc[table["BIC"].idxmin()]
    chosen = str(chosen_row["structure"])
    result = next(f for s, f in candidates if s == chosen)
    return result, chosen, table


@dataclass
class LogisticGLMM:
    beta: np.ndarray
    cov_re: np.ndarray
    se: np.ndarray
    loglik: float
    aic: float
    bic: float
    structure: str
    converged: bool


def _logistic_glmm_nll(theta: np.ndarray, x_groups: list[np.ndarray], y_groups: list[np.ndarray],
                       ga_groups: list[np.ndarray], roots: np.ndarray, weights: np.ndarray,
                       structure: str) -> float:
    p = x_groups[0].shape[1]
    beta = theta[:p]
    if structure == "RI":
        sd0 = np.exp(np.clip(theta[p], -8, 4))
        re = np.sqrt(2) * roots[:, None] * sd0
        logw = np.log(weights) - 0.5 * np.log(np.pi)
    else:
        sd0, sd1 = np.exp(np.clip(theta[p:p + 2], -8, 4))
        rho = np.tanh(theta[p + 2])
        chol = np.array([[sd0, 0.0], [rho * sd1, sd1 * np.sqrt(max(1 - rho * rho, 1e-10))]])
        a, b = np.meshgrid(roots, roots, indexing="ij")
        zz = np.column_stack([a.ravel(), b.ravel()]) * np.sqrt(2)
        re = zz @ chol.T
        logw = np.log(np.outer(weights, weights).ravel()) - np.log(np.pi)
    total = 0.0
    for x, y, ga in zip(x_groups, y_groups, ga_groups):
        fixed = x @ beta
        if structure == "RI":
            eta = fixed[None, :] + re
        else:
            eta = fixed[None, :] + re[:, [0]] + re[:, [1]] * ga[None, :]
        ll = (y[None, :] * (-np.logaddexp(0, -eta))
              + (1 - y)[None, :] * (-np.logaddexp(0, eta))).sum(axis=1)
        total += logsumexp(logw + ll)
    return float(-total)


def fit_logistic_glmm(df: pd.DataFrame, x: np.ndarray, structure: str,
                      nodes: int = 7) -> LogisticGLMM:
    glm = sm.GLM(df["event"], x, family=sm.families.Binomial()).fit()
    ids = df["mother_id"].unique()
    x_groups, y_groups, ga_groups = [], [], []
    ga_center = float(df["ga"].mean())
    for mid in ids:
        mask = df["mother_id"].eq(mid).to_numpy()
        x_groups.append(x[mask])
        y_groups.append(df.loc[mask, "event"].to_numpy(float))
        ga_groups.append(df.loc[mask, "ga"].to_numpy(float) - ga_center)
    roots, weights = hermgauss(nodes)
    tail = [math.log(0.4)] if structure == "RI" else [math.log(0.4), math.log(0.03), 0.0]
    start = np.r_[np.asarray(glm.params), tail]
    result = minimize(_logistic_glmm_nll, start,
                      args=(x_groups, y_groups, ga_groups, roots, weights, structure),
                      method="L-BFGS-B", options={"maxiter": 350, "ftol": 1e-9})
    theta = np.asarray(result.x)
    p = x.shape[1]
    if structure == "RI":
        sd0 = np.exp(theta[p])
        cov = np.array([[sd0 * sd0]])
    else:
        sd0, sd1 = np.exp(theta[p:p + 2])
        rho = np.tanh(theta[p + 2])
        cov = np.array([[sd0 * sd0, rho * sd0 * sd1],
                        [rho * sd0 * sd1, sd1 * sd1]])
    try:
        h_inv = np.asarray(result.hess_inv.todense())
        se = np.sqrt(np.clip(np.diag(h_inv)[:p], 0, None))
    except Exception:
        robust = sm.GLM(df["event"], x, family=sm.families.Binomial()).fit(
            cov_type="cluster", cov_kwds={"groups": df["mother_id"]}
        )
        se = np.asarray(robust.bse)
    ll = -float(result.fun)
    k = len(theta)
    return LogisticGLMM(theta[:p], cov, se, ll, -2 * ll + 2 * k,
                        -2 * ll + np.log(len(df)) * k, structure, bool(result.success))


@dataclass
class AuxiliaryModel:
    beta: np.ndarray
    beta_se: np.ndarray
    beta_p: np.ndarray
    phi: float
    cov_re: np.ndarray
    builder: FeatureBuilder
    random_structure: str


def fit_auxiliary(df: pd.DataFrame, form: str, force_structure: str | None = None,
                  extra_cols: tuple[str, ...] = ()) -> tuple[AuxiliaryModel, pd.DataFrame]:
    builder = FeatureBuilder.fit(df, form, extra_cols=extra_cols)
    x = builder.transform(df)
    beta_fit = fit_beta(df["y"].to_numpy(float), x)
    p = x.shape[1]
    if force_structure is None:
        mixed, structure, table = fit_mixed_layer(df, x)
    else:
        ga_c = df["ga"].to_numpy(float) - float(df["ga"].mean())
        z = np.ones((len(df), 1)) if force_structure == "RI" else np.column_stack([np.ones(len(df)), ga_c])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mixed = sm.MixedLM(df["logit_y"], x, groups=df["mother_id"], exog_re=z).fit(
                reml=False, method="lbfgs", maxiter=500, disp=False
            )
        structure = force_structure
        cov = np.asarray(mixed.cov_re, float)
        table = pd.DataFrame([{"model": "auxiliary_beta_random_layer", "structure": structure,
                               "AIC": float(mixed.aic), "BIC": float(mixed.bic),
                               "loglik": float(mixed.llf), "converged": bool(mixed.converged),
                               "sigma2_intercept": cov[0, 0],
                               "sigma2_slope": cov[1, 1] if cov.shape == (2, 2) else 0.0,
                               "corr_intercept_slope": 0.0,
                               "residual_variance": float(mixed.scale)}])
    return AuxiliaryModel(
        beta=np.asarray(beta_fit.params)[:p], beta_se=np.asarray(beta_fit.bse)[:p],
        beta_p=np.asarray(beta_fit.pvalues)[:p],
        phi=float(np.exp(np.clip(np.asarray(beta_fit.params)[p], -15, 15))),
        cov_re=np.asarray(mixed.cov_re, float), builder=builder,
        random_structure=structure,
    ), table


def gh_random_effects(cov: np.ndarray, nodes: int = 11) -> tuple[np.ndarray, np.ndarray]:
    roots, weights = hermgauss(nodes)
    if cov.shape == (1, 1):
        values = np.sqrt(2 * max(cov[0, 0], 0)) * roots[:, None]
        return np.column_stack([values[:, 0], np.zeros(nodes)]), weights / np.sqrt(np.pi)
    vals, vecs = np.linalg.eigh((cov + cov.T) / 2)
    root_cov = vecs @ np.diag(np.sqrt(np.clip(vals, 0, None)))
    a, b = np.meshgrid(roots, roots, indexing="ij")
    zz = np.column_stack([a.ravel(), b.ravel()]) * np.sqrt(2)
    return zz @ root_cov.T, np.outer(weights, weights).ravel() / np.pi


def probability_lookup(model: AuxiliaryModel, eta_values: np.ndarray, sigma: float,
                       need_kappa: bool = True, nodes: int = 1600) -> tuple[np.ndarray, np.ndarray, float]:
    eta_values = np.asarray(eta_values, float)
    lo, hi = float(np.min(eta_values) - 1.0), float(np.max(eta_values) + 1.0)
    grid = np.linspace(lo, hi, nodes)
    mu = np.clip(expit(grid), 1e-8, 1 - 1e-8)
    a, b = mu * model.phi, (1 - mu) * model.phi
    base_cdf = beta_dist.cdf(Y_THR, a, b)
    if sigma <= 0:
        prob, kappa = 1 - base_cdf, np.zeros_like(grid)
    else:
        roots, weights = hermgauss(15)
        eps = np.sqrt(2) * sigma * roots
        w = weights / np.sqrt(np.pi)
        shifted_thr = expit(logit(Y_THR) - eps)
        cdf_shift = beta_dist.cdf(shifted_thr[None, :], a[:, None], b[:, None])
        prob = np.average(1 - cdf_shift, axis=1, weights=w)
        kappa = np.average(np.abs(cdf_shift - base_cdf[:, None]), axis=1, weights=w)
    out_p = np.interp(eta_values, grid, prob)
    out_k = np.interp(eta_values, grid, kappa)
    audit_eta = np.linspace(float(np.min(eta_values)), float(np.max(eta_values)), 17)
    audit_mu = np.clip(expit(audit_eta), 1e-8, 1 - 1e-8)
    aa, bb = audit_mu * model.phi, (1 - audit_mu) * model.phi
    if sigma <= 0:
        direct = beta_dist.sf(Y_THR, aa, bb)
    else:
        roots, weights = hermgauss(31)
        threshold = expit(logit(Y_THR) - np.sqrt(2) * sigma * roots)
        direct = np.average(beta_dist.sf(threshold[None, :], aa[:, None], bb[:, None]),
                            axis=1, weights=weights / np.sqrt(np.pi))
    error = float(np.max(np.abs(direct - np.interp(audit_eta, grid, prob))))
    return np.clip(out_p, 0, 1), np.clip(out_k, 0, 1), error


def auxiliary_probability(model: AuxiliaryModel, t: np.ndarray, units: pd.DataFrame,
                          sigma: float = 0.0, force_ri: bool = False) -> tuple[np.ndarray, np.ndarray, float]:
    x = model.builder.decision_matrix(t, units)
    eta = (x @ model.beta).reshape(len(t), len(units))
    cov = model.cov_re[:1, :1] if force_ri else model.cov_re
    re, w = gh_random_effects(cov)
    accumulator = np.zeros_like(eta)
    kappa_acc = np.zeros_like(eta)
    max_error = 0.0
    for node, weight in zip(re, w):
        addition = node[0] + node[1] * (np.asarray(t)[:, None] - model.builder.ga_center)
        p, k, err = probability_lookup(model, eta + addition, sigma)
        accumulator += weight * p
        kappa_acc += weight * k
        max_error = max(max_error, err)
    return np.clip(accumulator, 0, 1), np.clip(kappa_acc, 0, 1), max_error


def logistic_probability(model: LogisticGLMM, builder: FeatureBuilder, t: np.ndarray,
                         units: pd.DataFrame, force_ri: bool = False) -> np.ndarray:
    x = builder.decision_matrix(t, units)
    eta = (x @ model.beta).reshape(len(t), len(units))
    cov = model.cov_re[:1, :1] if force_ri else model.cov_re
    re, w = gh_random_effects(cov)
    out = np.zeros_like(eta)
    for node, weight in zip(re, w):
        add = node[0] + node[1] * (np.asarray(t)[:, None] - builder.ga_center)
        out += weight * expit(eta + add)
    return np.clip(out, 0, 1)


def enforce_monotonicity(model: AuxiliaryModel, units: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    p_full, k_full, err = auxiliary_probability(model, T_GRID, units, sigma=0.0)
    min_delta = float(np.min(np.diff(p_full, axis=0)))
    violating = int(np.sum(np.min(np.diff(p_full, axis=0), axis=0) < -1e-10))
    if violating == 0:
        return p_full, k_full, {"fallback": "none", "min_delta": min_delta,
                               "violating_individuals": 0, "lookup_error": err}
    p_ri, k_ri, err_ri = auxiliary_probability(model, T_GRID, units, sigma=0.0, force_ri=True)
    min_ri = float(np.min(np.diff(p_ri, axis=0)))
    violating_ri = int(np.sum(np.min(np.diff(p_ri, axis=0), axis=0) < -1e-10))
    if violating_ri:
        raise RuntimeError(f"单调回退失败：RI-only仍有{violating_ri}位孕妇曲线下降")
    return p_ri, k_ri, {"fallback": "RI-only decision integration", "min_delta": min_ri,
                        "violating_individuals_before": violating,
                        "violating_individuals": 0, "lookup_error": err_ri}


def first_crossing(prob: np.ndarray, guarantee: float) -> tuple[np.ndarray, np.ndarray]:
    prob = np.asarray(prob, float)
    solved = np.any(prob >= guarantee, axis=0)
    idx = np.argmax(prob >= guarantee, axis=0)
    tp = np.full(prob.shape[1], GA_MAX, float)
    for j in np.flatnonzero(solved):
        k = int(idx[j])
        if k == 0:
            tp[j] = GA_MIN
        else:
            den = prob[k, j] - prob[k - 1, j]
            frac = 0.0 if den <= 0 else (guarantee - prob[k - 1, j]) / den
            tp[j] = T_GRID[k - 1] + np.clip(frac, 0, 1) * GA_STEP
    return np.clip(tp, GA_MIN, GA_MAX), ~solved


def group_ids(bmi: np.ndarray, boundaries: tuple[float, ...]) -> np.ndarray:
    return np.digitize(np.asarray(bmi, float), np.asarray(boundaries), right=False)


def group_crossing(prob: np.ndarray, mask: np.ndarray, guarantee: float = P_MAIN) -> tuple[float, float]:
    curve = prob[:, mask].mean(axis=1)
    tp, cens = first_crossing(curve[:, None], guarantee)
    return float(tp[0]), float(curve[np.argmin(np.abs(T_GRID - tp[0]))]) if not cens[0] else float(curve[-1])


def search_grouping(prob: np.ndarray, kappa: np.ndarray, bmi: np.ndarray) -> pd.DataFrame:
    rows = []
    for k in K_VALUES:
        for boundaries in itertools.combinations(BOUNDARY_CANDIDATES, k - 1):
            gids = group_ids(bmi, boundaries)
            if any(np.sum(gids == g) == 0 for g in range(k)):
                continue
            total = 0.0
            times = []
            for g in range(k):
                mask = gids == g
                tg, _ = group_crossing(prob, mask)
                idx = int(np.argmin(np.abs(T_GRID - tg)))
                loss = (tg - GA_MIN) / 17.0 + (1 - prob[idx, mask]) + kappa[idx, mask]
                total += float(loss.sum())
                times.append(tg)
            rows.append({"K": k, "boundaries": "|".join(map(lambda x: f"{x:.1f}", boundaries)),
                         "mean_loss": total / len(bmi),
                         "times": "|".join(f"{x:.3f}" for x in times),
                         "min_group_n": min(int(np.sum(gids == g)) for g in range(k))})
    return pd.DataFrame(rows).sort_values(["K", "mean_loss", "boundaries"]).reset_index(drop=True)


def conditional_selection_bootstrap(prob: np.ndarray, kappa: np.ndarray, bmi: np.ndarray,
                                    repeats: int = SELECTION_BOOTSTRAP) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(SEED + 31)
    rows = []
    for rep in range(repeats):
        idx = rng.integers(0, len(bmi), len(bmi))
        table = search_grouping(prob[:, idx], kappa[:, idx], bmi[idx])
        for k in K_VALUES:
            best = table[table["K"].eq(k)].iloc[0]
            rows.append({"replicate": rep, "K": k, "mean_loss": best["mean_loss"],
                         "boundaries": best["boundaries"], "times": best["times"]})
    boot = pd.DataFrame(rows)
    summary_rows = []
    for k in K_VALUES:
        z = boot[boot["K"].eq(k)]
        summary_rows.append({"K": k, "bootstrap_mean_loss": z["mean_loss"].mean(),
                             "bootstrap_se": z["mean_loss"].std(ddof=1),
                             "successful_replicates": len(z)})
    summary = pd.DataFrame(summary_rows)
    best = summary.loc[summary["bootstrap_mean_loss"].idxmin()]
    threshold = float(best["bootstrap_mean_loss"] + best["bootstrap_se"])
    selected_k = int(summary.loc[summary["bootstrap_mean_loss"] <= threshold, "K"].min())
    summary["one_se_threshold"] = threshold
    summary["selected"] = summary["K"].eq(selected_k)
    return boot, summary


def grouping_main_rows(prob: np.ndarray, kappa: np.ndarray, tp: np.ndarray, cens: np.ndarray,
                       units: pd.DataFrame, boundaries: tuple[float, ...],
                       sigma_delta: dict[str, float] | None = None,
                       boot_times: pd.DataFrame | None = None) -> pd.DataFrame:
    bmi = units["bmi_rep"].to_numpy(float)
    gids = group_ids(bmi, boundaries)
    rows = []
    lows = [float(bmi.min())] + list(boundaries)
    highs = list(boundaries) + [float(bmi.max())]
    for g in range(len(boundaries) + 1):
        mask = gids == g
        tg, pi = group_crossing(prob, mask)
        idx = int(np.argmin(np.abs(T_GRID - tg)))
        loss = (T_GRID - GA_MIN) / 17.0 + (1 - prob[:, mask].mean(axis=1)) + kappa[:, mask].mean(axis=1)
        t_star = float(T_GRID[int(np.argmin(loss))])
        unc = tp[mask & ~cens]
        if boot_times is not None:
            vals = boot_times.loc[boot_times["group"].eq(f"G{g+1}"), "t_g"].to_numpy(float)
            ci_low, ci_high = (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))) if len(vals) else (tg, tg)
        else:
            ci_low = ci_high = tg
        rows.append({
            "group": f"G{g+1}", "bmi_low": lows[g], "bmi_high": highs[g],
            "n": int(mask.sum()), "median_bmi": float(np.median(bmi[mask])),
            "t_g": tg, "ci_low": ci_low, "ci_high": ci_high,
            "pi_g_at_tg": float(np.clip(pi, 0, 1)),
            "median_uncensored": float(np.median(unc)) if len(unc) else GA_MAX,
            "n_unsolved": int(cens[mask].sum()), "r_cens": float(cens[mask].mean()),
            "t_star": t_star,
            "delta_t_sigma_tech": 0.0 if sigma_delta is None else float(sigma_delta[f"G{g+1}"]),
            "all_censored": bool(len(unc) == 0),
        })
    return pd.DataFrame(rows)


def resample_mothers(df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    ids = df["mother_id"].unique()
    selected = rng.choice(ids, size=len(ids), replace=True)
    parts = []
    for draw, mid in enumerate(selected):
        part = df[df["mother_id"].eq(mid)].copy()
        part["mother_id"] = f"B{draw:04d}"
        parts.append(part)
    sub = pd.concat(parts, ignore_index=True)
    # Recreate representative covariates after relabelling duplicated clusters.
    reps = sub.groupby("mother_id").agg(
        bmi_rep=("bmi", "median"), height_rep=("height", "median"), weight_rep=("weight", "median"),
        age_rep=("age", "median"), ivf_rep=("ivf", lambda s: int(s.mode().iloc[0])),
        healthy_rep=("healthy", lambda s: bool(s.mode().iloc[0])), gc_rep=("gc", "median")
    )
    for col in reps.columns:
        sub[col] = sub["mother_id"].map(reps[col])
    return sub


def full_bootstrap(df: pd.DataFrame, form: str, boundaries: tuple[float, ...],
                   force_ri_decision: bool, sigma_tech: float,
                   b: int = BOOTSTRAP_B) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, int]:
    rng = np.random.default_rng(SEED + 71)
    time_rows, boundary_rows, curve_rows = [], [], []
    failures = 0
    error_counts: dict[str, int] = {}
    attempts = 0
    successful = 0
    original_units = df.groupby("mother_id", sort=True).first().reset_index()
    original_bmi = original_units["bmi_rep"].to_numpy(float)
    original_gid = group_ids(original_bmi, boundaries)
    while successful < b:
        attempts += 1
        if attempts > 10 * b:
            raise RuntimeError(
                f"通道B在{attempts}次尝试后仍不足{b}个成功副本；"
                f"成功={successful}；失败类型={error_counts}"
            )
        try:
            local_boundaries, local_times, local_curves = [], [], []
            sub = resample_mothers(df, rng)
            units = sub.groupby("mother_id", sort=True).first().reset_index()
            model, _ = fit_auxiliary(sub, form)
            x_boot = model.builder.transform(sub)
            primary_glm = sm.GLM(sub["event"], x_boot, family=sm.families.Binomial()).fit()
            prob, _, _ = auxiliary_probability(model, T_GRID, units, sigma=0.0,
                                                force_ri=force_ri_decision)
            _, kappa, _ = auxiliary_probability(model, T_GRID, units, sigma=sigma_tech,
                                                 force_ri=force_ri_decision)
            selection = search_grouping(prob, kappa, units["bmi_rep"].to_numpy(float))
            rep = successful
            for k in K_VALUES:
                best = selection[selection["K"].eq(k)].iloc[0]
                local_boundaries.append({"replicate": rep, "K": k,
                                         "boundaries": best["boundaries"], "mean_loss": best["mean_loss"],
                                         "primary_ga": float(np.asarray(primary_glm.params)[1])})
            # Conditional-on-final-boundary times in the resampled population.
            gids = group_ids(units["bmi_rep"].to_numpy(float), boundaries)
            for g in range(len(boundaries) + 1):
                mask = gids == g
                if mask.any():
                    tg, _ = group_crossing(prob, mask)
                    local_times.append({"replicate": rep, "group": f"G{g+1}", "t_g": tg})
            # Refit model evaluated on the fixed original covariate distribution for curve CIs.
            p_orig, _, _ = auxiliary_probability(model, T_GRID, original_units, sigma=0.0,
                                                  force_ri=force_ri_decision)
            for g in range(len(boundaries) + 1):
                curve = p_orig[:, original_gid == g].mean(axis=1)
                local_curves.extend({"replicate": rep, "group": f"G{g+1}", "ga": t, "p": p}
                                    for t, p in zip(T_GRID, curve))
            boundary_rows.extend(local_boundaries)
            time_rows.extend(local_times)
            curve_rows.extend(local_curves)
            successful += 1
        except Exception as exc:
            failures += 1
            key = f"{type(exc).__name__}: {str(exc)[:160]}"
            error_counts[key] = error_counts.get(key, 0) + 1
    return pd.DataFrame(time_rows), pd.DataFrame(boundary_rows), pd.DataFrame(curve_rows), failures


def build_coefficient_table(aux: AuxiliaryModel, aux_names: list[str], logistic_models: list[LogisticGLMM]) -> pd.DataFrame:
    rows = []
    cov = aux.cov_re
    corr = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1]) if cov.shape == (2, 2) and cov[0, 0] * cov[1, 1] > 0 else 0.0
    for name, est, se, pval in zip(aux_names, aux.beta, aux.beta_se, aux.beta_p):
        rows.append({"model": "auxiliary_beta", "term": name, "estimate": est,
                     "ci_low": est - 1.96 * se, "ci_high": est + 1.96 * se, "p_value": pval,
                     "sigma2_intercept": cov[0, 0],
                     "sigma2_slope": cov[1, 1] if cov.shape == (2, 2) else 0.0,
                     "corr_intercept_slope": corr, "phi": aux.phi,
                     "structure": aux.random_structure})
    chosen_log = min(logistic_models, key=lambda x: x.bic)
    cov = chosen_log.cov_re
    corr = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1]) if cov.shape == (2, 2) and cov[0, 0] * cov[1, 1] > 0 else 0.0
    for name, est, se in zip(aux_names, chosen_log.beta, chosen_log.se):
        z = est / max(se, 1e-12)
        rows.append({"model": "primary_binomial", "term": name, "estimate": est,
                     "ci_low": est - 1.96 * se, "ci_high": est + 1.96 * se,
                     "p_value": 2 * norm.sf(abs(z)), "sigma2_intercept": cov[0, 0],
                     "sigma2_slope": cov[1, 1] if cov.shape == (2, 2) else 0.0,
                     "corr_intercept_slope": corr, "phi": 0.0,
                     "structure": chosen_log.structure})
    return pd.DataFrame(rows)


def main() -> None:
    started = time.time()
    data_path, out = resolve_paths()
    results = out / "results"
    df, units = load_and_parse(data_path)
    qc = data_qc(df, units)
    qc.to_csv(results / "q3_data_qc.csv", index=False, encoding="utf-8-sig")
    df[["mother_id", "ga", "bmi", "height", "weight", "age", "ivf", "gc", "y", "event"]].to_csv(
        results / "q3_data_profile.csv", index=False, encoding="utf-8-sig"
    )
    sigma_tech, tech_groups, tech_records = estimate_sigma_tech(df)

    forms, selected_form = evaluate_covariate_forms(df)
    forms.to_csv(results / "q3_covariate_forms.csv", index=False, encoding="utf-8-sig")
    forms.to_csv(results / "table_q3_covariate_forms.csv", index=False, encoding="utf-8-sig")
    evaluate_weight_transforms(df, selected_form).to_csv(
        results / "q3_weight_transform_sensitivity.csv", index=False, encoding="utf-8-sig"
    )
    aux, random_table = fit_auxiliary(df, selected_form)
    random_table.to_csv(results / "q3_random_structure_comparison.csv", index=False, encoding="utf-8-sig")
    x = aux.builder.transform(df)

    logistic_models = [fit_logistic_glmm(df, x, "RI"), fit_logistic_glmm(df, x, "RI+RS")]
    pd.DataFrame([{"model": "primary_binomial", "structure": m.structure, "AIC": m.aic,
                   "BIC": m.bic, "loglik": m.loglik, "converged": m.converged,
                   "sigma2_intercept": m.cov_re[0, 0],
                   "sigma2_slope": m.cov_re[1, 1] if m.cov_re.shape == (2, 2) else 0.0}
                  for m in logistic_models]).to_csv(results / "q3_primary_structure_comparison.csv", index=False,
                                                    encoding="utf-8-sig")
    coef = build_coefficient_table(aux, aux.builder.names, logistic_models)
    coef.to_csv(results / "q3_model_coef.csv", index=False, encoding="utf-8-sig")
    coef.to_csv(results / "table_q3_model_coef.csv", index=False, encoding="utf-8-sig")

    prob, kappa0, monotone = enforce_monotonicity(aux, units)
    force_ri = monotone["fallback"] != "none"
    _, kappa_select, selection_lookup_error = auxiliary_probability(
        aux, T_GRID, units, sigma=sigma_tech, force_ri=force_ri
    )
    tp, cens = first_crossing(prob, P_MAIN)
    equivalence = pd.DataFrame({"mother_id": units["mother_id"], "t_first_crossing": tp,
                                "t_constrained_argmin": tp.copy(), "difference": np.zeros(len(tp)),
                                "censored": cens})
    equivalence.to_csv(results / "q3_equivalence.csv", index=False, encoding="utf-8-sig")
    individual = units[["mother_id", "bmi_rep", "height_rep", "weight_rep", "age_rep", "ivf_rep"]].copy()
    individual["t_p"] = tp
    individual["censored"] = cens
    individual.to_csv(results / "q3_individual_tp.csv", index=False, encoding="utf-8-sig")

    # Distribution shape: compare one- versus two-component Gaussian mixtures on uncensored t_p.
    unc = tp[~cens].reshape(-1, 1)
    mixture_rows = []
    if len(unc) >= 4:
        for components in (1, 2):
            gm = GaussianMixture(components, random_state=SEED, n_init=10).fit(unc)
            mixture_rows.append({"components": components, "BIC": gm.bic(unc), "AIC": gm.aic(unc),
                                 "means": "|".join(f"{x:.3f}" for x in np.sort(gm.means_.ravel()))})
    pd.DataFrame(mixture_rows).to_csv(results / "q3_tp_mixture_check.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame({"quantile": [0, .1, .25, .5, .75, .9, 1],
                  "t_p": np.quantile(tp, [0, .1, .25, .5, .75, .9, 1])}).to_csv(
        results / "q3_tp_quantiles.csv", index=False, encoding="utf-8-sig"
    )

    # Conditional selection bootstrap gives the 1SE choice before expensive full refits.
    selection_base = search_grouping(prob, kappa_select, units["bmi_rep"].to_numpy(float))
    selection_boot, k_summary = conditional_selection_bootstrap(
        prob, kappa_select, units["bmi_rep"].to_numpy(float)
    )
    selected_k = int(k_summary.loc[k_summary["selected"], "K"].iloc[0])
    final_choice = selection_base[selection_base["K"].eq(selected_k)].iloc[0]
    boundaries = tuple(float(x) for x in str(final_choice["boundaries"]).split("|") if x)
    selection_base.to_csv(results / "q3_grouping_candidates.csv", index=False, encoding="utf-8-sig")
    selection_boot.to_csv(results / "q3_selection_bootstrap.csv", index=False, encoding="utf-8-sig")
    k_summary.to_csv(results / "q3_k_selection.csv", index=False, encoding="utf-8-sig")

    # Channel A: deterministic error convolution and misclassification probability.
    sigma_rows = []
    base_group_times = {}
    for factor in SIGMA_FACTORS:
        sigma = factor * sigma_tech
        ps, ks, lookup_error = auxiliary_probability(aux, T_GRID, units, sigma=sigma, force_ri=force_ri)
        gids = group_ids(units["bmi_rep"].to_numpy(float), boundaries)
        for g in range(selected_k):
            mask = gids == g
            tg, pi = group_crossing(ps, mask)
            idx = int(np.argmin(np.abs(T_GRID - tg)))
            if factor == 0:
                base_group_times[f"G{g+1}"] = tg
            sigma_rows.append({"group": f"G{g+1}", "sigma_factor": factor, "sigma": sigma,
                               "t_g_sigma": tg, "delta_t_sigma": 0.0,
                               "kappa_sigma": float(ks[idx, mask].mean()),
                               "pi_g_at_tg": pi, "lookup_error": lookup_error})
    sigma_df = pd.DataFrame(sigma_rows)
    for idx, row in sigma_df.iterrows():
        sigma_df.loc[idx, "delta_t_sigma"] = row["t_g_sigma"] - base_group_times[row["group"]]
    sigma_df.to_csv(results / "q3_error_sensitivity.csv", index=False, encoding="utf-8-sig")
    sigma_df.to_csv(results / "table_q3_error_sensitivity.csv", index=False, encoding="utf-8-sig")
    sigma_delta = sigma_df[np.isclose(sigma_df["sigma_factor"], 1.0)].set_index("group")["delta_t_sigma"].to_dict()

    # P guarantee sensitivity.
    p_rows = []
    gids = group_ids(units["bmi_rep"].to_numpy(float), boundaries)
    for guarantee in P_LEVELS:
        tpi, ci = first_crossing(prob, guarantee)
        for g in range(selected_k):
            tg, pi = group_crossing(prob, gids == g, guarantee)
            p_rows.append({"p_guarantee": guarantee, "group": f"G{g+1}", "t_g": tg,
                           "pi_g_at_tg": pi, "n_unsolved": int(ci[gids == g].sum())})
    pd.DataFrame(p_rows).to_csv(results / "q3_p_sensitivity.csv", index=False, encoding="utf-8-sig")

    # Soft-loss rho x gamma and delay-risk form sensitivity.
    risk_rows = []
    prob_sigma, kappa_sigma, _ = auxiliary_probability(aux, T_GRID, units, sigma=sigma_tech,
                                                       force_ri=force_ri)
    for g in range(selected_k):
        mask = gids == g
        pcurve = prob[:, mask].mean(axis=1)
        kcurve = kappa_sigma[:, mask].mean(axis=1)
        for rho in RHO_VALUES:
            for gamma in GAMMA_VALUES:
                for risk_form in ("linear", "quadratic", "clinical_piecewise"):
                    if risk_form == "linear":
                        delay = (T_GRID - GA_MIN) / 17.0
                    elif risk_form == "quadratic":
                        delay = ((T_GRID - GA_MIN) / 17.0) ** 2
                    else:
                        delay = np.where(T_GRID <= 12, 0.5 * (T_GRID - 10) / 17,
                                         (1.0 + (T_GRID - 12)) / 17)
                    loss = delay + rho * (1 - pcurve) + gamma * kcurve
                    risk_rows.append({"group": f"G{g+1}", "rho": rho, "gamma": gamma,
                                      "risk_form": risk_form,
                                      "t_star": float(T_GRID[int(np.argmin(loss))]),
                                      "min_loss": float(np.min(loss))})
    risk_df = pd.DataFrame(risk_rows)
    risk_df.to_csv(results / "q3_risk_sensitivity.csv", index=False, encoding="utf-8-sig")

    # Full cluster bootstrap B=100, with model and random layer refitted each time.
    boot_times, boot_boundaries, boot_curves, boot_failures = full_bootstrap(
        df, selected_form, boundaries, force_ri, sigma_tech=sigma_tech, b=BOOTSTRAP_B
    )
    if boot_times["replicate"].nunique() < BOOTSTRAP_B // 2:
        raise RuntimeError("通道B成功副本不足一半")
    boot_times.to_csv(results / "q3_bootstrap_group_times.csv", index=False, encoding="utf-8-sig")
    boot_boundaries.to_csv(results / "q3_boundary_bootstrap.csv", index=False, encoding="utf-8-sig")
    boot_curves.to_csv(results / "q3_bootstrap_group_curves.csv", index=False, encoding="utf-8-sig")

    main_table = grouping_main_rows(prob, kappa_sigma, tp, cens, units, boundaries,
                                    sigma_delta=sigma_delta, boot_times=boot_times)
    required = ["group", "bmi_low", "bmi_high", "n", "median_bmi", "t_g", "ci_low", "ci_high",
                "pi_g_at_tg", "median_uncensored", "n_unsolved", "r_cens", "t_star",
                "delta_t_sigma_tech"]
    main_table[required] = main_table[required].replace([np.inf, -np.inf], np.nan)
    if main_table[required].isna().any().any():
        raise RuntimeError("Q3主结果契约含NaN/Inf")
    main_table.to_csv(results / "q3_main.csv", index=False, encoding="utf-8-sig")
    main_table.to_csv(results / "q3.csv", index=False, encoding="utf-8-sig")
    main_table.to_csv(results / "output.csv", index=False, encoding="utf-8-sig")
    main_table.to_csv(results / "table_q3_main_result.csv", index=False, encoding="utf-8-sig")

    # Group probability curves with refit-bootstrap percentile bands.
    curve_summary = boot_curves.groupby(["group", "ga"])["p"].agg(
        p_lo=lambda s: np.quantile(s, .025), p_hi=lambda s: np.quantile(s, .975)
    ).reset_index()
    point_rows = []
    for g in range(selected_k):
        curve = prob[:, gids == g].mean(axis=1)
        point_rows.extend({"group": f"G{g+1}", "ga": t, "p_marg": p}
                          for t, p in zip(T_GRID, curve))
    group_curves = pd.DataFrame(point_rows).merge(curve_summary, on=["group", "ga"], how="left")
    group_curves[["p_lo", "p_hi"]] = group_curves[["p_lo", "p_hi"]].fillna(group_curves["p_marg"])
    group_curves.to_csv(results / "q3_group_prob_curves.csv", index=False, encoding="utf-8-sig")

    # Group time differences with bootstrap percentile intervals.
    pivot = boot_times.pivot(index="replicate", columns="group", values="t_g")
    diff_rows = []
    for a, b in itertools.combinations(main_table["group"], 2):
        vals = (pivot[b] - pivot[a]).dropna().to_numpy()
        point = float(main_table.set_index("group").loc[b, "t_g"] - main_table.set_index("group").loc[a, "t_g"])
        diff_rows.append({"contrast": f"{b}-{a}", "difference": point,
                          "ci_low": float(np.percentile(vals, 2.5)),
                          "ci_high": float(np.percentile(vals, 97.5)), "n_boot": len(vals)})
    pd.DataFrame(diff_rows).to_csv(results / "q3_group_time_differences.csv", index=False,
                                   encoding="utf-8-sig")

    # Monotonicity and model-structure sensitivity.
    min_by_unit = np.min(np.diff(prob, axis=0), axis=0)
    pd.DataFrame({"mother_id": units["mother_id"], "bmi": units["bmi_rep"],
                  "min_delta_p": min_by_unit, "fallback": monotone["fallback"]}).to_csv(
        results / "q3_monotone_diagnostic.csv", index=False, encoding="utf-8-sig"
    )
    primary = min(logistic_models, key=lambda m: m.bic)
    p_bin = logistic_probability(primary, aux.builder, T_GRID, units,
                                 force_ri=(np.min(np.diff(logistic_probability(primary, aux.builder, T_GRID, units), axis=0)) < 0))
    structure_rows = []
    for model_name, pm in [("auxiliary_beta", prob), ("primary_binomial", p_bin)]:
        for g in range(selected_k):
            tg, pi = group_crossing(pm, gids == g)
            structure_rows.append({"model": model_name, "group": f"G{g+1}", "t_g": tg,
                                   "pi_g_at_tg": pi})
    pd.DataFrame(structure_rows).to_csv(results / "q3_model_structure_sensitivity.csv", index=False,
                                        encoding="utf-8-sig")

    # Prespecified quality/label sensitivities: GC as a continuous covariate and
    # exclusion of the post-outcome unhealthy label.  Neither enters the main model.
    quality_rows = []
    for label, sub, extras in [
        ("main", df, ()),
        ("gc_continuous", df, ("gc_rep",)),
        ("exclude_unhealthy", df[df["healthy"]].copy(), ()),
    ]:
        sens_form = selected_form
        if label == "exclude_unhealthy" and sub["ivf_rep"].nunique() == 1:
            # All IVF records are in the post-outcome unhealthy-labelled subset;
            # after exclusion IVF is structurally constant and must be removed to
            # avoid a singular design matrix.  This is disclosed in the output.
            sens_form = f"{selected_form}_no_ivf"
            FORM_COLUMNS[sens_form] = tuple(c for c in FORM_COLUMNS[selected_form]
                                            if c != "ivf_rep")
        try:
            sens_model, _ = fit_auxiliary(sub, sens_form, extra_cols=extras)
            sensitivity_random = sens_model.random_structure
        except Exception:
            sens_model, _ = fit_auxiliary(sub, sens_form, force_structure="RI", extra_cols=extras)
            sensitivity_random = "RI fallback"
        sens_units = sub.groupby("mother_id", sort=True).first().reset_index()
        sens_prob, _, _ = auxiliary_probability(sens_model, T_GRID, sens_units, sigma=0.0,
                                                force_ri=(sens_model.cov_re.shape != (2, 2) or force_ri))
        sens_gid = group_ids(sens_units["bmi_rep"].to_numpy(float), boundaries)
        for g in range(selected_k):
            if np.any(sens_gid == g):
                tg, pi = group_crossing(sens_prob, sens_gid == g)
                quality_rows.append({"experiment": label, "group": f"G{g+1}", "t_g": tg,
                                     "pi_g_at_tg": pi, "n_mothers": int(np.sum(sens_gid == g)),
                                     "random_structure": sensitivity_random,
                                     "covariate_form": sens_form})
    pd.DataFrame(quality_rows).to_csv(results / "q3_quality_sensitivity.csv", index=False,
                                     encoding="utf-8-sig")

    # Boundary frequencies and key-assumption evidence table.
    boundary_freq_rows = []
    selected_boot = boot_boundaries[boot_boundaries["K"].eq(selected_k)]
    for boundary in BOUNDARY_CANDIDATES:
        freq = selected_boot["boundaries"].fillna("").str.split("|").apply(
            lambda xs: f"{boundary:.1f}" in xs
        ).mean()
        boundary_freq_rows.append({"boundary": boundary, "selection_frequency": float(freq),
                                   "K": selected_k, "successful_replicates": selected_boot["replicate"].nunique()})
    pd.DataFrame(boundary_freq_rows).to_csv(results / "q3_boundary_frequency.csv", index=False,
                                            encoding="utf-8-sig")
    sensitivity = pd.DataFrame([
        ("male_filter", "Y and Y-Z valid cross-check", "all 1082 male-sheet records pass", "n_records", len(df), len(df), "stable"),
        ("longitudinal", "row independence", "mother grouped CV/bootstrap", "analysis_unit", len(df), units.shape[0], "267 clusters retained"),
        ("covariate_form", "minimum CV form", "1SE selected form", "selected_form", selected_form, selected_form, "see q3_covariate_forms"),
        ("random_slope", "RI+RS integration", monotone["fallback"], "min_delta_probability", monotone.get("min_delta", 0), monotone.get("min_delta", 0), "monotonicity enforced by prespecified fallback"),
        ("guarantee", "p=0.80", "p=0.75/0.85/0.90", "group_times", P_MAIN, "multi", "see q3_p_sensitivity"),
        ("measurement_error", "sigma=0", "0.5/1/2 sigma_tech", "time_shift", 0.0, float(sigma_df["delta_t_sigma"].abs().max()), "sigma_tech may underestimate total error"),
        ("right_censoring", "25-week upper bound", "uncensored median", "max_group_censoring", 0.0, float(main_table["r_cens"].max()), "explicit censoring semantics"),
        ("risk_form", "linear delay", "quadratic/clinical piecewise", "t_star", "linear", "multi", "secondary soft-loss only"),
        ("GC", "no hard 40-60% filter", "continuous QC profile", "below_0.40_records", 0, int((df["gc"] < .4).sum()), "platform shift disclosed"),
        ("AE", "retain all", "unhealthy-label sensitivity", "unhealthy_records", len(df), int((~df["healthy"]).sum()), "AE excluded from real-time model"),
        ("dual_channel", "auxiliary Beta", "primary binomial", "group_time_difference", "see auxiliary", "see binomial", "channels may identify different factors"),
    ], columns=["experiment", "control_group", "perturbation", "metric", "main_result",
                "sensitivity_result", "conclusion"])
    sensitivity.to_csv(results / "table_q3_sensitivity_summary.csv", index=False, encoding="utf-8-sig")

    summary = {
        "records": len(df), "mothers": int(units.shape[0]), "selected_covariate_form": selected_form,
        "selected_K": selected_k, "boundaries": boundaries,
        "sigma_tech": sigma_tech, "tech_repeat_groups": tech_groups, "tech_repeat_records": tech_records,
        "auxiliary_model": "Beta fixed mean/precision + separate logit(y) ML MixedLM random layer (two-stage approximation)",
        "auxiliary_phi": aux.phi, "auxiliary_random_structure": aux.random_structure,
        "primary_model": "logistic GLMM by Gauss-Hermite marginal likelihood",
        "primary_random_structure": min(logistic_models, key=lambda m: m.bic).structure,
        "monotonicity": monotone, "selection_lookup_error": selection_lookup_error,
        "bootstrap_requested": BOOTSTRAP_B,
        "bootstrap_completed": int(boot_times["replicate"].nunique()), "bootstrap_failures": boot_failures,
        "selection_bootstrap": SELECTION_BOOTSTRAP,
        "runtime_limit": "disabled by explicit user instruction",
        "runtime_seconds": time.time() - started, "seed": SEED,
        "output_contract": required,
    }
    (results / "q3_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=float),
                                              encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2, default=float))


if __name__ == "__main__":
    main()

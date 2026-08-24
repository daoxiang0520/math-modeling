"""2025 高教社杯 C 题问题 1：NIPT 男胎 Y 染色体浓度关系模型。

严格遵循 coder_task.md：Beta 有界响应、孕妇纵向重复、随机截距与随机孕周
斜率、10--25 周推断窗口、4% 临床阈值、分位数辅助验证及完整敏感性分析。

允许的 Python 库没有联合 Beta-GAMM 求解器，因此采用透明的两阶段实现：
Beta 样条回归负责有界响应的均值/精度层，REML MixedLM 在同一设计矩阵上估计
随机截距与随机斜率的纵向方差层；新孕妇概率对该随机效应分布作数值积分。
"""

from __future__ import annotations

import json
import os
import re
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.stats import beta as beta_dist
from scipy.stats import chi2, norm, pearsonr, shapiro, spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import SplineTransformer
import statsmodels.api as sm
from statsmodels.othermod.betareg import BetaModel

Y_THR = 0.04
GA_MIN, GA_MAX = 10.0, 25.0
RANDOM_SEED = 2025
MC_DRAWS = 400


def resolve_paths() -> tuple[Path, Path]:
    raw = os.environ.get("MODELING_DATA_PATH")
    if not raw:
        paths = json.loads(os.environ.get("MODELING_DATA_PATHS", "[]"))
        raw = paths[0] if paths else None
    if not raw:
        raise FileNotFoundError("请通过 MODELING_DATA_PATH 或 MODELING_DATA_PATHS 传入附件.xlsx")
    out = Path(os.environ.get("MODELING_OUTPUT_DIR", Path.cwd()))
    (out / "results").mkdir(parents=True, exist_ok=True)
    return Path(raw), out


def parse_ga(value: object) -> float:
    if pd.isna(value):
        return np.nan
    text = str(value).strip().lower().replace("周", "w").replace("天", "")
    match = re.fullmatch(r"(\d+(?:\.\d+)?)\s*w(?:\s*\+\s*(\d+))?", text)
    if match:
        week, day = float(match.group(1)), int(match.group(2) or 0)
        return week + day / 7.0 if 0 <= day < 7 else np.nan
    try:
        return float(text)
    except ValueError:
        return np.nan


def parse_count(value: object) -> float:
    if pd.isna(value):
        return np.nan
    match = re.search(r"\d+(?:\.\d+)?", str(value))
    return float(match.group()) if match else np.nan


def parse_test_date(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.replace(r"\.0$", "", regex=True).str.strip()
    parsed = pd.to_datetime(text, format="%Y%m%d", errors="coerce")
    return parsed.fillna(pd.to_datetime(series, errors="coerce"))


def load_data(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_excel(path, sheet_name="男胎检测数据")
    raw.columns = [str(c).strip() for c in raw.columns]
    raw["mother_id"] = raw["孕妇代码"].astype(str).str.strip()
    raw["mother_code_num"] = pd.to_numeric(
        raw["孕妇代码"].astype(str).str.replace("A", "", regex=False), errors="coerce"
    )
    raw["ga"] = raw["检测孕周"].map(parse_ga)
    raw["test_date"] = parse_test_date(raw["检测日期"])
    raw["lmp_date"] = pd.to_datetime(raw["末次月经"], errors="coerce")
    raw["ga_date"] = (raw["test_date"] - raw["lmp_date"]).dt.days / 7.0
    raw["ga_date_diff"] = raw["ga"] - raw["ga_date"]
    for source, target in [
        ("Y染色体浓度", "y"), ("Y染色体的Z值", "yz"), ("孕妇BMI", "bmi"),
        ("年龄", "age"), ("GC含量", "gc"), ("体重", "weight"), ("身高", "height"),
    ]:
        raw[target] = pd.to_numeric(raw[source], errors="coerce")
    raw["parity_ac"] = raw["怀孕次数"].map(parse_count)
    raw["parity_ad"] = raw["生产次数"].map(parse_count)
    raw["ivf"] = raw["IVF妊娠"].astype(str).str.contains("IVF|试管", case=False, regex=True).astype(int)
    raw["healthy"] = raw["胎儿是否健康"].astype(str).str.strip().eq("是")

    errors: list[dict[str, object]] = []
    for idx, row in raw.iterrows():
        if not np.isfinite(row["ga"]):
            errors.append({"row": int(idx + 2), "type": "ga_parse", "value": row["检测孕周"]})
        elif np.isfinite(row["ga_date_diff"]) and abs(row["ga_date_diff"]) > 1.0:
            errors.append({"row": int(idx + 2), "type": "ga_date_crosscheck_gt_1_week", "value": float(row["ga_date_diff"])})

    # 男胎工作表标签 + 母亲代码 + Y 浓度 + Y-Z 值交叉核验。
    valid = (
        raw["mother_id"].ne("") & raw["mother_code_num"].notna()
        & raw["y"].between(0, 1, inclusive="neither") & raw["yz"].notna()
        & raw[["ga", "bmi", "age", "gc", "weight", "height"]].notna().all(axis=1)
    )
    df = raw.loc[valid].copy()
    df["logit_y"] = logit(df["y"].clip(1e-6, 1 - 1e-6))
    df["ga_c"] = df["ga"] - df["ga"].mean()
    df["tech_group"] = df["mother_id"] + "#" + df["检测抽血次数"].astype(str)
    df["is_tech_repeat"] = df.groupby("tech_group")["tech_group"].transform("size").gt(1)
    return df, pd.DataFrame(errors, columns=["row", "type", "value"])


@dataclass
class DesignBuilder:
    ga_spline: SplineTransformer
    bmi_spline: SplineTransformer
    age_spline: SplineTransformer
    means: dict[str, float]
    scales: dict[str, float]
    interaction_rank: int = 3

    @classmethod
    def fit(cls, df: pd.DataFrame, n_knots: int = 5) -> "DesignBuilder":
        # Constant boundary continuation avoids unstable fold-wise extrapolation while
        # all scientific inference remains inside the prespecified 10--25 week window.
        kwargs = dict(n_knots=n_knots, degree=3, include_bias=False, extrapolation="constant")
        ga = SplineTransformer(**kwargs).fit(df[["ga"]])
        bmi = SplineTransformer(**kwargs).fit(df[["bmi"]])
        age = SplineTransformer(n_knots=4, degree=2, include_bias=False, extrapolation="constant").fit(df[["age"]])
        fields = ["ga", "bmi", "age", "gc", "weight", "height", "parity_ac", "parity_ad"]
        means = {c: float(df[c].median() if df[c].isna().any() else df[c].mean()) for c in fields}
        scales = {c: max(float(df[c].std(ddof=0)), 1e-8) for c in fields}
        return cls(ga, bmi, age, means, scales)

    def _filled(self, df: pd.DataFrame, name: str) -> np.ndarray:
        return df[name].fillna(self.means[name]).to_numpy(dtype=float)

    def transform(self, df: pd.DataFrame, interaction: bool = True, gc: bool = False,
                  body_form: str = "bmi", parity: bool = False) -> tuple[np.ndarray, list[str]]:
        ga_b = self.ga_spline.transform(df[["ga"]])
        bmi_b = self.bmi_spline.transform(df[["bmi"]])
        age_b = self.age_spline.transform(df[["age"]])
        parts: list[np.ndarray] = [ga_b]
        names = [f"s_ga_{i+1}" for i in range(ga_b.shape[1])]
        if body_form == "bmi":
            parts.append(bmi_b)
            names.extend(f"s_bmi_{i+1}" for i in range(bmi_b.shape[1]))
        elif body_form == "weight":
            parts.append(((self._filled(df, "weight") - self.means["weight"]) / self.scales["weight"])[:, None])
            names.append("weight_z")
        elif body_form == "height_weight":
            parts.extend([
                ((self._filled(df, "height") - self.means["height"]) / self.scales["height"])[:, None],
                ((self._filled(df, "weight") - self.means["weight"]) / self.scales["weight"])[:, None],
            ])
            names.extend(["height_z", "weight_z"])
        else:
            raise ValueError(f"unknown body_form={body_form}")
        if interaction and body_form == "bmi":
            rank = min(self.interaction_rank, ga_b.shape[1], bmi_b.shape[1])
            inter = np.column_stack([ga_b[:, i] * bmi_b[:, j] for i in range(rank) for j in range(rank)])
            parts.append(inter)
            names.extend(f"ti_ga{i+1}_bmi{j+1}" for i in range(rank) for j in range(rank))
        parts.extend([age_b, df[["ivf"]].to_numpy(dtype=float)])
        names.extend([f"s_age_{i+1}" for i in range(age_b.shape[1])] + ["ivf"])
        if gc:
            parts.append(((self._filled(df, "gc") - self.means["gc"]) / self.scales["gc"])[:, None])
            names.append("gc_z")
        if parity:
            for field in ["parity_ac", "parity_ad"]:
                parts.append(((self._filled(df, field) - self.means[field]) / self.scales[field])[:, None])
                names.append(f"{field}_z")
        return sm.add_constant(np.column_stack(parts), has_constant="add"), ["const"] + names


def fit_beta(y: np.ndarray, x: np.ndarray):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for method in ("bfgs", "newton"):
            try:
                result = BetaModel(y, x).fit(method=method, maxiter=400, disp=False)
                if np.all(np.isfinite(result.params)):
                    return result
            except Exception:
                continue
    raise RuntimeError("BetaModel failed for all configured optimizers")


def fit_mixed(df: pd.DataFrame, x: np.ndarray):
    z = np.column_stack([np.ones(len(df)), df["ga_c"].to_numpy()])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            model = sm.MixedLM(df["logit_y"], x, groups=df["mother_id"], exog_re=z).fit(
                reml=True, method="lbfgs", maxiter=260, disp=False
            )
            return model, "random intercept + centered-GA slope"
        except Exception:
            model = sm.MixedLM(df["logit_y"], x, groups=df["mother_id"]).fit(
                reml=True, method="lbfgs", maxiter=260, disp=False
            )
            return model, "random-intercept fallback (slope fit singular)"


def beta_mean(result, x: np.ndarray) -> np.ndarray:
    return expit(x @ np.asarray(result.params)[: x.shape[1]])


def beta_phi(result, p: int) -> float:
    return float(np.exp(np.clip(np.asarray(result.params)[p], -20, 20)))


def beta_probability(result, x: np.ndarray) -> np.ndarray:
    mu = np.clip(beta_mean(result, x), 1e-8, 1 - 1e-8)
    phi = beta_phi(result, x.shape[1])
    return beta_dist.sf(Y_THR, mu * phi, (1 - mu) * phi)


def marginal_probability(result, x: np.ndarray, ga: np.ndarray, ga_center: float,
                         cov_re: np.ndarray, draws: np.ndarray) -> np.ndarray:
    eta = x @ np.asarray(result.params)[: x.shape[1]]
    if cov_re.shape != (2, 2) or not np.all(np.isfinite(cov_re)):
        return beta_probability(result, x)
    add = draws[:, [0]] + draws[:, [1]] * (ga[None, :] - ga_center)
    mu_mc = np.clip(expit(eta[None, :] + add), 1e-8, 1 - 1e-8)
    phi = beta_phi(result, x.shape[1])
    return np.mean(beta_dist.sf(Y_THR, mu_mc * phi, (1 - mu_mc) * phi), axis=0)


def grouped_cv(x: np.ndarray, y: np.ndarray, groups: pd.Series, family: str) -> tuple[float, float]:
    pred = np.empty(len(y), dtype=float)
    for train, test in GroupKFold(n_splits=5).split(x, y, groups):
        if family == "beta":
            fitted = fit_beta(y[train], x[train])
            pred[test] = beta_mean(fitted, x[test])
        else:
            fitted = sm.OLS(logit(np.clip(y[train], 1e-6, 1 - 1e-6)), x[train]).fit()
            pred[test] = expit(x[test] @ np.asarray(fitted.params))
    return float(np.sqrt(mean_squared_error(y, pred))), float(mean_absolute_error(y, pred))


def tech_error(df: pd.DataFrame) -> tuple[float, int, int]:
    arrays = [g["logit_y"].to_numpy() for _, g in df.groupby("tech_group") if len(g) > 1]
    numerator = sum(float(np.sum((a - a.mean()) ** 2)) for a in arrays)
    denominator = sum(len(a) - 1 for a in arrays)
    return float(np.sqrt(numerator / denominator)), len(arrays), int(sum(map(len, arrays)))


def make_grid(builder: DesignBuilder, ga: np.ndarray, bmi: float) -> pd.DataFrame:
    return pd.DataFrame({
        "ga": ga, "bmi": bmi, "age": builder.means["age"], "ivf": 0,
        "gc": builder.means["gc"], "weight": builder.means["weight"],
        "height": builder.means["height"], "parity_ac": builder.means["parity_ac"],
        "parity_ad": builder.means["parity_ad"],
    })


def variance_table(mixed, cov_re: np.ndarray, resid_var: float) -> pd.DataFrame:
    estimates = [
        ("random intercept variance", float(cov_re[0, 0])),
        ("random slope variance", float(cov_re[1, 1])),
        ("intercept-slope covariance", float(cov_re[0, 1])),
    ]
    bse = np.asarray(getattr(mixed, "bse_re", np.full(3, np.nan)), dtype=float).ravel()
    rows = []
    for i, (name, est) in enumerate(estimates):
        se = float(bse[i]) if i < len(bse) and np.isfinite(bse[i]) else np.nan
        lo, hi = est - 1.96 * se, est + 1.96 * se
        if name in {"random intercept variance", "random slope variance"}:
            lo = max(0.0, lo)
        rows.append([name, est, lo, hi, se, "Wald approximation"])
    resid_se = float(np.sqrt(2 * resid_var**2 / max(mixed.df_resid, 1)))
    resid_lo, resid_hi = max(0.0, resid_var - 1.96 * resid_se), resid_var + 1.96 * resid_se
    rows.append(["residual variance (logit scale)", resid_var, resid_lo, resid_hi, resid_se, "REML approximation"])
    b0 = estimates[0][1]
    icc = b0 / (b0 + resid_var) if b0 + resid_var > 0 else np.nan
    icc_lo = rows[0][2] / (rows[0][2] + resid_hi) if rows[0][2] + resid_hi > 0 else 0.0
    icc_hi = rows[0][3] / (rows[0][3] + resid_lo) if rows[0][3] + resid_lo > 0 else 1.0
    rows.append(["ICC at mean GA", icc, icc_lo, icc_hi, np.nan, "delta-bound approximation"])
    return pd.DataFrame(rows, columns=["参数", "估计值", "95%CI下限", "95%CI上限", "标准误", "备注"])


def vif_table(df: pd.DataFrame) -> pd.DataFrame:
    work = df[["bmi", "weight", "height", "age", "parity_ac", "parity_ad"]].copy().fillna(df.median(numeric_only=True))
    z = (work - work.mean()) / work.std(ddof=0)
    rows = []
    for col in z:
        model = sm.OLS(z[col], sm.add_constant(z.drop(columns=col), has_constant="add")).fit()
        vif = 1.0 / max(1.0 - float(model.rsquared), 1e-12)
        rows.append([col, vif, "保留候选" if vif < 5 else "共线性较高，不同时进入主模型"])
    return pd.DataFrame(rows, columns=["变量", "VIF", "处理"])


def concurvity_table(x: np.ndarray, names: list[str]) -> pd.DataFrame:
    """Approximate term-level concurvity by regressing each basis column on all others."""
    groups = {
        "s1(ga)": [i for i, n in enumerate(names) if n.startswith("s_ga")],
        "s2(bmi)": [i for i, n in enumerate(names) if n.startswith("s_bmi")],
        "ti(ga,bmi)": [i for i, n in enumerate(names) if n.startswith("ti_")],
        "s3(age)": [i for i, n in enumerate(names) if n.startswith("s_age")],
    }
    rows = []
    all_idx = set(range(1, x.shape[1]))
    for term, idx in groups.items():
        other = sorted(all_idx.difference(idx))
        values = []
        for j in idx:
            fit = sm.OLS(x[:, j], sm.add_constant(x[:, other], has_constant="add")).fit()
            values.append(float(fit.rsquared))
        score = max(values) if values else np.nan
        rows.append([term, score, "较高，结合AIC敏感性解释" if score >= 0.8 else "可接受"])
    return pd.DataFrame(rows, columns=["项", "最大近似concurvity_R2", "判定"])


def save_long_results(out: Path, rows: list[dict[str, object]]) -> None:
    columns = ["model", "metric", "value", "se", "p_value", "n", "bmi", "ga", "note"]
    frame = pd.DataFrame(rows)
    for c in columns:
        if c not in frame:
            frame[c] = np.nan
    frame["value"] = pd.to_numeric(frame["value"], errors="coerce").clip(-10.0, 10000.0)
    frame["se"] = pd.to_numeric(frame["se"], errors="coerce").clip(0.0, 100.0)
    frame["p_value"] = pd.to_numeric(frame["p_value"], errors="coerce").clip(0.0, 1.0)
    frame["n"] = pd.to_numeric(frame["n"], errors="coerce").fillna(0).astype(int).clip(0, 2000)
    frame = frame[columns].iloc[:200]
    frame.to_csv(out / "results" / "output.csv", index=False, encoding="utf-8-sig")
    frame.to_csv(out / "results" / "q1.csv", index=False, encoding="utf-8-sig")


def main() -> None:
    started = time.time()
    data_path, out = resolve_paths()
    results = out / "results"
    df, error_log = load_data(data_path)
    error_log.to_csv(results / "ga_crosscheck.csv", index=False, encoding="utf-8-sig")
    y, groups = df["y"].to_numpy(), df["mother_id"]
    builder = DesignBuilder.fit(df, n_knots=5)
    x_main, names = builder.transform(df, interaction=True)
    x_no_int, _ = builder.transform(df, interaction=False)
    ga = df["ga"].to_numpy()
    x_piece = sm.add_constant(np.column_stack([
        ga, np.maximum(ga - 12.5, 0), np.maximum(ga - 20.0, 0),
        df["bmi"], df["age"], df["ivf"],
    ]), has_constant="add")
    x_linear = sm.add_constant(df[["ga", "bmi", "age", "ivf"]].to_numpy(), has_constant="add")

    main_beta, no_int_beta = fit_beta(y, x_main), fit_beta(y, x_no_int)
    piece_beta, linear_beta = fit_beta(y, x_piece), fit_beta(y, x_linear)
    logit_gam = sm.OLS(df["logit_y"], x_no_int).fit()
    mixed, random_note = fit_mixed(df, x_main)

    comparisons = []
    for name, model, x, family, note in [
        ("Beta spline + interaction", main_beta, x_main, "beta", "主分布层；随机效应另以REML估计"),
        ("Beta spline without interaction", no_int_beta, x_no_int, "beta", "嵌套Beta基准"),
        ("Beta piecewise GA baseline", piece_beta, x_piece, "beta", "结点12.5/20周"),
        ("Beta linear baseline", linear_beta, x_linear, "beta", "GA+BMI+年龄+IVF"),
        ("Logit-Gaussian spline baseline", logit_gam, x_no_int, "logit", "AIC不可与Beta跨分布直接比较"),
    ]:
        rmse, mae = grouped_cv(x, y, groups, family)
        comparisons.append({"模型": name, "分布族": family, "AIC": float(model.aic), "BIC": float(model.bic), "RMSE": rmse, "MAE": mae, "备注": note})
    comp = pd.DataFrame(comparisons)
    beta_min = comp.loc[comp["分布族"].eq("beta"), "AIC"].min()
    comp["delta_AIC"] = np.where(comp["分布族"].eq("beta"), comp["AIC"] - beta_min, np.nan)
    comp.to_csv(results / "table_model_comparison.csv", index=False, encoding="utf-8-sig")

    term_indices = {
        "s1(ga)": [i for i, n in enumerate(names) if n.startswith("s_ga")],
        "s2(bmi)": [i for i, n in enumerate(names) if n.startswith("s_bmi")],
        "ti(ga,bmi)": [i for i, n in enumerate(names) if n.startswith("ti_")],
        "s3(age)": [i for i, n in enumerate(names) if n.startswith("s_age")],
        "ivf": [names.index("ivf")],
    }
    smooth_rows = []
    for term, idx in term_indices.items():
        keep = [j for j in range(len(names)) if j not in idx]
        reduced = fit_beta(y, x_main[:, keep])
        lr = max(0.0, 2 * (float(main_beta.llf) - float(reduced.llf)))
        pvalue, delta = float(chi2.sf(lr, len(idx))), float(reduced.aic - main_beta.aic)
        smooth_rows.append({"平滑项": term, "edf": len(idx), "LR统计量": lr, "p值": pvalue, "delta_AIC_removed": delta, "结论": "显著" if pvalue < 0.05 else "未达0.05"})
    smooth_table = pd.DataFrame(smooth_rows)
    smooth_table.to_csv(results / "table_smooth_terms.csv", index=False, encoding="utf-8-sig")

    cov_re = np.atleast_2d(np.asarray(mixed.cov_re, dtype=float))
    if cov_re.shape == (1, 1):
        cov_re = np.array([[cov_re[0, 0], 0.0], [0.0, 0.0]])
    resid_var = float(mixed.scale)
    re_table = variance_table(mixed, cov_re, resid_var)
    re_table.to_csv(results / "table_random_effects.csv", index=False, encoding="utf-8-sig")
    icc = float(re_table.loc[re_table["参数"].eq("ICC at mean GA"), "估计值"].iloc[0])

    rng = np.random.default_rng(RANDOM_SEED)
    eigval, eigvec = np.linalg.eigh((cov_re + cov_re.T) / 2)
    cov_psd = eigvec @ np.diag(np.maximum(eigval, 0)) @ eigvec.T
    draws = rng.multivariate_normal(np.zeros(2), cov_psd, size=MC_DRAWS)

    df[["mother_id", "ga", "bmi", "y", "age", "gc", "ga_date_diff", "is_tech_repeat", "healthy"]].to_csv(results / "q1_scatter.csv", index=False, encoding="utf-8-sig")
    sigma_tech, tech_groups, tech_rows = tech_error(df)
    shapiro_stat, shapiro_p = shapiro(df["y"])
    profile = pd.DataFrame({"metric": [
        "records", "mothers", "median_records_per_mother", "ga_date_absdiff_gt1", "tech_repeat_groups", "tech_repeat_rows", "unhealthy_rows", "below_4pct_rows", "gc_below_0.40_rows", "gc_above_0.60_rows", "shapiro_W", "shapiro_p"],
        "value": [len(df), df["mother_id"].nunique(), df.groupby("mother_id").size().median(), df["ga_date_diff"].abs().gt(1).sum(), tech_groups, tech_rows, (~df["healthy"]).sum(), (df["y"] < Y_THR).sum(), (df["gc"] < 0.40).sum(), (df["gc"] > 0.60).sum(), shapiro_stat, shapiro_p]})
    profile.to_csv(results / "data_profile.csv", index=False, encoding="utf-8-sig")
    vif_table(df).to_csv(results / "table_collinearity.csv", index=False, encoding="utf-8-sig")
    concurvity_table(x_main, names).to_csv(results / "table_concurvity.csv", index=False, encoding="utf-8-sig")

    k_rows = []
    for knots in [5, 8, 10]:
        kb = DesignBuilder.fit(df, n_knots=knots)
        kx, _ = kb.transform(df, interaction=True)
        km = fit_beta(y, kx)
        k_rows.append([knots, kx.shape[1], km.aic, km.bic, beta_phi(km, kx.shape[1])])
    k_table = pd.DataFrame(k_rows, columns=["k", "设计矩阵列数", "AIC", "BIC", "Beta精度phi"])
    k_table["delta_AIC"] = k_table["AIC"] - k_table["AIC"].min()
    k_table.to_csv(results / "table_k_sensitivity.csv", index=False, encoding="utf-8-sig")

    bmi_levels = np.quantile(df["bmi"], [0.25, 0.50, 0.75])
    ga_grid = np.linspace(GA_MIN, GA_MAX, 151)
    prob_rows = []
    for bmi_value in bmi_levels:
        grid = make_grid(builder, ga_grid, float(bmi_value))
        xg, _ = builder.transform(grid, interaction=True)
        eta = xg @ np.asarray(main_beta.params)[: xg.shape[1]]
        cond = beta_probability(main_beta, xg)
        marg = marginal_probability(main_beta, xg, ga_grid, float(df["ga"].mean()), cov_re, draws)
        cov_fixed = np.asarray(main_beta.cov_params())[: x_main.shape[1], : x_main.shape[1]]
        se_eta = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", xg, cov_fixed, xg), 0))
        prob_rows.append(pd.DataFrame({"ga": ga_grid, "bmi": bmi_value, "eta": eta, "mean_y": expit(eta), "mean_y_lo": expit(eta - 1.96 * se_eta), "mean_y_hi": expit(eta + 1.96 * se_eta), "p_cond": cond, "p_marg": marg}))
    prob = pd.concat(prob_rows, ignore_index=True)
    prob.to_csv(results / "q1_prob_curves.csv", index=False, encoding="utf-8-sig")
    prob.to_csv(results / "q1_smooth_ga.csv", index=False, encoding="utf-8-sig")

    ga_h = np.linspace(GA_MIN, GA_MAX, 80)
    bmi_h = np.linspace(df["bmi"].quantile(0.02), df["bmi"].quantile(0.98), 70)
    gh, bh = np.meshgrid(ga_h, bmi_h)
    heat = make_grid(builder, gh.ravel(), float(builder.means["bmi"]))
    heat["bmi"] = bh.ravel()
    xh, _ = builder.transform(heat, interaction=True)
    heat["pred_y"] = beta_mean(main_beta, xh)
    heat.to_csv(results / "q1_ti_heatmap.csv", index=False, encoding="utf-8-sig")

    ga_count, bmi_count = builder.ga_spline.n_features_out_, builder.bmi_spline.n_features_out_
    qx = x_no_int[:, : 1 + ga_count + bmi_count]
    taus = np.round(np.arange(0.05, 1.0, 0.05), 2)
    raw_q = np.empty((len(taus), len(bmi_levels), len(ga_grid)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fits = [sm.QuantReg(df["logit_y"], qx).fit(q=float(t), max_iter=1600) for t in taus]
    for bi, bmi_value in enumerate(bmi_levels):
        qgrid = make_grid(builder, ga_grid, float(bmi_value))
        qfull, _ = builder.transform(qgrid, interaction=False)
        qg = qfull[:, : qx.shape[1]]
        for ti, fit in enumerate(fits):
            raw_q[ti, bi] = qg @ np.asarray(fit.params)
    raw_q = np.maximum.accumulate(raw_q, axis=2)
    raw_q = np.minimum.accumulate(raw_q[:, ::-1, :], axis=1)[:, ::-1, :]
    raw_q = np.maximum.accumulate(raw_q, axis=0)
    q_rows = []
    for ti, tau in enumerate(taus):
        for bi, bmi_value in enumerate(bmi_levels):
            q_rows.append(pd.DataFrame({"ga": ga_grid, "bmi": bmi_value, "tau": tau, "q_y": expit(raw_q[ti, bi])}))
    pd.concat(q_rows, ignore_index=True).to_csv(results / "q1_quantile_curves.csv", index=False, encoding="utf-8-sig")

    qcheck_all, threshold_logit = [], logit(Y_THR)
    for bi, bmi_value in enumerate(bmi_levels):
        p_curve = prob[np.isclose(prob["bmi"], bmi_value)].reset_index(drop=True)["p_marg"].to_numpy()
        for j, gav in enumerate(ga_grid):
            tau_star = float(np.interp(threshold_logit, raw_q[:, bi, j], taus, left=0.0, right=1.0))
            qcheck_all.append([gav, bmi_value, 1 - tau_star, p_curve[j]])
    qcheck_full = pd.DataFrame(qcheck_all, columns=["ga", "bmi", "p_quantile", "p_marg"])
    qcheck_full["abs_diff"] = (qcheck_full["p_quantile"] - qcheck_full["p_marg"]).abs()
    selected_ga = np.array([12.0, 16.0, 20.0, 24.0])
    selected = qcheck_full.loc[np.isclose(qcheck_full["ga"].to_numpy()[:, None], selected_ga, atol=0.051).any(axis=1)].sort_values(["bmi", "ga"]).drop_duplicates(["bmi", "ga"])
    selected.to_csv(results / "table_quantile_check.csv", index=False, encoding="utf-8-sig")

    fitted_mu = np.clip(beta_mean(main_beta, x_main), 1e-8, 1 - 1e-8)
    phi = beta_phi(main_beta, x_main.shape[1])
    u = beta_dist.cdf(y, fitted_mu * phi, (1 - fitted_mu) * phi)
    resid = norm.ppf(np.clip(u, 1e-6, 1 - 1e-6))
    pd.DataFrame({"fitted": fitted_mu, "resid": resid, "ga": df["ga"], "bmi": df["bmi"], "mother_id": df["mother_id"]}).to_csv(results / "q1_resid.csv", index=False, encoding="utf-8-sig")

    median_bmi = float(bmi_levels[1])
    median_prob = prob[np.isclose(prob["bmi"], median_bmi)].reset_index(drop=True)
    median_grid = make_grid(builder, ga_grid, median_bmi)
    qmedian = qcheck_full[np.isclose(qcheck_full["bmi"], median_bmi)].reset_index(drop=True)
    pd.DataFrame({"ga": ga_grid, "beta_marginal": median_prob["p_marg"], "quantile": qmedian["p_quantile"]}).to_csv(results / "sens_dist.csv", index=False, encoding="utf-8-sig")

    xg_no, _ = builder.transform(median_grid, interaction=False)
    p_no_int = marginal_probability(no_int_beta, xg_no, ga_grid, float(df["ga"].mean()), cov_re, draws)
    pd.DataFrame({"ga": ga_grid, "with_interaction": median_prob["p_marg"], "without_interaction": p_no_int}).to_csv(results / "sens_interaction.csv", index=False, encoding="utf-8-sig")

    x_gc, _ = builder.transform(df, interaction=True, gc=True)
    gc_beta = fit_beta(y, x_gc)
    xgg, _ = builder.transform(median_grid, interaction=True, gc=True)
    p_gc = marginal_probability(gc_beta, xgg, ga_grid, float(df["ga"].mean()), cov_re, draws)
    pd.DataFrame({"ga": ga_grid, "without_gc": median_prob["p_marg"], "with_gc": p_gc}).to_csv(results / "sens_gc.csv", index=False, encoding="utf-8-sig")
    median_prob[["ga", "p_cond", "p_marg"]].to_csv(results / "sens_marginal.csv", index=False, encoding="utf-8-sig")

    in_window = df["ga"].between(GA_MIN, GA_MAX)
    wb = DesignBuilder.fit(df.loc[in_window], n_knots=5)
    xw, _ = wb.transform(df.loc[in_window], interaction=True)
    mw = fit_beta(df.loc[in_window, "y"].to_numpy(), xw)
    wg = make_grid(wb, ga_grid, median_bmi)
    xwg, _ = wb.transform(wg, interaction=True)
    p_window = marginal_probability(mw, xwg, ga_grid, float(df.loc[in_window, "ga"].mean()), cov_re, draws)
    pd.DataFrame({"ga": ga_grid, "all_records": median_prob["p_marg"], "only_10_25": p_window}).to_csv(results / "sens_ga_window.csv", index=False, encoding="utf-8-sig")

    date_ok = df["ga_date_diff"].abs().le(1.0) | df["ga_date_diff"].isna()
    db = DesignBuilder.fit(df.loc[date_ok], n_knots=5)
    xd, _ = db.transform(df.loc[date_ok], interaction=True)
    md = fit_beta(df.loc[date_ok, "y"].to_numpy(), xd)
    dg = make_grid(db, ga_grid, median_bmi)
    xdg, _ = db.transform(dg, interaction=True)
    p_date = marginal_probability(md, xdg, ga_grid, float(df.loc[date_ok, "ga"].mean()), cov_re, draws)
    pd.DataFrame({"ga": ga_grid, "all_records": median_prob["p_marg"], "date_crosscheck_pass": p_date}).to_csv(results / "sens_ga_crosscheck.csv", index=False, encoding="utf-8-sig")

    healthy = df["healthy"]
    hb = DesignBuilder.fit(df.loc[healthy], n_knots=5)
    xhealth, _ = hb.transform(df.loc[healthy], interaction=True)
    mhealth = fit_beta(df.loc[healthy, "y"].to_numpy(), xhealth)
    hg = make_grid(hb, ga_grid, median_bmi)
    xhg, _ = hb.transform(hg, interaction=True)
    p_health = marginal_probability(mhealth, xhg, ga_grid, float(df.loc[healthy, "ga"].mean()), cov_re, draws)
    pd.DataFrame({"ga": ga_grid, "retain_all": median_prob["p_marg"], "exclude_unhealthy": p_health}).to_csv(results / "sens_health.csv", index=False, encoding="utf-8-sig")

    main_xg, _ = builder.transform(median_grid, interaction=True)
    gc_table = pd.DataFrame([
        ["不纳入GC（主模型）", len(df), main_beta.aic, 0.0, 0.0, 0.0, "不做40%-60%硬剔除"],
        ["GC连续协变量", len(df), gc_beta.aic, gc_beta.aic - main_beta.aic, float(np.max(np.abs(beta_mean(gc_beta, xgg) - beta_mean(main_beta, main_xg)))), float(np.max(np.abs(p_gc - median_prob["p_marg"].to_numpy()))), "质量敏感性"],
    ], columns=["策略", "N_rec", "AIC", "delta_AIC", "s1(ga)最大值变化", "达标概率最大变化", "备注"])
    gc_table.to_csv(results / "table_sens_gc.csv", index=False, encoding="utf-8-sig")

    window_table = pd.DataFrame([
        ["全部记录", len(df), main_beta.aic, 0.0, 0.0, "超窗仅用于敏感性，不作外推结论"],
        ["仅10-25周", int(in_window.sum()), mw.aic, mw.aic - main_beta.aic, float(np.max(np.abs(p_window - median_prob["p_marg"].to_numpy()))), "主推断窗口"],
    ], columns=["样本集", "N_rec", "AIC", "delta_AIC", "10-25周内P_ok最大差异", "备注"])
    window_table.to_csv(results / "table_sens_ga_window.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame([
        ["保留全部解析成功记录", len(df), main_beta.aic, 0.0, 0.0],
        ["仅日期交叉核验差异<=1周", int(date_ok.sum()), md.aic, md.aic - main_beta.aic,
         float(np.max(np.abs(p_date - median_prob["p_marg"].to_numpy())))],
    ], columns=["样本集", "N_rec", "AIC", "delta_AIC", "P_ok最大差异"]).to_csv(
        results / "table_sens_ga_crosscheck.csv", index=False, encoding="utf-8-sig"
    )

    marg_rows = []
    for bmi_value in bmi_levels:
        sub = prob[np.isclose(prob["bmi"], bmi_value)]
        for gav in selected_ga:
            row = sub.iloc[int(np.argmin(np.abs(sub["ga"].to_numpy() - gav)))]
            marg_rows.append([bmi_value, row["ga"], row["p_marg"], row["p_cond"], row["p_cond"] - row["p_marg"]])
    pd.DataFrame(marg_rows, columns=["BMI", "孕周", "P_marg", "P_cond", "delta_P_cond_minus_marg"]).to_csv(results / "table_sens_marginal.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame([
        ["保留不健康记录", len(df), main_beta.aic, 0.0, 0.0],
        ["剔除不健康记录", int(healthy.sum()), mhealth.aic, mhealth.aic - main_beta.aic, float(np.max(np.abs(p_health - median_prob["p_marg"].to_numpy())))],
    ], columns=["策略", "N_rec", "AIC", "delta_AIC", "P_ok最大差异"]).to_csv(results / "table_sens_health.csv", index=False, encoding="utf-8-sig")

    body_rows = []
    for body in ["bmi", "weight", "height_weight"]:
        xb, _ = builder.transform(df, interaction=False, body_form=body)
        mb = fit_beta(y, xb)
        low = make_grid(builder, np.array([builder.means["ga"]]), builder.means["bmi"])
        high = low.copy()
        if body == "bmi":
            low["bmi"], high["bmi"] = df["bmi"].quantile(0.25), df["bmi"].quantile(0.75)
        else:
            low["weight"], high["weight"] = df["weight"].quantile(0.25), df["weight"].quantile(0.75)
        xl, _ = builder.transform(low, interaction=False, body_form=body)
        xh_body, _ = builder.transform(high, interaction=False, body_form=body)
        effect = float(beta_mean(mb, xh_body)[0] - beta_mean(mb, xl)[0])
        direction = ("正向" if effect > 0 else "负向") + f"（Q3-Q1预测差={effect:.5f}）"
        body_rows.append([body, mb.aic, mb.bic, direction])
    body_table = pd.DataFrame(body_rows, columns=["协变量形态", "AIC", "BIC", "主效应方向"])
    body_table["delta_AIC"] = body_table["AIC"] - body_table["AIC"].min()
    body_table.to_csv(results / "table_covariate_forms.csv", index=False, encoding="utf-8-sig")
    xp, _ = builder.transform(df, interaction=False, parity=True)
    parity_model = fit_beta(y, xp)
    pd.DataFrame([["年龄+IVF+BMI", no_int_beta.aic, no_int_beta.bic], ["再加孕次与产次", parity_model.aic, parity_model.bic]], columns=["候选模型", "AIC", "BIC"]).to_csv(results / "table_parity_candidates.csv", index=False, encoding="utf-8-sig")

    df[["ga", "bmi", "y", "mother_id"]].assign(y_thr=Y_THR).to_csv(results / "q1_threshold_anchor.csv", index=False, encoding="utf-8-sig")

    pear_ga, spear_ga = pearsonr(df["ga"], df["y"]), spearmanr(df["ga"], df["y"])
    pear_bmi, spear_bmi = pearsonr(df["bmi"], df["y"]), spearmanr(df["bmi"], df["y"])
    long_rows: list[dict[str, object]] = []
    for _, row in comp.iterrows():
        for metric in ["AIC", "BIC", "RMSE", "MAE", "delta_AIC"]:
            if pd.notna(row[metric]):
                long_rows.append({"model": row["模型"], "metric": metric, "value": row[metric], "n": len(df), "note": row["备注"] + ("；长表按结果契约下限截断" if metric in ["AIC", "BIC"] and row[metric] < -10 else "")})
    for row in smooth_rows:
        long_rows.extend([
            {"model": row["平滑项"], "metric": "edf", "value": row["edf"], "n": len(df)},
            {"model": row["平滑项"], "metric": "p_value", "value": row["p值"], "p_value": row["p值"], "n": len(df)},
            {"model": row["平滑项"], "metric": "delta_AIC_removed", "value": row["delta_AIC_removed"], "n": len(df)},
        ])
    for _, row in re_table.iterrows():
        long_rows.append({"model": "random_effects", "metric": row["参数"], "value": row["估计值"], "se": row["标准误"], "n": len(df), "note": row["备注"]})
    for metric, value, pvalue in [
        ("sigma_tech", sigma_tech, np.nan), ("pearson_ga", pear_ga.statistic, pear_ga.pvalue),
        ("spearman_ga", spear_ga.statistic, spear_ga.pvalue), ("pearson_bmi", pear_bmi.statistic, pear_bmi.pvalue),
        ("spearman_bmi", spear_bmi.statistic, spear_bmi.pvalue), ("beta_precision_phi", phi, np.nan),
    ]:
        long_rows.append({"model": "main_model", "metric": metric, "value": value, "p_value": pvalue, "n": len(df), "note": f"tech groups={tech_groups}, rows={tech_rows}" if metric == "sigma_tech" else ""})
    for _, row in gc_table.iterrows():
        long_rows.append({"model": "sens_gc", "metric": "P_ok_diff", "value": row["达标概率最大变化"], "n": row["N_rec"], "note": row["策略"]})
    for _, row in window_table.iterrows():
        long_rows.append({"model": "sens_ga_window", "metric": "P_ok_diff", "value": row["10-25周内P_ok最大差异"], "n": row["N_rec"], "note": row["样本集"]})
    save_long_results(out, long_rows)

    summary = {
        "records": len(df), "mothers": int(df["mother_id"].nunique()),
        "ga_range": [float(df["ga"].min()), float(df["ga"].max())],
        "bmi_range": [float(df["bmi"].min()), float(df["bmi"].max())],
        "y_range": [float(df["y"].min()), float(df["y"].max())],
        "below_threshold_rate": float((df["y"] < Y_THR).mean()),
        "beta_precision_phi": phi, "sigma_tech_logit": sigma_tech, "icc": icc,
        "random_effect_structure": random_note,
        "correlations": {"ga_pearson": [float(pear_ga.statistic), float(pear_ga.pvalue)], "ga_spearman": [float(spear_ga.statistic), float(spear_ga.pvalue)], "bmi_pearson": [float(pear_bmi.statistic), float(pear_bmi.pvalue)], "bmi_spearman": [float(spear_bmi.statistic), float(spear_bmi.pvalue)]},
        "max_sensitivity_changes": {"distribution": float(np.max(np.abs(qmedian["p_quantile"] - median_prob["p_marg"]))), "interaction": float(np.max(np.abs(p_no_int - median_prob["p_marg"]))), "gc": float(np.max(np.abs(p_gc - median_prob["p_marg"]))), "ga_window": float(np.max(np.abs(p_window - median_prob["p_marg"]))), "ga_date_crosscheck": float(np.max(np.abs(p_date - median_prob["p_marg"]))), "health": float(np.max(np.abs(p_health - median_prob["p_marg"])))},
        "elapsed_seconds": time.time() - started,
        "model_note": "two-stage Beta-GAMM: Beta likelihood + REML random intercept/slope variance layer",
        "clinical_threshold": Y_THR, "confidence_level_for_intervals": 0.95,
    }
    (results / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

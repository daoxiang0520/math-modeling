"""2025 C题问题1：NIPT Y染色体浓度的重复测量关系模型。

只使用任务允许的科学计算库。主模型为 logit(Y) 尺度的样条混合效应模型，
它是当前 Python 生态中对 Beta-GAMM 的可执行近似；概率层显式加入残差与随机
效应方差，从而区分条件预测和面向新孕妇的边缘预测。
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
from scipy.stats import chi2, norm, pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import SplineTransformer
import statsmodels.api as sm


Y_THR = 0.04
GA_MIN, GA_MAX = 10.0, 25.0
RANDOM_SEED = 2025


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
        week = float(match.group(1))
        day = int(match.group(2) or 0)
        return week + day / 7.0 if 0 <= day < 7 else np.nan
    try:
        return float(text)
    except ValueError:
        return np.nan


def parse_test_date(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.replace(r"\.0$", "", regex=True).str.strip()
    parsed = pd.to_datetime(text, format="%Y%m%d", errors="coerce")
    return parsed.fillna(pd.to_datetime(series, errors="coerce"))


def load_data(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="男胎检测数据")
    df.columns = [str(c).strip() for c in df.columns]
    df["mother_id"] = df["孕妇代码"].astype(str).str.strip()
    df["ga"] = df["检测孕周"].map(parse_ga)
    df["test_date"] = parse_test_date(df["检测日期"])
    df["lmp_date"] = pd.to_datetime(df["末次月经"], errors="coerce")
    df["ga_date"] = (df["test_date"] - df["lmp_date"]).dt.days / 7.0
    df["ga_date_diff"] = df["ga"] - df["ga_date"]
    df["y"] = pd.to_numeric(df["Y染色体浓度"], errors="coerce")
    df["yz"] = pd.to_numeric(df["Y染色体的Z值"], errors="coerce")
    df["bmi"] = pd.to_numeric(df["孕妇BMI"], errors="coerce")
    df["age"] = pd.to_numeric(df["年龄"], errors="coerce")
    df["gc"] = pd.to_numeric(df["GC含量"], errors="coerce")
    df["weight"] = pd.to_numeric(df["体重"], errors="coerce")
    df["height"] = pd.to_numeric(df["身高"], errors="coerce")
    df["ivf"] = df["IVF妊娠"].astype(str).str.contains("IVF|试管", case=False, regex=True).astype(int)
    df["healthy"] = df["胎儿是否健康"].astype(str).str.strip().eq("是")
    valid = (
        df["mother_id"].ne("")
        & df["y"].between(0, 1, inclusive="neither")
        & df["yz"].notna()
        & df[["ga", "bmi", "age", "gc"]].notna().all(axis=1)
    )
    df = df.loc[valid].copy()
    df["logit_y"] = logit(df["y"].clip(1e-6, 1 - 1e-6))
    df["ga_c"] = df["ga"] - df["ga"].mean()
    df["tech_group"] = df["mother_id"] + "#" + df["检测抽血次数"].astype(str)
    df["is_tech_repeat"] = df.groupby("tech_group")["tech_group"].transform("size").gt(1)
    return df


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
        kwargs = dict(n_knots=n_knots, degree=3, include_bias=False, extrapolation="linear")
        ga = SplineTransformer(**kwargs).fit(df[["ga"]])
        bmi = SplineTransformer(**kwargs).fit(df[["bmi"]])
        age = SplineTransformer(n_knots=4, degree=2, include_bias=False, extrapolation="linear").fit(df[["age"]])
        means = {c: float(df[c].mean()) for c in ["ga", "bmi", "age", "gc", "weight", "height"]}
        scales = {c: max(float(df[c].std(ddof=0)), 1e-8) for c in means}
        return cls(ga, bmi, age, means, scales)

    def transform(self, df: pd.DataFrame, interaction: bool = True, gc: bool = False,
                  body_form: str = "bmi") -> tuple[np.ndarray, list[str]]:
        ga_b = self.ga_spline.transform(df[["ga"]])
        bmi_b = self.bmi_spline.transform(df[["bmi"]])
        age_b = self.age_spline.transform(df[["age"]])
        parts = [ga_b]
        names = [f"s_ga_{i+1}" for i in range(ga_b.shape[1])]
        if body_form == "bmi":
            parts.append(bmi_b)
            names.extend(f"s_bmi_{i+1}" for i in range(bmi_b.shape[1]))
        elif body_form == "weight":
            parts.append(((df[["weight"]].to_numpy() - self.means["weight"]) / self.scales["weight"]))
            names.append("weight_z")
        elif body_form == "height_weight":
            parts.extend([
                (df[["height"]].to_numpy() - self.means["height"]) / self.scales["height"],
                (df[["weight"]].to_numpy() - self.means["weight"]) / self.scales["weight"],
            ])
            names.extend(["height_z", "weight_z"])
        if interaction and body_form == "bmi":
            r = min(self.interaction_rank, ga_b.shape[1], bmi_b.shape[1])
            inter = np.column_stack([ga_b[:, i] * bmi_b[:, j] for i in range(r) for j in range(r)])
            parts.append(inter)
            names.extend(f"ti_ga{i+1}_bmi{j+1}" for i in range(r) for j in range(r))
        parts.extend([age_b, df[["ivf"]].to_numpy()])
        names.extend([f"s_age_{i+1}" for i in range(age_b.shape[1])] + ["ivf"])
        if gc:
            parts.append((df[["gc"]].to_numpy() - self.means["gc"]) / self.scales["gc"])
            names.append("gc_z")
        return sm.add_constant(np.column_stack(parts), has_constant="add"), ["const"] + names


def fit_ols(y: np.ndarray, x: np.ndarray, groups: pd.Series | None = None):
    model = sm.OLS(y, x).fit()
    if groups is not None:
        return model.get_robustcov_results(cov_type="cluster", groups=np.asarray(groups))
    return model


def fit_mixed(df: pd.DataFrame, x: np.ndarray):
    z = np.column_stack([np.ones(len(df)), df["ga_c"].to_numpy()])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            return sm.MixedLM(df["logit_y"], x, groups=df["mother_id"], exog_re=z).fit(
                reml=True, method="lbfgs", maxiter=220, disp=False
            ), "random intercept + GA slope"
        except Exception:
            return sm.MixedLM(df["logit_y"], x, groups=df["mother_id"]).fit(
                reml=True, method="lbfgs", maxiter=220, disp=False
            ), "random intercept fallback"


def grouped_cv(x: np.ndarray, y: np.ndarray, groups: pd.Series) -> tuple[float, float]:
    pred = np.empty(len(y))
    splitter = GroupKFold(n_splits=5)
    for train, test in splitter.split(x, y, groups):
        model = LinearRegression(fit_intercept=False).fit(x[train], y[train])
        pred[test] = model.predict(x[test])
    return float(np.sqrt(mean_squared_error(y, pred))), float(mean_absolute_error(y, pred))


def model_metrics(name: str, model, x: np.ndarray, y: np.ndarray, groups: pd.Series,
                  note: str) -> dict[str, object]:
    k = x.shape[1]
    llf = float(model.llf)
    return {
        "模型": name,
        "AIC": 2 * k - 2 * llf,
        "BIC": np.log(len(y)) * k - 2 * llf,
        "RMSE": grouped_cv(x, y, groups)[0],
        "MAE": grouped_cv(x, y, groups)[1],
        "备注": note,
    }


def predict_fixed(model, x: np.ndarray) -> np.ndarray:
    params = np.asarray(model.fe_params if hasattr(model, "fe_params") else model.params)
    return x @ params[: x.shape[1]]


def tech_error(df: pd.DataFrame) -> tuple[float, int, int]:
    groups = [g["logit_y"].to_numpy() for _, g in df.groupby("tech_group") if len(g) > 1]
    numerator = sum(float(np.sum((arr - arr.mean()) ** 2)) for arr in groups)
    denominator = sum(len(arr) - 1 for arr in groups)
    return float(np.sqrt(numerator / denominator)), len(groups), int(sum(map(len, groups)))


def save_long_results(out: Path, rows: list[dict[str, object]]) -> None:
    columns = ["model", "metric", "value", "se", "p_value", "n", "bmi", "ga", "note"]
    frame = pd.DataFrame(rows)
    for c in columns:
        if c not in frame:
            frame[c] = np.nan
    frame[columns].to_csv(out / "results" / "output.csv", index=False, encoding="utf-8-sig")
    frame[columns].to_csv(out / "results" / "q1.csv", index=False, encoding="utf-8-sig")


def main() -> None:
    started = time.time()
    data_path, out = resolve_paths()
    df = load_data(data_path)
    y = df["logit_y"].to_numpy()
    builder = DesignBuilder.fit(df, n_knots=5)
    x_main, names = builder.transform(df, interaction=True)
    x_no_int, _ = builder.transform(df, interaction=False)
    main_model, random_note = fit_mixed(df, x_main)
    main_ols = fit_ols(y, x_main, df["mother_id"])
    no_int = fit_ols(y, x_no_int, df["mother_id"])
    x_linear = sm.add_constant(df[["ga", "bmi"]].to_numpy(), has_constant="add")
    linear = fit_ols(y, x_linear, df["mother_id"])
    ga = df["ga"].to_numpy()
    x_piece = sm.add_constant(np.column_stack([
        ga, np.maximum(ga - 12.5, 0), np.maximum(ga - 20.0, 0), df["bmi"].to_numpy()
    ]), has_constant="add")
    piece = fit_ols(y, x_piece, df["mother_id"])

    comparisons = [
        model_metrics("Spline mixed model", main_model, x_main, y, df["mother_id"], random_note),
        model_metrics("Spline fixed model", main_ols, x_main, y, df["mother_id"], "cluster-robust SE"),
        model_metrics("Spline without interaction", no_int, x_no_int, y, df["mother_id"], "nested baseline"),
        model_metrics("Piecewise GA baseline", piece, x_piece, y, df["mother_id"], "knots 12.5/20 weeks"),
        model_metrics("Linear baseline", linear, x_linear, y, df["mother_id"], "GA + BMI"),
    ]
    comp = pd.DataFrame(comparisons)
    comp["delta_AIC"] = comp["AIC"] - comp["AIC"].min()
    comp.to_csv(out / "results" / "table_model_comparison.csv", index=False, encoding="utf-8-sig")

    term_indices = {
        "s1(ga)": [i for i, n in enumerate(names) if n.startswith("s_ga")],
        "s2(bmi)": [i for i, n in enumerate(names) if n.startswith("s_bmi")],
        "ti(ga,bmi)": [i for i, n in enumerate(names) if n.startswith("ti_")],
        "s3(age)": [i for i, n in enumerate(names) if n.startswith("s_age")],
        "ivf": [names.index("ivf")],
    }
    smooth_rows = []
    for term, idx in term_indices.items():
        r = np.zeros((len(idx), len(names)))
        r[np.arange(len(idx)), idx] = 1
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            test = main_ols.wald_test(r, scalar=True)
        keep = [j for j in range(len(names)) if j not in idx]
        reduced = sm.OLS(y, x_main[:, keep]).fit()
        smooth_rows.append({
            "平滑项": term,
            "edf": len(idx),
            "p值": float(test.pvalue),
            "delta_AIC_removed": float(reduced.aic - sm.OLS(y, x_main).fit().aic),
            "结论": "显著" if float(test.pvalue) < 0.05 else "未达0.05",
        })
    pd.DataFrame(smooth_rows).to_csv(out / "results" / "table_smooth_terms.csv", index=False, encoding="utf-8-sig")

    cov_re = np.atleast_2d(np.asarray(main_model.cov_re))
    resid_var = float(main_model.scale)
    b0 = float(cov_re[0, 0]) if cov_re.size else 0.0
    b1 = float(cov_re[1, 1]) if cov_re.shape[0] > 1 else 0.0
    b01 = float(cov_re[0, 1]) if cov_re.shape[0] > 1 else 0.0
    icc = b0 / (b0 + resid_var) if b0 + resid_var > 0 else np.nan
    re_table = pd.DataFrame([
        ["random intercept variance", b0], ["random slope variance", b1],
        ["intercept-slope covariance", b01], ["residual variance", resid_var], ["ICC at mean GA", icc],
    ], columns=["参数", "估计值"])
    re_table.to_csv(out / "results" / "table_random_effects.csv", index=False, encoding="utf-8-sig")

    # 原始散点、重复轨迹与质量核查数据
    df[["mother_id", "ga", "bmi", "y", "age", "gc", "ga_date_diff", "is_tech_repeat"]].to_csv(
        out / "results" / "q1_scatter.csv", index=False, encoding="utf-8-sig"
    )
    profile = pd.DataFrame({
        "metric": ["records", "mothers", "median_records_per_mother", "ga_date_absdiff_gt1",
                   "tech_repeat_groups", "tech_repeat_rows", "unhealthy_rows", "below_4pct_rows"],
        "value": [len(df), df["mother_id"].nunique(), df.groupby("mother_id").size().median(),
                  df["ga_date_diff"].abs().gt(1).sum(), df.loc[df["is_tech_repeat"], "tech_group"].nunique(),
                  df["is_tech_repeat"].sum(), (~df["healthy"]).sum(), (df["y"] < Y_THR).sum()],
    })
    profile.to_csv(out / "results" / "data_profile.csv", index=False, encoding="utf-8-sig")

    bmi_levels = np.quantile(df["bmi"], [0.25, 0.5, 0.75])
    ga_grid = np.linspace(GA_MIN, GA_MAX, 151)
    grid_rows = []
    for bmi_value in bmi_levels:
        grid = pd.DataFrame({
            "ga": ga_grid, "bmi": bmi_value, "age": builder.means["age"], "ivf": 0,
            "gc": builder.means["gc"], "weight": builder.means["weight"], "height": builder.means["height"]
        })
        xg, _ = builder.transform(grid, interaction=True)
        eta = predict_fixed(main_model, xg)
        ga_c = ga_grid - df["ga"].mean()
        if cov_re.shape[0] > 1:
            random_var = b0 + 2 * ga_c * b01 + ga_c**2 * b1
        else:
            random_var = np.full_like(ga_grid, b0)
        sd_cond = np.sqrt(max(resid_var, 1e-12))
        sd_marg = np.sqrt(np.maximum(resid_var + random_var, 1e-12))
        p_cond = 1 - norm.cdf((logit(Y_THR) - eta) / sd_cond)
        p_marg = 1 - norm.cdf((logit(Y_THR) - eta) / sd_marg)
        # 固定效应近似置信带
        cov_fixed = np.asarray(main_model.cov_params())[: x_main.shape[1], : x_main.shape[1]]
        se_eta = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", xg, cov_fixed, xg), 0))
        grid_rows.append(pd.DataFrame({
            "ga": ga_grid, "bmi": bmi_value, "eta": eta, "mean_y": expit(eta),
            "mean_y_lo": expit(eta - 1.96 * se_eta), "mean_y_hi": expit(eta + 1.96 * se_eta),
            "p_cond": p_cond, "p_marg": p_marg,
        }))
    prob = pd.concat(grid_rows, ignore_index=True)
    prob.to_csv(out / "results" / "q1_prob_curves.csv", index=False, encoding="utf-8-sig")
    prob.to_csv(out / "results" / "q1_smooth_ga.csv", index=False, encoding="utf-8-sig")

    # 二维预测面：中心年龄、自然受孕，随机效应为0。
    ga_h = np.linspace(GA_MIN, GA_MAX, 80)
    bmi_h = np.linspace(df["bmi"].quantile(0.02), df["bmi"].quantile(0.98), 70)
    gh, bh = np.meshgrid(ga_h, bmi_h)
    heat = pd.DataFrame({
        "ga": gh.ravel(), "bmi": bh.ravel(), "age": builder.means["age"], "ivf": 0,
        "gc": builder.means["gc"], "weight": builder.means["weight"], "height": builder.means["height"]
    })
    xh, _ = builder.transform(heat, interaction=True)
    heat["pred_y"] = expit(predict_fixed(main_model, xh))
    heat.to_csv(out / "results" / "q1_ti_heatmap.csv", index=False, encoding="utf-8-sig")

    # 分位数回归辅助模型，使用无交互样条并保存中心BMI的完整分位数族。
    qx = x_no_int[:, : 1 + (builder.ga_spline.n_features_out_ + builder.bmi_spline.n_features_out_)]
    q_grid = pd.DataFrame({
        "ga": ga_grid, "bmi": builder.means["bmi"], "age": builder.means["age"], "ivf": 0,
        "gc": builder.means["gc"], "weight": builder.means["weight"], "height": builder.means["height"]
    })
    qxg_full, _ = builder.transform(q_grid, interaction=False)
    qxg = qxg_full[:, : qx.shape[1]]
    q_rows = []
    q_pred_matrix = []
    taus = np.round(np.arange(0.05, 1.0, 0.05), 2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for tau in taus:
            qr = sm.QuantReg(y, qx).fit(q=float(tau), max_iter=1200)
            pred = qxg @ np.asarray(qr.params)
            q_pred_matrix.append(pred)
            q_rows.append(pd.DataFrame({"ga": ga_grid, "bmi": builder.means["bmi"], "tau": tau, "q_y": expit(pred)}))
    qdf = pd.concat(q_rows, ignore_index=True)
    qdf.to_csv(out / "results" / "q1_quantile_curves.csv", index=False, encoding="utf-8-sig")

    # 分位数阈值反演，与混合模型边缘概率对照。
    qmat = np.asarray(q_pred_matrix)
    threshold_logit = logit(Y_THR)
    inv_rows = []
    median_prob = prob[np.isclose(prob["bmi"], bmi_levels[1])].reset_index(drop=True)
    for j, gav in enumerate(ga_grid):
        vals = np.maximum.accumulate(qmat[:, j])
        tau_star = float(np.interp(threshold_logit, vals, taus, left=0.0, right=1.0))
        inv_rows.append([gav, builder.means["bmi"], 1 - tau_star, median_prob.loc[j, "p_marg"]])
    qcheck = pd.DataFrame(inv_rows, columns=["ga", "bmi", "p_quantile", "p_marg"])
    qcheck["abs_diff"] = (qcheck["p_quantile"] - qcheck["p_marg"]).abs()
    qcheck.to_csv(out / "results" / "table_quantile_check.csv", index=False, encoding="utf-8-sig")

    fitted = predict_fixed(main_model, x_main)
    resid = y - fitted
    pd.DataFrame({"fitted": fitted, "resid": resid, "ga": df["ga"], "bmi": df["bmi"],
                  "mother_id": df["mother_id"]}).to_csv(out / "results" / "q1_resid.csv", index=False, encoding="utf-8-sig")

    sigma_tech, tech_groups, tech_rows = tech_error(df)
    # 敏感性：分布/交互/GC/条件边缘/孕周窗口。
    sens_dist = qcheck.rename(columns={"p_quantile": "quantile", "p_marg": "mixed_normal"})
    sens_dist.to_csv(out / "results" / "sens_dist.csv", index=False, encoding="utf-8-sig")

    xg_no, _ = builder.transform(pd.DataFrame({
        "ga": ga_grid, "bmi": builder.means["bmi"], "age": builder.means["age"], "ivf": 0,
        "gc": builder.means["gc"], "weight": builder.means["weight"], "height": builder.means["height"]
    }), interaction=False)
    no_int_plain = sm.OLS(y, x_no_int).fit()
    sens_int = median_prob[["ga", "p_marg"]].copy()
    eta_no = xg_no @ np.asarray(no_int_plain.params)
    sens_int["with_interaction"] = sens_int.pop("p_marg")
    sens_int["without_interaction"] = 1 - norm.cdf((logit(Y_THR) - eta_no) / np.sqrt(max(no_int_plain.scale, 1e-12)))
    sens_int.to_csv(out / "results" / "sens_interaction.csv", index=False, encoding="utf-8-sig")

    x_gc, _ = builder.transform(df, interaction=True, gc=True)
    gc_model = sm.OLS(y, x_gc).fit()
    gc_grid = pd.DataFrame({
        "ga": ga_grid, "bmi": builder.means["bmi"], "age": builder.means["age"], "ivf": 0,
        "gc": builder.means["gc"], "weight": builder.means["weight"], "height": builder.means["height"]
    })
    xgg, _ = builder.transform(gc_grid, interaction=True, gc=True)
    p_gc = 1 - norm.cdf((logit(Y_THR) - xgg @ np.asarray(gc_model.params)) / np.sqrt(max(gc_model.scale, 1e-12)))
    pd.DataFrame({"ga": ga_grid, "without_gc": median_prob["p_marg"], "with_gc": p_gc}).to_csv(
        out / "results" / "sens_gc.csv", index=False, encoding="utf-8-sig"
    )
    median_prob[["ga", "p_cond", "p_marg"]].to_csv(out / "results" / "sens_marginal.csv", index=False, encoding="utf-8-sig")

    in_window = df["ga"].between(GA_MIN, GA_MAX)
    window_builder = DesignBuilder.fit(df.loc[in_window], n_knots=5)
    xw, _ = window_builder.transform(df.loc[in_window], interaction=True)
    mw = sm.OLS(df.loc[in_window, "logit_y"], xw).fit()
    wg = pd.DataFrame({
        "ga": ga_grid, "bmi": builder.means["bmi"], "age": builder.means["age"], "ivf": 0,
        "gc": builder.means["gc"], "weight": builder.means["weight"], "height": builder.means["height"]
    })
    xwg, _ = window_builder.transform(wg, interaction=True)
    p_window = 1 - norm.cdf((logit(Y_THR) - xwg @ np.asarray(mw.params)) / np.sqrt(max(mw.scale, 1e-12)))
    pd.DataFrame({"ga": ga_grid, "all_records": median_prob["p_marg"], "only_10_25": p_window}).to_csv(
        out / "results" / "sens_ga_window.csv", index=False, encoding="utf-8-sig"
    )

    # 阈值锚定数据不做额外聚合，保留每条记录。
    df[["ga", "bmi", "y", "mother_id"]].assign(y_thr=Y_THR).to_csv(
        out / "results" / "q1_threshold_anchor.csv", index=False, encoding="utf-8-sig"
    )

    body_rows = []
    for body in ["bmi", "weight", "height_weight"]:
        xb, _ = builder.transform(df, interaction=False, body_form=body)
        mb = sm.OLS(y, xb).fit()
        body_rows.append([body, mb.aic, mb.bic, "fixed-effect comparison"])
    body_table = pd.DataFrame(body_rows, columns=["协变量形态", "AIC", "BIC", "备注"])
    body_table["delta_AIC"] = body_table["AIC"] - body_table["AIC"].min()
    body_table.to_csv(out / "results" / "table_covariate_forms.csv", index=False, encoding="utf-8-sig")

    pear_ga = pearsonr(df["ga"], df["y"])
    spear_ga = spearmanr(df["ga"], df["y"])
    pear_bmi = pearsonr(df["bmi"], df["y"])
    spear_bmi = spearmanr(df["bmi"], df["y"])
    long_rows: list[dict[str, object]] = []
    for _, row in comp.iterrows():
        for metric in ["AIC", "BIC", "RMSE", "MAE", "delta_AIC"]:
            long_rows.append({"model": row["模型"], "metric": metric, "value": row[metric], "n": len(df), "note": row["备注"]})
    for row in smooth_rows:
        long_rows.append({"model": row["平滑项"], "metric": "p_value", "value": row["p值"],
                          "p_value": row["p值"], "n": len(df), "note": row["结论"]})
        long_rows.append({"model": row["平滑项"], "metric": "delta_AIC_removed", "value": row["delta_AIC_removed"], "n": len(df)})
    for metric, value in [("sigma_b0_sq", b0), ("sigma_b1_sq", b1), ("ICC", icc),
                          ("sigma_tech", sigma_tech), ("pearson_ga", pear_ga.statistic),
                          ("spearman_ga", spear_ga.statistic), ("pearson_bmi", pear_bmi.statistic),
                          ("spearman_bmi", spear_bmi.statistic)]:
        long_rows.append({"model": "main_model", "metric": metric, "value": value, "n": len(df),
                          "note": f"tech groups={tech_groups}, rows={tech_rows}" if metric == "sigma_tech" else ""})
    save_long_results(out, long_rows)

    summary = {
        "records": len(df), "mothers": int(df["mother_id"].nunique()),
        "ga_range": [float(df["ga"].min()), float(df["ga"].max())],
        "bmi_range": [float(df["bmi"].min()), float(df["bmi"].max())],
        "y_range": [float(df["y"].min()), float(df["y"].max())],
        "below_threshold_rate": float((df["y"] < Y_THR).mean()),
        "sigma_tech_logit": sigma_tech, "icc": icc,
        "correlations": {
            "ga_pearson": [float(pear_ga.statistic), float(pear_ga.pvalue)],
            "ga_spearman": [float(spear_ga.statistic), float(spear_ga.pvalue)],
            "bmi_pearson": [float(pear_bmi.statistic), float(pear_bmi.pvalue)],
            "bmi_spearman": [float(spear_bmi.statistic), float(spear_bmi.pvalue)],
        },
        "elapsed_seconds": time.time() - started,
        "model_note": "logit-normal spline mixed model; executable approximation to Beta-GAMM",
    }
    (out / "results" / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

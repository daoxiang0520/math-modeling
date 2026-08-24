"""2025 高教社杯 C 题问题 2：BMI 分组与最佳 NIPT 时点。

严格遵循 coder_task.md 的动态 LTM：男胎孕妇为分析单位，重新拟合简约
分段线性 Beta 均值模型，以 REML 混合模型估计随机截距/随机孕周斜率方差，
对新孕妇随机效应积分后构造达标概率；用一维动态规划监督分箱，以综合风险
损失直接选择组时点；用测量误差 MC 与孕妇层 bootstrap 双通道传播不确定性。

实现口径：statsmodels 没有联合 Beta 混合模型求解器，故主拟合采用可审计的
两阶段近似（Beta 固定效应 + MixedLM 方差层）。bootstrap 为满足 90 秒契约，
重拟合 Beta 固定效应，并用孕妇残差截距/斜率的经验协方差快速更新方差层；
实际重复次数与任何运行时降级均写入 q2_summary.json 和报告数据表。
"""

from __future__ import annotations

import json
import os
import re
import time
import warnings
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.stats import beta as beta_dist
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import SplineTransformer
import statsmodels.api as sm
from statsmodels.othermod.betareg import BetaModel


Y_THR = 0.04
GA_MIN, GA_MAX, GA_STEP = 10.0, 25.0, 0.1
KNOTS = (12.5, 20.0)
RHO_DEFAULT = 1.0
N_MIN_DEFAULT = 20
SIGMA_TECH = 0.133
RANDOM_SEED = 2025
MC_DRAWS = 1000
BOOTSTRAP_TARGET = 100
BOOTSTRAP_ACTUAL = 100
BOOTSTRAP_MC_DRAWS = 1000


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


def parse_test_date(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.replace(r"\.0$", "", regex=True).str.strip()
    parsed = pd.to_datetime(text, format="%Y%m%d", errors="coerce")
    return parsed.fillna(pd.to_datetime(series, errors="coerce"))


def load_data(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name="男胎检测数据")
    raw.columns = [str(c).strip() for c in raw.columns]
    raw["mother_code_num"] = pd.to_numeric(
        raw["孕妇代码"].astype(str).str.replace("A", "", regex=False), errors="coerce"
    )
    raw["mother_id"] = raw["孕妇代码"].astype(str).str.strip()
    raw["lmp_date"] = pd.to_datetime(raw["末次月经"], errors="coerce")
    raw["test_date"] = parse_test_date(raw["检测日期"])
    raw["ga"] = raw["检测孕周"].map(parse_ga)
    raw["ga_date"] = (raw["test_date"] - raw["lmp_date"]).dt.days / 7.0
    raw["ga_date_diff"] = raw["ga"] - raw["ga_date"]
    for source, target in [
        ("Y染色体浓度", "y"), ("Y染色体的Z值", "yz"), ("孕妇BMI", "bmi"),
        ("年龄", "age"), ("GC含量", "gc"),
    ]:
        raw[target] = pd.to_numeric(raw[source], errors="coerce")
    raw["ivf"] = raw["IVF妊娠"].astype(str).str.contains("IVF|试管", case=False, regex=True).astype(int)
    raw["healthy"] = raw["胎儿是否健康"].astype(str).str.strip().eq("是")
    # 按 parse_hints 解析非整倍体；组合文本只用于完成解析核验，不进入 Q2 决策。
    aneup = raw["染色体的非整倍体"].astype("string")
    raw["aneuploid_numeric_parse"] = pd.to_numeric(
        aneup.str.replace("T", "", regex=False), errors="coerce"
    )
    valid = (
        raw["mother_code_num"].notna() & raw["mother_id"].ne("")
        & raw["y"].between(0, 1, inclusive="neither") & raw["yz"].notna()
        & raw[["ga", "bmi", "age", "gc"]].notna().all(axis=1)
    )
    df = raw.loc[valid].copy()
    df["logit_y"] = logit(df["y"].clip(1e-7, 1 - 1e-7))
    bmi_unit = df.groupby("mother_id")["bmi"].median().rename("b_i")
    df = df.join(bmi_unit, on="mother_id")
    return df


def design(df: pd.DataFrame, form: str = "piecewise", gc: bool = False) -> tuple[np.ndarray, list[str]]:
    ga = df["ga"].to_numpy(float)
    cols = [ga]
    names = ["ga"]
    if form == "piecewise":
        cols += [np.maximum(ga - KNOTS[0], 0), np.maximum(ga - KNOTS[1], 0)]
        names += ["hinge_12_5", "hinge_20"]
    elif form != "linear":
        raise ValueError(form)
    cols += [df["b_i"].to_numpy(float), df["age"].to_numpy(float), df["ivf"].to_numpy(float)]
    names += ["bmi", "age", "ivf"]
    if gc:
        cols.append(df["gc"].to_numpy(float))
        names.append("gc")
    return sm.add_constant(np.column_stack(cols), has_constant="add"), ["const"] + names


def decision_design(t: np.ndarray, b: np.ndarray, age_mean: float, form: str = "piecewise",
                    gc_value: float | None = None) -> np.ndarray:
    tt, bb = np.meshgrid(np.asarray(t, float), np.asarray(b, float), indexing="ij")
    cols = [tt.ravel()]
    if form == "piecewise":
        cols += [np.maximum(tt.ravel() - KNOTS[0], 0), np.maximum(tt.ravel() - KNOTS[1], 0)]
    cols += [bb.ravel(), np.full(tt.size, age_mean), np.zeros(tt.size)]
    if gc_value is not None:
        cols.append(np.full(tt.size, gc_value))
    return sm.add_constant(np.column_stack(cols), has_constant="add")


def fit_beta(y: np.ndarray, x: np.ndarray):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for method in ("bfgs", "newton", "nm"):
            try:
                result = BetaModel(y, x).fit(method=method, maxiter=350, disp=False)
                if np.all(np.isfinite(result.params)):
                    return result
            except Exception:
                continue
    raise RuntimeError("BetaModel failed")


def phi_of(result, p: int) -> float:
    return float(np.exp(np.clip(np.asarray(result.params)[p], -15, 15)))


def psd_cov(cov: np.ndarray) -> np.ndarray:
    cov = np.atleast_2d(np.asarray(cov, float))
    if cov.shape == (1, 1):
        cov = np.array([[cov[0, 0], 0.0], [0.0, 0.0]])
    if cov.shape != (2, 2) or not np.all(np.isfinite(cov)):
        return np.zeros((2, 2))
    val, vec = np.linalg.eigh((cov + cov.T) / 2)
    return vec @ np.diag(np.maximum(val, 0)) @ vec.T


def quick_random_cov(df: pd.DataFrame, x: np.ndarray, beta_params: np.ndarray, ga_mean: float) -> np.ndarray:
    resid = df["logit_y"].to_numpy() - x @ beta_params
    coefs = []
    for _, idx in df.groupby("boot_id" if "boot_id" in df else "mother_id").groups.items():
        ii = np.asarray(list(idx), dtype=int)
        z = np.column_stack([np.ones(len(ii)), df.iloc[ii]["ga"].to_numpy() - ga_mean])
        coefs.append(np.linalg.lstsq(z, resid[ii], rcond=None)[0])
    cov = np.cov(np.asarray(coefs).T, ddof=1) if len(coefs) > 2 else np.zeros((2, 2))
    return psd_cov(cov)


@dataclass
class Q2Model:
    beta: object
    cov_re: np.ndarray
    form: str
    ga_mean: float
    age_mean: float
    gc_mean: float
    random_note: str

    @property
    def params(self) -> np.ndarray:
        p = self.beta.model.exog.shape[1]
        return np.asarray(self.beta.params)[:p]

    @property
    def phi(self) -> float:
        return phi_of(self.beta, self.beta.model.exog.shape[1])


def fit_model(df: pd.DataFrame, form: str = "piecewise", gc: bool = False,
              mixed: bool = True) -> Q2Model:
    x, _ = design(df, form, gc)
    beta = fit_beta(df["y"].to_numpy(), x)
    ga_mean = float(df["ga"].mean())
    cov = None
    note = "empirical residual intercept/slope covariance"
    if mixed:
        z = np.column_stack([np.ones(len(df)), df["ga"].to_numpy() - ga_mean])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                mlm = sm.MixedLM(df["logit_y"], x, groups=df["mother_id"], exog_re=z).fit(
                    reml=True, method="lbfgs", maxiter=220, disp=False
                )
                cov = psd_cov(np.asarray(mlm.cov_re))
                note = "REML random intercept + centered-GA slope"
            except Exception:
                cov = None
    if cov is None:
        cov = quick_random_cov(df.reset_index(drop=True), x, np.asarray(beta.params)[:x.shape[1]], ga_mean)
    return Q2Model(beta, cov, form, ga_mean, float(df["age"].mean()), float(df["gc"].mean()), note)


def grouped_cv(df: pd.DataFrame, form: str) -> tuple[float, float]:
    pred = np.empty(len(df))
    y = df["y"].to_numpy()
    groups = df["mother_id"].to_numpy()
    for train, test in GroupKFold(n_splits=5).split(df, y, groups):
        xtr, _ = design(df.iloc[train], form)
        xte, _ = design(df.iloc[test], form)
        fit = fit_beta(y[train], xtr)
        pred[test] = expit(xte @ np.asarray(fit.params)[:xtr.shape[1]])
    return float(np.sqrt(mean_squared_error(y, pred))), float(mean_absolute_error(y, pred))


def random_draws(cov: np.ndarray, n: int, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).multivariate_normal(np.zeros(2), psd_cov(cov), size=n)


def probability_from_eta(eta: np.ndarray, t: np.ndarray, model: Q2Model, sigma: float = 0.0,
                         n_draws: int = MC_DRAWS, seed: int = RANDOM_SEED) -> np.ndarray:
    """eta shape=(T,N); MC over random effects and optional measurement error.

    Beta survival probability is a smooth one-dimensional function of the logit mean.
    A 768-node adaptive lookup per GA point avoids tens of millions of expensive special-
    function calls while retaining every requested MC draw.  The interpolation error is
    audited against direct evaluation in the report tables.
    """
    rng = np.random.default_rng(seed)
    re = rng.multivariate_normal(np.zeros(2), psd_cov(model.cov_re), size=n_draws)
    u = rng.normal(0, sigma, size=n_draws) if sigma > 0 else np.zeros(n_draws)
    acc = np.zeros_like(eta, dtype=float)
    for j, week in enumerate(t):
        add = re[:, 0] + re[:, 1] * (week - model.ga_mean) + u
        z = eta[j][None, :] + add[:, None]
        lo, hi = float(z.min()) - 0.02, float(z.max()) + 0.02
        lookup_z = np.linspace(lo, hi, 768)
        lookup_mu = np.clip(expit(lookup_z), 1e-8, 1 - 1e-8)
        lookup_p = beta_dist.sf(Y_THR, lookup_mu * model.phi, (1 - lookup_mu) * model.phi)
        acc[j] = np.interp(z.ravel(), lookup_z, lookup_p).reshape(z.shape).mean(axis=0)
    return np.clip(acc, 0, 1)


def probability_matrix(model: Q2Model, t: np.ndarray, b: np.ndarray, sigma: float = 0.0,
                       n_draws: int = MC_DRAWS, seed: int = RANDOM_SEED,
                       gc: bool = False) -> np.ndarray:
    x = decision_design(t, b, model.age_mean, model.form, model.gc_mean if gc else None)
    eta = (x @ model.params).reshape(len(t), len(b))
    return probability_from_eta(eta, t, model, sigma, n_draws, seed)


def interpolation_audit(model: Q2Model, seed: int = RANDOM_SEED) -> pd.DataFrame:
    """Compare lookup-MC with direct Beta survival evaluation on a small fixed grid."""
    tt = np.array([10.0, 12.0, 16.0, 20.0, 25.0])
    bb = np.array([28.0, 32.0, 36.0])
    x = decision_design(tt, bb, model.age_mean, model.form)
    eta = (x @ model.params).reshape(len(tt), len(bb))
    approx = probability_from_eta(eta, tt, model, 0, 400, seed)
    rng = np.random.default_rng(seed)
    re = rng.multivariate_normal(np.zeros(2), psd_cov(model.cov_re), size=400)
    direct = np.empty_like(approx)
    for j, week in enumerate(tt):
        z = eta[j][None, :] + re[:, 0, None] + re[:, 1, None] * (week - model.ga_mean)
        mu = np.clip(expit(z), 1e-8, 1 - 1e-8)
        direct[j] = beta_dist.sf(Y_THR, mu * model.phi, (1 - mu) * model.phi).mean(axis=0)
    rows = []
    for i, week in enumerate(tt):
        for j, bmi in enumerate(bb):
            rows.append({"ga": week, "bmi": bmi, "p_lookup": approx[i, j],
                         "p_direct": direct[i, j], "absolute_error": abs(approx[i, j] - direct[i, j])})
    return pd.DataFrame(rows)


def delay_risk(t: np.ndarray, kind: str = "slope", threshold: float = 12.0) -> np.ndarray:
    t = np.asarray(t, float)
    if kind == "step":
        return np.where(t <= threshold, 0.0, 1.0)
    multiplier = 2.0 if kind == "steep" else 1.0
    return np.where(t <= threshold, 0.0, multiplier * (t - threshold) / (27.0 - threshold))


def refine_minimum(t: np.ndarray, loss: np.ndarray) -> float:
    j = int(np.argmin(loss))
    if j == 0 or j == len(t) - 1:
        return float(t[j])
    x = t[j - 1:j + 2]
    y = loss[j - 1:j + 2]
    coef = np.polyfit(x, y, 2)
    if coef[0] <= 0 or not np.all(np.isfinite(coef)):
        return float(t[j])
    vertex = -coef[1] / (2 * coef[0])
    return float(np.clip(vertex, x[0], x[-1]))


def group_optimum(p: np.ndarray, t: np.ndarray, rho: float = RHO_DEFAULT,
                  risk_kind: str = "slope", threshold: float = 12.0) -> tuple[float, float, float, np.ndarray]:
    pbar = np.mean(p, axis=1)
    loss = rho * (1 - pbar) + delay_risk(t, risk_kind, threshold)
    tt = refine_minimum(t, loss)
    return tt, float(np.interp(tt, t, pbar)), float(np.interp(tt, t, loss)), loss


def pava(values: np.ndarray) -> np.ndarray:
    x = np.arange(len(values))
    return IsotonicRegression(increasing=True, out_of_bounds="clip").fit_transform(x, values)


def segment_tables(p_sorted: np.ndarray, t: np.ndarray, b_sorted: np.ndarray,
                   rho: float, n_min: int, risk_kind: str = "slope",
                   threshold: float = 12.0) -> tuple[np.ndarray, np.ndarray]:
    n = p_sorted.shape[1]
    cumulative = np.column_stack([np.zeros(len(t)), np.cumsum(p_sorted, axis=1)])
    costs = np.full((n + 1, n + 1), np.inf)
    times = np.full((n + 1, n + 1), np.nan)
    risk = delay_risk(t, risk_kind, threshold)[:, None]
    for i in range(n):
        js = np.arange(i + n_min, n + 1)
        if len(js) == 0:
            continue
        means = (cumulative[:, js] - cumulative[:, [i]]) / (js - i)[None, :]
        losses = rho * (1 - means) + risk
        idx = np.argmin(losses, axis=0)
        valid_end = (js == n) | (b_sorted[js - 1] < b_sorted[np.minimum(js, n - 1)])
        for col, j in enumerate(js[valid_end]):
            jj = np.flatnonzero(js == j)[0]
            costs[i, j] = (j - i) * losses[idx[jj], jj]
            times[i, j] = t[idx[jj]]
    return costs, times


def solve_bins(p: np.ndarray, t: np.ndarray, b: np.ndarray, n_min: int = N_MIN_DEFAULT,
               rho: float = RHO_DEFAULT, risk_kind: str = "slope", threshold: float = 12.0,
               max_k: int = 6) -> dict[int, dict[str, object]]:
    order = np.argsort(b, kind="mergesort")
    bs, ps = np.asarray(b)[order], p[:, order]
    n = len(bs)
    costs, seg_t = segment_tables(ps, t, bs, rho, n_min, risk_kind, threshold)
    dp = np.full((max_k + 1, n + 1), np.inf)
    prev = np.full((max_k + 1, n + 1), -1, int)
    last_t = np.full((max_k + 1, n + 1), np.nan)
    dp[0, 0] = 0.0
    for k in range(1, max_k + 1):
        for j in range(k * n_min, n + 1):
            for i in range((k - 1) * n_min, j - n_min + 1):
                if not np.isfinite(dp[k - 1, i]) or not np.isfinite(costs[i, j]):
                    continue
                if k > 1 and np.isfinite(last_t[k - 1, i]) and seg_t[i, j] + 1e-10 < last_t[k - 1, i]:
                    continue
                value = dp[k - 1, i] + costs[i, j]
                if value < dp[k, j]:
                    dp[k, j], prev[k, j], last_t[k, j] = value, i, seg_t[i, j]
    solutions: dict[int, dict[str, object]] = {}
    for k in range(2, max_k + 1):
        if not np.isfinite(dp[k, n]):
            continue
        cuts, j = [n], n
        for kk in range(k, 0, -1):
            j = int(prev[kk, j])
            cuts.append(j)
        cuts = sorted(cuts)
        raw = [(bs[c - 1] + bs[c]) / 2 for c in cuts[1:-1]]
        boundaries = round_boundaries(raw, bs, n_min)
        gid = assign_groups(bs, boundaries)
        group_times = [group_optimum(ps[:, gid == g], t, rho, risk_kind, threshold)[0]
                       for g in range(k)]
        solutions[k] = {
            "J": float(sum((gid == g).sum() * group_optimum(ps[:, gid == g], t, rho, risk_kind, threshold)[2]
                           for g in range(k))),
            "boundaries": boundaries,
            "times": group_times,
            "order": order,
        }
    return solutions


def round_boundaries(raw: list[float], b_sorted: np.ndarray, n_min: int) -> list[float]:
    if not raw:
        return []
    n, chosen, prev_idx = len(b_sorted), [], 0
    for slot, value in enumerate(raw):
        remain = len(raw) - slot
        candidates = np.arange(15.0, 50.01, 0.5)
        candidates = sorted(candidates, key=lambda q: (abs(q - value), q))
        pick = None
        for q in candidates:
            idx = int(np.searchsorted(b_sorted, q, side="left"))
            if idx - prev_idx >= n_min and n - idx >= remain * n_min and (not chosen or q > chosen[-1]):
                pick = float(q)
                prev_idx = idx
                break
        if pick is None:
            idx = int(np.searchsorted(b_sorted, value, side="right"))
            pick = float(np.ceil(((b_sorted[max(idx - 1, 0)] + b_sorted[min(idx, n - 1)]) / 2) * 2) / 2)
            prev_idx = int(np.searchsorted(b_sorted, pick, side="left"))
        chosen.append(pick)
    return chosen


def assign_groups(b: np.ndarray, boundaries: list[float]) -> np.ndarray:
    return np.digitize(np.asarray(b, float), np.asarray(boundaries, float), right=False)


def elbow_k(solutions: dict[int, dict[str, object]]) -> int:
    ks = np.array(sorted(solutions), float)
    js = np.array([solutions[int(k)]["J"] for k in ks], float)
    if len(ks) <= 2 or np.ptp(js) < 1e-10:
        return int(ks[0])
    x = (ks - ks.min()) / np.ptp(ks)
    y = (js - js.min()) / np.ptp(js)
    distance = np.abs((y[-1] - y[0]) * x - (x[-1] - x[0]) * y + x[-1] * y[0] - y[-1] * x[0])
    return int(ks[int(np.argmax(distance))])


def summarize_groups(p: np.ndarray, t: np.ndarray, b: np.ndarray, boundaries: list[float],
                     rho: float = RHO_DEFAULT, risk_kind: str = "slope",
                     threshold: float = 12.0) -> tuple[pd.DataFrame, pd.DataFrame]:
    gid = assign_groups(b, boundaries)
    rows, curves = [], []
    k = len(boundaries) + 1
    lows = [float(np.min(b))] + boundaries
    highs = boundaries + [float(np.max(b))]
    for g in range(k):
        mask = gid == g
        tt, pp, ll, loss = group_optimum(p[:, mask], t, rho, risk_kind, threshold)
        rows.append({"group": f"G{g+1}", "bmi_low": lows[g], "bmi_high": highs[g],
                     "n": int(mask.sum()), "median_bmi": float(np.median(b[mask])),
                     "optimal_week": tt, "pbar_at_opt": pp, "expected_loss": ll})
        pbar = p[:, mask].mean(axis=1)
        curves.extend({"group": f"G{g+1}", "ga": float(x), "pbar": float(y),
                       "loss": float(z), "optimal_week": tt} for x, y, z in zip(t, pbar, loss))
    return pd.DataFrame(rows), pd.DataFrame(curves)


def misclassification(model: Q2Model, t0: float, b: np.ndarray, sigma: float,
                      n_draws: int = 800, seed: int = RANDOM_SEED) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    re = rng.multivariate_normal(np.zeros(2), psd_cov(model.cov_re), size=n_draws)
    u = rng.normal(0, sigma, size=n_draws) if sigma > 0 else np.zeros(n_draws)
    x = decision_design(np.array([t0]), b, model.age_mean, model.form)
    eta = (x @ model.params).reshape(1, len(b))[0]
    add = re[:, 0, None] + re[:, 1, None] * (t0 - model.ga_mean)
    mu_true = np.clip(expit(eta[None, :] + add), 1e-8, 1 - 1e-8)
    mu_obs = np.clip(expit(eta[None, :] + add + u[:, None]), 1e-8, 1 - 1e-8)
    p_true = beta_dist.sf(Y_THR, mu_true * model.phi, (1 - mu_true) * model.phi)
    p_obs = beta_dist.sf(Y_THR, mu_obs * model.phi, (1 - mu_obs) * model.phi)
    # 共用潜在分位数耦合：sigma=0 时同一次检测不会凭空错分。
    return float(np.maximum(p_true - p_obs, 0).mean()), float(np.maximum(p_obs - p_true, 0).mean())


def spline_probability(df: pd.DataFrame, t: np.ndarray, b: np.ndarray, interaction: bool,
                       cov: np.ndarray, seed: int) -> np.ndarray:
    kw = dict(n_knots=5, degree=3, include_bias=False, extrapolation="constant")
    ga_sp = SplineTransformer(**kw).fit(df[["ga"]])
    bmi_sp = SplineTransformer(**kw).fit(df[["b_i"]])
    age_sp = SplineTransformer(n_knots=4, degree=2, include_bias=False,
                               extrapolation="constant").fit(df[["age"]])

    def xmake(frame: pd.DataFrame) -> np.ndarray:
        a, c, d = ga_sp.transform(frame[["ga"]]), bmi_sp.transform(frame[["b_i"]]), age_sp.transform(frame[["age"]])
        parts = [a, c]
        if interaction:
            rank = min(3, a.shape[1], c.shape[1])
            parts.append(np.column_stack([a[:, i] * c[:, j] for i in range(rank) for j in range(rank)]))
        parts += [d, frame[["ivf"]].to_numpy(float)]
        return sm.add_constant(np.column_stack(parts), has_constant="add")

    x = xmake(df)
    fit = fit_beta(df["y"].to_numpy(), x)
    tt, bb = np.meshgrid(t, b, indexing="ij")
    grid = pd.DataFrame({"ga": tt.ravel(), "b_i": bb.ravel(),
                         "age": float(df["age"].mean()), "ivf": 0})
    xg = xmake(grid)
    dummy = Q2Model(fit, cov, "linear", float(df["ga"].mean()), float(df["age"].mean()),
                    float(df["gc"].mean()), "spline sensitivity")
    eta = (xg @ np.asarray(fit.params)[:xg.shape[1]]).reshape(len(t), len(b))
    return probability_from_eta(eta, t, dummy, 0, 300, seed)


def quantile_probability(df: pd.DataFrame, t: np.ndarray, b: np.ndarray) -> np.ndarray:
    x, _ = design(df, "piecewise")
    taus = np.arange(0.05, 1.0, 0.05)
    params = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for tau in taus:
            params.append(sm.QuantReg(df["logit_y"], x).fit(q=tau, max_iter=900).params)
    xg = decision_design(t, b, float(df["age"].mean()), "piecewise")
    q = (xg @ np.asarray(params).T).reshape(len(t), len(b), len(taus))
    q = np.maximum.accumulate(q, axis=2)
    q = np.maximum.accumulate(q, axis=0)
    q = np.minimum.accumulate(q[:, ::-1, :], axis=1)[:, ::-1, :]
    target = logit(Y_THR)
    out = np.empty((len(t), len(b)))
    for i in range(len(t)):
        for j in range(len(b)):
            out[i, j] = 1 - np.interp(target, q[i, j], taus, left=0.0, right=1.0)
    return np.clip(out, 0, 1)


def resample_mothers(df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    ids = df["mother_id"].unique()
    sampled = rng.choice(ids, len(ids), replace=True)
    blocks = []
    for k, mid in enumerate(sampled):
        block = df[df["mother_id"].eq(mid)].copy()
        block["boot_id"] = f"B{k:03d}"
        blocks.append(block)
    return pd.concat(blocks, ignore_index=True)


def main() -> None:
    started = time.time()
    data_path, out = resolve_paths()
    results = out / "results"
    df = load_data(data_path).reset_index(drop=True)
    if df["mother_id"].nunique() != 267 or len(df) != 1082:
        warnings.warn(f"男胎有效样本为 {df['mother_id'].nunique()} 位/{len(df)} 条，与预期267/1082不同")

    units = df.groupby("mother_id", as_index=False).agg(
        b_i=("b_i", "first"), age=("age", "mean"), first_ga=("ga", "min")
    ).sort_values("b_i").reset_index(drop=True)
    b = units["b_i"].to_numpy()
    t = np.round(np.arange(GA_MIN, GA_MAX + 0.001, GA_STEP), 1)

    cv_rows = []
    for form in ("piecewise", "linear"):
        rmse, mae = grouped_cv(df, form)
        cv_rows.append({"model_form": form, "cv_rmse": rmse, "mae": mae})
    cv = pd.DataFrame(cv_rows)
    # coder_task 指定默认分段线性；同时完整披露若 CV 选择不同。
    selected_cv = str(cv.loc[cv["cv_rmse"].idxmin(), "model_form"])
    cv["selected"] = cv["model_form"].eq("piecewise")
    cv["cv_winner"] = cv["model_form"].eq(selected_cv)
    cv.to_csv(results / "tab_q2_model_selection.csv", index=False, encoding="utf-8-sig")

    model = fit_model(df, "piecewise", mixed=True)
    p = probability_matrix(model, t, b, 0, MC_DRAWS, RANDOM_SEED)
    if not np.isfinite(p).all():
        raise FloatingPointError("p_marg contains nan/inf")
    interpolation_check = interpolation_audit(model)
    interpolation_check.to_csv(results / "q2_probability_interpolation_audit.csv", index=False,
                               encoding="utf-8-sig")

    # 个体时点、保序与概率反演锚。
    indiv_t = np.array([group_optimum(p[:, [i]], t)[0] for i in range(len(b))])
    indiv_mono = pava(indiv_t)
    anchor = {}
    for target in (0.80, 0.85):
        vals = []
        for i in range(len(b)):
            hit = np.flatnonzero(p[:, i] >= target)
            vals.append(float(t[hit[0]]) if len(hit) else np.nan)
        anchor[target] = vals
    pd.DataFrame({"mother_id": units["mother_id"], "bmi": b, "tstar_raw": indiv_t,
                  "tstar_pava": indiv_mono, "t_p80": anchor[0.80], "t_p85": anchor[0.85]}).to_csv(
        results / "q2_individual_tstar.csv", index=False, encoding="utf-8-sig"
    )

    full_solutions = solve_bins(p, t, b)
    preliminary_k = elbow_k(full_solutions)

    # 用户要求忽略90秒限制：完整执行100次cluster bootstrap，副本内同样采用
    # 0.1周网格与1000次固定种子随机效应/测量误差MC。
    boot_t = np.round(np.arange(GA_MIN, GA_MAX + 0.001, GA_STEP), 1)
    rng = np.random.default_rng(RANDOM_SEED + 17)
    boot_boundaries, boot_times = [], []
    boot_completed = 0
    for it in range(BOOTSTRAP_ACTUAL):
        sub = resample_mothers(df, rng)
        boot_units = sub.groupby("boot_id", as_index=False).agg(b_i=("b_i", "first"))
        bm = fit_model(sub.reset_index(drop=True), "piecewise", mixed=False)
        pb = probability_matrix(bm, boot_t, boot_units["b_i"].to_numpy(), SIGMA_TECH,
                                n_draws=BOOTSTRAP_MC_DRAWS, seed=RANDOM_SEED + 1000 + it)
        sols = solve_bins(pb, boot_t, boot_units["b_i"].to_numpy(), max_k=preliminary_k)
        for kk, sol in sols.items():
            for slot, boundary in enumerate(sol["boundaries"], 1):
                boot_boundaries.append({"replicate": it + 1, "K": kk,
                                        "boundary_index": slot, "boundary": boundary})
            for g, week in enumerate(sol["times"], 1):
                boot_times.append({"replicate": it + 1, "K": kk, "group": f"G{g}",
                                   "optimal_week_total_error": week})
        boot_completed = it + 1

    boot_b = pd.DataFrame(boot_boundaries)
    boot_timing = pd.DataFrame(boot_times)
    boot_b.to_csv(results / "q2_bootstrap_boundaries.csv", index=False, encoding="utf-8-sig")
    boot_timing.to_csv(results / "q2_bootstrap_times.csv", index=False, encoding="utf-8-sig")

    stability = {}
    for kk in sorted(full_solutions):
        sub = boot_b[boot_b["K"].eq(kk)]
        slot_scores = []
        for _, ss in sub.groupby("boundary_index"):
            slot_scores.append(float(ss["boundary"].value_counts(normalize=True).iloc[0]))
        stability[kk] = float(np.mean(slot_scores)) if slot_scores else np.nan
    selected_k = preliminary_k
    while selected_k > 2 and stability.get(selected_k, 0) < 0.5:
        selected_k -= 1
    # 若所有候选分箱的组时点跨度均不足0.5周，额外分组没有可报告决策价值；
    # 按结果契约保留最小允许K=2，并在主表统一标记 distinct_required=False。
    if max(float(np.ptp(sol["times"])) for sol in full_solutions.values()) < 0.5:
        selected_k = 2
    boundaries = list(full_solutions[selected_k]["boundaries"])
    main_table, loss_curves = summarize_groups(p, t, b, boundaries)

    # 通道A：固定主分组，只替换观察概率，以隔离测量误差作用。
    sigma_values = [0.0, 0.5 * SIGMA_TECH, SIGMA_TECH, 2 * SIGMA_TECH]
    sigma_rows, error_rows = [], []
    gid = assign_groups(b, boundaries)
    for si, sigma in enumerate(sigma_values):
        ps = p if sigma == 0 else probability_matrix(model, t, b, sigma, MC_DRAWS,
                                                      RANDOM_SEED + 50 + si)
        for g in range(selected_k):
            tt, _, _, _ = group_optimum(ps[:, gid == g], t)
            base_t = float(main_table.loc[g, "optimal_week"])
            fnr, fpr = misclassification(model, base_t, b[gid == g], sigma,
                                         seed=RANDOM_SEED + 80 + si * 10 + g)
            row = {"sigma": sigma, "sigma_multiple": sigma / SIGMA_TECH,
                   "group": f"G{g+1}", "optimal_week_sigma": tt,
                   "delta_t": tt - base_t, "FNR": fnr, "FPR": fpr}
            sigma_rows.append(row)
            error_rows.append(row.copy())
    sigma_table = pd.DataFrame(sigma_rows)
    sigma_table.to_csv(results / "tab_q2_sensitivity_sigma.csv", index=False, encoding="utf-8-sig")
    sigma_table[["sigma", "sigma_multiple", "group", "optimal_week_sigma", "delta_t"]].to_csv(
        results / "q2_sigma_shift.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(error_rows).to_csv(results / "q2_error_classification.csv", index=False, encoding="utf-8-sig")

    # bootstrap 总误差区间与主结果契约。
    for g in range(selected_k):
        vals = boot_timing[(boot_timing["K"].eq(selected_k)) &
                           (boot_timing["group"].eq(f"G{g+1}"))]["optimal_week_total_error"]
        main_table.loc[g, "ci_low"] = float(vals.quantile(0.025))
        main_table.loc[g, "ci_high"] = float(vals.quantile(0.975))
        tech = sigma_table[(sigma_table["group"].eq(f"G{g+1}")) &
                           np.isclose(sigma_table["sigma"], SIGMA_TECH)]
        main_table.loc[g, "delta_t_sigma_tech"] = float(tech["delta_t"].iloc[0])
    distinct = bool(main_table["optimal_week"].max() - main_table["optimal_week"].min() >= 0.5)
    main_table["distinct_required"] = distinct
    contract_cols = ["group", "bmi_low", "bmi_high", "n", "median_bmi", "optimal_week",
                     "ci_low", "ci_high", "pbar_at_opt", "expected_loss",
                     "delta_t_sigma_tech", "distinct_required"]
    main_table[contract_cols].to_csv(results / "output.csv", index=False, encoding="utf-8-sig")
    main_table[contract_cols].to_csv(results / "tab_q2_main_results.csv", index=False, encoding="utf-8-sig")
    loss_curves.to_csv(results / "q2_group_loss_curves.csv", index=False, encoding="utf-8-sig")

    # 达标概率曲线：选择组中位 BMI，连续轴用折线，颜色+线型冗余。
    prob_rows = []
    for _, row in main_table.iterrows():
        bm = float(row["median_bmi"])
        pp = probability_matrix(model, t, np.array([bm]), 0, MC_DRAWS,
                                RANDOM_SEED + int(row.name) + 200)[:, 0]
        prob_rows.extend({"group": row["group"], "median_bmi": bm, "ga": float(tt),
                          "p_marg": float(v)} for tt, v in zip(t, pp))
    pd.DataFrame(prob_rows).to_csv(results / "q2_prob_curves.csv", index=False, encoding="utf-8-sig")

    # rho敏感性。
    rho_rows = []
    for rho in (0.25, 0.5, 1.0, 2.0, 4.0):
        for g in range(selected_k):
            tt, _, _, _ = group_optimum(p[:, gid == g], t, rho=rho)
            rho_rows.append({"rho": rho, "group": f"G{g+1}", "optimal_week": tt})
    rho_table = pd.DataFrame(rho_rows)
    rho_table.to_csv(results / "tab_q2_sensitivity_rho.csv", index=False, encoding="utf-8-sig")
    rho_table.to_csv(results / "q2_rho_sensitivity.csv", index=False, encoding="utf-8-sig")

    # n_min × K敏感性。
    bin_rows = []
    for nmin in (10, 15, 20, 30):
        sols = solve_bins(p, t, b, n_min=nmin)
        for kk, sol in sols.items():
            ts = np.asarray(sol["times"])
            bin_rows.append({"n_min": nmin, "K": kk,
                             "boundaries": "|".join(f"{q:.1f}" for q in sol["boundaries"]),
                             "total_loss_J": sol["J"], "max_tstar_gap": float(np.ptp(ts))})
    bins_table = pd.DataFrame(bin_rows)
    bins_table.to_csv(results / "tab_q2_sensitivity_bins.csv", index=False, encoding="utf-8-sig")

    # 内部校准：首次观测按最终组及预设GA带。
    first = df.sort_values(["mother_id", "ga"]).groupby("mother_id", as_index=False).first()
    first["group_id"] = assign_groups(first["b_i"].to_numpy(), boundaries)
    calibration = []
    for g in range(selected_k):
        sub = first[first["group_id"].eq(g)]
        ga_ref = float(sub["ga"].median())
        pm = probability_matrix(model, np.array([ga_ref]), sub["b_i"].to_numpy(), 0,
                                MC_DRAWS, RANDOM_SEED + 300 + g).mean()
        calibration.append({"type": "BMI组", "label": f"G{g+1}", "n": len(sub),
                            "ga_reference": ga_ref, "observed_rate": float((sub["y"] >= Y_THR).mean()),
                            "model_probability": float(pm)})
    bands = [(11, 12), (12, 13), (13, 15), (15, 20)]
    for lo, hi in bands:
        sub = first[first["ga"].ge(lo) & first["ga"].lt(hi)]
        if len(sub) == 0:
            continue
        ga_ref = float(sub["ga"].median())
        pm = probability_matrix(model, np.array([ga_ref]), sub["b_i"].to_numpy(), 0,
                                MC_DRAWS, RANDOM_SEED + 350 + lo).mean()
        calibration.append({"type": "GA带", "label": f"{lo}-{hi}周", "n": len(sub),
                            "ga_reference": ga_ref, "observed_rate": float((sub["y"] >= Y_THR).mean()),
                            "model_probability": float(pm)})
    calibration_df = pd.DataFrame(calibration)
    calibration_df["calibration_error"] = calibration_df["model_probability"] - calibration_df["observed_rate"]
    calibration_df.to_csv(results / "q2_calibration.csv", index=False, encoding="utf-8-sig")

    # s2：q1样条与分位数反演；s8：含交互仅作敏感性。
    p_spline = spline_probability(df, t, b, False, model.cov_re, RANDOM_SEED + 400)
    p_inter = spline_probability(df, t, b, True, model.cov_re, RANDOM_SEED + 401)
    p_quant = quantile_probability(df, t, b)
    source_rows = []
    for name, pp in [("简约分段Beta混合", p), ("q1无交互样条", p_spline),
                     ("分位数反演", p_quant), ("含ti交互样条", p_inter)]:
        for g in range(selected_k):
            tt, _, _, _ = group_optimum(pp[:, gid == g], t)
            source_rows.append({"source": name, "group": f"G{g+1}", "optimal_week": tt})
    source_df = pd.DataFrame(source_rows)
    source_df.to_csv(results / "q2_probability_source_sensitivity.csv", index=False, encoding="utf-8-sig")

    # s3 风险函数、阈值；s6 数据质量；s7 GC连续协变量。
    assumption_rows = []
    base_times = main_table["optimal_week"].to_numpy()
    base_boundary = "|".join(map(lambda x: f"{x:.1f}", boundaries))

    def add_assumption(aid: str, base_spec: str, perturbed: str, times_alt: np.ndarray,
                       boundary_alt: str, conclusion: str) -> None:
        assumption_rows.append({"assumption_id": aid, "base_spec": base_spec,
                                "perturbed_spec": perturbed,
                                "max_delta_t_k": float(np.max(np.abs(times_alt - base_times))),
                                "boundary_change": boundary_alt != base_boundary,
                                "conclusion": conclusion})

    linear_model = fit_model(df, "linear", mixed=False)
    p_linear = probability_matrix(linear_model, t, b, 0, 400, RANDOM_SEED + 500)
    alt, _ = summarize_groups(p_linear, t, b, boundaries)
    add_assumption("A3", "分段线性", "纯线性", alt["optimal_week"].to_numpy(), base_boundary,
                   "固定效应形式敏感性")
    p_cond = probability_matrix(Q2Model(model.beta, np.zeros((2, 2)), model.form, model.ga_mean,
                                        model.age_mean, model.gc_mean, "conditional"), t, b, 0, 1,
                                RANDOM_SEED + 501)
    alt, _ = summarize_groups(p_cond, t, b, boundaries)
    add_assumption("A4", "随机效应边缘概率", "随机效应=0条件概率", alt["optimal_week"].to_numpy(),
                   base_boundary, "条件概率可能推早时点")
    for kind, thr in [("step", 12.0), ("steep", 12.0), ("slope", 13.0)]:
        alt, _ = summarize_groups(p, t, b, boundaries, risk_kind=kind, threshold=thr)
        add_assumption("A5", "12周阈值线性斜坡", f"{kind},阈值{thr:g}",
                       alt["optimal_week"].to_numpy(), base_boundary, "延迟风险形式敏感性")
    alt_rho = rho_table[rho_table["rho"].eq(4.0)]["optimal_week"].to_numpy()
    add_assumption("A6", "rho=1", "rho=4", alt_rho, base_boundary, "主观损失比是主要决策不确定性")
    sol_n30 = solve_bins(p, t, b, n_min=30)[selected_k]
    alt, _ = summarize_groups(p, t, b, sol_n30["boundaries"])
    add_assumption("A7", "n_min=20", "n_min=30", alt["optimal_week"].to_numpy(),
                   "|".join(f"{x:.1f}" for x in sol_n30["boundaries"]), "分箱约束敏感性")
    tech_times = sigma_table[np.isclose(sigma_table["sigma"], SIGMA_TECH)]["optimal_week_sigma"].to_numpy()
    add_assumption("A8", "sigma=0", "sigma=sigma_tech", tech_times, base_boundary, "测量误差通道A")
    cal_err = float(calibration_df["calibration_error"].abs().max())
    assumption_rows.append({"assumption_id": "A9", "base_spec": "概率层内部校准",
                            "perturbed_spec": "观察首次达标比例", "max_delta_t_k": np.nan,
                            "boundary_change": False,
                            "conclusion": f"最大绝对校准差={cal_err:.3f}；非外部验证"})

    # 数据质量三类剔除，均重新拟合、重分箱。
    lmp_drop = set(df["ga_date_diff"].abs().nlargest(20).index)
    jumps = []
    for mid, sub in df.sort_values("ga").groupby("mother_id"):
        if sub["y"].diff().abs().max() > 0.05:
            jumps.append(mid)
    quality_specs = [
        ("s6_LMP", df.drop(index=list(lmp_drop)), "剔除日期孕周差异最大的20条"),
        ("s6_jump", df[~df["mother_id"].isin(jumps)], f"剔除组内Y跳变>0.05的{len(jumps)}人"),
        ("s6_health", df[df["healthy"]], "剔除胎儿不健康记录"),
    ]
    for aid, sub, label in quality_specs:
        qm = fit_model(sub.reset_index(drop=True), "piecewise", mixed=False)
        qb = sub.groupby("mother_id")["b_i"].first().to_numpy()
        qp = probability_matrix(qm, t, qb, 0, 300, RANDOM_SEED + 600 + len(assumption_rows))
        qs = solve_bins(qp, t, qb, n_min=min(N_MIN_DEFAULT, max(10, len(qb) // (selected_k + 1))))[selected_k]
        qt, _ = summarize_groups(qp, t, qb, qs["boundaries"])
        add_assumption(aid, "全样本", label, qt["optimal_week"].to_numpy(),
                       "|".join(f"{x:.1f}" for x in qs["boundaries"]), "数据质量敏感性")

    gc_model = fit_model(df, "piecewise", gc=True, mixed=False)
    pgc = probability_matrix(gc_model, t, b, 0, 300, RANDOM_SEED + 700, gc=True)
    gc_alt, _ = summarize_groups(pgc, t, b, boundaries)
    add_assumption("s7_GC", "不含GC", "GC连续协变量", gc_alt["optimal_week"].to_numpy(),
                   base_boundary, "不采用40%-60%硬剔除")
    src_inter = source_df[source_df["source"].eq("含ti交互样条")]["optimal_week"].to_numpy()
    add_assumption("s8_ti", "无交互", "含ti(GA,BMI)交互", src_inter, base_boundary,
                   "交互仅作敏感性，不进入主决策")
    assumptions = pd.DataFrame(assumption_rows)
    assumptions.to_csv(results / "tab_q2_sensitivity_key_assumptions.csv", index=False, encoding="utf-8-sig")

    # 参数、随机效应与分箱选择审计表。
    x_main, names = design(df, "piecewise")
    coef = pd.DataFrame({"term": names, "estimate": model.params,
                         "std_error": np.asarray(model.beta.bse)[:len(names)],
                         "p_value": np.asarray(model.beta.pvalues)[:len(names)]})
    coef.to_csv(results / "tab_q2_model_coefficients.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame([
        {"component": "random_intercept_variance", "estimate": model.cov_re[0, 0]},
        {"component": "random_slope_variance", "estimate": model.cov_re[1, 1]},
        {"component": "intercept_slope_covariance", "estimate": model.cov_re[0, 1]},
        {"component": "beta_precision_phi", "estimate": model.phi},
    ]).to_csv(results / "tab_q2_random_effects.csv", index=False, encoding="utf-8-sig")
    selection_rows = []
    for kk, sol in full_solutions.items():
        selection_rows.append({"K": kk, "total_loss_J": sol["J"],
                               "boundaries": "|".join(f"{x:.1f}" for x in sol["boundaries"]),
                               "bootstrap_boundary_stability": stability.get(kk, np.nan),
                               "preliminary_elbow": kk == preliminary_k,
                               "selected": kk == selected_k})
    pd.DataFrame(selection_rows).to_csv(results / "tab_q2_k_selection.csv", index=False, encoding="utf-8-sig")

    summary = {
        "records": int(len(df)), "mothers": int(df["mother_id"].nunique()),
        "clinical_threshold": Y_THR, "optimization_window": [GA_MIN, GA_MAX],
        "rho": RHO_DEFAULT, "n_min": N_MIN_DEFAULT, "selected_K": selected_k,
        "boundaries": boundaries, "distinct_required": distinct,
        "beta_precision_phi": model.phi, "random_effect_method": model.random_note,
        "cv_winner": selected_cv, "prespecified_model": "piecewise",
        "mc_draws_main": MC_DRAWS, "bootstrap_target": BOOTSTRAP_TARGET,
        "bootstrap_completed": boot_completed, "bootstrap_grid_step": GA_STEP,
        "bootstrap_mc_draws": BOOTSTRAP_MC_DRAWS, "sigma_tech": SIGMA_TECH,
        "first_observation_pass_rate": float((first["y"] >= Y_THR).mean()),
        "max_internal_calibration_error": cal_err,
        "max_lookup_interpolation_error": float(interpolation_check["absolute_error"].max()),
        "runtime_seconds": time.time() - started,
        "method_limit": "two-stage approximate Beta mixed model; bootstrap covariance uses empirical residual slopes",
    }
    (results / "q2_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

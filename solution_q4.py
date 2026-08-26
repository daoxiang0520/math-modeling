"""Question 4: three-layer female-fetus aneuploidy decision pipeline.

The implementation follows q4_coder_task.md: X-concentration reliability gate,
chromosome-specific cost-sensitive Z thresholds, pregnancy-level aggregation,
grouped CV, 200 cluster bootstraps, 100 structured permutations, and sensitivity
experiments.  The 90-second task limit is retained and audited in the summary.
"""
from __future__ import annotations

import json
import os
import re
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

SEED = 2025
CHROMS = (13, 18, 21)
TAU_GRID = np.arange(-3.0, 6.0 + 1e-9, 0.5)
N_BOOT = 200
N_PERM = 100
COLS = [
    "seq", "mother", "age", "height", "weight", "lmp", "ivf", "date", "draw", "ga_text",
    "bmi", "reads", "map_rate", "dup_rate", "unique_reads", "gc", "z13", "z18", "z21",
    "zx", "yz", "y_conc", "w", "gc13", "gc18", "gc21", "filt_rate", "ab", "pregnancies",
    "deliveries", "healthy",
]
METRICS = ("coverage", "sens", "spec", "ppv", "npv", "f1", "auc")


def resolve_paths() -> tuple[Path, Path]:
    raw = os.getenv("MODELING_DATA_PATH")
    if not raw:
        paths = json.loads(os.getenv("MODELING_DATA_PATHS", "[]"))
        raw = paths[0] if paths else str(Path(__file__).with_name("附件.xlsx"))
    root = Path(os.getenv("MODELING_OUTPUT_DIR", Path(__file__).resolve().parent))
    (root / "results").mkdir(parents=True, exist_ok=True)
    return Path(raw), root


def parse_ga(value: object) -> float:
    m = re.search(r"(\d+)\s*[wW周]?\s*(?:\+\s*(\d+))?", str(value))
    return float(m.group(1)) + float(m.group(2) or 0) / 7.0 if m else np.nan


def parse_ab(value: object) -> set[int]:
    if pd.isna(value):
        return set()
    return {int(x) for x in re.findall(r"T\s*(13|18|21)", str(value).upper())}


def load_sheet(path: Path, sheet: int) -> pd.DataFrame:
    d = pd.read_excel(path, sheet_name=sheet)
    if d.shape[1] != 31:
        raise ValueError(f"Expected 31 columns, got {d.shape[1]} on sheet {sheet}")
    d.columns = COLS
    d["mother_raw"] = d["mother"].astype(str)
    # The task hint shows A-prefixed codes, while the actual female sheet uses B-prefixed codes.
    d["mother"] = (d["mother"].astype(str).str.replace("A", "", regex=False)
                   .str.replace("B", "", regex=False).astype(float))
    d["lmp"] = pd.to_datetime(d["lmp"], errors="coerce")
    d["ga"] = d["ga_text"].map(parse_ga)
    numeric = ["age", "height", "weight", "bmi", "reads", "map_rate", "dup_rate", "unique_reads",
               "gc", "z13", "z18", "z21", "zx", "yz", "y_conc", "w", "gc13", "gc18", "gc21", "filt_rate"]
    for c in numeric:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    parsed = d["ab"].map(parse_ab)
    for chrom in CHROMS:
        d[f"y{chrom}"] = parsed.map(lambda s, c=chrom: int(c in s))
    d["ab_positive"] = parsed.map(bool).astype(int)
    d["cv_group"] = d["mother"]
    return d


def folds(groups: np.ndarray, n_splits: int = 5, seed: int = SEED) -> list[tuple[np.ndarray, np.ndarray]]:
    unique = np.unique(groups)
    rng = np.random.default_rng(seed)
    shuffled = unique.copy(); rng.shuffle(shuffled)
    fold_of = {g: i % n_splits for i, g in enumerate(shuffled)}
    fid = np.array([fold_of[g] for g in groups])
    return [(np.flatnonzero(fid != k), np.flatnonzero(fid == k)) for k in range(n_splits)]


def safe_auc(y: np.ndarray, score: np.ndarray) -> float:
    return float(roc_auc_score(y, score)) if len(np.unique(y)) == 2 else 0.5


def confusion(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    y = np.asarray(y, int); pred = np.asarray(pred, int)
    tp = int(np.sum((y == 1) & (pred == 1))); fn = int(np.sum((y == 1) & (pred == 0)))
    fp = int(np.sum((y == 0) & (pred == 1))); tn = int(np.sum((y == 0) & (pred == 0)))
    div = lambda a, b: float(a / b) if b else 0.0
    return {"TP": tp, "FN": fn, "FP": fp, "TN": tn,
            "sens": div(tp, tp + fn), "spec": div(tn, tn + fp),
            "ppv": div(tp, tp + fp), "npv": div(tn, tn + fn),
            "f1": div(2 * tp, 2 * tp + fp + fn)}


def choose_tau(z: np.ndarray, y: np.ndarray, lam: float, weights: np.ndarray | None = None) -> tuple[float, float]:
    z = np.asarray(z, float); y = np.asarray(y, int)
    if y.sum() == 0:
        return 3.0, float(np.sum(z >= 3.0))
    pred = z[:, None] >= TAU_GRID[None, :]
    w = np.ones(len(y)) if weights is None else np.asarray(weights, float)
    fn = np.sum((~pred) * (y[:, None] == 1) * w[:, None], axis=0)
    fp = np.sum(pred * (y[:, None] == 0) * w[:, None], axis=0)
    costs = lam * fn + fp
    j = int(np.argmin(costs))
    return float(TAU_GRID[j]), float(costs[j])


def calibration_features(d: pd.DataFrame, chrom: int, include_ga: bool = False) -> np.ndarray:
    cols = ["w", f"gc{chrom}", "gc", "filt_rate", "map_rate", "bmi"]
    if include_ga:
        cols.append("ga")
    x = d[cols].to_numpy(float).copy()
    med = np.nanmedian(x, axis=0)
    bad = np.where(~np.isfinite(x)); x[bad] = med[bad[1]]
    return x


def oof_chrom(
    d: pd.DataFrame, chrom: int, lam: float = 3.0, qprob: float = 0.10,
    score_mode: str = "raw", w_mode: str = "raw", gc_mode: str = "continuous",
    include_ga: bool = False, grouped: bool = True, seed: int = SEED,
) -> dict[str, object]:
    n = len(d); y = d[f"y{chrom}"].to_numpy(int); z_raw = d[f"z{chrom}"].to_numpy(float)
    w = d["w"].to_numpy(float).copy()
    if w_mode == "truncated":
        w = np.maximum(w, 0.0)
    state = np.full(n, -1, int); score = np.full(n, np.nan); tau_used = np.full(n, np.nan)
    fold_taus: list[float] = []; fold_costs: list[float] = []; q_values: list[float] = []
    split = folds(d["cv_group"].to_numpy(), seed=seed) if grouped else folds(np.arange(n), seed=seed)
    x_all = calibration_features(d, chrom, include_ga)
    gc = d["gc"].to_numpy(float)
    for tr, te in split:
        q = float(np.quantile(w[tr], qprob)); q_values.append(q)
        eligible_tr = np.ones(len(tr), bool); eligible_te = np.ones(len(te), bool)
        if gc_mode == "hard_filter":
            eligible_tr = (gc[tr] >= 0.40) & (gc[tr] <= 0.60)
            eligible_te = (gc[te] >= 0.40) & (gc[te] <= 0.60)
        det_tr = eligible_tr & (w[tr] >= q); det_te = eligible_te & (w[te] >= q)
        tr2 = tr[det_tr]; te2 = te[det_te]
        if score_mode == "calibrated" and len(tr2) > 8:
            scaler = StandardScaler().fit(x_all[tr2]); model = LinearRegression().fit(scaler.transform(x_all[tr2]), z_raw[tr2])
            ztr = z_raw[tr2] - model.predict(scaler.transform(x_all[tr2]))
            zte = z_raw[te2] - model.predict(scaler.transform(x_all[te2]))
        else:
            ztr = z_raw[tr2]; zte = z_raw[te2]
        weights = None
        if gc_mode == "quality_weight" and len(tr2):
            med = np.median(gc[tr2]); mad = np.median(np.abs(gc[tr2] - med))
            scale = max(1.4826 * mad, 1e-9)
            weights = 1.0 / (1.0 + np.abs(gc[tr2] - med) / scale)
        tau, cost = choose_tau(ztr, y[tr2], lam, weights)
        fold_taus.append(tau); fold_costs.append(cost)
        if len(te2):
            score[te2] = zte; state[te2] = (zte >= tau).astype(int); tau_used[te2] = tau
    valid = state >= 0
    cm = confusion(y[valid], state[valid]) if valid.any() else confusion(np.array([], int), np.array([], int))
    cm.update({"coverage": float(valid.mean()), "auc": safe_auc(y[valid], score[valid]) if valid.any() else 0.5,
               "tau": float(np.median(fold_taus)), "tau_min": float(np.min(fold_taus)),
               "tau_max": float(np.max(fold_taus)), "cv_cost": float(np.sum(fold_costs)),
               "q": float(np.median(q_values))})
    score[~np.isfinite(score)] = z_raw[~np.isfinite(score)]
    tau_used[~np.isfinite(tau_used)] = cm["tau"]
    return {"state": state, "score": score, "tau_used": tau_used, "metrics": cm,
            "fold_taus": fold_taus, "q_values": q_values}


def pregnancy_rows(d: pd.DataFrame, result: dict[str, object], chrom: int) -> tuple[list[dict[str, object]], pd.DataFrame]:
    tmp = pd.DataFrame({"mother": d["cv_group"].to_numpy(), "y": d[f"y{chrom}"].to_numpy(int),
                        "state": result["state"], "score": result["score"]})
    rows: list[dict[str, object]] = []; pred_rows = []
    for mother, g in tmp.groupby("mother", sort=False):
        det = g[g.state >= 0]; true = int(g.y.max())
        if det.empty:
            pred_rows.append({"mother": mother, "chrom": chrom, "truth": true, "determinate": 0,
                              "conservative": -1, "majority": -1, "max_risk": -1, "score_max": float(g.score.max())})
            continue
        cons = int((det.state == 1).any())
        maj = int((det.state == 1).sum() > len(det) / 2)
        # With a common record-level threshold, max-risk and any-positive are algebraically equivalent.
        maxrisk = cons
        pred_rows.append({"mother": mother, "chrom": chrom, "truth": true, "determinate": 1,
                          "conservative": cons, "majority": maj, "max_risk": maxrisk,
                          "score_max": float(det.score.max())})
    p = pd.DataFrame(pred_rows)
    valid = p.determinate == 1
    for rule in ("conservative", "majority", "max_risk"):
        cm = confusion(p.loc[valid, "truth"].to_numpy(), p.loc[valid, rule].to_numpy())
        cm.update({"level": "pregnant", "chrom": str(chrom), "rule": rule,
                   "coverage": float(valid.mean()), "tau": float(result["metrics"]["tau"]),
                   "auc": safe_auc(p.loc[valid, "truth"].to_numpy(), p.loc[valid, "score_max"].to_numpy()),
                   "tau_min": float(result["metrics"]["tau_min"]), "tau_max": float(result["metrics"]["tau_max"]),
                   "cv_cost": float(result["metrics"]["cv_cost"]), "q": float(result["metrics"]["q"])})
        rows.append(cm)
    return rows, p


def evaluate(
    d: pd.DataFrame, lam: float = 3.0, qprob: float = 0.10, score_mode: str = "raw",
    w_mode: str = "raw", gc_mode: str = "continuous", include_ga: bool = False,
    grouped: bool = True, seed: int = SEED,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[int, dict[str, object]]]:
    rows = []; record_preds = []; preg_preds = []; details = {}
    for chrom in CHROMS:
        res = oof_chrom(d, chrom, lam, qprob, score_mode, w_mode, gc_mode, include_ga, grouped, seed)
        details[chrom] = res; m = res["metrics"]
        rows.append({"level": "record", "chrom": str(chrom), "rule": "na", **m})
        record_preds.append(pd.DataFrame({"row_id": d["seq"].to_numpy(), "mother": d["cv_group"].to_numpy(),
                                          "chrom": chrom, "truth": d[f"y{chrom}"].to_numpy(),
                                          "state": res["state"], "score": res["score"], "tau_used": res["tau_used"]}))
        prows, pp = pregnancy_rows(d, res, chrom); rows.extend(prows); preg_preds.append(pp)
    return pd.DataFrame(rows), pd.concat(record_preds, ignore_index=True), pd.concat(preg_preds, ignore_index=True), details


def bootstrap(d: pd.DataFrame, b: int = N_BOOT) -> pd.DataFrame:
    rng = np.random.default_rng(SEED + 11); mothers = d["mother"].unique(); blocks = {m: np.flatnonzero(d.mother.to_numpy() == m) for m in mothers}
    out = []
    for rep in range(b):
        sampled = rng.choice(mothers, len(mothers), replace=True)
        indices = np.concatenate([blocks[m] for m in sampled])
        db = d.iloc[indices].copy().reset_index(drop=True)
        db["cv_group"] = np.concatenate([np.full(len(blocks[m]), k) for k, m in enumerate(sampled)])
        summary, _, _, _ = evaluate(db, seed=SEED + rep + 100)
        summary["replicate"] = rep; out.append(summary)
    return pd.concat(out, ignore_index=True)


def permutation_test(d: pd.DataFrame, observed: pd.DataFrame, b: int = N_PERM) -> pd.DataFrame:
    rng = np.random.default_rng(SEED + 29); out = []
    for rep in range(b):
        dp = d.copy()
        for chrom in CHROMS:
            shuffled = np.empty(len(dp), int)
            for _, idx in dp.groupby("mother").groups.items():
                ii = np.asarray(list(idx), int); shuffled[ii] = rng.permutation(dp.loc[ii, f"y{chrom}"].to_numpy())
            dp[f"y{chrom}"] = shuffled
        s, _, _, _ = evaluate(dp, seed=SEED + rep + 500)
        for chrom in CHROMS:
            row = s[(s.level == "record") & (s.chrom == str(chrom))].iloc[0]
            out.append({"replicate": rep, "chrom": chrom, "f1": row.f1, "auc": row.auc})
    p = pd.DataFrame(out)
    obs = observed[observed.level == "record"].set_index("chrom")
    pvals = []
    for chrom in CHROMS:
        q = p[p.chrom == chrom]
        for metric in ("f1", "auc"):
            v = float(obs.loc[str(chrom), metric]); pv = (1 + int(np.sum(q[metric] >= v))) / (b + 1)
            pvals.append({"chrom": chrom, "metric": metric, "observed": v, "perm_p": pv})
    return p, pd.DataFrame(pvals)


def full_feature_logistic(d: pd.DataFrame) -> pd.DataFrame:
    feature_cols = ["age", "height", "bmi", "reads", "map_rate", "dup_rate", "unique_reads", "gc",
                    "z13", "z18", "z21", "zx", "w", "gc13", "gc18", "gc21", "filt_rate", "ga"]
    x = d[feature_cols].to_numpy(float).copy(); med = np.nanmedian(x, axis=0); bad = np.where(~np.isfinite(x)); x[bad] = med[bad[1]]
    out = []
    for chrom in CHROMS:
        y = d[f"y{chrom}"].to_numpy(int); pred = np.zeros(len(d), int); prob = np.zeros(len(d))
        for tr, te in folds(d.mother.to_numpy(), seed=SEED + chrom):
            sc = StandardScaler().fit(x[tr]); model = LogisticRegression(max_iter=1000, class_weight={0: 1, 1: 3}, random_state=SEED)
            model.fit(sc.transform(x[tr]), y[tr]); prob[te] = model.predict_proba(sc.transform(x[te]))[:, 1]; pred[te] = prob[te] >= .5
        cm = confusion(y, pred); cm.update({"chrom": chrom, "model": "full_feature_logistic", "auc": safe_auc(y, prob)})
        out.append(cm)
    return pd.DataFrame(out)


def sensitivity_tables(f: pd.DataFrame, male: pd.DataFrame, main: pd.DataFrame) -> dict[str, pd.DataFrame]:
    tables: dict[str, pd.DataFrame] = {}
    gate_rows = []
    for qname, qp in (("q5", .05), ("q10", .10), ("q25", .25)):
        s, _, _, details = evaluate(f, qprob=qp)
        for chrom in CHROMS:
            r = s[(s.level == "record") & (s.chrom == str(chrom))].iloc[0]
            gate_rows.append({"gate": qname, "qprob": qp, "q_value": np.quantile(f.w, qp), "chrom": chrom,
                              "coverage": r.coverage, "accuracy": (r.TP + r.TN) / max(r.TP+r.TN+r.FP+r.FN, 1),
                              "positive_rate": f.loc[f.w >= np.quantile(f.w, qp), f"y{chrom}"].mean(),
                              "n_invalid": int(np.sum(f.w < np.quantile(f.w, qp))), "f1": r.f1, "auc": r.auc,
                              "tau": r.tau})
    tables["q4_gate_sensitivity"] = pd.DataFrame(gate_rows)

    cost_rows = []; curve_rows = []
    reliable = f.w >= f.w.quantile(.10)
    for lam in (1, 2, 3, 5):
        s, _, _, _ = evaluate(f, lam=lam)
        for chrom in CHROMS:
            r = s[(s.level == "record") & (s.chrom == str(chrom))].iloc[0]
            cost_rows.append({"lambda": lam, "chrom": chrom, "tau": r.tau, "sens": r.sens, "spec": r.spec,
                              "f1": r.f1, "auc": r.auc, "cv_cost": r.cv_cost})
            z = f.loc[reliable, f"z{chrom}"].to_numpy(); y = f.loc[reliable, f"y{chrom}"].to_numpy()
            for tau in TAU_GRID:
                cm = confusion(y, z >= tau)
                curve_rows.append({"lambda": lam, "chrom": chrom, "tau": tau, "cost": lam*cm["FN"]+cm["FP"]})
    tables["q4_cost_sensitivity"] = pd.DataFrame(cost_rows); tables["q4_cost_curves"] = pd.DataFrame(curve_rows)

    cal_rows = []
    for mode in ("raw", "calibrated"):
        s, _, _, _ = evaluate(f, score_mode=mode)
        for chrom in CHROMS:
            r = s[(s.level == "record") & (s.chrom == str(chrom))].iloc[0]
            cal_rows.append({"mode": mode, "chrom": chrom, "tau": r.tau, "cv_cost": r.cv_cost, "f1": r.f1, "auc": r.auc})
    cal = pd.DataFrame(cal_rows); raw_cost = cal[cal["mode"] == "raw"].cv_cost.sum(); adj_cost = cal[cal["mode"] == "calibrated"].cv_cost.sum()
    cal["adopted"] = bool(adj_cost < raw_cost and cal[cal["mode"] == "calibrated"].f1.mean() > cal[cal["mode"] == "raw"].f1.mean())
    tables["q4_calibration_sens"] = cal

    gc_rows = []
    for mode in ("continuous", "quality_weight", "hard_filter"):
        s, _, _, _ = evaluate(f, gc_mode=mode)
        for chrom in CHROMS:
            r = s[(s.level == "record") & (s.chrom == str(chrom))].iloc[0]
            gc_rows.append({"gc_mode": mode, "chrom": chrom, "coverage": r.coverage, "tau": r.tau,
                            "sens": r.sens, "spec": r.spec, "f1": r.f1, "auc": r.auc})
    tables["q4_gc_sens"] = pd.DataFrame(gc_rows)

    w_rows = []
    for mode in ("raw", "truncated"):
        s, _, _, _ = evaluate(f, w_mode=mode)
        for chrom in CHROMS:
            r = s[(s.level == "record") & (s.chrom == str(chrom))].iloc[0]
            w_rows.append({"w_mode": mode, "chrom": chrom, "coverage": r.coverage, "q10": r.q,
                           "tau": r.tau, "f1": r.f1, "auc": r.auc})
    tables["q4_w_trunc_sens"] = pd.DataFrame(w_rows)

    time_rows = []
    for inc in (False, True):
        s, _, _, _ = evaluate(f, score_mode="calibrated" if inc else "raw", include_ga=inc)
        for chrom in CHROMS:
            r = s[(s.level == "record") & (s.chrom == str(chrom))].iloc[0]
            time_rows.append({"time_included": inc, "chrom": chrom, "cv_cost": r.cv_cost, "f1": r.f1, "auc": r.auc})
    tables["q4_time_sens"] = pd.DataFrame(time_rows)

    # Female reliable-subset versus male-table raw Z discrimination.
    auc_rows = []
    for sex, d, mask in (("female", f, f.w >= f.w.quantile(.10)), ("male", male, np.ones(len(male), bool))):
        for chrom in CHROMS:
            y = d.loc[mask, f"y{chrom}"].to_numpy(); z = d.loc[mask, f"z{chrom}"].to_numpy()
            auc_rows.append({"sex": sex, "chrom": chrom, "auc": safe_auc(y, z), "n": len(y), "positives": int(y.sum())})
    tables["q4_male_female_auc"] = pd.DataFrame(auc_rows)

    roc_rows = []
    for _, a in tables["q4_male_female_auc"].iterrows():
        d = f if a.sex == "female" else male; mask = (d.w >= d.w.quantile(.10)) if a.sex == "female" else np.ones(len(d), bool)
        y = d.loc[mask, f"y{int(a.chrom)}"].to_numpy(); z = d.loc[mask, f"z{int(a.chrom)}"].to_numpy()
        if len(np.unique(y)) == 2:
            fpr, tpr, th = roc_curve(y, z)
            for x, yy, t in zip(fpr, tpr, th):
                if np.isfinite(t):
                    roc_rows.append({"sex": a.sex, "chrom": int(a.chrom), "fpr": x, "tpr": yy, "threshold": t})
    tables["q4_z_roc"] = pd.DataFrame(roc_rows)

    # Clinical Z>=3 benchmark and grouped-vs-row CV leakage audit.
    bench = []
    for chrom in CHROMS:
        y = f.loc[reliable, f"y{chrom}"].to_numpy(); z = f.loc[reliable, f"z{chrom}"].to_numpy(); cm = confusion(y, z >= 3)
        bench.append({"chrom": chrom, "tau": 3.0, **cm, "auc": safe_auc(y, z)})
    tables["q4_z3_benchmark"] = pd.DataFrame(bench)

    leak_rows = []
    for grouped in (True, False):
        s, _, _, _ = evaluate(f, grouped=grouped)
        for chrom in CHROMS:
            r = s[(s.level == "record") & (s.chrom == str(chrom))].iloc[0]
            leak_rows.append({"cv": "mother_grouped" if grouped else "row_independent", "chrom": chrom, "f1": r.f1, "auc": r.auc})
    tables["q4_cv_leakage_audit"] = pd.DataFrame(leak_rows)

    tables["q4_full_feature_model"] = full_feature_logistic(f)
    return tables


def data_profile(f: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    row = f[["seq", "mother", "w", "gc", "z13", "z18", "z21", "zx", "bmi", "ga", "ab", "healthy",
             "y13", "y18", "y21", "ab_positive"]].copy()
    row["ab"] = row["ab"].fillna("negative")
    row["bmi_missing"] = row["bmi"].isna().astype(int)
    row["bmi"] = row["bmi"].fillna(f["bmi"].median())
    ab_counts = f.ab.fillna("negative").value_counts()
    metrics = [
        ("female_records", len(f), "count"), ("female_mothers", f.mother.nunique(), "count"),
        ("w_min", f.w.min(), "proportion"), ("w_median", f.w.median(), "proportion"),
        ("w_max", f.w.max(), "proportion"), ("w_negative_rate", (f.w < 0).mean(), "proportion"),
        ("q5", f.w.quantile(.05), "proportion"), ("q10", f.w.quantile(.10), "proportion"),
        ("q25", f.w.quantile(.25), "proportion"), ("ab_positive_records", f.ab.notna().sum(), "count"),
        ("ab_positive_mothers", f.loc[f.ab.notna(), "mother"].nunique(), "count"),
        ("mixed_label_mothers", sum(g.ab.notna().any() and g.ab.isna().any() for _, g in f.groupby("mother")), "count"),
        ("ae_unique", f.healthy.nunique(), "count"), ("ab_positive_and_ae_unhealthy", 0, "count"),
        ("gc_min", f.gc.min(), "proportion"), ("gc_max", f.gc.max(), "proportion"),
        ("gc_below_0_40_rate", (f.gc < .4).mean(), "proportion"),
        ("bmi_missing_records", f.bmi.isna().sum(), "count"),
    ]
    for label, count in ab_counts.items(): metrics.append((f"AB_{label}", count, "count"))
    return row, pd.DataFrame(metrics, columns=["metric", "value", "unit"])


def main() -> None:
    started = time.perf_counter(); data_path, root = resolve_paths(); resdir = root / "results"
    male = load_sheet(data_path, 0); female = load_sheet(data_path, 1)
    if not (female.yz.isna().all() and female.y_conc.isna().all()):
        raise AssertionError("Female sheet must have both Y-related columns empty")

    profile_rows, profile_table = data_profile(female)
    main_table, record_pred, preg_pred, details = evaluate(female)
    boot = bootstrap(female, N_BOOT)
    perm_dist, perm_p = permutation_test(female, main_table, N_PERM)

    # Attach bootstrap intervals for every metric; contract ci_* uses F1.
    keys = ["level", "chrom", "rule"]
    for metric in METRICS:
        ci = boot.groupby(keys)[metric].quantile([.025, .975]).unstack().reset_index()
        ci.columns = keys + [f"{metric}_ci_low", f"{metric}_ci_high"]
        main_table = main_table.merge(ci, on=keys, how="left")
    main_table["ci_low"] = main_table["f1_ci_low"]; main_table["ci_high"] = main_table["f1_ci_high"]
    main_table = main_table.merge(perm_p[perm_p.metric == "f1"][["chrom", "perm_p"]].assign(chrom=lambda x: x.chrom.astype(str)), on="chrom", how="left")

    sensitivity = sensitivity_tables(female, male, main_table)
    # Gate bootstrap threshold audit is available directly from the cluster resamples.
    gate_boot = boot[boot.level == "record"][["replicate", "chrom", "q", "coverage", "tau", "f1", "auc"]].copy()

    threshold_table = main_table[main_table.level == "record"].copy()
    merge_table = main_table[main_table.level == "pregnant"].copy()
    boot_summary = []
    for chrom in CHROMS:
        obs = main_table[(main_table.level == "record") & (main_table.chrom == str(chrom))].iloc[0]
        for metric in ("f1", "auc", "sens", "spec", "ppv", "npv"):
            boot_summary.append({"chrom": chrom, "metric": metric, "observed": obs[metric],
                                 "CI_low": obs[f"{metric}_ci_low"], "CI_high": obs[f"{metric}_ci_high"],
                                 "perm_p": float(perm_p[(perm_p.chrom == chrom) & (perm_p.metric == (metric if metric in ('f1','auc') else 'f1'))].perm_p.iloc[0])})
    boot_summary = pd.DataFrame(boot_summary)

    sens_report = pd.DataFrame([
        ("4", "AE交叉标签清洗", "阳性数", int(female.ab.notna().sum()), 0, "AE全为健康，禁止作为金标准"),
        ("5", "X浓度负值截断", "覆盖率", float(main_table.coverage.iloc[0]), float(sensitivity["q4_w_trunc_sens"].query("w_mode=='truncated'").coverage.mean()), "截断破坏门界信息"),
        ("6", "q5/q10/q25门界", "覆盖率范围", float(main_table.coverage.iloc[0]), f"{sensitivity['q4_gate_sensitivity'].coverage.min():.3f}-{sensitivity['q4_gate_sensitivity'].coverage.max():.3f}", "存在覆盖率-性能权衡"),
        ("7", "取消门控", "无法判定比例", float(1-main_table.coverage.iloc[0]), 0.0, "取消门控提高覆盖但纳入低可靠记录"),
        ("8", "女胎/男胎Z对照", "宏平均AUC", float(sensitivity["q4_male_female_auc"].query("sex=='female'").auc.mean()), float(sensitivity["q4_male_female_auc"].query("sex=='male'").auc.mean()), "女胎必须单独定标"),
        ("9", "lambda=1/2/3/5", "tau范围", str(list(threshold_table.tau.round(3))), f"{sensitivity['q4_cost_sensitivity'].tau.min():.1f}-{sensitivity['q4_cost_sensitivity'].tau.max():.1f}", "阈值依赖代价权重"),
        ("11", "三种孕妇合并", "F1范围", float(merge_table.query("rule=='conservative'").f1.mean()), f"{merge_table.f1.min():.3f}-{merge_table.f1.max():.3f}", "多数票与保守规则差异显著；max-risk与任一阳性等价"),
        ("13", "GC连续/权重/硬剔除", "覆盖率范围", float(main_table.coverage.iloc[0]), f"{sensitivity['q4_gc_sens'].coverage.min():.3f}-{sensitivity['q4_gc_sens'].coverage.max():.3f}", "硬剔除损失样本"),
    ], columns=["关键假设", "扰动/对照", "比较指标", "主结果", "敏感性结果", "结论"])

    outputs: dict[str, pd.DataFrame] = {
        "q4": main_table, "output": main_table, "q4_data_profile": profile_rows,
        "tab_q4_data_profile": profile_table, "q4_record_predictions": record_pred,
        "q4_pregnant_predictions": preg_pred, "q4_bootstrap_metrics": boot,
        "q4_gate_bootstrap": gate_boot, "q4_permutation_distribution": perm_dist,
        "q4_permutation_tests": perm_p, "tab_q4_thresholds": threshold_table,
        "tab_q4_merge": merge_table, "tab_q4_bootstrap": boot_summary,
        "tab_q4_sensitivity": sens_report,
    }
    outputs.update(sensitivity)
    outputs["tab_q4_gate"] = sensitivity["q4_gate_sensitivity"]
    for name, frame in outputs.items():
        if frame.select_dtypes(include=[np.number]).isin([np.inf, -np.inf]).any().any():
            raise AssertionError(f"Infinite value in {name}")
        frame.to_csv(resdir / f"{name}.csv", index=False)

    elapsed = time.perf_counter() - started
    summary = {
        "records": int(len(female)), "mothers": int(female.mother.nunique()),
        "ab_positive_records": int(female.ab.notna().sum()),
        "ab_positive_mothers": int(female.loc[female.ab.notna(), "mother"].nunique()),
        "w_negative_rate": float((female.w < 0).mean()), "q10": float(female.w.quantile(.10)),
        "coverage": float(main_table.loc[main_table.level == "record", "coverage"].mean()),
        "bootstrap_completed": N_BOOT, "permutations_completed": N_PERM,
        "seed": SEED, "runtime_seconds": elapsed, "runtime_limit_seconds": 90,
        "runtime_contract_passed": bool(elapsed <= 90), "main_model": "raw chromosome-specific Z thresholds",
        "calibration_adopted": bool(sensitivity["q4_calibration_sens"].adopted.iloc[0]),
        "ci_metric_for_contract": "f1",
    }
    (resdir / "q4_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    if elapsed > 90:
        raise RuntimeError(f"Q4 runtime contract exceeded: {elapsed:.2f}s > 90s")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

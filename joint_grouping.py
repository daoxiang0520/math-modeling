"""Joint BMI grouping and group-time selection for Q2.

The individual target ``tp`` is the earliest week at which the marginal
probability of exceeding 4% reaches ``p_guarantee``.  A common group time is
chosen with asymmetric absolute loss.  With early:late cost ratio 4:1, the
minimizer is the 0.80 quantile, so the group rule targets 80% member coverage
without confusing that coverage with the individual 0.80 probability level.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class JointSolution:
    K: int
    boundaries: tuple[float, ...]
    times: tuple[float, ...]
    counts: tuple[int, ...]
    coverages: tuple[float, ...]
    total_loss: float
    mean_loss: float
    valid_time_separation: bool


def assign_groups(bmi: np.ndarray, boundaries: tuple[float, ...] | list[float]) -> np.ndarray:
    return np.digitize(np.asarray(bmi, float), np.asarray(boundaries, float), right=False)


def _group_time_loss(tp: np.ndarray, coverage: float, early_weight: float,
                     late_weight: float, step: float) -> tuple[float, float, float]:
    values = np.asarray(tp, float)
    try:
        raw = float(np.quantile(values, coverage, method="higher"))
    except TypeError:  # numpy < 1.22
        raw = float(np.quantile(values, coverage, interpolation="higher"))
    week = float(np.clip(np.ceil((raw - 1e-12) / step) * step, 10.0, 25.0))
    loss = early_weight * np.maximum(values - week, 0.0) + late_weight * np.maximum(week - values, 0.0)
    return week, float(np.sum(loss)), float(np.mean(values <= week + 1e-12))


def _rounded_boundaries(raw: list[float], bmi_sorted: np.ndarray, n_min: int) -> tuple[float, ...]:
    if not raw:
        return ()
    candidates = np.arange(15.0, 50.01, 0.5)
    chosen: list[float] = []
    previous_index = 0
    for slot, value in enumerate(raw):
        remaining = len(raw) - slot
        pick = None
        for q in sorted(candidates, key=lambda x: (abs(x - value), x)):
            idx = int(np.searchsorted(bmi_sorted, q, side="left"))
            if (idx - previous_index >= n_min and
                    len(bmi_sorted) - idx >= remaining * n_min and
                    (not chosen or q > chosen[-1])):
                pick = float(q)
                previous_index = idx
                break
        if pick is None:
            return ()
        chosen.append(pick)
    return tuple(chosen)


def solve_joint(tp: np.ndarray, bmi: np.ndarray, k_max: int = 5, n_min: int = 30,
                coverage: float = 0.80, early_weight: float = 4.0,
                late_weight: float = 1.0, min_time_gap: float = 0.5,
                step: float = 0.1) -> dict[int, JointSolution]:
    """Dynamic program over contiguous BMI intervals for K=1..k_max."""
    tp = np.asarray(tp, float)
    bmi = np.asarray(bmi, float)
    order = np.argsort(bmi, kind="mergesort")
    bs, ts = bmi[order], tp[order]
    n = len(bs)
    # Clinically executable half-BMI cut grid.  Collapsing grid values that
    # induce the same split reduces a 267x267 segment table to roughly 30x30.
    split_to_boundary: dict[int, float] = {}
    for q in np.arange(15.0, 50.01, 0.5):
        idx = int(np.searchsorted(bs, q, side="left"))
        if idx <= 0 or idx >= n or bs[idx - 1] >= bs[idx]:
            continue
        midpoint = (bs[idx - 1] + bs[idx]) / 2
        old = split_to_boundary.get(idx)
        if old is None or abs(q - midpoint) < abs(old - midpoint):
            split_to_boundary[idx] = float(q)
    positions = np.array([0] + sorted(split_to_boundary) + [n], int)
    m = len(positions)
    cost = np.full((m, m), np.inf)
    opt_time = np.full((m, m), np.nan)
    for a in range(m - 1):
        i = int(positions[a])
        for c in range(a + 1, m):
            j = int(positions[c])
            if j - i < n_min:
                continue
            week, value, _ = _group_time_loss(ts[i:j], coverage, early_weight, late_weight, step)
            cost[a, c] = value
            opt_time[a, c] = week

    dp = np.full((k_max + 1, m), np.inf)
    previous = np.full((k_max + 1, m), -1, int)
    dp[0, 0] = 0.0
    for k in range(1, k_max + 1):
        for c in range(1, m):
            if positions[c] < k * n_min:
                continue
            for a in range(c):
                value = dp[k - 1, a] + cost[a, c]
                if value < dp[k, c]:
                    dp[k, c] = value
                    previous[k, c] = a

    solutions: dict[int, JointSolution] = {}
    for k in range(1, k_max + 1):
        if not np.isfinite(dp[k, m - 1]):
            continue
        cut_slots = [m - 1]
        cursor = m - 1
        for kk in range(k, 0, -1):
            cursor = int(previous[kk, cursor])
            if cursor < 0:
                break
            cut_slots.append(cursor)
        if cursor != 0:
            continue
        cut_positions = [int(positions[c]) for c in sorted(cut_slots)]
        boundaries = tuple(split_to_boundary[c] for c in cut_positions[1:-1])
        if len(boundaries) != k - 1:
            continue
        gid = assign_groups(bs, boundaries)
        times: list[float] = []
        counts: list[int] = []
        coverages: list[float] = []
        total = 0.0
        for g in range(k):
            mask = gid == g
            week, value, group_coverage = _group_time_loss(
                ts[mask], coverage, early_weight, late_weight, step
            )
            times.append(week)
            counts.append(int(mask.sum()))
            coverages.append(group_coverage)
            total += value
        separated = k == 1 or bool(np.all(np.diff(times) >= min_time_gap - 1e-12))
        solutions[k] = JointSolution(
            K=k, boundaries=boundaries, times=tuple(times), counts=tuple(counts),
            coverages=tuple(coverages), total_loss=float(total), mean_loss=float(total / n),
            valid_time_separation=separated,
        )
    return solutions


def bootstrap_selection(tp: np.ndarray, bmi: np.ndarray, repeats: int = 300,
                        seed: int = 2026, **kwargs) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    n = len(tp)
    for replicate in range(repeats):
        idx = rng.integers(0, n, n)
        sols = solve_joint(np.asarray(tp)[idx], np.asarray(bmi)[idx], **kwargs)
        for k in range(1, int(kwargs.get("k_max", 5)) + 1):
            sol = sols.get(k)
            rows.append({
                "replicate": replicate,
                "K": k,
                "mean_loss": np.nan if sol is None else sol.mean_loss,
                "boundaries": "" if sol is None else "|".join(f"{x:.1f}" for x in sol.boundaries),
                "times": "" if sol is None else "|".join(f"{x:.1f}" for x in sol.times),
                "valid_time_separation": False if sol is None else sol.valid_time_separation,
            })
    return pd.DataFrame(rows)


def select_one_se(full: dict[int, JointSolution], boot: pd.DataFrame,
                  min_valid_rate: float = 0.80) -> tuple[int, pd.DataFrame]:
    rows = []
    for k, z in boot.groupby("K"):
        finite = z[np.isfinite(z["mean_loss"])]
        rows.append({
            "K": int(k),
            "bootstrap_mean_loss": float(finite["mean_loss"].mean()),
            "bootstrap_se": float(finite["mean_loss"].std(ddof=1)),
            "valid_rate": float(finite["valid_time_separation"].mean()),
            "successful_replicates": int(len(finite)),
        })
    summary = pd.DataFrame(rows).sort_values("K").reset_index(drop=True)
    eligible = summary[
        summary["K"].isin([k for k, sol in full.items() if sol.valid_time_separation]) &
        (summary["valid_rate"] >= min_valid_rate)
    ]
    if eligible.empty:
        selected = 1
        threshold = float("nan")
    else:
        best = eligible.loc[eligible["bootstrap_mean_loss"].idxmin()]
        threshold = float(best["bootstrap_mean_loss"] + best["bootstrap_se"])
        selected = int(eligible.loc[eligible["bootstrap_mean_loss"] <= threshold, "K"].min())
    summary["one_se_threshold"] = threshold
    summary["selected"] = summary["K"].eq(selected)
    return selected, summary

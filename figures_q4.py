"""Generate all registered Question 4 scientific figures."""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(os.getenv("MODELING_OUTPUT_DIR", Path(__file__).resolve().parent))
RES = ROOT / "results"; OUT = ROOT / "figures"; OUT.mkdir(parents=True, exist_ok=True)
COLORS = {13: "#0072B2", 18: "#D55E00", 21: "#009E73"}
MARKERS = {13: "o", 18: "s", 21: "^"}
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
    "axes.unicode_minus": False, "font.size": 9, "axes.titlesize": 11, "axes.labelsize": 9.5,
    "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 130, "savefig.dpi": 320,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def save(fig: plt.Figure, name: str, left: float = .11, bottom: float = .15) -> None:
    fig.subplots_adjust(left=left, right=.97, bottom=bottom, top=.89, wspace=.30, hspace=.34)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{name}.{ext}", facecolor="white")
    plt.close(fig)


def panel(ax, label: str) -> None:
    ax.text(-.12, 1.03, label, transform=ax.transAxes, fontweight="bold", va="bottom")


def data_profile() -> None:
    d = pd.read_csv(RES / "q4_data_profile.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.55))
    ax = axes[0]
    groups = [d.loc[d.ab_positive == k, "w"].to_numpy() for k in (0, 1)]
    bp = ax.boxplot(groups, tick_labels=["AB negative\n(n=538)", "AB positive\n(n=67)"], patch_artist=True,
                    showfliers=False, widths=.52)
    for b, c in zip(bp["boxes"], ["#56B4E9", "#E69F00"]): b.set_facecolor(c); b.set_alpha(.65)
    rng = np.random.default_rng(2025)
    for i, vals in enumerate(groups, 1):
        ax.scatter(rng.normal(i, .055, len(vals)), vals, s=8, alpha=.22, color="0.25", linewidths=0)
    ax.axhline(0, color="0.4", ls=":", lw=1)
    ax.set(ylabel="X-chromosome concentration", title="X concentration retains negative values")
    ax.text(.02, .96, f"Negative: {(d.w < 0).mean():.1%}\nAE healthy: {len(d)}/{len(d)}", transform=ax.transAxes, va="top")
    panel(ax, "a")
    ax = axes[1]
    counts = d.ab.value_counts().sort_values()
    ax.barh(np.arange(len(counts)), counts.values, color=["#999999" if x == "negative" else "#CC79A7" for x in counts.index])
    ax.set_yticks(np.arange(len(counts)), counts.index); ax.set(xlabel="Number of records", title="AB multilabel composition")
    for i, v in enumerate(counts.values): ax.text(v + 3, i, str(v), va="center", fontsize=7.5)
    panel(ax, "b")
    save(fig, "fig_q4_data_profile")


def gate_curve() -> None:
    d = pd.read_csv(RES / "q4_gate_sensitivity.csv")
    g = d.groupby(["gate", "qprob"], as_index=False).agg(coverage=("coverage", "mean"), accuracy=("accuracy", "mean"),
                                                            f1=("f1", "mean"))
    g = g.sort_values("qprob")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
    axes[0].plot(g.qprob * 100, g.coverage, "o-", color="#0072B2", lw=1.8)
    axes[1].plot(g.qprob * 100, g.accuracy, "s-", color="#D55E00", lw=1.8)
    for ax, ylab, title in zip(axes, ["Determinate-record coverage", "Macro record accuracy"],
                               ["Coverage decreases with stricter gate", "Accuracy among determinate records"]):
        ax.set(xticks=g.qprob*100, ylim=(0, 1.03), xlabel="Lower-tail gate quantile (%)", ylabel=ylab, title=title)
    panel(axes[0], "a"); panel(axes[1], "b")
    save(fig, "fig_q4_gate_curve")


def z_roc() -> None:
    d = pd.read_csv(RES / "q4_z_roc.csv"); auc = pd.read_csv(RES / "q4_male_female_auc.csv")
    fig, ax = plt.subplots(figsize=(5.4, 4.35))
    for chrom in (13, 18, 21):
        q = d[(d.sex == "female") & (d.chrom == chrom)]
        av = auc[(auc.sex == "female") & (auc.chrom == chrom)].auc.iloc[0]
        ax.plot(q.fpr, q.tpr, color=COLORS[chrom], marker=MARKERS[chrom], markevery=max(1, len(q)//7),
                ms=3, lw=1.6, label=f"T{chrom}, AUC={av:.3f}")
    ax.plot([0, 1], [0, 1], "--", color="0.45", lw=1, label="Random")
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="False-positive rate", ylabel="True-positive rate",
           title="Female-fetus Z scores have weak discrimination")
    ax.legend(frameon=False, loc="lower right")
    save(fig, "fig_q4_z_roc")


def cost_sensitivity() -> None:
    d = pd.read_csv(RES / "q4_cost_sensitivity.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.45))
    for chrom, q in d.groupby("chrom"):
        c = COLORS[int(chrom)]; mk = MARKERS[int(chrom)]
        axes[0].plot(q["lambda"], q.tau, marker=mk, color=c, lw=1.6, label=f"T{int(chrom)}")
        axes[1].plot(q["lambda"], q.sens, marker=mk, color=c, lw=1.6, label=f"T{int(chrom)}")
    axes[0].set(xlabel="False-negative cost weight λ", ylabel="Selected Z threshold", title="Threshold sensitivity")
    axes[1].set(xlabel="False-negative cost weight λ", ylabel="Out-of-fold sensitivity", title="Operating-point sensitivity", ylim=(-.02, 1.02))
    axes[0].legend(frameon=False); panel(axes[0], "a"); panel(axes[1], "b")
    save(fig, "fig_q4_cost_sens")


def merge_compare() -> None:
    d = pd.read_csv(RES / "q4.csv"); d = d[d.level == "pregnant"]
    rules = ["conservative", "majority", "max_risk"]
    metrics = ["sens", "spec", "npv", "f1"]
    mat = []
    labels = []
    for chrom in (13, 18, 21):
        for rule in rules:
            r = d[(d.chrom == chrom) & (d.rule == rule)].iloc[0]
            mat.append([r[m] for m in metrics]); labels.append(f"T{chrom} · {rule}")
    mat = np.asarray(mat)
    fig, ax = plt.subplots(figsize=(7.2, 4.25))
    im = ax.imshow(mat, aspect="auto", cmap="cividis", vmin=0, vmax=1)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]): ax.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center", color="white" if mat[i,j] < .35 else "black", fontsize=7.5)
    ax.set_xticks(range(len(metrics)), ["Sensitivity", "Specificity", "NPV", "F1"])
    ax.set_yticks(range(len(labels)), labels); ax.set_title("Pregnancy-level aggregation rules")
    cb = fig.colorbar(im, ax=ax, pad=.02); cb.set_label("Metric value")
    save(fig, "fig_q4_merge_compare", left=.20)


def calibration() -> None:
    d = pd.read_csv(RES / "q4_calibration_sens.csv")
    wide = d.pivot(index="chrom", columns="mode", values=["auc", "cv_cost"])
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
    for chrom in wide.index:
        c = COLORS[int(chrom)]; mk = MARKERS[int(chrom)]
        axes[0].plot([0, 1], [wide.loc[chrom, ("auc", "raw")], wide.loc[chrom, ("auc", "calibrated")]], marker=mk, color=c, label=f"T{int(chrom)}")
        axes[1].plot([0, 1], [wide.loc[chrom, ("cv_cost", "raw")], wide.loc[chrom, ("cv_cost", "calibrated")]], marker=mk, color=c)
    for ax in axes: ax.set_xticks([0, 1], ["Raw Z", "Residualized Z"])
    axes[0].set(ylabel="Out-of-fold AUC", title="Discrimination", ylim=(0, 1)); axes[0].legend(frameon=False)
    axes[1].set(ylabel="Cross-validated weighted cost", title="Decision cost")
    panel(axes[0], "a"); panel(axes[1], "b")
    save(fig, "fig_q4_calibration")


def gc_sensitivity() -> None:
    d = pd.read_csv(RES / "q4_gc_sens.csv"); modes = ["continuous", "quality_weight", "hard_filter"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.45))
    for chrom, q in d.groupby("chrom"):
        q = q.set_index("gc_mode").loc[modes].reset_index(); c = COLORS[int(chrom)]; mk = MARKERS[int(chrom)]
        axes[0].scatter(range(3), q.coverage, color=c, marker=mk, s=48, label=f"T{int(chrom)}")
        axes[1].scatter(range(3), q.auc, color=c, marker=mk, s=48)
    for ax in axes: ax.set_xticks(range(3), ["Continuous", "Robust weight", "Hard 40–60%"], rotation=10)
    axes[0].set(ylabel="Coverage", title="Sample retention", ylim=(0, 1.03)); axes[0].legend(frameon=False)
    axes[1].set(ylabel="Out-of-fold AUC", title="Discrimination", ylim=(0, 1.03))
    panel(axes[0], "a"); panel(axes[1], "b")
    save(fig, "fig_q4_gc_sens")


def w_truncation() -> None:
    d = pd.read_csv(RES / "q4_w_trunc_sens.csv"); modes = ["raw", "truncated"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
    for chrom, q in d.groupby("chrom"):
        q = q.set_index("w_mode").loc[modes].reset_index(); c = COLORS[int(chrom)]
        axes[0].plot([0, 1], q.coverage, color=c, marker=MARKERS[int(chrom)], label=f"T{int(chrom)}")
        axes[1].plot([0, 1], q.auc, color=c, marker=MARKERS[int(chrom)])
    for ax in axes: ax.set_xticks([0, 1], ["Retain negatives", "Truncate at 0"])
    axes[0].set(ylabel="Coverage", title="Gate behavior", ylim=(0, 1.03)); axes[0].legend(frameon=False)
    axes[1].set(ylabel="Out-of-fold AUC", title="Z-layer discrimination", ylim=(0, 1.03))
    panel(axes[0], "a"); panel(axes[1], "b")
    save(fig, "fig_q4_w_trunc")


def cv_bootstrap() -> None:
    d = pd.read_csv(RES / "q4_bootstrap_metrics.csv"); d = d[d.level == "record"]
    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    vals = [d.loc[d.chrom == chrom, "auc"].to_numpy() for chrom in (13, 18, 21)]
    bp = ax.boxplot(vals, tick_labels=["T13", "T18", "T21"], patch_artist=True, showfliers=False)
    for b, chrom in zip(bp["boxes"], (13, 18, 21)): b.set_facecolor(COLORS[chrom]); b.set_alpha(.55)
    rng = np.random.default_rng(2025)
    for i, v in enumerate(vals, 1): ax.scatter(rng.normal(i, .05, len(v)), v, s=7, alpha=.14, color="0.2", linewidths=0)
    ax.axhline(.5, color="#D55E00", ls="--", lw=1.1, label="Random AUC = 0.5")
    ax.set(ylim=(0, 1), xlabel="Chromosome", ylabel="Bootstrap out-of-fold AUC",
           title="Mother-cluster bootstrap uncertainty (B=200)")
    ax.legend(frameon=False)
    save(fig, "fig_q4_cv_bootstrap")


def male_female_auc() -> None:
    d = pd.read_csv(RES / "q4_male_female_auc.csv")
    fig, ax = plt.subplots(figsize=(6.2, 3.9))
    for chrom, q in d.groupby("chrom"):
        q = q.set_index("sex").loc[["female", "male"]]
        ax.plot([0, 1], q.auc, color=COLORS[int(chrom)], marker=MARKERS[int(chrom)], lw=1.8, ms=7, label=f"T{int(chrom)}")
    ax.axhline(.5, color="0.45", ls="--", lw=1)
    ax.set_xticks([0, 1], ["Female fetus", "Male fetus"]); ax.set_ylim(0, 1)
    ax.set(ylabel="Raw Z-score AUC", title="Sex-specific discrimination requires separate calibration")
    ax.legend(frameon=False)
    save(fig, "fig_male_female_auc")


def time_sensitivity() -> None:
    d = pd.read_csv(RES / "q4_time_sens.csv")
    fig, ax = plt.subplots(figsize=(6.2, 3.9))
    for chrom, q in d.groupby("chrom"):
        q = q.sort_values("time_included")
        ax.plot([0, 1], q.auc.to_numpy(), color=COLORS[int(chrom)], marker=MARKERS[int(chrom)], lw=1.8, ms=7, label=f"T{int(chrom)}")
    ax.axhline(.5, color="0.45", ls="--", lw=1)
    ax.set_xticks([0, 1], ["No gestational age", "+ Gestational age"])
    ax.set(ylim=(0, 1), ylabel="Out-of-fold AUC", title="Gestational-age sensitivity")
    ax.legend(frameon=False)
    save(fig, "fig_q4_time_sens", bottom=.17)


def main() -> None:
    data_profile(); gate_curve(); z_roc(); cost_sensitivity(); merge_compare(); calibration()
    gc_sensitivity(); w_truncation(); cv_bootstrap(); male_female_auc(); time_sensitivity()
    print(f"Generated 11 Q4 figures in {OUT}")


if __name__ == "__main__":
    main()

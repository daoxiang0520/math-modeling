"""Publication-style figures required by revised q2_plan1_coder_task.md."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figures import (RESULTS, BLUE, ORANGE, GREEN, PURPLE, SKY, VERMILLION,
                     GRAY_DARK, GRAY_LIGHT, set_style, polish, panel_labels, save)

COLORS = {"G1": BLUE, "G2": ORANGE, "G3": GREEN, "G4": PURPLE, "G5": SKY}
MARKERS = {"G1": "o", "G2": "s", "G3": "^", "G4": "D", "G5": "P"}
LINES = {"G1": "-", "G2": "--", "G3": "-.", "G4": ":", "G5": (0, (5, 2))}


def probability_curves():
    d = pd.read_csv(RESULTS / "q2_plan1_prob_curves.csv")
    fig, ax = plt.subplots(figsize=(6.7, 3.7), constrained_layout=True)
    for group, z in d.groupby("group"):
        z = z.sort_values("ga")
        ax.fill_between(z.ga, z.ci_low, z.ci_high, color=COLORS[group], alpha=.16,
                        label=f"{group} bootstrap 95%区间")
        ax.plot(z.ga, z.p_marg, color=COLORS[group], linestyle=LINES[group],
                marker=MARKERS[group], markevery=15,
                label=f"{group}（中位BMI={z.median_bmi.iloc[0]:.1f}）")
    ax.axhline(.80, color=VERMILLION, linestyle=":", label="80%保证水平")
    ax.set(xlabel="检测孕周（周）", ylabel="Y浓度达到4%的边缘概率",
           title="单调主模型的孕周—达标概率曲线", xlim=(10, 25), ylim=(0, 1))
    ax.legend(frameon=False, ncol=2, loc="lower right")
    polish(ax, "y")
    save(fig, "fig_q2_ga_bmi_prob_curves")


def individual_t80():
    d = pd.read_csv(RESULTS / "q2_plan1_individual_t80.csv")
    main = pd.read_csv(RESULTS / "q2.csv")
    fig, ax = plt.subplots(figsize=(6.7, 3.75), constrained_layout=True)
    for group, z in d.groupby("group"):
        unc = z[~z.censored]
        cen = z[z.censored]
        ax.scatter(unc.bmi, unc["t_p0.80"], s=16, alpha=.52, color=COLORS[group],
                   marker=MARKERS[group], label=f"{group} 已穿越（n={len(unc)}）")
        if len(cen):
            ax.scatter(cen.bmi, cen["t_p0.80"], s=27, facecolors="none",
                       edgecolors=COLORS[group], marker="^", label=f"{group} ≥25周（n={len(cen)}）")
    boundaries = sorted(set(main.bmi_high.iloc[:-1].round(6)))
    for j, boundary in enumerate(boundaries):
        ax.axvline(boundary, color=VERMILLION, linestyle="--", linewidth=1.0,
                   label="联合选择BMI边界" if j == 0 else None)
    for _, row in main.iterrows():
        lo = row.bmi_low; hi = row.bmi_high
        ax.hlines(row["t_p0.80_median"], lo, hi, color=COLORS[row.group], linewidth=3.5)
    ax.set(xlabel="孕妇BMI（kg/m²）", ylabel="最早达到80%保证的孕周（周）",
           title="联合选择的BMI分组、组时点与个体右删失", ylim=(9.5, 25.7))
    ax.legend(frameon=False, ncol=2, loc="upper left")
    polish(ax, "y")
    save(fig, "fig_q2_bmi_bins_tstar")


def loss_curves():
    d = pd.read_csv(RESULTS / "q2_plan1_loss_curves.csv")
    fig, ax = plt.subplots(figsize=(6.7, 3.6), constrained_layout=True)
    for group, z in d.groupby("group"):
        ax.plot(z.ga, z.loss, color=COLORS[group], linestyle=LINES[group],
                marker=MARKERS[group], markevery=15, label=f"{group} 综合损失")
        star = z.t_star.iloc[0]
        ax.scatter([star], [z.loc[(z.ga-star).abs().idxmin(), "loss"]],
                   color=COLORS[group], marker=MARKERS[group], s=45, zorder=3)
    ax.axvline(10, color=GRAY_DARK, linestyle=":", label="优化域下界")
    ax.set(xlabel="检测孕周（周）", ylabel="综合损失（无量纲）",
           title="联合目标：过早风险权重4、延迟风险权重1", xlim=(10, 25))
    ax.legend(frameon=False, loc="upper left")
    polish(ax, "y")
    save(fig, "fig_q2_loss_curves_optimal")


def sigma_shift():
    d = pd.read_csv(RESULTS / "q2_plan1_sigma_sensitivity.csv")
    fig, ax = plt.subplots(figsize=(6.7, 3.55), constrained_layout=True)
    for group, z in d.groupby("group"):
        ax.plot(z.sigma_factor, z.delta_t, color=COLORS[group], linestyle=LINES[group],
                marker=MARKERS[group], label=group)
    ax.axhline(0, color=GRAY_DARK, linewidth=.8)
    ax.set(xlabel=r"技术误差尺度（$\sigma/\sigma_{tech}$）", ylabel="组推荐时点偏移（周）",
           title=r"技术测量误差对 $t_{0.80}$ 的影响（$\sigma_{tech}=0.133$）",
           xticks=[0, .5, 1, 2])
    ax.legend(frameon=False)
    polish(ax, "y")
    save(fig, "fig_q2_error_shift_sigma")


def boundary_heatmap():
    d = pd.read_csv(RESULTS / "q2_plan1_boundary_bootstrap.csv")
    candidates = np.sort(d.selected_boundary.unique())
    slots = np.sort(d.boundary_slot.unique())
    freq = np.vstack([
        d.loc[d.boundary_slot.eq(slot), "selected_boundary"].value_counts(normalize=True)
        .reindex(candidates, fill_value=0).to_numpy() for slot in slots
    ])
    fig, ax = plt.subplots(figsize=(6.7, 3.15), constrained_layout=True)
    im = ax.imshow(freq, cmap="viridis", vmin=0, vmax=max(.5, freq.max()), aspect="auto")
    ax.set_xticks(np.arange(len(candidates)), [f"{x:g}" for x in candidates])
    ax.set_yticks(np.arange(len(slots)), [f"切点{int(x)}" for x in slots])
    ax.set(xlabel="联合选择的BMI边界（kg/m²）", title="300次选择Bootstrap的边界稳定性")
    for i in range(len(slots)):
        for j in range(len(candidates)):
            value = freq[i, j]
            if value >= .01:
                ax.text(j, i, f"{value:.0%}", ha="center", va="center",
                        color="white" if value > .24 else "black", fontsize=7)
    cb = fig.colorbar(im, ax=ax, pad=.02)
    cb.set_label("重现频率")
    save(fig, "fig_q2_bootstrap_boundary_heatmap")


def calibration():
    d = pd.read_csv(RESULTS / "q2_plan1_calibration.csv")
    fig, axes = plt.subplots(1, 2, figsize=(6.7, 3.15), constrained_layout=True)
    for ax, kind, title in zip(axes, ["BMI组", "孕周带"], ["BMI组校准", "首次观测孕周带校准"]):
        z = d[d.calibration_type == kind]
        ax.plot([0, 1], [0, 1], color=GRAY_DARK, linestyle="--", linewidth=1, label="理想校准")
        ax.scatter(z.p_model, z.p_observed, s=18 + 1.1*z.n, color=BLUE if kind == "BMI组" else GREEN,
                   alpha=.75, edgecolor="white", linewidth=.5)
        listing = "\n".join(
            f"{r['group']}: n={int(r.n)}, 预测={r.p_model:.2f}, 观察={r.p_observed:.2f}"
            for _, r in z.iterrows()
        )
        ax.text(.03, .96, listing, transform=ax.transAxes, ha="left", va="top", fontsize=6.3,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": .82, "pad": 2})
        ax.set(xlabel="模型预测达标概率", ylabel="观察达标比例", title=title,
               xlim=(0, 1), ylim=(0, 1))
        polish(ax, "both")
    panel_labels(fig, axes)
    save(fig, "fig_q2_calibration")


def rho_sensitivity():
    d = pd.read_csv(RESULTS / "q2_plan1_rho_sensitivity.csv")
    fig, ax = plt.subplots(figsize=(6.7, 3.25), constrained_layout=True)
    for group, z in d.groupby("group"):
        ax.plot(z.rho, z.t_star, color=COLORS[group], linestyle=LINES[group],
                marker=MARKERS[group], label=group)
    ax.axhline(10, color=GRAY_DARK, linestyle=":", label="优化域下界")
    ax.set(xlabel=r"不达标损失权重 $\rho$", ylabel="二级风险最优时点（周）",
           title="过早检测风险权重敏感性",
           xticks=[.5, 1, 2], ylim=(9.5, 25.5))
    ax.legend(frameon=False, ncol=3)
    polish(ax, "y")
    save(fig, "fig_q2_rho_sensitivity")


def error_rates():
    d = pd.read_csv(RESULTS / "q2_plan1_error_classification.csv")
    fig, axes = plt.subplots(1, 2, figsize=(6.7, 3.05), sharey=True, constrained_layout=True)
    for ax, metric, title in zip(axes, ["FNR", "FPR"], ["假阴性联合概率", "假阳性联合概率"]):
        for group, z in d.groupby("group"):
            ax.plot(z.sigma_factor, z[metric], color=COLORS[group], linestyle=LINES[group],
                    marker=MARKERS[group], label=group)
        ax.set(xlabel=r"$\sigma/\sigma_{tech}$", ylabel="联合概率", title=title,
               xticks=[0, .5, 1, 2], ylim=(0, max(.001, d[["FNR", "FPR"]].to_numpy().max()*1.25)))
        polish(ax, "y")
    axes[0].legend(frameon=False)
    panel_labels(fig, axes)
    save(fig, "fig_q2_fnr_fpr_sigma")


def monotone_diagnostic():
    d = pd.read_csv(RESULTS / "q2_plan1_monotone_diagnostic.csv")
    fig, axes = plt.subplots(1, 2, figsize=(6.7, 3.15), sharey=True, constrained_layout=True)
    specs = [("linear_monotone", "主模型：纯线性，单调"),
             ("piecewise_diagnostic", "诊断模型：分段线性")]
    for ax, (form, title) in zip(axes, specs):
        z = d[d.model_form == form]
        for group, q in z.groupby("group"):
            ax.plot(q.ga, q.p_marg, color=COLORS[group], linestyle=LINES[group],
                    marker=MARKERS[group], markevery=20, label=group)
        ax.axhline(.8, color=VERMILLION, linestyle=":")
        ax.set(xlabel="检测孕周（周）", ylabel="边缘达标概率", title=title,
               xlim=(10, 25), ylim=(0, 1))
        polish(ax, "y")
    axes[0].legend(frameon=False)
    panel_labels(fig, axes)
    save(fig, "fig_q2_monotone_diagnostic")


def joint_k_selection():
    d = pd.read_csv(RESULTS / "tab_q2_joint_selection.csv")
    fig, ax = plt.subplots(figsize=(6.7, 3.45), constrained_layout=True)
    ax.errorbar(d.K, d.bootstrap_mean_loss, yerr=d.bootstrap_se, color=BLUE,
                marker="o", linestyle="-", capsize=3, label="Bootstrap均值 ± 1 SE")
    chosen = d[d.selected].iloc[0]
    ax.scatter([chosen.K], [chosen.bootstrap_mean_loss], s=85, marker="*",
               color=VERMILLION, zorder=4, label=f"选择 K={int(chosen.K)}")
    ax.axhline(chosen.one_se_threshold, color=GRAY_DARK, linestyle="--",
               label="一标准误差阈值")
    for _, row in d.iterrows():
        ax.annotate(str(row.boundaries) if pd.notna(row.boundaries) else "统一策略",
                    (row.K, row.bootstrap_mean_loss), xytext=(0, 8),
                    textcoords="offset points", ha="center", fontsize=6.5)
    ax.set(xlabel="BMI组数 K", ylabel="人均非对称决策损失（周）",
           title="组数、BMI切点与组时点的联合选择", xticks=d.K)
    ax.legend(frameon=False)
    polish(ax, "y")
    save(fig, "fig_q2_joint_k_selection")


def main():
    set_style()
    probability_curves(); individual_t80(); loss_curves(); sigma_shift()
    boundary_heatmap(); calibration(); rho_sensitivity(); error_rates(); monotone_diagnostic()
    joint_k_selection()
    print("generated 10 joint Plan 1 figures")


if __name__ == "__main__":
    main()

"""为 C题问题1生成数据、逻辑、模型原理与稳健性可视化。

输入均来自 solution.py 生成的 results/*.csv；输出 PNG/PDF/SVG，PNG 为 300 DPI。
"""

from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path.cwd() / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch
import networkx as nx
import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.stats import norm


ROOT = Path(os.environ.get("MODELING_OUTPUT_DIR", Path.cwd()))
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
GRAY = FIGURES / "grayscale"
FIGURES.mkdir(parents=True, exist_ok=True)
GRAY.mkdir(parents=True, exist_ok=True)

Y_THR = 0.04
BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
PURPLE = "#CC79A7"
SKY = "#56B4E9"
YELLOW = "#F0E442"
VERMILLION = "#D55E00"
GRAY_DARK = "#4D4D4D"
GRAY_LIGHT = "#D8D8D8"
PALETTE = [BLUE, ORANGE, GREEN, PURPLE, SKY, VERMILLION]


def set_style() -> None:
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Microsoft YaHei", "SimHei", "Arial Unicode MS", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.7,
        "lines.linewidth": 1.4,
        "savefig.dpi": 300,
        "figure.dpi": 150,
    })


def polish(ax: plt.Axes, grid_axis: str | None = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRAY_LIGHT, linewidth=0.5, alpha=0.7)
        ax.set_axisbelow(True)


def panel_labels(fig: plt.Figure, axes, labels=None) -> None:
    flat = np.ravel(np.asarray(axes, dtype=object))
    labels = labels or [chr(ord("a") + i) for i in range(len(flat))]
    for ax, label in zip(flat, labels):
        ax.text(-0.11, 1.04, label, transform=ax.transAxes, fontsize=9,
                fontweight="bold", ha="left", va="bottom")


def save(fig: plt.Figure, name: str, grayscale: bool = True) -> None:
    fig.savefig(FIGURES / f"{name}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / f"{name}.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    if grayscale:
        image = mpimg.imread(FIGURES / f"{name}.png")
        rgb = image[..., :3]
        lum = rgb @ np.array([0.2126, 0.7152, 0.0722])
        plt.imsave(GRAY / f"{name}_grayscale.png", lum, cmap="gray", vmin=0, vmax=1)


def figure_roadmap() -> None:
    fig, ax = plt.subplots(figsize=(6.7, 4.35))
    ax.axis("off")
    nodes = {
        "原始附件\n1082条记录": (0.07, 0.52),
        "解析与核验\n孕周/日期/男胎": (0.26, 0.52),
        "重复结构\n267位孕妇": (0.45, 0.76),
        "主模型\nBeta样条+混合方差层": (0.65, 0.76),
        "边缘达标概率\nP(Y>=4%)": (0.87, 0.76),
        "技术重复\n测量误差": (0.45, 0.24),
        "分位数回归\n阈值反演": (0.65, 0.24),
        "敏感性验证\n交互/GC/窗口": (0.87, 0.24),
    }
    graph = nx.DiGraph()
    graph.add_nodes_from(nodes)
    edges = [
        ("原始附件\n1082条记录", "解析与核验\n孕周/日期/男胎"),
        ("解析与核验\n孕周/日期/男胎", "重复结构\n267位孕妇"),
        ("解析与核验\n孕周/日期/男胎", "技术重复\n测量误差"),
        ("重复结构\n267位孕妇", "主模型\nBeta样条+混合方差层"),
        ("主模型\nBeta样条+混合方差层", "边缘达标概率\nP(Y>=4%)"),
        ("技术重复\n测量误差", "分位数回归\n阈值反演"),
        ("分位数回归\n阈值反演", "边缘达标概率\nP(Y>=4%)"),
        ("边缘达标概率\nP(Y>=4%)", "敏感性验证\n交互/GC/窗口"),
    ]
    graph.add_edges_from(edges)
    node_colors = [SKY, SKY, ORANGE, BLUE, GREEN, ORANGE, PURPLE, GREEN]
    nx.draw_networkx_edges(graph, pos=nodes, ax=ax, edge_color="#777777", width=1.2,
                           arrows=True, arrowsize=14, connectionstyle="arc3,rad=0.03")
    nx.draw_networkx_nodes(graph, pos=nodes, ax=ax, node_color=node_colors,
                           node_size=2350, edgecolors="white", linewidths=1.5)
    nx.draw_networkx_labels(graph, pos=nodes, ax=ax, font_family="Microsoft YaHei", font_size=8)
    ax.set_title("问题1研究逻辑：从重复测量数据到达标概率", pad=8, fontweight="bold")
    save(fig, "fig_roadmap")


def figure_scatter() -> None:
    df = pd.read_csv(RESULTS / "q1_scatter.csv")
    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.75), constrained_layout=True)
    hb1 = axes[0].hexbin(df["ga"], 100 * df["y"], C=df["bmi"], gridsize=32,
                         reduce_C_function=np.mean, cmap="viridis", mincnt=1)
    cb1 = fig.colorbar(hb1, ax=axes[0], pad=0.02)
    cb1.set_label("六边形内平均BMI (kg/m2)")
    axes[0].axhline(4, color=VERMILLION, linestyle="--", linewidth=1.2, label="4%达标阈值")
    axes[0].set(xlabel="检测孕周 (周)", ylabel="Y染色体浓度 (%)", title="孕周-Y浓度关系")
    axes[0].legend(frameon=False, loc="upper left")

    hb2 = axes[1].hexbin(df["bmi"], 100 * df["y"], C=df["ga"], gridsize=30,
                         reduce_C_function=np.mean, cmap="cividis", mincnt=1)
    cb2 = fig.colorbar(hb2, ax=axes[1], pad=0.02)
    cb2.set_label("六边形内平均孕周 (周)")
    axes[1].axhline(4, color=VERMILLION, linestyle="--", linewidth=1.2)
    axes[1].set(xlabel="孕妇BMI (kg/m2)", ylabel="Y染色体浓度 (%)", title="BMI-Y浓度关系")
    for ax in axes:
        polish(ax)
        ax.set_ylim(0, max(24, 100 * df["y"].max() * 1.03))
    panel_labels(fig, axes)
    save(fig, "fig_q1_scatter")


def figure_smooth_ga() -> None:
    df = pd.read_csv(RESULTS / "q1_smooth_ga.csv")
    levels = sorted(df["bmi"].unique())
    fig, ax = plt.subplots(figsize=(6.7, 3.3), constrained_layout=True)
    markers = ["o", "s", "^"]
    for i, level in enumerate(levels):
        sub = df[np.isclose(df["bmi"], level)]
        ax.plot(sub["ga"], 100 * sub["mean_y"], color=PALETTE[i], marker=markers[i],
                markevery=20, markersize=3.4, label=f"BMI={level:.1f}")
        ax.fill_between(sub["ga"], 100 * sub["mean_y_lo"], 100 * sub["mean_y_hi"],
                        color=PALETTE[i], alpha=0.13, linewidth=0)
    ax.axhline(4, color=VERMILLION, linestyle="--", linewidth=1.2, label="4%达标阈值")
    ax.set(xlabel="检测孕周 (周)", ylabel="模型预测Y浓度 (%)",
           title="不同BMI水平下的孕周平滑效应（阴影为95%固定效应区间）")
    ax.set_xlim(10, 25)
    polish(ax, "y")
    ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    save(fig, "fig_q1_smooth_ga")


def figure_interaction_heatmap() -> None:
    df = pd.read_csv(RESULTS / "q1_ti_heatmap.csv")
    ga = np.sort(df["ga"].unique())
    bmi = np.sort(df["bmi"].unique())
    z = df.pivot(index="bmi", columns="ga", values="pred_y").loc[bmi, ga].to_numpy() * 100
    fig, ax = plt.subplots(figsize=(6.7, 3.6), constrained_layout=True)
    mesh = ax.pcolormesh(ga, bmi, z, shading="auto", cmap="viridis")
    contour = ax.contour(ga, bmi, z, levels=[4, 6, 8, 10], colors="white", linewidths=0.8)
    ax.clabel(contour, inline=True, fmt="%g%%", fontsize=7)
    cb = fig.colorbar(mesh, ax=ax, pad=0.02)
    cb.set_label("条件均值Y浓度 (%)")
    ax.set(xlabel="检测孕周 (周)", ylabel="孕妇BMI (kg/m2)",
           title="孕周-BMI交互预测面（预设交互模型，未获显著支持，仅描述性；中心年龄、自然受孕）")
    save(fig, "fig_q1_smooth_bmi_int")


def figure_3d_relationship() -> None:
    """Supplementary 3-D view with raw observations, model surface and 4% plane."""
    heat = pd.read_csv(RESULTS / "q1_ti_heatmap.csv")
    obs = pd.read_csv(RESULTS / "q1_scatter.csv")
    obs = obs[obs["ga"].between(10, 25)].copy()
    ga = np.sort(heat["ga"].unique())
    bmi = np.sort(heat["bmi"].unique())
    ga_mesh, bmi_mesh = np.meshgrid(ga, bmi)
    z = heat.pivot(index="bmi", columns="ga", values="pred_y").loc[bmi, ga].to_numpy() * 100

    fig = plt.figure(figsize=(7.2, 5.2), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    surface = ax.plot_surface(
        ga_mesh, bmi_mesh, z, cmap="viridis", alpha=0.72,
        linewidth=0, antialiased=True, rcount=45, ccount=55,
    )
    # Preserve the original data cloud instead of showing a model-only decorative surface.
    ax.scatter(
        obs["ga"], obs["bmi"], 100 * obs["y"], s=7, c=GRAY_DARK,
        alpha=0.22, depthshade=False, edgecolors="none",
    )
    threshold = np.full_like(ga_mesh, 100 * Y_THR)
    ax.plot_surface(
        ga_mesh, bmi_mesh, threshold, color=VERMILLION, alpha=0.10,
        linewidth=0, shade=False,
    )
    ax.contour(
        ga_mesh, bmi_mesh, z, zdir="z", offset=0,
        levels=[4, 6, 8, 10], cmap="cividis", linewidths=0.9,
    )
    ax.set(
        xlabel="检测孕周 (周)", ylabel="孕妇BMI (kg/m2)",
        zlabel="Y染色体浓度 (%)", xlim=(10, 25),
        ylim=(float(obs["bmi"].min()), float(obs["bmi"].max())),
        zlim=(0, max(24.0, float(100 * obs["y"].max()) * 1.03)),
        title="Y染色体浓度与孕周、BMI的三维关系\n原始记录 + 含预设交互的Beta模型描述性曲面",
    )
    ax.view_init(elev=27, azim=-125)
    ax.set_box_aspect((1.45, 1.0, 0.85))
    ax.grid(True, linewidth=0.35, alpha=0.35)
    cb = fig.colorbar(surface, ax=ax, shrink=0.66, pad=0.08, aspect=22)
    cb.set_label("Beta模型条件均值 (%)")
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY_DARK,
                   markersize=4.5, label=f"原始记录 (n={len(obs)})"),
            Patch(facecolor=VERMILLION, alpha=0.18, label="4%临床阈值平面"),
        ],
        loc="upper left", bbox_to_anchor=(0.02, 0.96), frameon=False,
    )
    fig.text(
        0.5, 0.015,
        "注：三维透视仅作关系展示；定量比较优先参考二维交互热力图。交互项在似然比检验中未显著。",
        ha="center", va="bottom", fontsize=7, color=GRAY_DARK,
    )
    save(fig, "fig_q1_3d_relationship")


def figure_quantiles() -> None:
    df = pd.read_csv(RESULTS / "q1_quantile_curves.csv")
    median_bmi = float(np.median(np.sort(df["bmi"].unique())))
    df = df[np.isclose(df["bmi"], median_bmi)]
    chosen = [0.10, 0.25, 0.50, 0.75, 0.90]
    fig, ax = plt.subplots(figsize=(6.7, 3.35), constrained_layout=True)
    for i, tau in enumerate(chosen):
        sub = df[np.isclose(df["tau"], tau)]
        ax.plot(sub["ga"], 100 * sub["q_y"], color=PALETTE[i], marker=["v", "s", "o", "D", "^"][i],
                markevery=22, markersize=3, label=f"tau={tau:.2f}")
    ax.axhline(4, color=VERMILLION, linestyle="--", linewidth=1.2, label="4%阈值")
    ax.set(xlabel="检测孕周 (周)", ylabel="logit尺度回译的Y浓度分位数 (%)",
           title=f"中心BMI={median_bmi:.1f}下的单调条件分位数族")
    ax.set_xlim(10, 25)
    polish(ax, "y")
    ax.legend(frameon=False, ncol=6, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    save(fig, "fig_q1_quantile_curves")


def figure_prob_curves() -> None:
    df = pd.read_csv(RESULTS / "q1_prob_curves.csv")
    levels = sorted(df["bmi"].unique())
    fig, ax = plt.subplots(figsize=(6.7, 3.45), constrained_layout=True)
    markers = ["o", "s", "^"]
    for i, level in enumerate(levels):
        sub = df[np.isclose(df["bmi"], level)]
        if "p_marg_se" in sub.columns:
            ax.fill_between(sub["ga"], (sub["p_marg"] - 1.96 * sub["p_marg_se"]).clip(0, 1),
                            (sub["p_marg"] + 1.96 * sub["p_marg_se"]).clip(0, 1),
                            color=PALETTE[i], alpha=0.12, linewidth=0)
        ax.plot(sub["ga"], sub["p_marg"], color=PALETTE[i], linestyle="-",
                marker=markers[i], markevery=22, markersize=3, label=f"BMI={level:.1f}")
        ax.plot(sub["ga"], sub["p_cond"], color=PALETTE[i], linestyle="--", alpha=0.8)
    try:
        bs = pd.read_csv(RESULTS / "table_bootstrap_robust.csv")
        if len(bs) and np.isclose(bs["bmi"].iloc[0], levels[len(levels) // 2], atol=0.05):
            ax.fill_between(bs["ga"], bs["p_marg_lo_bs"], bs["p_marg_hi_bs"],
                            color=GRAY_DARK, alpha=0.14, linewidth=0, label="孕妇bootstrap 95%")
    except FileNotFoundError:
        pass
    ax.axhline(0.90, color=GRAY_DARK, linewidth=0.8, linestyle=":")
    ax.text(24.9, 0.905, "90%", ha="right", va="bottom", color=GRAY_DARK)
    ax.set(xlabel="检测孕周 (周)", ylabel="P(Y>=4%)", title="达标概率：边缘 vs 条件（浅色=MC误差带，灰色=孕妇bootstrap 95%）",
           xlim=(10, 25), ylim=(0, 1.02))
    polish(ax, "y")
    legend1 = ax.legend(frameon=False, title="BMI分位水平", loc="lower right")
    ax.add_artist(legend1)
    ax.legend(handles=[Line2D([0], [0], color=GRAY_DARK, linestyle="-", label="边缘（新孕妇）"),
                       Line2D([0], [0], color=GRAY_DARK, linestyle="--", label="条件（随机效应=0）")],
              frameon=False, loc="lower left")
    save(fig, "fig_q1_prob_curves")


def figure_diagnostics() -> None:
    df = pd.read_csv(RESULTS / "q1_resid.csv")
    fig, axes = plt.subplots(2, 2, figsize=(6.7, 5.0), constrained_layout=True)
    axes[0, 0].scatter(df["fitted"], df["resid"], s=7, alpha=0.28, color=BLUE, edgecolors="none")
    axes[0, 0].axhline(0, color=VERMILLION, linestyle="--", linewidth=1)
    axes[0, 0].set(xlabel="Beta模型拟合均值", ylabel="随机分位数残差", title="残差-拟合值")
    q = np.linspace(0.5 / len(df), 1 - 0.5 / len(df), len(df))
    theo = norm.ppf(q)
    sample = np.sort((df["resid"] - df["resid"].mean()) / df["resid"].std())
    axes[0, 1].scatter(theo, sample, s=7, alpha=0.45, color=ORANGE, edgecolors="none")
    lim = max(abs(theo).max(), abs(sample).max())
    axes[0, 1].plot([-lim, lim], [-lim, lim], color=GRAY_DARK, linestyle="--", linewidth=0.9)
    axes[0, 1].set(xlabel="理论正态分位数", ylabel="标准化残差分位数", title="Q-Q图")
    axes[1, 0].hist(df["resid"], bins=35, density=True, color=SKY, alpha=0.8, edgecolor="white")
    xx = np.linspace(df["resid"].min(), df["resid"].max(), 250)
    axes[1, 0].plot(xx, norm.pdf(xx, df["resid"].mean(), df["resid"].std()), color=VERMILLION)
    axes[1, 0].set(xlabel="随机分位数残差", ylabel="密度", title="Beta分位数残差分布")
    group_mean = df.groupby("mother_id", as_index=False)["resid"].mean()
    axes[1, 1].hist(group_mean["resid"], bins=25, color=GREEN, alpha=0.8, edgecolor="white")
    axes[1, 1].axvline(0, color=GRAY_DARK, linestyle="--", linewidth=0.9)
    axes[1, 1].set(xlabel="孕妇内平均残差", ylabel="孕妇数", title="组水平偏差")
    for ax in axes.flat:
        polish(ax, "y")
    panel_labels(fig, axes)
    save(fig, "fig_diag_resid")


def two_curve_sensitivity(filename: str, name: str, cols: tuple[str, str], labels: tuple[str, str],
                          title: str, ylabel: str = "P(Y>=4%)") -> None:
    df = pd.read_csv(RESULTS / filename)
    fig, ax = plt.subplots(figsize=(6.7, 3.15), constrained_layout=True)
    for i, (col, label) in enumerate(zip(cols, labels)):
        ax.plot(df["ga"], df[col], color=PALETTE[i], linestyle=["-", "--"][i],
                marker=["o", "s"][i], markevery=22, markersize=3, label=label)
    ax.set(xlabel="检测孕周 (周)", ylabel=ylabel, title=title, xlim=(10, 25))
    if ylabel.startswith("P"):
        ax.set_ylim(0, 1.02)
    polish(ax, "y")
    ax.legend(frameon=False)
    save(fig, name)


def figure_sensitivities() -> None:
    two_curve_sensitivity("sens_dist.csv", "fig_sens_dist", ("beta_marginal", "quantile"),
                          ("Beta边缘概率（主模型）", "分位数阈值反演"), "分布假设敏感性")
    two_curve_sensitivity("sens_interaction.csv", "fig_sens_interaction",
                          ("with_interaction", "without_interaction"), ("含孕周-BMI交互", "无交互"),
                          "交互项敏感性")
    two_curve_sensitivity("sens_gc.csv", "fig_sens_gc", ("without_gc", "with_gc"),
                          ("主模型", "加入GC连续协变量"), "GC处理敏感性（不采用硬阈值剔除）")
    two_curve_sensitivity("sens_marginal.csv", "fig_sens_marginal", ("p_marg", "p_cond"),
                          ("边缘预测", "条件预测"), "个体随机效应积分的影响")
    two_curve_sensitivity("sens_ga_window.csv", "fig_sens_ga_window", ("all_records", "only_10_25"),
                          ("全部记录拟合", "仅10-25周拟合"), "孕周建模窗口敏感性")
    two_curve_sensitivity("sens_ga_crosscheck.csv", "fig_sens_ga_crosscheck",
                          ("all_records", "date_crosscheck_pass"),
                          ("全部解析成功记录", "日期核验差异≤1周"), "孕周日期交叉核验敏感性")
    two_curve_sensitivity("sens_health.csv", "fig_sens_health", ("retain_all", "exclude_unhealthy"),
                          ("保留不健康记录", "剔除不健康记录"), "胎儿健康记录处理敏感性")


def figure_threshold_anchor() -> None:
    df = pd.read_csv(RESULTS / "q1_threshold_anchor.csv")
    bins = np.arange(np.floor(df["ga"].min()), np.ceil(df["ga"].max()) + 1)
    df["ga_bin"] = pd.cut(df["ga"], bins=bins, include_lowest=True)
    agg = df.groupby("ga_bin", observed=True).agg(ga=("ga", "mean"), rate=("y", lambda s: (s >= Y_THR).mean()), n=("y", "size"))
    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.85), constrained_layout=True)
    hb = axes[0].hexbin(df["ga"], 100 * df["y"], gridsize=34, bins="log", cmap="cividis", mincnt=1)
    cb = fig.colorbar(hb, ax=axes[0], pad=0.02)
    cb.set_label("log10(记录密度)")
    axes[0].axhline(4, color=VERMILLION, linestyle="--", linewidth=1.2)
    axes[0].set(xlabel="检测孕周 (周)", ylabel="Y染色体浓度 (%)", title="4%阈值锚定于原始分布")
    axes[1].plot(agg["ga"], agg["rate"], color=BLUE, marker="o", markersize=4)
    sizes = 18 + 70 * agg["n"] / agg["n"].max()
    axes[1].scatter(agg["ga"], agg["rate"], s=sizes, facecolor="white", edgecolor=BLUE, zorder=3)
    axes[1].set(xlabel="孕周分箱中心 (周)", ylabel="箱内达标比例", title="观察达标比例（点大小表示n）", ylim=(0, 1.02))
    for ax in axes:
        polish(ax, "y")
    panel_labels(fig, axes)
    save(fig, "fig_anchor_threshold")


def figure_data_quality() -> None:
    df = pd.read_csv(RESULTS / "q1_scatter.csv")
    counts = df.groupby("mother_id").size()
    fig, axes = plt.subplots(1, 3, figsize=(6.7, 2.65), constrained_layout=True)
    axes[0].hist(counts, bins=np.arange(0.5, counts.max() + 1.5), color=BLUE, edgecolor="white")
    axes[0].set(xlabel="每位孕妇记录数", ylabel="孕妇数", title="纵向重复结构")
    diff = df["ga_date_diff"].dropna()
    axes[1].hist(diff, bins=30, color=ORANGE, edgecolor="white")
    axes[1].axvspan(-1, 1, color=GREEN, alpha=0.14, label="容差+/-1周")
    axes[1].set(xlabel="文本孕周-日期推算孕周 (周)", ylabel="记录数", title="孕周交叉核验")
    axes[1].legend(frameon=False)
    axes[2].hist(df["bmi"], bins=30, color=GREEN, edgecolor="white")
    axes[2].axvline(df["bmi"].median(), color=GRAY_DARK, linestyle="--", label=f"中位数 {df['bmi'].median():.1f}")
    axes[2].set(xlabel="BMI (kg/m2)", ylabel="记录数", title="高BMI样本分布")
    axes[2].legend(frameon=False)
    for ax in axes:
        polish(ax, "y")
    panel_labels(fig, axes)
    save(fig, "fig_data_quality")


def figure_model_principle() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(6.7, 3.25), constrained_layout=True)
    axes[0].axis("off")
    components = [
        (0.03, 0.78, "s1(孕周)", SKY), (0.03, 0.56, "s2(BMI)", ORANGE),
        (0.03, 0.34, "s3(年龄)+IVF", GREEN), (0.37, 0.62, "随机截距+斜率", YELLOW),
        (0.73, 0.52, "logit均值 eta", BLUE),
    ]
    for x, y, text, color in components:
        patch = FancyBboxPatch((x, y), 0.24, 0.14, boxstyle="round,pad=0.02,rounding_size=0.02",
                               facecolor=color, alpha=0.78, edgecolor="white", transform=axes[0].transAxes)
        axes[0].add_patch(patch)
        axes[0].text(x + 0.12, y + 0.07, text, transform=axes[0].transAxes, ha="center", va="center")
    for start in [(0.27, 0.85), (0.27, 0.63), (0.27, 0.41), (0.61, 0.69)]:
        axes[0].annotate("", xy=(0.73, 0.59), xytext=start, xycoords=axes[0].transAxes,
                         arrowprops=dict(arrowstyle="->", color=GRAY_DARK, lw=1))
    axes[0].text(0.5, 0.04, "两阶段估计：Beta 均值层 + REML 随机效应方差层；主模型不含孕周-BMI交互",
                 transform=axes[0].transAxes, ha="center", color=GRAY_DARK)
    axes[0].set_title("关系模型的组成")

    x = np.linspace(-5.2, -1.2, 500)
    eta = logit(0.055)
    cond_sd, marg_sd = 0.42, 0.78
    axes[1].plot(x, norm.pdf(x, eta, cond_sd), color=ORANGE, label="条件分布")
    axes[1].plot(x, norm.pdf(x, eta, marg_sd), color=BLUE, linestyle="--", label="边缘分布")
    threshold = logit(Y_THR)
    axes[1].axvline(threshold, color=VERMILLION, linestyle=":", linewidth=1.4)
    axes[1].fill_between(x, 0, norm.pdf(x, eta, marg_sd), where=x >= threshold, color=BLUE, alpha=0.16)
    axes[1].text(threshold, axes[1].get_ylim()[1] * 0.92, "logit(4%)", ha="center", va="top", color=VERMILLION)
    axes[1].set(xlabel="logit(Y浓度)", ylabel="概率密度", title="边缘化：把个体差异积分进达标概率")
    polish(axes[1])
    axes[1].legend(frameon=False)
    panel_labels(fig, axes)
    save(fig, "fig_model_principle")


def figure_model_comparison() -> None:
    df = pd.read_csv(RESULTS / "table_model_comparison.csv")
    order = df[df["分布族"].eq("beta")].sort_values("delta_AIC", ascending=True)
    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.8), constrained_layout=True)
    y = np.arange(len(order))
    axes[0].barh(y, order["delta_AIC"], color=[BLUE] + [GRAY_LIGHT] * (len(order) - 1))
    axes[0].set_yticks(y, order["模型"])
    axes[0].invert_yaxis()
    axes[0].set(xlabel="delta AIC（越小越好）", title="拟合优度")
    rm = df.sort_values("RMSE")
    yy = np.arange(len(rm))
    axes[1].scatter(rm["RMSE"], yy, s=38, color=ORANGE)
    axes[1].set_yticks(yy, rm["模型"])
    axes[1].invert_yaxis()
    axes[1].set(xlabel="按孕妇分组5折CV RMSE（Y浓度比例）", title="泛化误差")
    for ax in axes:
        polish(ax, "x")
    panel_labels(fig, axes)
    save(fig, "fig_model_comparison")


def contact_sheet(names: list[str]) -> None:
    ncols = 3
    nrows = int(np.ceil(len(names) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 3.3 * nrows), constrained_layout=True)
    for ax, name in zip(axes.flat, names):
        ax.imshow(mpimg.imread(FIGURES / f"{name}.png"))
        ax.set_title(name, fontsize=8)
        ax.axis("off")
    for ax in axes.flat[len(names):]:
        ax.axis("off")
    fig.savefig(FIGURES / "all_figures_contact_sheet.png", dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    set_style()
    figure_roadmap()
    figure_scatter()
    figure_smooth_ga()
    figure_interaction_heatmap()
    figure_3d_relationship()
    figure_quantiles()
    figure_prob_curves()
    figure_diagnostics()
    figure_sensitivities()
    figure_threshold_anchor()
    figure_data_quality()
    figure_model_principle()
    figure_model_comparison()
    names = [
        "fig_roadmap", "fig_data_quality", "fig_q1_scatter", "fig_q1_smooth_ga",
        "fig_q1_smooth_bmi_int", "fig_q1_3d_relationship", "fig_q1_quantile_curves", "fig_q1_prob_curves",
        "fig_model_principle", "fig_model_comparison", "fig_diag_resid", "fig_sens_dist",
        "fig_sens_interaction", "fig_sens_gc", "fig_sens_marginal", "fig_sens_ga_window",
        "fig_sens_ga_crosscheck", "fig_sens_health", "fig_anchor_threshold",
    ]
    contact_sheet(names)
    print(f"generated {len(names)} figures in {FIGURES}")


if __name__ == "__main__":
    main()

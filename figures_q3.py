"""Generate the 13 publication-style figures required by q3_coder_task.md."""
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
RES = ROOT / "results"
OUT = ROOT / "figures"
OUT.mkdir(parents=True, exist_ok=True)

COLORS = {"G1": "#0072B2", "G2": "#D55E00", "auxiliary_beta": "#009E73", "primary_binomial": "#CC79A7"}
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
    "axes.unicode_minus": False,
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 9.5,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 130,
    "savefig.dpi": 320,
})


def save(fig: plt.Figure, name: str) -> None:
    fig.subplots_adjust(left=0.11, right=0.97, bottom=0.14, top=0.90, wspace=0.30, hspace=0.34)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{name}.{ext}", facecolor="white")
    plt.close(fig)


def panel_label(ax, label: str) -> None:
    ax.text(-0.12, 1.04, label, transform=ax.transAxes, fontweight="bold", va="bottom")


def data_profile() -> None:
    d = pd.read_csv(RES / "q3_data_profile.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35))
    ax = axes[0]
    sc = ax.scatter(d.ga, d.y * 100, c=d.bmi, s=11, alpha=.55, cmap="viridis", linewidths=0)
    cb = fig.colorbar(sc, ax=ax, pad=.02); cb.set_label("BMI (kg/m²)")
    ax.axhline(4, color="#D55E00", ls="--", lw=1, label="4% threshold")
    ax.set(xlabel="Gestational age (weeks)", ylabel="Y-chromosome concentration (%)", title="Repeated observations")
    ax.legend(frameon=False, loc="upper left"); panel_label(ax, "a")
    ax = axes[1]
    m = d.groupby("mother_id", as_index=False).agg(bmi=("bmi", "median"), repeats=("ga", "size"))
    bins = np.arange(np.floor(m.bmi.min()), np.ceil(m.bmi.max()) + 1.5, 1.5)
    ax.hist(m.bmi, bins=bins, color="#56B4E9", edgecolor="white", linewidth=.5)
    ax.axvline(36, color="#D55E00", ls="--", lw=1.2, label="Selected boundary = 36")
    ax.set(xlabel="Mother-level median BMI (kg/m²)", ylabel="Number of mothers", title=f"Mother-level BMI distribution (n={len(m)})")
    ax.legend(frameon=False); panel_label(ax, "b")
    save(fig, "fig_q3_data_profile")


def covariate_selection() -> None:
    d = pd.read_csv(RES / "q3_covariate_forms.csv")
    labels = {"bmi": "BMI", "weight": "Weight", "height_weight": "Height + weight", "all": "All variables"}
    x = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(7.2, 3.7))
    ax.errorbar(x, d.CV_RMSE, yerr=d.CV_RMSE_SE, fmt="o", ms=7, capsize=4, color="#0072B2", ecolor="#0072B2")
    threshold = float(d.one_se_threshold.iloc[0])
    ax.axhline(threshold, color="#D55E00", ls="--", lw=1.2, label=f"1-SE threshold = {threshold:.4f}")
    chosen = d.index[d.chosen_by_1se.astype(bool)][0]
    ax.scatter([chosen], [d.loc[chosen, "CV_RMSE"]], s=125, facecolors="none", edgecolors="#D55E00", lw=1.8, label="Selected")
    for i, row in d.iterrows():
        ax.text(i, row.CV_RMSE + row.CV_RMSE_SE + .00012, f"VIF={row.VIF_max:.1f}", ha="center", fontsize=7.4)
    ax.set_xticks(x, [labels[v] for v in d.form], rotation=8)
    ax.set(xlabel="Covariate form", ylabel="Grouped-CV RMSE", title="Covariate-form selection by mother-grouped cross-validation")
    ax.set_ylim(d.CV_RMSE.min() - .0013, max(threshold + .0018, (d.CV_RMSE + d.CV_RMSE_SE).max() + .0012))
    ax.legend(frameon=False, ncol=2, loc="lower right")
    save(fig, "fig_q3_covariate_selection")


def prob_curves() -> None:
    d = pd.read_csv(RES / "q3_group_prob_curves.csv")
    main = pd.read_csv(RES / "q3_main.csv")
    fig, ax = plt.subplots(figsize=(7.2, 4.1))
    for g, q in d.groupby("group"):
        c = COLORS[g]
        ax.fill_between(q.ga, q.p_lo, q.p_hi, color=c, alpha=.16, linewidth=0)
        ax.plot(q.ga, q.p_marg, color=c, lw=2, label=f"{g} marginal probability")
        t = float(main.loc[main.group == g, "t_g"].iloc[0])
        p = float(main.loc[main.group == g, "pi_g_at_tg"].iloc[0])
        ax.plot(t, p, "o", color=c, ms=6)
        ax.annotate(f"{t:.1f} weeks", (t, p), xytext=((-52, 8) if g == "G2" else (4, -15)), textcoords="offset points", color=c)
    ax.axhline(.80, color="black", ls="--", lw=1, label="Guarantee level 0.80")
    ax.set(xlim=(10, 25), ylim=(0, 1.01), xlabel="Gestational age (weeks)", ylabel="Marginal P(Y ≥ 4%)", title="Group-level threshold probabilities with 95% bootstrap bands")
    ax.legend(frameon=False, loc="lower right")
    save(fig, "fig_q3_prob_curves")


def individual_hist() -> None:
    d = pd.read_csv(RES / "q3_individual_tp.csv")
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    unc = d.loc[~d.censored.astype(bool), "t_p"]
    cen = d.censored.astype(bool).sum()
    ax.hist(unc, bins=np.arange(10, 25.51, .75), color="#0072B2", alpha=.82, edgecolor="white")
    ax.axvline(25, color="#D55E00", ls="--", lw=1.4)
    ax.text(24.7, ax.get_ylim()[1] * .82, f"Right-censored at 25 weeks: {cen}", ha="right", color="#D55E00")
    ax.set(xlabel="Individual first crossing time $t_p$ (weeks)", ylabel="Number of mothers", title="Distribution of individual 0.80 crossing times")
    save(fig, "fig_q3_individual_tp_hist")


def group_tradeoff() -> None:
    d = pd.read_csv(RES / "q3_main.csv")
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    x = np.arange(len(d))
    for i, row in d.iterrows():
        c = COLORS[row.group]
        ax.errorbar(i, row.t_g, yerr=[[row.t_g-row.ci_low], [row.ci_high-row.t_g]], fmt="o", ms=8, capsize=5, color=c)
        ax.text(i, row.t_g + (1.0 if row.t_g < 24 else -2.0), f"n={int(row.n)}\ncensored={row.r_cens:.1%}", ha="center", va="bottom" if row.t_g < 24 else "top", color=c)
    labels = [f"{r.group}\nBMI [{r.bmi_low:g}, {r.bmi_high:g}{')' if i == 0 else ']'}" for i, (_, r) in enumerate(d.iterrows())]
    ax.set_xticks(x, labels)
    ax.set(ylim=(9.5, 26), xlabel="Selected BMI group", ylabel="Recommended time (weeks)", title="Recommended timing and mother-cluster bootstrap interval")
    ax.axhline(25, color="0.5", ls=":", lw=1)
    save(fig, "fig_q3_group_t_tradeoff")


def equivalence() -> None:
    d = pd.read_csv(RES / "q3_equivalence.csv")
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    ax.scatter(d.t_first_crossing, d.t_constrained_argmin, s=18, alpha=.45, color="#0072B2", edgecolors="none")
    ax.plot([10, 25], [10, 25], color="#D55E00", ls="--", lw=1.2, label="Identity")
    maxdiff = d.difference.abs().max()
    ax.text(.03, .94, f"max |difference| = {maxdiff:.2g} weeks", transform=ax.transAxes, va="top")
    ax.set(xlim=(9.7, 25.3), ylim=(9.7, 25.3), xlabel="First-crossing solution (weeks)", ylabel="Constrained-argmin solution (weeks)", title="Equivalence of the two timing definitions")
    ax.legend(frameon=False)
    save(fig, "fig_q3_equivalence_check")


def monotone() -> None:
    d = pd.read_csv(RES / "q3_monotone_diagnostic.csv")
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.scatter(d.bmi, d.min_delta_p, s=18, alpha=.62, color="#009E73", edgecolors="none")
    ax.axhline(0, color="#D55E00", ls="--", lw=1.1)
    ax.text(.02, .95, f"Violations after fallback: {(d.min_delta_p < -1e-10).sum()} / {len(d)}", transform=ax.transAxes, va="top")
    ax.set(xlabel="BMI (kg/m²)", ylabel="Minimum adjacent probability increment", title="Monotonicity diagnostic after prespecified random-intercept fallback")
    save(fig, "fig_q3_monotone_diagnostic")


def error_shift() -> None:
    d = pd.read_csv(RES / "q3_error_sensitivity.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35))
    for g, q in d.groupby("group"):
        axes[0].plot(q.sigma_factor, q.t_g_sigma, marker="o", lw=1.8, color=COLORS[g], label=g)
        axes[1].plot(q.sigma_factor, q.kappa_sigma, marker="o", lw=1.8, color=COLORS[g], label=g)
    axes[0].set(xlabel=r"Measurement error factor $\sigma/\sigma_{tech}$", ylabel="Recommended time (weeks)", title="Timing shift")
    axes[1].set(xlabel=r"Measurement error factor $\sigma/\sigma_{tech}$", ylabel=r"Robustness penalty $\kappa_\sigma$", title="Probability degradation")
    axes[0].legend(frameon=False); panel_label(axes[0], "a"); panel_label(axes[1], "b")
    save(fig, "fig_q3_error_shift_sigma")


def boundary_bootstrap() -> None:
    d = pd.read_csv(RES / "q3_boundary_bootstrap.csv")
    rows = []
    for _, r in d.iterrows():
        for b in str(r.boundaries).split("|"):
            rows.append((int(r.K), float(b)))
    z = pd.DataFrame(rows, columns=["K", "boundary"])
    tab = pd.crosstab(z.K, z.boundary).reindex(index=[2, 3, 4], columns=[24, 26, 28, 30, 32, 34, 36], fill_value=0)
    denom = d.groupby("K").size().reindex(tab.index).to_numpy()[:, None]
    freq = tab.to_numpy() / denom
    fig, ax = plt.subplots(figsize=(7.2, 3.5))
    im = ax.imshow(freq, aspect="auto", cmap="Blues", vmin=0, vmax=max(.01, freq.max()))
    for i in range(freq.shape[0]):
        for j in range(freq.shape[1]):
            ax.text(j, i, f"{freq[i,j]:.0%}", ha="center", va="center", color="white" if freq[i,j] > freq.max()*.55 else "black", fontsize=7.5)
    ax.set_xticks(range(len(tab.columns)), [f"{v:g}" for v in tab.columns]); ax.set_yticks(range(3), ["K=2", "K=3", "K=4"])
    ax.set(xlabel="Candidate BMI boundary (kg/m²)", ylabel="Number of groups", title="Bootstrap frequency of selected BMI boundaries")
    cb = fig.colorbar(im, ax=ax, pad=.02); cb.set_label("Selection frequency")
    save(fig, "fig_q3_bmi_boundary_bootstrap")


def p_sensitivity() -> None:
    d = pd.read_csv(RES / "q3_p_sensitivity.csv")
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    for g, q in d.groupby("group"):
        ax.plot(q.p_guarantee, q.t_g, marker="o", lw=1.8, color=COLORS[g], label=g)
    ax.axhline(25, color="0.5", ls=":", lw=1)
    ax.set(xticks=sorted(d.p_guarantee.unique()), ylim=(9.5, 25.5), xlabel="Required guarantee probability", ylabel="Recommended time (weeks)", title="Sensitivity to the target guarantee probability")
    ax.legend(frameon=False)
    save(fig, "fig_q3_p_sensitivity")


def risk_sensitivity() -> None:
    d = pd.read_csv(RES / "q3_risk_sensitivity.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), sharey=True)
    markers = {"linear": "o", "quadratic": "s", "clinical_piecewise": "^"}
    form_colors = {"linear": "#0072B2", "quadratic": "#D55E00", "clinical_piecewise": "#009E73"}
    gamma_styles = {0.0: ":", 0.5: "--", 1.0: "-", 2.0: "-."}
    for ax, (g, qg) in zip(axes, d.groupby("group")):
        for form, q in qg.groupby("risk_form"):
            for gamma, qq in q.groupby("gamma"):
                ax.plot(qq.rho, qq.t_star, marker=markers[form], color=form_colors[form], ls=gamma_styles[float(gamma)], lw=1.25, alpha=.85)
        ax.set_title(f"{g}: loss-optimal time"); ax.set_xlabel(r"Late-risk weight $\rho$"); panel_label(ax, "a" if g == "G1" else "b")
    axes[0].set_ylabel("Loss-minimizing time (weeks)")
    handles = [Line2D([0], [0], color=form_colors[f], marker=markers[f], label=f.replace("_", " ")) for f in form_colors]
    handles += [Line2D([0], [0], color="0.35", ls=gamma_styles[g], label=f"γ={g:g}") for g in gamma_styles]
    fig.legend(handles=handles, frameon=False, fontsize=7, ncol=4, loc="lower center", bbox_to_anchor=(.54, -.005))
    save(fig, "fig_q3_risk_sensitivity")


def structure_sensitivity() -> None:
    d = pd.read_csv(RES / "q3_model_structure_sensitivity.csv")
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    offsets = {"auxiliary_beta": -.12, "primary_binomial": .12}
    labels = {"auxiliary_beta": "Auxiliary Beta channel", "primary_binomial": "Primary binary channel"}
    for model, q in d.groupby("model"):
        x = np.array([0 if g == "G1" else 1 for g in q.group]) + offsets[model]
        ax.scatter(x, q.t_g, s=70, color=COLORS[model], marker="o" if model == "auxiliary_beta" else "s", label=labels[model])
    ax.set_xticks([0, 1], ["G1", "G2"]); ax.set_ylim(9.5, 25.5)
    ax.set(xlabel="BMI group", ylabel="Recommended time (weeks)", title="Sensitivity to probability-model channel")
    ax.legend(frameon=False)
    save(fig, "fig_q3_model_structure_sensitivity")


def coefficient_dualchannel() -> None:
    d = pd.read_csv(RES / "q3_model_coef.csv")
    d = d[d.term != "const"].copy()
    term_map = {"ga_c": "Gestational age", "weight_z": "Weight (z)", "age_z": "Age (z)", "ivf": "IVF"}
    order = ["ivf", "age_z", "weight_z", "ga_c"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.55), sharey=True)
    models = ["auxiliary_beta", "primary_binomial"]
    titles = ["Auxiliary Beta mean model", "Primary binary GLMM"]
    for ax, model, title in zip(axes, models, titles):
        q = d[d.model == model].set_index("term").reindex(order).reset_index()
        y = np.arange(len(q))
        xerr = np.vstack([q.estimate-q.ci_low, q.ci_high-q.estimate])
        ax.errorbar(q.estimate, y, xerr=xerr, fmt="o", ms=6, capsize=3, color=COLORS[model])
        ax.axvline(0, color="0.45", ls="--", lw=1)
        ax.set_yticks(y, [term_map[t] for t in q.term]); ax.set_xlabel("Coefficient (95% CI)"); ax.set_title(title)
    panel_label(axes[0], "a"); panel_label(axes[1], "b")
    save(fig, "fig_q3_model_coef_dualchannel")


def main() -> None:
    data_profile(); covariate_selection(); prob_curves(); individual_hist(); group_tradeoff()
    equivalence(); monotone(); error_shift(); boundary_bootstrap(); p_sensitivity()
    risk_sensitivity(); structure_sensitivity(); coefficient_dualchannel()
    print(f"Generated 13 figures in {OUT}")


if __name__ == "__main__":
    main()

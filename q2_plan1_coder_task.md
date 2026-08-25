# Q2 方案一（plan_1，联合优化版）实现架构说明书——GA 单调 Beta 随机截距模型 + BMI分组/时点联合选择

> 本任务包为精编版（框架 Architect 因 token 上限输出兜底，已由人工按 R4 LTM 补全）。
> 并列方案见 `outputs/tasks/q2_plan2/coder_task.md`（贝叶斯），两案共享数据层与决策语义。

## 0. 当前小题
问题2：临床证明男胎孕妇 BMI 是影响 Y 染色体浓度最早达标时间（≥4%）的主要因素。
试对男胎孕妇 BMI 合理分组，给出每组 BMI 区间和最佳 NIPT 时点，使潜在风险最小，
并分析检测误差对结果的影响。

## 1. 建模目标（R4 LTM，不得修改）
以 GA 纯线性 Beta 随机截距模型为决策层（β_ga>0 断言 → P_marg 单调 → t_p 适定），
在 K∈{1,...,5} 中联合选择组数、0.5 kg/m²刻度的连续BMI切点与各组统一检测时点。
个体保证概率固定为0.80，组覆盖率主值也为0.80；用过早:延迟=4:1的非对称绝对损失选择组时点，
并以孕妇层bootstrap区间量化模型不确定性。不允许降级为方向性结论；双通道误差语义分离。

## 2. 关键建模设定（15 条假设摘要，完整见动态 LTM）
1. 男胎筛选：孕妇代码 + Y浓度(V列) + Y-Z值(U列) 有效值交叉核查。
2. 孕妇层纵向重复测量；主决策模型采用孕妇随机截距 b0i，不以行为独立样本；随机斜率仅作单调性诊断，不进入主决策概率。
3. 孕周 ga=w+d/7，纯周 d=0；末次月经+检测日期交叉核对（±1 周）。
4. 【关键】决策层 GA 纯线性 Beta 随机截距模型：η=β0+β_ga·ga+β_bmi·bmi+β_age·age+β_ivf·ivf
   +b0i；无 GA-BMI 交互/平滑项和随机斜率；断言 β_ga>0。
5. 【关键】y_thr=0.04；p_guarantee=0.80 为决策保证水平（与 y_thr 区分），敏感性 0.75/0.85/0.90。
6. 【关键】P_marg(t,b) 仅作 BMI 函数：age_ref=男胎个体年龄中位数、ivf_ref=众数；b_i=median_j(BMI)。
7. 【关键】数值断言 ΔP_marg≥−1e−10（0.1 周网格）。主模型唯一固定为随机截距模型，不允许在删除随机斜率、截断随机斜率等路径间事后选择；含随机斜率模型只进入诊断图。
8. 【关键】主交付始终给具体时点；t_p0.80(b)=inf{t∈[10,25]:P_marg≥0.80}；无解右删失
   t_p=25.0、c_i=1，报告 n_unsolved 与 r_cens；联合损失中将其视为晚于25周，组覆盖率不计为已达标。
9. 【关键】联合选择 K∈{1,...,5}、连续BMI切点与组时点；切点限0.5 kg/m²刻度，每组n≥30，
   相邻组推荐时点差≥0.5周。所有候选均包含K=1统一策略，不得预设K=2或BMI=30。
10. 【关键】对个体t_i=t_p0.80(b_i)，组时点T_g最小化
    Σ_i[4·max(t_i−T_g,0)+max(T_g−t_i,0)]，等价于组内t_i的80%分位数，目标组覆盖率q=0.80。
    分组用连续区间动态规划求解；K用300次条件选择bootstrap的一标准误差规则选取最简单合格模型。
11. 【关键】联合决策损失为过早检测风险4·max(t_i−T_g,0)与延迟风险max(T_g−t_i,0)之和；不再使用会把主结果压到下界的连续延迟概率损失。
12. 为兼容既有12列契约，t_star与t_p0.80_median均写入联合优化组时点T_g；rho敏感性解释为过早风险权重乘数ρ∈{0.5,1,2}。
13. 【关键】双通道分离：通道 B=孕妇层 cluster bootstrap（B=100）主区间；通道 A=σ_tech 卷积
    （σ∈{0,0.5σ,σ,2σ}）独立敏感性列（Δt 与错分率），不叠加。
14. GC/测序质量列仅作敏感性；不按 40–60% 硬剔除（实测 0.386–0.421，平台系统偏差）。
15. 推断域 t∈[10,25]（主交付）；风险视图至 27 周；parse_hints 按规范解析字符串列。

## 3. 算法与求解（伪代码步骤）
1. 加载附件.xlsx（男胎 sheet）；parse_hints：孕妇代码去 A 转 float、末次月经 to_datetime、
   非整倍体去 T 转 float、孕周 w+d/7；男胎筛选（U/V 列有效值交叉核查）；报告筛选前后样本量。
2. 聚合孕妇单元：b_i=median_j(BMI_ij)（267 人）；age_ref=中位数、ivf_ref=众数；ḡ=孕周均值。
3. 拟合决策层模型（两阶段，复用 Q1 结构）：
   a) statsmodels.othermod.betareg.BetaModel：y~Beta(μ,φ)，logit(μ)=β0+β_ga·ga+β_bmi·bmi
      +β_age·age+β_ivf·ivf（分段线性仅作 S2 敏感性）；
   b) MixedLM(REML)：同一固定效应设计矩阵上估计随机截距方差 σ_b0²；随机斜率模型仅作诊断；
   c) 断言 β_ga>0；对 0.1 周网格断言 ΔP_marg≥−1e−10；不满足即停止并报错，不事后切换模型。
4. 边缘概率矩阵：P_marg(t,b_i) = E_{b_i}[1−F_Beta(0.04; logit^{-1}(η(t,b_i,b_i)), φ)]，
   t∈{10.0,...,25.0}（151 点）× 267 孕妇，MC=1000（固定种子），向量化；clip 到 [0,1]。
5. t_p0.80 反演：首次穿越 inf{t:P_marg≥0.80}；无解→25.0、c_i=1；输出个体 t_p0.80 与
   n_unsolved、r_cens；删失>20% 组以未删失中位数为主口径并标注。
6. 联合优化：按BMI排序，为所有满足n≥30的连续区间预计算非对称损失和最优T_g；动态规划求
   K=1,...,5的最优切点/时点；用300次条件选择bootstrap和一标准误差规则选K，并报告切点稳定性。
7. 通道 B（主区间）：B=100 次孕妇层 cluster bootstrap（按孕妇有放回重抽样→重拟合模型→
   在最终分区重算组推荐），取 2.5/97.5% 分位数；90秒预案允许B=60、0.2周、500 MC。
8. 通道 A（独立敏感性）：σ∈{0,0.5σ_tech,σ_tech,2σ_tech} logit 加性卷积重算 P_marg 与
   个体t_p和联合组时点，输出组级 Δt_sigma_tech 与错分率（FNR/FPR）。
9. 权重敏感性：将过早风险权重4乘以ρ∈{0.5,1,2}，重新计算各组非对称损失最优时点。
10. 敏感性：p∈{0.75,0.85,0.90}；过早损失权重×{0.5,1,2}；K=1,...,5完整比较；
    S2分段模型仅作单调性诊断；GC/数据质量不进入主决策层。
11. 内部校准：按 GA 带（11–12/12–13/13–15/15–20）观察达标比例 vs 模型 P_marg（中位首次观测孕周）。
12. 输出：results/q2.csv（主表）+ 图表数据源 CSV + summary。

## 4. 结果契约（q2.csv 主表，1–5 行）
列：group, bmi_low, bmi_high, n, median_bmi, t_p0.80_median, ci_low, ci_high, t_star,
distinct_required, delta_t_sigma_tech, n_unsolved
- 为兼容既有12列契约，`t_p0.80_median`字段保存联合优化得到的组统一时点（即组内t_p0.80的80%分位数）；
- t_p0.80_median/ci_low/ci_high/t_star ∈ [10,25]；bmi ∈ [15,50]；n ∈ [30,267]；
- distinct_required 为 True/False；False 时附统一时点建议值；
- 允许扩展列；禁止 NaN/Inf。

## 5. 图表计划（按 plan_id 命名到 figures/）
- fig_q2_ga_bmi_prob_curves：所选各组中位 BMI 的 P_marg 曲线（含95%带）
- fig_q2_bmi_bins_tstar：个体 t_p0.80 vs BMI 散点 + 分组中位数线 + 无解标注
- fig_q2_loss_curves_optimal：各组非对称决策损失曲线与联合最优时点
- fig_q2_error_shift_sigma：通道 A σ 偏移与错分率（含"≈3 周"换算标注）
- fig_q2_bootstrap_boundary_heatmap：所选K各切点的300次选择bootstrap稳定性热图
- fig_q2_calibration：内部校准散点
- fig_q2_rho_sensitivity：ρ 敏感性折线
- fig_q2_fnr_fpr_sigma：σ 下 FNR/FPR
- fig_q2_monotone_diagnostic：单调 vs 分段 P_marg 诊断（展示 12.5–20 凹陷伪影）
- fig_q2_joint_k_selection：K=1,...,5的bootstrap损失、1 SE误差棒、选择阈值与最终K

## 6. 表格计划
- tab_q2_main_results：分组主结果表（第 4 节列）
- tab_q2_joint_selection / tab_q2_k_selection：K=1,...,5的全数据损失、bootstrap损失、一标准误差阈值与选择标记
- tab_q2_model_selection：线性 vs 分段 CV 比较（S2）
- tab_q2_sensitivity_p：p∈{0.75/0.85/0.90} 敏感性
- tab_q2_sensitivity_rho / sigma / bins：ρ、σ、边界/删失阈值敏感性
- tab_q2_calibration：内部校准表

## 7. 实现约束
- 允许库：numpy/pandas/scipy/sklearn/statsmodels/matplotlib/networkx/pulp；总时长 ≤90s；
- 结果写入 MODELING_OUTPUT_DIR/results/q2.csv（+图表数据源）；
- 数据经 MODELING_DATA_PATH / MODELING_DATA_PATHS 传入；
- 主程序固定为 `solution_q2_plan1.py`，绘图程序固定为 `figures_q2_plan1.py`；
- 图表数据源固定写入 `results/q2_plan1_*.csv`，汇总写入 `results/q2_plan1_summary.json`；
- 代码结构可复用 Q1 solution.py：数据层/边缘化/bootstrap/σ 卷积；
- 超时预案：bootstrap B 降 60、GA 网格 0.2 周、MC 降 500，并在报告中披露。

## 8. 报告结论规则（不得预写数值结果）
- 主结论必须从 `results/q2.csv` 读取，不得预先指定组数、切点或推荐周数；
- 若 `distinct_required=False`，报告统一时点；若为 True，分别报告所选各组时点；
- 必须同时披露10周下界人数、25周窗口内无解人数、p=0.75/0.80/0.85/0.90敏感性及边界诊断；
- 明确 `p_guarantee=0.80` 是决策保证水平而非题目常量，结论对其敏感时不得只报告单一阈值。

# Q2 方案一（plan_1）实现架构说明书——GA 单调 Beta 混合 + t_p0.80（频率派）

> 本任务包为精编版（框架 Architect 因 token 上限输出兜底，已由人工按 R4 LTM 补全）。
> 并列方案见 `outputs/tasks/q2_plan2/coder_task.md`（贝叶斯），两案共享数据层与决策语义。

## 0. 当前小题
问题2：临床证明男胎孕妇 BMI 是影响 Y 染色体浓度最早达标时间（≥4%）的主要因素。
试对男胎孕妇 BMI 合理分组，给出每组 BMI 区间和最佳 NIPT 时点，使潜在风险最小，
并分析检测误差对结果的影响。

## 1. 建模目标（R4 LTM，不得修改）
以 GA 纯线性 Beta 混合模型为决策层（β_ga>0 断言 → P_marg 单调 → t_p 适定），
对男胎孕妇按 BMI 固定二分组（K=2、边界 30.0），主交付 = 组内个体 t_p0.80 中位数
+ cluster bootstrap 95% 区间；不允许降级为方向性结论；双通道误差语义分离。

## 2. 关键建模设定（15 条假设摘要，完整见动态 LTM）
1. 男胎筛选：孕妇代码 + Y浓度(V列) + Y-Z值(U列) 有效值交叉核查。
2. 孕妇层纵向重复测量；随机截距 b0i + 中心化 GA 随机斜率 b1i；不以行为独立样本。
3. 孕周 ga=w+d/7，纯周 d=0；末次月经+检测日期交叉核对（±1 周）。
4. 【关键】决策层 GA 纯线性 Beta 混合：η=β0+β_ga·ga+β_bmi·bmi+β_age·age+β_ivf·ivf
   +b0i+b1i(ga−ḡ)；无 GA-BMI 交互/平滑项；断言 β_ga>0。
5. 【关键】y_thr=0.04；p_guarantee=0.80 为决策保证水平（与 y_thr 区分），敏感性 0.75/0.85/0.90。
6. 【关键】P_marg(t,b) 仅作 BMI 函数：age_ref=男胎个体年龄中位数、ivf_ref=众数；b_i=median_j(BMI)。
7. 【关键】数值断言 ΔP_marg≥0（0.1 周网格）；不满足回退固定效应 GA 线性+BMI 或截断 b1i 使 β_ga+b1i>0。
8. 【关键】主交付始终给具体时点；t_p0.80(b)=inf{t∈[10,25]:P_marg≥0.80}；无解右删失
   t_p=25.0、c_i=1，报告 n_unsolved 与 r_cens；r_cens>20% 以未删失中位数为主口径并标注。
9. 【关键】分组固定 K=2、τ_BMI=30.0（G1=[20,30)、G2=[30,∞)）；候选边界 C={24,26,28,30,32,34,36}
   的 bootstrap 重现频率作诊断；30.0 扰动 28/32 敏感性。
10. 【关键】T_g=median_{i∈g} t_p0.80(b_i)；ΔT=|T_G2−T_G1|；distinct_required=(ΔT≥0.5)；
    False 时给出统一时点=全体个体 t_p0.80 合并中位数及其 CI + 统计证据，不得只报方向。
11. 【关键】r_delay(t)=(t−10)/17，t∈[10,27]，无 12 周零平台；与原文三档风险兼容（12 内低/13–27 高/28+ 极高）。
12. t* 为二级参考列：L=ρ[1−P_marg]+r_delay，ρ=1、t0=10，0.1 周网格+三点抛物线细分；ρ∈{0.5,1,2} 敏感性。
13. 【关键】双通道分离：通道 B=孕妇层 cluster bootstrap（B=100）主区间；通道 A=σ_tech 卷积
    （σ∈{0,0.5σ,σ,2σ}）独立敏感性列（Δt 与错分率），不叠加。
14. GC/测序质量列仅作敏感性；不按 40–60% 硬剔除（实测 0.386–0.421，平台系统偏差）。
15. 推断域 t∈[10,25]（主交付）；风险视图至 27 周；parse_hints 按规范解析字符串列。

## 3. 算法与求解（伪代码步骤）
1. 加载附件.xlsx（男胎 sheet）；parse_hints：孕妇代码去 A 转 float、末次月经 to_datetime、
   非整倍体去 T 转 float、孕周 w+d/7；男胎筛选（U/V 列有效值交叉核查）；报告筛选前后样本量。
2. 聚合孕妇单元：b_i=median_j(BMI_ij)（267 人）；age_ref=中位数、ivf_ref=众数；ḡ=孕周均值。
3. 拟合决策层模型（两阶段，复用 Q1 结构）：
   a) statsmodels.othermod.betareg.BetaModel：y~Beta(μ,φ)，logit(μ)=β0+β_ga(ga−10)+β_bmi·bmi
      +β_age·age+β_ivf·ivf（分段线性仅作 S2 敏感性）；
   b) MixedLM(REML)：同一固定效应设计矩阵上估计随机截距+中心化 GA 随机斜率 Σ_b；
   c) 断言 β_ga>0；对 0.1 周网格断言 ΔP_marg≥0；不满足按假设7 回退/截断并披露。
4. 边缘概率矩阵：P_marg(t,b_i) = E_{b_i}[1−F_Beta(0.04; logit^{-1}(η(t,b_i,b_i)), φ)]，
   t∈{10.0,...,25.0}（151 点）× 267 孕妇，MC=1000（固定种子），向量化；clip 到 [0,1]。
5. t_p0.80 反演：首次穿越 inf{t:P_marg≥0.80}；无解→25.0、c_i=1；输出个体 t_p0.80 与
   n_unsolved、r_cens；删失>20% 组以未删失中位数为主口径并标注。
6. 组推荐：T_g=median_{i∈g} t_p0.80(b_i)；ΔT；distinct_required；False→统一时点+差异 CI。
7. 通道 B（主区间）：B=100 次孕妇层 cluster bootstrap（按孕妇有放回重抽样→重拟合模型→
   重算组推荐），取 2.5/97.5% 分位数；同批给出边界候选 C 的重现频率。
8. 通道 A（独立敏感性）：σ∈{0,0.5σ_tech,σ_tech,2σ_tech} logit 加性卷积重算 P_marg 与
   t_p/t*，输出组级 Δt_sigma_tech 与错分率（FNR/FPR）。
9. t* 二级参考：r_delay=(t−10)/17，L=ρ[1−P_marg]+r_delay（ρ=1），0.1 周网格+抛物线细分。
10. 敏感性：p∈{0.75,0.85,0.90}；边界 28/32；ρ∈{0.5,1,2}；S2 分段模型；s6 数据质量剔除
    （LMP 不一致/组内 Y 跳变/不健康）；s7 GC 连续协变量。
11. 内部校准：按 GA 带（11–12/12–13/13–15/15–20）观察达标比例 vs 模型 P_marg（中位首次观测孕周）。
12. 输出：results/q2.csv（主表）+ 图表数据源 CSV + summary。

## 4. 结果契约（q2.csv 主表，2–6 行）
列：group, bmi_low, bmi_high, n, median_bmi, t_p0.80_median, ci_low, ci_high, t_star,
distinct_required, delta_t_sigma_tech, n_unsolved
- t_p0.80_median/ci_low/ci_high/t_star ∈ [10,25]；bmi ∈ [15,50]；n ∈ [10,267]；
- distinct_required 为 True/False；False 时附统一时点建议值；
- 允许扩展列；禁止 NaN/Inf。

## 5. 图表计划（按 plan_id 命名到 figures/）
- fig_q2_ga_bmi_prob_curves：组中位 P_marg 曲线（G1/G2，含 95% 带）
- fig_q2_bmi_bins_tstar：个体 t_p0.80 vs BMI 散点 + 分组中位数线 + 无解标注
- fig_q2_loss_curves_optimal：新损失曲线（斜坡 t0=10）与 t*
- fig_q2_error_shift_sigma：通道 A σ 偏移与错分率（含"≈3 周"换算标注）
- fig_q2_bootstrap_boundary_heatmap：边界候选重现频率热图
- fig_q2_calibration：内部校准散点
- fig_q2_rho_sensitivity：ρ 敏感性折线
- fig_q2_fnr_fpr_sigma：σ 下 FNR/FPR
- fig_q2_monotone_diagnostic：单调 vs 分段 P_marg 诊断（展示 12.5–20 凹陷伪影）

## 6. 表格计划
- tab_q2_main_results：分组主结果表（第 4 节列）
- tab_q2_model_selection：线性 vs 分段 CV 比较（S2）
- tab_q2_sensitivity_p：p∈{0.75/0.85/0.90} 敏感性
- tab_q2_sensitivity_rho / sigma / bins：ρ、σ、边界/删失阈值敏感性
- tab_q2_calibration：内部校准表

## 7. 实现约束
- 允许库：numpy/pandas/scipy/sklearn/statsmodels/matplotlib/networkx/pulp；总时长 ≤90s；
- 结果写入 MODELING_OUTPUT_DIR/results/q2.csv（+图表数据源）；
- 数据经 MODELING_DATA_PATH / MODELING_DATA_PATHS 传入；
- 代码结构可复用 Q1 solution.py：数据层/边缘化/bootstrap/σ 卷积；
- 超时预案：bootstrap B 降 60、GA 网格 0.2 周、MC 降 500，并在报告中披露。

## 8. 报告结论（必须改写为）
- 删除"数据支持统一约 12 周检测"；
- 改为：低 BMI 组建议约 11.5–12 周、高 BMI 组约 12.3 周（80% 达标保证）；
  旧"统一 12 周"结论系损失函数零风险平台所致，非数据结论；
- 诊断依据（12 周零平台压平、四种概率源钉在 12.0–12.5、t_p80 随 BMI 分层 r≈0.82）写入报告。

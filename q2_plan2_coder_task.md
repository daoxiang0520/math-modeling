# Q2 方案二（plan_2）实现架构说明书——贝叶斯单调 Beta 混合后验决策

> 本任务包与 `outputs/tasks/q2/coder_task.md`（plan_1，频率派）并列，供编程手分别实现、
> 对比效果。数据层与决策语义完全一致（男胎 267 孕妇、b_i=BMI 中位数、t_p0.80 主交付、
> R4 不允许降级、右删失处理）；差异仅在「概率模型的推断方式」：plan_1 用频率派点估计 +
> cluster bootstrap，plan_2 用贝叶斯后验（先验约束单调）+ 后验可信区间。

## 0. 与 plan_1 共享的设定（不得修改）

- 男胎筛选：孕妇代码 + Y染色体浓度（V列）+ Y染色体Z值（U列）有效数值交叉核查。
- 孕周解析 ga = w + d/7，纯周格式 d=0；末次月经+检测日期交叉核对（容差 ±1 周）。
- 分析单元：267 位孕妇，b_i = median_j(BMI_ij)；age_ref=男胎个体年龄中位数，ivf_ref=众数。
- 阈值：y_thr = 0.04（唯一题目阈值）；p_guarantee = 0.80 主值，敏感性 0.75/0.85/0.90。
- 分组：固定 K=2、边界 30.0 kg/m²（G1=[20,30)，G2=[30,∞)∩样本范围）。
- 主交付始终给具体时点（R4 不允许降级）；个体窗口内无解按右删失：t_p=25.0、c_i=1，
  报告 n_unsolved 与删失比例；删失比例>20% 以未删失个体中位数为主口径并显式标注。
- 推断域 t∈[10.0,25.0]；风险视图 r_delay(t)=(t−10)/17 定义于 [10,27]。
- GC/测序质量列仅作敏感性；不用 40–60% 硬阈值。

## 1. 贝叶斯模型（plan_2 的核心差异）

似然：y_ij ~ Beta(μ_ij, φ)，logit(μ_ij) = β0 + β_ga(ga_ij−10) + β_bmi·b_i + β_age·age_i
+ β_ivf·ivf_i + b0i + b1i(ga_ij−10)；b_i=(b0i,b1i)' ~ N(0, Σ_b)。

先验（保证 GA 单调 + 正则化）：
- β_ga ~ 截断正态 N(μ0, 1.0) 限制在 (0, +∞)（μ0 可取 plan_1 的 β_ga 点估计，如 0.05），
  从先验结构保证 P_marg(t,b) 对 t 单调不减；
- β0, β_bmi, β_age, β_ivf ~ N(0, 2.5)；
- φ ~ LogNormal(ln(70), 0.5)；
- Σ_b：标准差 σ0, σ1 ~ HalfNormal(0.5)，相关 ρ ~ LKJ(2)（或等价先验）；
- 个体单调性：对 b1i 施加截断使 β_ga + b1i > 0（数值断言 ΔP_marg≥0，同 plan_1 假设7）。

后验推断（允许库：numpy/scipy/statsmodels，**无 PyMC**）：
- 首选：Laplace 近似——scipy.optimize 求 MAP（负对数后验），Hessian 求后验正态近似，
  参数后验协方差 = inv(Hessian)；给出参数后验区间与预测区间。
- 备选：轻量 Metropolis-Hastings（纯 numpy，链长 ≤2000、burn-in 500、thin 2，多链 2 条）；
  若 90s 预算内不收敛，降为 Laplace 并在报告中披露。
- 后验预测：从后验抽样（Laplace：多元正态抽样；MH：后验样本）→ 对随机效应数值积分
  → 得到 P_marg(t,b) 的后验分布 → t_p0.80 的后验中位数与 95% 等尾可信区间。

## 2. 主交付（与 plan_1 同表结构，CI 语义为后验可信区间）

主表列：group, bmi_low, bmi_high, n, median_bmi, t_p0.80_posterior_median, ci_low,
ci_high, t_star, distinct_required, delta_t_sigma_tech, n_unsolved。

- 组推荐 t_g(0.80) = median_{i∈g} t_p0.80^{post}(b_i)（个体 t_p 取后验中位数）；
- ci_low/ci_high：组推荐的后验可信区间（参数不确定性来自后验，抽样层面不另做 bootstrap
  ——通道 B 由后验区间取代）；
- distinct_required = |组间后验中位差| ≥ 0.5 周；同时报告 P(|Δ_t_p| ≥ 0.5) 后验概率；
  False 时给统一时点 = 全体个体 t_p0.80 后验中位数 + 可信区间（不得只报方向）；
- t_star：二级风险视图参考列（r_delay + ρ[1−P_marg] 最小化，P_marg 取后验中位，ρ=1）；
- delta_t_sigma_tech：通道 A σ_tech 卷积（σ∈{0,0.5σ,σ,2σ}）对组推荐的偏移，独立敏感性列。

## 3. 对比交付（plan_1 vs plan_2）

输出对比表：方案、组推荐（G1/G2）、CI 宽度、distinct_required、删失比例、计算时长，
并给出结论句：两方案方向是否一致、差异来源（先验/后验 vs bootstrap）。

## 4. 图（按 plan_1 的 plan_id 命名，另加后验诊断图）

- fig_q2_ga_bmi_prob_curves：后验中位 P_marg 曲线（叠加 95% 后验带）
- fig_q2_bmi_bins_tstar：个体 t_p0.80 后验中位数 vs BMI 散点 + 组中位数 + 无解标注
- fig_q2_loss_curves_optimal：新损失曲线（斜坡 t0=10）与 t*（后验中位）
- fig_q2_error_shift_sigma：通道 A σ 偏移与错分率
- fig_q2_calibration：内部校准（沿用）
- fig_q2_rho_sensitivity：ρ 敏感性（可选）
- fig_q2_posterior_diag（新增）：β_ga 后验分布；G1/G2 组推荐后验分布（或 MH 轨迹摘要）

## 5. 约束与输出

- 7 库限定；总时长 ≤90s（超时预案：Laplace 优先、MH 链长减半，报告中披露）；
- 结果写入 results/q2_plan2.csv（主表）+ 图表数据源 CSV；
- 环境变量 MODELING_DATA_PATH / MODELING_DATA_PATHS / MODELING_OUTPUT_DIR 与 plan_1 一致；
- 代码结构可复用 plan_1 的 solution_q2.py：数据加载/男胎筛选/孕周解析/边缘化/σ 卷积，
  仅替换"点估计+重拟合 bootstrap"为"后验推断"。

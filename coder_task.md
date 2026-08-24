# 建模方案与实现架构说明书

## 0. 当前小题
**当前小题（2/4）**：问题2  临床证明，男胎孕妇的BMI 是影响胎儿Y 染色体浓度的最早达标时间（即浓度达到或超
过4%的最早时间）的主要因素。试对男胎孕妇的BMI 进行合理分组，给出每组的BMI 区间和最佳NIPT
时点，使得孕妇可能的潜在风险最小，并分析检测误差对结果的影响。

## 1. 题目与目标
**题目**：2025 年高教社杯全国大学生数学建模竞赛题目 
（请先阅读“全国大学生数学建模竞赛论文格式规范”） 
 
C 题  NIPT 的时点选择与胎儿的异常判定 
 
NIPT（Non-invasive Prenatal Test，即无创产前检测）是一种通过采集母体血液、检测胎儿的游离
DNA 片段、分析胎儿染色体是否存在异常的产前检测技术，目的是通过早期检测确定胎儿的健康状况。
根据临床经验，畸型胎儿主要有唐氏综合征、爱德华氏综合征和帕陶氏综合征，这三种体征分别由胎儿
21 号、18 号和13 号“染色体游离DNA 片段的比例”（简称“染色体浓度”）是否异常决定。NIPT 的
准确性主要由胎儿性染色体（男胎XY，女胎XX）浓度判断。通常孕妇的孕期在10 周~25 周之间可以
检测胎儿性染色体浓度，且如果男胎的Y 染色体浓度达到或高于4%、女胎的X 染色体浓度没有异常，
则可认为NIPT 的结果是基本准确的，否则难以保证结果准确性要求。同时，实际中应尽早发现不健康
的胎儿，否则会带来治疗窗口期缩短的风险，早期发现（12 周以内）风险较低；中期发现（13－27 周）
风险高；晚期发现（28 周
**建模目标**：对男胎孕妇按 BMI 进行监督单调分箱，给出每组 BMI 区间与风险损失最小的最佳 NIPT 时点，并用双通道误差传播和内部校准验证量化检测误差与不确定性对时点的影响。

## 2. 建模设定（动态 LTM，编程手不得修改）
**假设**：
- 假设1：分析对象沿用问题1筛选出的男胎样本，以孕妇为单位进行分组与决策；默认男胎孕妇数 n_w=267（依据：问题1 LTM 与数据智能摘要提示 267 位孕妇、1082 条记录）；每位孕妇用其多次记录中孕妇BMI的中位数作为问题2分组变量 b_i（该取法未在原文中明确给出，为推断值；可验证：比较BMI取均值/首次/末次时的分箱边界变化）。风险：若误纳女胎或对重复测量处理错误，会污染达标概率。
- 假设2：数值孕周必须由文本“周数+天数”解析为 ga_ij = w_ij + d_ij/7，纯周格式 d_ij=0；Coder 必须按照 parse_hints 解析字符串列：孕妇代码用 df['孕妇代码'].str.replace('A','',regex=False).astype(float)，末次月经用 pd.to_datetime(df['末次月经'])，染色体的非整倍体用 df['染色体的非整倍体'].str.replace('T','',regex=False).astype(float)，检测孕周按 w+d/7 自行解析。依据：数据列解析建议。
- 假设3【关键假设】：问题2决策概率层不直接沿用问题1未惩罚样条主模型，而是在同一男胎数据上重新拟合简约 Beta 混合模型；固定效应采用 GA 分段线性（节点 12.5 周、20 周）+ BMI + 年龄 + IVF，随机效应按研究设计保留随机截距和中心化 GA 随机斜率；分段线性与纯线性两种形式用孕妇分组五折交叉验证择优，默认取分段线性，纯线性作为敏感性。依据：q1 分组 CV 实测分段线性/线性 RMSE 0.0328/0.0331，而样条模型 RMSE 0.078 且在 19–20 周出现疑似过拟合凹陷；风险：若线性形式偏差，p_marg 与最佳时点会有偏；可验证：孕妇分组五折 CV、残差诊断，并比较线性/分段线性下 t_k* 与边界变化。
- 假设4【关键假设】：对新孕妇的达标概率决策必须采用对随机效应积分的边缘概率 p_marg(t,b)，不使用随机效应=0 的条件概率 p_cond；决策输入按问题2含义固定 age=样本年龄均值 age_bar、ivf=0（BMI 是问题2唯一分组变量，其余因素留给问题3），且决策公式不含 ti(t,b) 交互项。依据：q1 实测 p_cond 比 p_marg 高 0.066–0.096；q1 主模型 LR 检验 p=0.766 剔除 ti(t,b)；风险：用条件概率会把时点推早，含交互会与衔接点矛盾；可验证：比较 p_cond 与 p_marg 的达标概率及 t_k*，交互项仅作敏感性 s8。
- 假设5【关键假设】：风险分级按题目三档连续化：早期发现（12.0 周以内）风险较低、中期发现（13–27 周）风险高、晚期发现（28.0 周以后）风险极高（原文：早期发现 12 周以内 风险较低；中期发现 13－27 周 风险高；晚期发现 28 周以后 风险极高）。为消除阶梯函数在带内损失平坦导致最优解全部钉在 12/13 周边界的问题，采用带内线性递增风险 r_delay(t)；r_ext=3.0 为无单位极高风险档尺度（该值未在原文中明确给出，为推断值，仅用于域外定义，不参与优化）。风险：若临床真实风险非连续，最优时点会偏差；可验证：阶梯/斜坡/更陡斜坡三形式敏感性 s3。
- 假设6【关键假设】：综合损失定义为 l(t;b)=ρ[1−p_marg(t,b)]+r_delay(t)，其中 ρ=c_under/c_delay 为检测过早不达标损失与延迟发现损失之比，默认 ρ=1；绝对尺度不影响 argmin，只保留比值，故不做量纲缩放。优化域取 t∈[10.0,25.0] 周（原文：通常孕妇的孕期在10 周~25 周之间可以检测胎儿性染色体浓度），28.0 周档只在域外定义、不参与优化；若最优解落在 25.0 周边界，结论写“不晚于 25 周”且不外推。依据：原文检测窗口与风险语义；风险：ρ 主观设定是结果主要不确定性来源；可验证：s1 中 ρ∈{0.25,0.5,1,2,4} 扫描。
- 假设7【关键假设】：BMI 分组以 267 位实际孕妇为单元、按 b_i 排序后进行一维动态规划分箱，不使用 bmi 均匀网格；分箱目标为直接最小化组内期望损失 J=Σ_k n_k L_k(t_k*)，其中组时点 t_k* 是组平均达标概率曲线的损失最小点，不是组内个体时点的简单平均。约束为 n_k≥n_min=20（敏感性 10/15/30）、K∈{2,...,6}、BMI 区间连续且边界圆整到 0.5 kg/m²；分箱后要求 t_1*≤...≤t_K*，若个体 t*(b_i) 出现非单调，则先做保序回归（PAVA）调整个体时点再分箱。依据：BMI 数据分布集中，5%–95% 分位为 28.27–37.18，<28 仅 5 人、≥36 仅 20 人，低尾/高尾组推荐属外推性结论；风险：均匀网格会受尾部稀疏影响，外推组时点不稳定；可验证：s5 扫描 n_min 与 K，bootstrap 边界重现频率检验边界稳定性。
- 假设8【关键假设】：检测误差传播采用双通道。通道 A：测量误差在 logit 浓度尺度作加性扰动 u~N(0,σ²)，σ∈{0,0.5σ_tech,σ_tech,2σ_tech}，σ_tech=0.133 为问题1技术重复估计值（q1 实测）；输出每组时点偏移 δt_k(σ)=t_k*(σ)−t_k*(0) 与错分率，并报告“σ_tech≈0.133 logit 约等于孕周均值效应 3 周的幅度，单点时点无临床意义，必须给时点窗口”的结论。通道 B：参数/抽样不确定性复用 q1 孕妇层 cluster bootstrap（默认 B=100），每次重抽样重拟合简约模型、重算分箱边界与 t_k*，用 0.025 与 0.975 经验分位数给出时点区间；每个 bootstrap 副本内再叠加 σ=σ_tech 或保守 σ=2σ_tech 的 MC，得到总区间作为主交付。个体异质性 ICC=0.809 已通过随机效应积分进 p_marg，σ_tech 只承担检测可重复性误差，二者不重复计入。依据：q1 技术重复与 bootstrap 结果；风险：σ 或 bootstrap 规模不足会低估总不确定性；可验证：比较 σ 四档与 B 变化后的时点区间宽度。
- 假设9【关键假设】：观察到首次达标时间左截尾严重（81.3% 首次观测已 ≥4%，83.8% 的首次达标等于首次观测），因此不以观察事件时间作生存分析主模型；将其改为概率层内部校准验证：每组用“首次观测达标比例” vs 模型在该组中位首次观测孕周处的 p̄_k，并按 GA 带（11–12、12–13、13–15、15–20 周）做观察比例与模型 p_marg 的校准表/图。依据：q1/数据观察；风险：同数据内部校准只能说明模型可信度，不是外部验证；可验证：校准表偏差大小与分箱结果对照。
- 假设10：敏感性分析清单固定为 s1 ρ∈{0.5,1,2}∪{0.25,4}；s2 p_marg 来源=简约模型 vs q1 样条 vs 分位数反演（q1 最大差 0.151，报告跨度或取保守更晚时点）；s3 r_delay 阶梯/斜坡/更陡斜坡及 12/13 周阈值；s4 σ 四档；s5 n_min∈{10,15,20,30}、K∈{2..6}；s6 数据质量剔除 LMP 不一致 20 行、组内 y 跳变>0.05 的 31 人、不健康 38 条，每类报最大 δt_k；s7 GC 作连续协变量（q1 δp≈0.0001，走流程）；s8 含 ti(t,b) 交互（q1 δp=0.028，文档化）。依据：人类架构审核反馈。
- 假设11：GC 含量不采用 40.0%–60.0% 硬阈值剔除（原文：正常 GC 含量范围为40% ~ 60%）；因实测总 GC 约在 0.386–0.421，低于 0.4 常见，判断为平台系统偏差，仅作连续敏感性 s7。依据：问题1 LTM 与题目常量；风险：硬剔除会损失大量记录并引入选择偏差；可验证：s7 纳入 GC 前后 t_k* 与边界变化。
- 假设12：主结果表必须包含 distinct_required 标识；若组间最佳时点差异 <0.5 周，必须给出统计证据说明“数据支持统一时点”，不得硬造差异。依据：人类架构审核交付约束。
**符号表**：
- n_w: 男胎孕妇数，默认 n_w=267
- n_rec: 男胎检测记录数
- i: 孕妇索引，i=1,...,n_w
- j: 孕妇内检测记录索引，j=1,...,n_i
- b_i: 第 i 位孕妇的分组 BMI，kg/m²，取该孕妇多次记录孕妇BMI的中位数
- ga_ij: 第 i 位孕妇第 j 次检测的数值孕周，单位：周
- w_ij: 检测孕周文本中的整数周部分
- d_ij: 检测孕周文本中的天数部分，0≤d_ij<7
- ḡ: 男胎样本孕周均值，用于随机斜率中心化，单位：周
- y_ij: 第 i 位孕妇第 j 次检测的 Y 染色体浓度，无量纲比例，0<y<1
- y_thr: 临床达标阈值，y_thr=0.04，即 4.0%
- µ_ij: Beta 分布条件均值，0<µ_ij<1
- φ: Beta 分布精度参数，φ>0
- β0: 固定截距
- β1: GA 线性固定效应
- β2: GA 分段线性效应，对应 (ga−12.5)_+
- β3: GA 分段线性效应，对应 (ga−20)_+
- β4: BMI 固定效应
- β5: 年龄固定效应
- β6: IVF 固定效应
- (x)_+: 正部函数，max(x,0)
- b0i: 第 i 位孕妇随机截距
- b1i: 第 i 位孕妇随机孕周斜率
- b_i_vec: 随机效应向量 (b0i,b1i)'
- Σ_b: 随机效应协方差矩阵
- σ_b0²: 随机截距方差
- σ_b1²: 随机斜率方差
- σ_b01: 随机截距与斜率协方差
- age_i: 第 i 位孕妇年龄，单位：岁
- age_bar: 男胎样本年龄均值，单位：岁
- ivf_i: 第 i 位孕妇 IVF 指示，0/1
- η_ij: logit 尺度线性预测器
- η(t,b,b_i): 固定 age=age_bar、ivf=0 时的决策线性预测器，含随机效应 b0i,b1i
- logit(u): logit 变换，logit(u)=ln(u/(1-u))
- F_Beta(y;µ,φ): Beta(µ,φ) 分布的累积分布函数在 y 处的值
- p_marg(t,b): 对随机效应积分的边缘达标概率 P(y≥y_thr | t,b)
- σ: logit 尺度加性测量误差标准差
- σ_tech: 问题1技术重复估计的测量误差标准差，σ_tech=0.133
- u: logit 尺度加性测量误差，u~N(0,σ²)
- p_obs(t,b;σ): 叠加测量误差后的观察到达标概率
- ρ: 过早不达标损失与延迟发现损失之比，ρ=c_under/c_delay，默认 ρ=1
- r_delay(t): 延迟发现风险函数，无单位
- r_ext: 28 周以后极高风险档尺度，r_ext=3.0，推断值
- T_grid: 优化时间网格，T_grid={10.0,10.1,...,25.0}
- t*(b): 个体最佳 NIPT 时点，单位：周
- t_i*: 第 i 位孕妇个体最佳 NIPT 时点，单位：周
- t_mono_i: 对个体时点序列保序回归后的单调时点，单位：周
- K: BMI 分组数，K∈{2,...,6}
- g_k: 第 k 个连续性 BMI 区间
- n_k: 第 k 组孕妇数
- pbar_k(t): 第 k 组平均达标概率曲线
- L_k(t): 第 k 组在时点 t 的组期望损失
- t_k*: 第 k 组最佳 NIPT 时点，单位：周
- J: 分箱总损失目标
- n_min: 每组最小孕妇数，默认 n_min=20
- B: 孕妇层 cluster bootstrap 重复次数，默认 B=100
- t_k*^{(b)}: 第 b 次 bootstrap 中的第 k 组最佳时点
- Q_α: 经验分位数函数
- δt_k(σ): 测量误差 σ 引起的第 k 组时点偏移，δt_k(σ)=t_k*(σ)−t_k*(0)
- FNR_k(σ): 第 k 组真实达标被误判未达标的模型错分率
- FPR_k(σ): 第 k 组真实未达标被误判达标的模型错分率
- t_p(b): 达标概率反演时点，t_p(b)=inf{t:p_marg(t,b)≥p}
- p: 达标概率反演目标值，p∈{0.80,0.85}
- p_obs_band: 校准带内的观察达标比例
- p_model_band: 校准带内模型达标概率均值
- distinct_required: 主结果表标识，指示组间时点差异是否达到可报告差异
**公式/方程**：
- 孕周数值解析：ga_ij = w_ij + d_ij/7，其中 0≤d_ij<7；纯周格式 d_ij=0。
- logit 变换：logit(u)=ln(u/(1-u))，定义域 u∈(0,1)。
- 简约 Beta 混合模型：y_ij ~ Beta(µ_ij, φ)，η_ij = logit(µ_ij) = β0 + β1·ga_ij + β2·(ga_ij−12.5)_+ + β3·(ga_ij−20)_+ + β4·b_i + β5·age_i + β6·ivf_i + b0i + b1i·(ga_ij−ḡ)。
- 随机效应分布：b_i_vec=(b0i,b1i)' ~ N(0,Σ_b)，Σ_b=[[σ_b0²,σ_b01],[σ_b01,σ_b1²]]。
- 决策线性预测器：η(t,b,b_i) = β0 + β1·t + β2·(t−12.5)_+ + β3·(t−20)_+ + β4·b + β5·age_bar + β6·0 + b0i + b1i·(t−ḡ)。
- 边缘达标概率：p_marg(t,b)=E_{b_i~N(0,Σ_b)}[1−F_Beta(y_thr; logit^{-1}(η(t,b,b_i)), φ)]。
- 延迟风险函数：r_delay(t)=0 当 t≤12.0；r_delay(t)=(t−12.0)/(27.0−12.0) 当 12.0<t≤27.0；r_delay(t)=r_ext=3.0 当 t≥28.0。
- 综合损失：l(t;b)=ρ[1−p_marg(t,b)]+r_delay(t)。
- 个体最优时点：t*(b)=argmin_{t∈T_grid} l(t;b)，T_grid={10.0,10.1,...,25.0}；在网格最小值附近用相邻三点的二次抛物线细分。
- 组平均达标概率：pbar_k(t)=1/n_k · Σ_{i∈g_k} p_marg(t,b_i)。
- 组最佳时点：t_k*=argmin_{t∈T_grid} L_k(t)，其中 L_k(t)=ρ[1−pbar_k(t)]+r_delay(t)。
- 分箱总目标：min_{K,{g_k}} J(g_1,...,g_K)=Σ_{k=1..K} n_k·L_k(t_k*)，约束 n_k≥n_min=20、K∈{2,...,6}、g_k 为连续 BMI 区间。
- 单调性约束：若个体最优时点序列 t_i*=t*(b_i) 在 BMI 排序下非单调，先以 PAVA 保序回归得到 t_mono_i；分箱后仍要求 t_1*≤...≤t_K*。
- 测量误差达标概率：p_obs(t,b;σ)=E_{b_i~N(0,Σ_b), u~N(0,σ²)}[1−F_Beta(y_thr; logit^{-1}(η(t,b,b_i)+u), φ)]。
- 测量误差时点偏移：δt_k(σ)=t_k*(σ)−t_k*(0)，其中 t_k*(σ) 用 p_obs(t,b;σ) 替代 p_marg 计算得到。
- 测量误差错分率：FNR_k(σ)=Pr(y_true≥y_thr 且 y_obs<y_thr | g_k,t_k*(0))；FPR_k(σ)=Pr(y_true<y_thr 且 y_obs≥y_thr | g_k,t_k*(0))。
- Bootstrap 时点区间：对 B=100 次孕妇层 cluster bootstrap 得到的 t_k*^{(b)} 取 [Q_0.025(t_k*^{(b)}), Q_0.975(t_k*^{(b)})]。
- 达标时间反演锚：t_p(b)=inf{t∈T_grid : p_marg(t,b)≥p}，p∈{0.80,0.85}。
- 内部校准：p_obs_band = n_ok,band / n_band；p_model_band = mean(p_marg(ga_ij,b_i))，其中 ga_ij 取该 GA 带中位首次观测孕周。
**解题思路**：步骤1：读取附件.xlsx；Coder 必须按照 parse_hints 解析字符串列：孕妇代码用 df['孕妇代码'].str.replace('A','',regex=False).astype(float)，末次月经用 pd.to_datetime，染色体的非整倍体用 df['染色体的非整倍体'].str.replace('T','',regex=False).astype(float)，检测孕周解析为 w+d/7。步骤2：复用问题1男胎筛选、孕周解析与重复测量结构；对每位男胎孕妇计算孕妇BMI中位数作为 b_i，形成 267 位孕妇分析单元。步骤3：在同一男胎数据上重拟合简约 Beta 混合模型，男生固定效应为 GA 分段线性（12.5、20 节点）+ BMI + 年龄 + IVF，随机效应为随机截距加中心化 GA 随机斜率；分段线性与纯线性用孕妇分组五折 CV 择优，默认分段线性，纯线性作敏感性。步骤4：向量化计算 p_marg(t,b_i) 查询矩阵，维度约为 151 个 GA 网格 × 267 位孕妇 × MC 1000，供后续损失、分箱与误差分析查表；MC 对随机效应 b_i 抽样并按 Beta(µ,φ) 计算 P(y≥y_thr)。步骤5：对每个 t∈{10.0,...,25.0} 计算损失 l(t;b_i)，取网格最小后用三点抛物线细分得到 t_i*；检查 t_i* 随 BMI 的单调性，若非单调先做 PAVA 保序。步骤6：按 b_i 排序做一维动态规划分箱，目标最小化 total expected loss J，约束 n_k≥n_min=20、K∈{2..6}、边界连续；K 用损失肘部 + bootstrap 边界重现频率选择，边界圆整到 0.5 kg/m²，并标注低尾/高尾外推风险。步骤7：通道 A 测量误差：对 σ∈{0,0.5σ_tech,σ_tech,2σ_tech} 计算 p_obs(t,b;σ)，重算 t_k*(σ)，输出 δt_k(σ)、FNR_k(σ)、FPR_k(σ) 和“σ_tech≈0.133 logit 相当于孕周均值效应约 3 周”的结论。步骤8：通道 B 参数不确定性：孕妇层 cluster bootstrap B=100，每次重抽样重拟合简约模型、重算 p_marg、分箱边界与 t_k*；若总时长超 90 秒，降为 B=60、GA 网格改 0.2 周并在报告中披露；组合误差区间作为主交付。步骤9：内部校准验证：按分组输出第一观测达标比例 vs 组模型达标概率，按 GA 带（11–12、12–13、13–15、15–20 周）做校准表/图。步骤10：执行敏感性 s1–s8：ρ、p_marg 来源、r_delay 形式、σ 四档、n_min 与 K、数据质量剔除、GC 连续、ti(t,b) 交互；决策公式默认不含交互，交互仅在 s8 出现。步骤11：写 results/output.csv，主结果表含组号、BMI 区间、n、组中位 BMI、t_k*、时点区间、p̄_k(t_k*)、期望损失、δt_k(σ_tech)，并设置 distinct_required 列；同时输出图数据源：损失曲线与最优时点、t*(b) vs t_p(b) 反演曲线、σ 时点偏移、bootstrap 边界重现热图、校准散点、ρ 敏感性折线。

## 3. 算法与求解
**算法摘要**：
**伪代码/实现步骤**：
1. 步骤1: 读取附件.xlsx，按parse_hints解析字符串列：孕妇代码=df['孕妇代码'].str.replace('A','',regex=False).astype(float)；末次月经=pd.to_datetime(df['末次月经'])；染色体的非整倍体=df['染色体的非整倍体'].str.replace('T','',regex=False).astype(float)；检测孕周解析为ga_ij=w_ij+d_ij/7，纯周格式d_ij=0。明确题目唯一达标阈值为y_thr=0.04；其他统计比例如81.3%、83.8%等仅作描述性统计，不得硬编码为题目阈值。
2. 步骤2: 复用问题1男胎筛选逻辑，得到男胎样本；以孕妇代码分组，对每位孕妇计算孕妇BMI中位数为b_i，形成267位孕妇分析单元；同时计算年龄均值age_bar、孕周均值ḡ。
3. 步骤3: 在同一男胎数据上重拟合简约Beta混合模型：y_ij~Beta(µ_ij,φ)，线性预测器η_ij=β0+β1·ga_ij+β2·(ga_ij-12.5)_+ +β3·(ga_ij-20)_+ +β4·b_i+β5·age_i+β6·ivf_i+b0i+b1i·(ga_ij-ḡ)。固定效应先拟合分段线性，备选纯线性；随机效应b_i~N(0,Σ_b)。用孕妇分组五折交叉验证的RMSE比较分段线性与纯线性，默认取分段线性，纯线性作敏感性记录。拟合使用scipy.optimize + 数值积分或近似降级方案；若总时长超90秒，降级为固定效应Beta回归+按孕妇经验贝叶斯随机效应近似。
4. 步骤4: 构建T_grid={10.0,10.1,...,25.0}，对每位孕妇b_i和每个t∈T_grid向量化计算η(t,b_i,b_i_vec)，用MC抽样随机效应b_i~N(0,Σ_b)计算边缘达标概率p_marg(t,b_i)=E_{b_i}[1-F_Beta(y_thr;logit^{-1}(η),φ)]。MC默认1000次；输出p_marg矩阵，维度约为151×267，并进行nan/inf检查和np.clip到[0,1]。
5. 步骤5: 对每个b_i计算综合损失l(t;b_i)=ρ[1-p_marg(t,b_i)]+r_delay(t)，ρ默认1；r_delay(t)按t≤12为0、12<t≤27为(t-12)/15、t≥28为r_ext=3.0，优化域仅到25周。取损失最小的网格点，在网格最小值附近用相邻三点二次抛物线细分得到个体最优时点t_i*；检查t_i*随b_i排序的单调性，若非单调先做PAVA保序得到t_mono_i。
6. 步骤6: 按b_i排序，用一维动态规划分箱：目标min_{K,{g_k}} J=Σ_k n_k L_k(t_k*)，约束n_k≥n_min=20、K∈{2,...,6}、区间连续；组损失L_k(t)=ρ[1-pbar_k(t)]+r_delay(t)，pbar_k(t)=mean_{i∈g_k} p_marg(t,b_i)。K用损失肘部+bootstrap边界重现频率选择，边界圆整到0.5 kg/m²。若分箱后组间t_k*差异<0.5周，标记distinct_required为False并输出统一时点统计说明；否则为True。
7. 步骤7: 通道A测量误差：对σ∈{0,0.5σ_tech,σ_tech,2σ_tech}（σ_tech=0.133）计算p_obs(t,b;σ)=E_{b_i,u}[1-F_Beta(y_thr;logit^{-1}(η+u),φ)]，u~N(0,σ²)。重算t_k*(σ)，输出δt_k(σ)=t_k*(σ)-t_k*(0)，并按t_k*(0)计算FNR_k(σ)与FPR_k(σ)；记录“σ_tech≈0.133 logit约等于孕周均值效应约3周”的结论。
8. 步骤8: 通道B参数不确定性：孕妇层cluster bootstrap B=100，每次重抽样孕妇ID并重拟合简约模型、重算p_marg、重算分箱边界与t_k*；每个bootstrap副本内叠加σ=σ_tech或保守σ=2σ_tech的MC误差，得到组合时点区间[Q_0.025(t_k*),Q_0.975(t_k*)]。若总时长超90秒，降为B=60、GA网格改0.2周并在报告披露。
9. 步骤9: 内部校准验证：对每组计算首次观测达标比例，与模型在组中位首次观测孕周处的p̄_k比较；按GA带11-12、12-13、13-15、15-20周输出观察比例与模型p_marg校准表/图。不得将内部校准解释为外部验证。
10. 步骤10: 执行敏感性s1-s8：s1 ρ∈{0.5,1,2}∪{0.25,4}；s2 p_marg来源=简约模型 vs q1样条 vs 分位数反演；s3 r_delay阶梯/斜坡/更陡斜坡；s4 σ四档；s5 n_min∈{10,15,20,30}、K∈{2..6}；s6 剔除LMP不一致20行、组内y跳变>0.05的31人、不健康38条；s7 GC作连续协变量；s8 含ti(t,b)交互。记录每项最大δt_k与边界变化。
11. 步骤11: 保存主结果到MODELING_OUTPUT_DIR/results/output.csv，必须包含group、bmi_low、bmi_high、n、median_bmi、optimal_week、ci_low、ci_high、pbar_at_opt、expected_loss、delta_t_sigma_tech、distinct_required；同时保存图数据源：损失曲线与最优时点、t*(b) vs t_p(b)反演曲线、σ时点偏移、bootstrap边界重现热图、校准散点、ρ敏感性折线。代码末尾显式df.to_csv(MODELING_OUTPUT_DIR/results/output.csv, index=False)。

## 4. 数据使用
**数据概要**（原始数据只由代码运行时读取）：
```json
{
  "total_rows": 1082,
  "total_cols": 31,
  "issues": [
    "以下列缺失率超过 50%：['染色体的非整倍体']"
  ],
  "files": [
    {
      "path": "C:\\Users\\lingi\\Desktop\\Research\\test1\\C题\\附件.xlsx",
      "rows": 1082,
      "cols": 31,
      "issues": [
        "以下列缺失率超过 50%：['染色体的非整倍体']"
      ],
      "columns": [
        {
          "name": "序号",
          "dtype": "int",
          "missing_rate": 0.0,
          "min": 1.0,
          "max": 1082.0,
          "mean": 541.5,
          "std": 312.490799864572,
          "unique_count": 1082,
          "sample_values": [
            "1",
            "2",
            "3"
          ],
          "parse_hint": ""
        },
        {
          "name": "孕妇代码",
          "dtype": "text",
          "missing_rate": 0.0,
          "min": null,
          "max": null,
          "mean": null,
          "std": null,
          "unique_count": 267,
          "sample_values": [
            "A001",
            "A001",
            "A001"
          ],
          "parse_hint": "df['孕妇代码'].str.replace('A', '', regex=False).astype(float)"
        },
        {
          "name": "年龄",
          "dtype": "int",
          "missing_rate": 0.0,
          "min": 21.0,
          "max": 43.0,
          "mean": 28.93992606284658,
          "std": 3.6562638956254556,
          "unique_count": 21,
          "sample_values": [
            "31",
            "31",
            "31"
          ],
          "parse_hint": ""
        },
        {
          "name": "身高",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": 144.0,
          "max": 175.0,
          "mean": 161.0637707948244,
          "std": 5.23217645959195,
          "unique_count": 29,
          "sample_values": [
            "160.0",
            "160.0",
            "160.0"
          ],
          "parse_hint": ""
        },
        {
          "name": "体重",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": 53.0,
          "max": 140.0,
          "mean": 83.8948336414048,
          "std": 9.917950284907176,
          "unique_count": 513,
          "sample_values": [
            "72.0",
            "73.0",
            "73.0"
          ],
          "parse_hint": ""
        },
        {
          "name": "末次月经",
          "dtype": "text",
          "missing_rate": 0.011090573012939002,
          "min": null,
          "max": null,
          "mean": null,
          "std": null,
          "unique_count": 189,
          "sample_values": [
            "2023-02-01 00:00:00",
            "2023-02-01 00:00:00",
            "2023-02-01 00:00:00"
          ],
          "parse_hint": "pd.to_datetime(df['末次月经'])"
        },
        {
          "name": "IVF妊娠",
          "dtype": "category",
          "missing_rate": 0.0,
          "min": null,
          "max": null,
          "mean": null,
          "std": null,
          "unique_count": 3,
          "sample_values": [
            "自然受孕",
            "自然受孕",
            "自然受孕"
          ],
          "parse_hint": ""
        },
        {
          "name": "检测日期",
          "dtype": "text",
          "missing_rate": 0.0,
          "min": null,
          "max": null,
          "mean": null,
          "std": null,
          "unique_count": 379,
          "sample_values": [
            "20230429",
            "20230531",
            "20230625"
          ],
          "parse_hint": ""
        },
        {
          "name": "检测抽血次数",
          "dtype": "int",
          "missing_rate": 0.0,
          "min": 1.0,
          "max": 5.0,
          "mean": 2.487985212569316,
          "std": 1.1275475042916057,
          "unique_count": 5,
          "sample_values": [
            "1",
            "2",
            "3"
          ],
          "parse_hint": ""
        },
        {
          "name": "检测孕周",
          "dtype": "category",
          "missing_rate": 0.0,
          "min": null,
          "max": null,
          "mean": null,
          "std": null,
          "unique_count": 108,
          "sample_values": [
            "11w+6",
            "15w+6",
            "20w+1"
          ],
          "parse_hint": ""
        },
        {
          "name": "孕妇BMI",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": 20.703125,
          "max": 46.875,
          "mean": 32.28879078882743,
          "std": 2.972432009058641,
          "unique_count": 783,
          "sample_values": [
            "28.125",
            "28.515625",
            "28.515625"
          ],
          "parse_hint": ""
        },
        {
          "name": "原始读段数",
          "dtype": "int",
          "missing_rate": 0.0,
          "min": 1342544.0,
          "max": 9895358.0,
          "mean": 4692190.2356746765,
          "std": 948737.9532073574,
          "unique_count": 1082,
          "sample_values": [
            "5040534",
            "3198810",
            "3848846"
          ],
          "parse_hint": ""
        },
        {
          "name": "在参考基因组上比对的比例",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": 0.5986381,
          "max": 0.846620403,
          "mean": 0.7974720650804067,
          "std": 0.014951434940374831,
          "unique_count": 1080,
          "sample_values": [
            "0.8067259",
            "0.8063927",
            "0.8038578"
          ],
          "parse_hint": ""
        },
        {
          "name": "重复读段的比例",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": 0.02115006,
          "max": 0.04651945,
          "mean": 0.030473940258780037,
          "std": 0.0027470825667931806,
          "unique_count": 1082,
          "sample_values": [
            "0.0276035",
            "0.02827083",
            "0.03259621"
          ],
          "parse_hint": ""
        },
        {
          "name": "唯一比对的读段数  ",
          "dtype": "int",
          "missing_rate": 0.0,
          "min": 980606.0,
          "max": 7342907.0,
          "mean": 3546560.647874307,
          "std": 712817.8684176581,
          "unique_count": 1082,
          "sample_values": [
            "3845411",
            "2457402",
            "2926292"
          ],
          "parse_hint": ""
        },
        {
          "name": "GC含量",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": 0.3862499,
          "max": 0.4213731,
          "mean": 0.4006764933197782,
          "std": 0.003307198507781393,
          "unique_count": 1080,
          "sample_values": [
            "0.3992619",
            "0.3932988",
            "0.3998897"
          ],
          "parse_hint": ""
        },
        {
          "name": "13号染色体的Z值",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": -3.527318775,
          "max": 5.676687906,
          "mean": 0.31214734582347503,
          "std": 1.211392468132073,
          "unique_count": 1082,
          "sample_values": [
            "0.782096634",
            "0.692855699",
            "-0.888701998"
          ],
          "parse_hint": ""
        },
        {
          "name": "18号染色体的Z值",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": -3.262149334,
          "max": 6.076343685,
          "mean": 0.5832121545295749,
          "std": 1.2885900722044448,
          "unique_count": 1082,
          "sample_values": [
            "-2.321211659",
            "1.168520758",
            "-1.01823645"
          ],
          "parse_hint": ""
        },
        {
          "name": "21号染色体的Z值",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": -3.289375849,
          "max": 3.137182924,
          "mean": -0.11570941281977819,
          "std": 1.0989187042609343,
          "unique_count": 1082,
          "sample_values": [
            "-1.026002604",
            "-2.595098987",
            "-1.308661706"
          ],
          "parse_hint": ""
        },
        {
          "name": "X染色体的Z值",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": -3.919147889,
          "max": 7.867669632,
          "mean": 0.32517843544177455,
          "std": 1.2928932169573666,
          "unique_count": 1082,
          "sample_values": [
            "-0.062103083",
            "0.582182673",
            "-0.342563969"
          ],
          "parse_hint": ""
        },
        {
          "name": "Y染色体的Z值",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": -4.005487139,
          "max": 7.000856057,
          "mean": 0.15614253445194085,
          "std": 1.3077288631066426,
          "unique_count": 1082,
          "sample_values": [
            "-1.035610255",
            "-0.363518671",
            "-0.734502556"
          ],
          "parse_hint": ""
        },
        {
          "name": "Y染色体浓度",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": 0.010003887,
          "max": 0.234217554,
          "mean": 0.07718697843807763,
          "std": 0.03351841045704469,
          "unique_count": 1082,
          "sample_values": [
            "0.02593584",
            "0.034886856",
            "0.066171003"
          ],
          "parse_hint": ""
        },
        {
          "name": "X染色体浓度",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": -0.076508263,
          "max": 0.223932573,
          "mean": 0.057024936271719034,
          "std": 0.041464976956572246,
          "unique_count": 1082,
          "sample_values": [
            "0.038061019",
            "0.059572251",
            "0.075994548"
          ],
          "parse_hint": ""
        },
        {
          "name": "13号染色体的GC含量",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": 0.366486967,
          "max": 0.402934581,
          "mean": 0.3786933978410351,
          "std": 0.0031643940914512212,
          "unique_count": 1081,
          "sample_values": [
            "0.377068639",
            "0.3715415",
            "0.377449453"
          ],
          "parse_hint": ""
        },
        {
          "name": "18号染色体的GC含量",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": 0.378464788,
          "max": 0.412192792,
          "mean": 0.39147874860905735,
          "std": 0.0030202677864658084,
          "unique_count": 1082,
          "sample_values": [
            "0.389803052",
            "0.384770662",
            "0.390582472"
          ],
          "parse_hint": ""
        },
        {
          "name": "21号染色体的GC含量",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": 0.385214418,
          "max": 0.425052136,
          "mean": 0.40085072870332716,
          "std": 0.0038095896395479436,
          "unique_count": 1082,
          "sample_values": [
            "0.399399221",
            "0.391706139",
            "0.399479687"
          ],
          "parse_hint": ""
        },
        {
          "name": "被过滤掉读段数的比例",
          "dtype": "float",
          "missing_rate": 0.0,
          "min": 0.011982671,
          "max": 0.037834786,
          "mean": 0.023043272059149723,
          "std": 0.0034032095124037457,
          "unique_count": 1082,
          "sample_values": [
            "0.027483794",
            "0.01961667",
            "0.022312402"
          ],
          "parse_hint": ""
        },
        {
          "name": "染色体的非整倍体",
          "dtype": "category",
          "missing_rate": 0.8835489833641405,
          "min": null,
          "max": null,
          "mean": null,
          "std": null,
          "unique_count": 5,
          "sample_values": [
            "T18",
            "T13T18",
            "T21"
          ],
          "parse_hint": "df['染色体的非整倍体'].str.replace('T', '', regex=False).astype(float)"
        },
        {
          "name": "怀孕次数",
          "dtype": "category",
          "missing_rate": 0.0,
          "min": null,
          "max": null,
          "mean": null,
          "std": null,
          "unique_count": 3,
          "sample_values": [
            "1",
            "1",
            "1"
          ],
          "parse_hint": ""
        },
        {
          "name": "生产次数",
          "dtype": "int",
          "missing_rate": 0.0,
          "min": 0.0,
          "max": 3.0,
          "mean": 0.3798521256931608,
          "std": 0.6371775009776978,
          "unique_count": 4,
          "sample_values": [
            "0",
            "0",
            "0"
          ],
          "parse_hint": ""
        },
        {
          "name": "胎儿是否健康",
          "dtype": "category",
          "missing_rate": 0.0,
          "min": null,
          "max": null,
          "mean": null,
          "std": null,
          "unique_count": 2,
          "sample_values": [
            "是",
            "是",
            "是"
          ],
          "parse_hint": ""
        }
      ]
    }
  ]
}
```
**数据智能摘要**：
- 附件.xlsx是NIPT检测记录明细表，每行对应一次检测/样本记录，而非一名孕妇一行；孕妇代码唯一数远小于总行数，且存在检测抽血次数字段，说明同一孕妇可能存在多次采血/多次检测，建模时应按孕妇代码识别重复测量并处理个体内相关，不能把每行当作独立样本。
- 核心分组键为孕妇代码，样本序号只标识记录；时间相关字段包括末次月经、检测日期和检测孕周。检测孕周为文本格式（如11w+6），需解析为数值孕周（周+天数换算）后才能用于回归与达标时间计算。
- 问题1至问题3应主要针对男胎样本；理论上女胎的Y染色体Z值和Y染色体浓度应为空白，但数据概要显示这两列缺失率为0且均为数值，因此需要先在代码中核查是否存在女胎样本或性别标识，否则男胎筛选和问题4的女胎分析都会不可靠。
- Y染色体浓度是核心响应变量，为0到1之间的比例值，普通线性回归可能不适用，可考虑logit/对数变换、Beta回归或混合效应模型；检测孕周、BMI、年龄、身高、体重是直接候选解释变量。
- BMI、身高和体重存在近似确定性关系，问题3若同时纳入身高、体重和BMI，应处理共线性，避免直接并列导致估计不稳定。
- 达标时间建模需要从重复检测记录中为每名男胎孕妇推导Y染色体浓度首次达到或超过4%的最早孕周；对始终未达标的孕妇属于右删失数据，适合用生存分析或达标时间模型，而不是仅对已达标记录做普通回归。
- 检测误差会影响达标时点和分组决策，建议利用同一孕妇重复测量或多次检测评估Y浓度的噪声水平，并据此对最佳NIPT时点进行稳健性分析。
- AB列（染色体非整倍体）缺失率约88%，但题目说明空白表示无异常，属于有意义缺失，应重新编码为正常/异常而不是删除；其取值如T18、T13T18、T21为文本组合，问题4需分别构造13、18、21号染色体异常标签。
- 问题4女胎异常判定以AB列为标签，但正常与异常高度不平衡，模型评估应使用分层训练、AUC或精确率-召回率等指标，避免只看总体准确率；X染色体浓度、X染色体Z值、21/18/13号染色体Z值和各染色体GC含量、读段数相关列是候选特征。
- 测序质量相关字段（总读段数、比对比例、重复读段比例、唯一比对读段数、GC含量、各染色体GC含量、被过滤比例）主要用于检测误差分析和问题4质量控制；正常GC含量范围约40%-60%，概要显示总GC含量范围约0.386-0.421，存在低于40%的记录，建模前需决定是否过滤。
- 末次月经和检测日期可用于推算或校验孕周，但末次月经有少量缺失，检测日期为文本格式，处理时需转换为日期并谨慎处理解析失败。
**已完成小题的结果（供本题复用）**：
- 小题 1（passed）：C:\Users\lingi\Desktop\Research\test1\outputs\results\q1.csv, C:\Users\lingi\Desktop\Research\test1\outputs\results\output.csv

## 5. 预期图表
- fig_roadmap [roadmap] 展示问题2从数据到分组时点决策的完整技术路线（数据来源：）
- fig_q1_scatter [scatter] （数据来源：）
- fig_q1_smooth_ga [line] （数据来源：）
- fig_q1_smooth_bmi_int [heatmap] （数据来源：）
- fig_q1_quantile_curves [line] （数据来源：）
- fig_q1_prob_curves [line] （数据来源：）
- fig_diag_resid [scatter] （数据来源：）
- fig_sens_dist [line] （数据来源：）
- fig_sens_interaction [line] （数据来源：）
- fig_sens_gc [line] （数据来源：）
- fig_sens_marginal [line] （数据来源：）
- fig_sens_ga_window [line] （数据来源：）
- fig_anchor_threshold [scatter] （数据来源：）
- fig_q2_ga_bmi_prob_curves [line] 回答BMI对Y染色体浓度达标概率的影响方向，为分组提供依据（数据来源：results/q2_prob_curves.csv）
- fig_q2_loss_curves_optimal [line] 展示分组后各组损失曲线及t_k*位置，支撑分组时点决策（数据来源：results/q2_group_loss_curves.csv）
- fig_q2_bmi_bins_tstar [scatter] 展示个体最优时点与BMI的关系、组边界与组最佳时点的单调性（数据来源：results/q2_individual_tstar.csv）
- fig_q2_error_shift_sigma [line] 量化logit尺度测量误差对最佳时点的影响，支撑误差传播结论（数据来源：results/q2_sigma_shift.csv）
- fig_q2_bootstrap_boundary_heatmap [heatmap] 评估BMI分箱边界的稳定性，支撑分组方案可靠性（数据来源：results/q2_bootstrap_boundaries.csv）
- fig_q2_calibration [scatter] 验证模型概率在GA带上的内部校准，说明概率层可信度（数据来源：results/q2_calibration.csv）
- fig_q2_rho_sensitivity [line] 分析主观参数ρ对分组时点决策的影响，识别主要不确定性来源（数据来源：results/q2_rho_sensitivity.csv）
- fig_q2_fnr_fpr_sigma [line] 展示检测误差导致的假阴性与假阳性率变化，支撑检测误差影响分析（数据来源：results/q2_error_classification.csv）

## 6. 预期表格
- table_model_comparison：模型比较结果（列：待定；）
- table_smooth_terms：平滑项显著性检验（列：待定；）
- table_random_effects：随机效应方差组分与ICC（列：待定；）
- table_covariate_forms：BMI/体重/身高+体重三种协变量形态比较（列：待定；）
- table_sens_gc：GC处理策略敏感性（列：待定；）
- table_sens_ga_window：孕周外推界限敏感性（列：待定；）
- table_sens_marginal：边缘与条件达标概率差异（列：待定；）
- table_quantile_check：分位数回归与Beta-GAMM达标概率一致性（列：待定；）
- tab_q2_main_results：男胎孕妇BMI分组与最佳NIPT时点主结果表（列：group, bmi_low, bmi_high, n, median_bmi, optimal_week, ci_low, ci_high, pbar_at_opt, expected_loss, delta_t_sigma_tech, distinct_required；给出每组BMI区间、样本量和最佳NIPT时点及其不确定性，作为问题2核心交付）
- tab_q2_model_selection：简约Beta混合模型固定效应形式选择（列：model_form, cv_rmse, mae, selected；用孕妇分组五折CV比较GA分段线性与纯线性，支撑模型选择）
- tab_q2_sensitivity_rho：ρ敏感性：不同损失比下的最佳时点（列：rho, group, optimal_week；扫描ρ∈{0.25,0.5,1,2,4}，展示主观参数对时点的影响）
- tab_q2_sensitivity_sigma：σ四档测量误差时点偏移与错分率（列：sigma, group, optimal_week_sigma, delta_t, FNR, FPR；量化σ∈{0,0.5σ_tech,σ_tech,2σ_tech}时各组时点偏移和错分率）
- tab_q2_sensitivity_bins：n_min与K敏感性：分箱边界与总损失（列：n_min, K, boundaries, total_loss_J, max_tstar_gap；扫描n_min∈{10,15,20,30}与K∈{2..6}，评估分箱方案稳定性）
- tab_q2_sensitivity_key_assumptions：关键假设扰动与对照实验汇总（列：assumption_id, base_spec, perturbed_spec, max_delta_t_k, boundary_change, conclusion；对关键假设3-9逐条做扰动/对照，回答假设不成立时结论是否仍稳）

## 7. 结果契约
```json
{
  "description": "问题2主结果表：每个BMI分组一行，给出BMI区间、最佳NIPT时点及不确定性",
  "allow_single_row": false,
  "min_rows": 2,
  "max_rows": 6,
  "columns": [
    {
      "name": "group",
      "dtype": "category",
      "min": null,
      "max": null,
      "distinct_required": false,
      "description": "BMI分组编号"
    },
    {
      "name": "bmi_low",
      "dtype": "float",
      "min": 15.0,
      "max": 50.0,
      "distinct_required": false,
      "description": "分组BMI下界"
    },
    {
      "name": "bmi_high",
      "dtype": "float",
      "min": 15.0,
      "max": 50.0,
      "distinct_required": false,
      "description": "分组BMI上界"
    },
    {
      "name": "n",
      "dtype": "int",
      "min": 10.0,
      "max": 267.0,
      "distinct_required": false,
      "description": "组内孕妇数"
    },
    {
      "name": "median_bmi",
      "dtype": "float",
      "min": 15.0,
      "max": 50.0,
      "distinct_required": false,
      "description": "组内b_i中位数"
    },
    {
      "name": "optimal_week",
      "dtype": "float",
      "min": 10.0,
      "max": 25.0,
      "distinct_required": true,
      "description": "组最佳NIPT时点，周"
    },
    {
      "name": "ci_low",
      "dtype": "float",
      "min": 10.0,
      "max": 25.0,
      "distinct_required": false,
      "description": "最佳时点95%区间下界"
    },
    {
      "name": "ci_high",
      "dtype": "float",
      "min": 10.0,
      "max": 25.0,
      "distinct_required": false,
      "description": "最佳时点95%区间上界"
    },
    {
      "name": "pbar_at_opt",
      "dtype": "float",
      "min": 0.0,
      "max": 1.0,
      "distinct_required": false,
      "description": "最佳时点组平均达标概率"
    },
    {
      "name": "expected_loss",
      "dtype": "float",
      "min": 0.0,
      "max": 10.0,
      "distinct_required": false,
      "description": "最佳时点组期望损失"
    },
    {
      "name": "delta_t_sigma_tech",
      "dtype": "float",
      "min": -3.0,
      "max": 3.0,
      "distinct_required": false,
      "description": "σ_tech引起的时点偏移，周"
    },
    {
      "name": "distinct_required",
      "dtype": "category",
      "min": null,
      "max": null,
      "distinct_required": false,
      "description": "是否达到可报告组间时点差异；False时需附统一时点统计说明"
    }
  ],
  "allow_extra_columns": true
}
```

## 8. 实现约束
- 只允许使用：numpy、pandas、scipy、sklearn、statsmodels、matplotlib、networkx、pulp
- 代码总执行时间不超过 90 秒
- 结果必须写入 MODELING_OUTPUT_DIR/results/output.csv
- 数值常量必须与题目 problem_facts 一致，不得自创物理参数

## 9. 最近一次失败（必须针对性修复，不得生成相同代码）
```
LTM 中出现 '60.0 %'，但题目常量为 [4.0] %（可能 LLM 记错了数值）
```

## 10. 编程手交付要求
- 只在本任务目录编写 `solution.py`，不要修改任何建模设定文件。
- `solution.py` 必须是完整可执行的 Python 代码，遵守第 8 节约束。
- 如需图表，另写 `figures.py`，把图片保存到 `figures/` 子目录（如 `figures/figure1.png`）。
- V17：`figures.py` 的图片必须按 figures_plan 的 `id` 命名（如 `figures/fig_q1_corr.png`，文件名 = plan.id）；
  系统按文件名把图登记到图表注册表，未按 plan_id 命名的图片不会被论文引用。
- 数据文件路径通过环境变量 `MODELING_DATA_PATHS`（JSON 数组）与 `MODELING_DATA_PATH`（第一个文件）传入；
  多附件时按第 4 节数据概要中的文件边界分别读取，不要假设已合并成一张表。
- 不要读取/写入任务包之外的原始数据之外的内容；数据路径见第 4 节。

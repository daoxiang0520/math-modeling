# 建模方案与实现架构说明书

## 0. 当前小题
**当前小题（3/4）**：问题3  男胎Y 染色体浓度达标时间受多种因素(身高、体重、年龄等)的影响，试综合考虑这些因
素、检测误差和胎儿的Y 染色体浓度达标比例（即浓度达到或超过4%的比例），根据男胎孕妇的BMI，
给出合理分组以及每组的最佳NIPT 时点，使得孕妇潜在风险最小，并分析检测误差对结果的影响。

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
**建模目标**：对男胎孕妇综合考虑身高、体重、年龄、检测误差和Y染色体浓度≥4%的达标比例，按BMI给出合理分组及每组最佳NIPT时点，使延迟风险最小，并分析检测误差对分组时点的影响。

## 2. 建模设定（动态 LTM，编程手不得修改）
**假设**：
- 假设1：分析对象限定为男胎样本；以孕妇代码、Y染色体浓度（V列）与Y染色体Z值（U列）同时具有有效数值交叉核查男胎，不能仅依据缺失率判定女胎。依据：数据概要显示Y染色体浓度列缺失率为0，与“女胎数据此列为空白”存在矛盾；风险：误纳女胎会引入Y≈0噪声；可验证：交叉表核对Y浓度/Y-Z值与孕妇代码、X染色体浓度及胎儿是否健康列的一致性。
- 假设2：数据结构为孕妇层纵向重复测量，每行是一次检测记录，同一孕妇可能多次检测；分组、达标时间与不确定性推断均以孕妇为聚类单位，不得把行当作独立样本。依据：数据智能摘要显示1082条检测记录、267位孕妇，且年龄列lag-1自相关约0.999（该值为数据自相关估计，非题目常量）；风险：行级独立会低估标准误并产生信息泄漏；可验证：行级独立模型与孕妇随机效应/聚类bootstrap的标准误差异。
- 假设3：孕周变量必须解析为数值周数，并分离天数符号与二元达标事件符号：ga_ij=w_ij+day_ij/7，day_ij为天数部分；D_ij=1(y_obs,ij≥y_thr)为二元达标事件指示，禁止再使用d_ij同时表示天数与达标事件。依据：人类架构审核指出旧LTM中d_ij符号重载；风险：符号重载会使核心公式不可推导；可验证：输出解析后孕周与该孕妇记录行核对，公式中不再出现歧义d_ij。
- 假设4【关键假设】：主交付采用约束风险最小化语义：对给定个体协变量x_i，最佳时点为最小化延迟风险r_delay(t)，同时满足边缘达标保证π_marg(t,x_i)≥p；当r_delay(t)在[10.0,25.0]上严格递增且π_marg(t,x_i)对t单调不减时，该约束最优解等价于达标保证首穿时点t_p,i=inf{t∈[10.0,25.0]:π_marg(t,x_i)≥p}，主值p=0.80。依据：题目要求“使得潜在风险最小”，人类审核要求解决主交付与风险最小化目标的矛盾；风险：若π_marg非单调或r_delay非严格递增，该等价性失效，需回退为显式约束优化求解；可验证：在0.1周网格上同时计算约束最优解与首穿时点，比较二者差异；并对单调度做数值检查。
- 假设5【关键假设】：个体边缘达标概率π_marg(t,x_i)在t∈[10.0,25.0]上对t单调不减；主模型不使用PAVA后处理，而是在0.1周网格上数值断言差分Δπ_marg≥0。若违反，则回退为固定效应孕周线性且不引入孕周随机斜率，或对孕周随机斜率截断使孕周效应保持为正。依据：t_p首穿时点适定性与约束风险最小化等价性均依赖单调性；风险：随机斜率为负可能破坏个体曲线单调性；可验证：网格差分检验、截断前后t_p分布比较，以及与分段线性模型S2的RMSE比较。
- 假设6【关键假设】：临床达标阈值y_thr=0.04（原文：如果男胎的Y染色体浓度达到或高于4%）；决策保证水平p_guarantee=0.80与浓度阈值严格区分。敏感性p∈{0.75,0.85,0.90}。依据：题目常量与人类架构反馈；风险：Beta尾部参数失配会影响达标概率；可验证：DHARMa类残差诊断、Beta与高斯/logit-GAM基准比较、不同p值下t_p偏移。
- 假设7【关键假设】：问题3主模型采用“主二项模型+辅助Beta混合模型”的二元-连续双通道结构。主二项模型直接对达标事件D_ij建模；辅助Beta模型对观测Y浓度y_obs,ij建模，用于生成边缘达标比例π_marg并传播测量误差。线性预测器默认形态为“体重+身高+年龄+IVF”（体重优先，身高校正），但必须显式比较BMI/体重/身高+体重/三者全入四种形态，以AIC/BIC/孕妇层CV和1SE规则确定最终报告形态；残差化身高体重仅作敏感性，且残差估计不确定不传播进主结果。依据：qiao 2023显示IVF与胎儿分数负相关，hou 2015提示体重优于BMI，mousavi 2022显示二元LFF与连续FF可能受不同因素支配；风险：BMI与体重/身高共线性会导致系数不稳定；可验证：VIF、concurvity、四种形态AIC/BIC、孕妇层五折CV RMSE。
- 假设8【关键假设】：纵向重复测量下随机效应按设计进入：主二项模型和辅助Beta模型均设置孕妇随机截距与中心化孕周随机斜率；报告方差组分、ICC及随机斜率是否必要的比较。依据：同一孕妇多次检测的纵向结构；风险：缺失随机斜率会低估个体间达标时间差异；可验证：带/不带随机斜率的AIC/BIC、方差组分及CI比较。
- 假设9【关键假设】：边缘达标概率π_marg(t,x_i)定义为对随机效应和logit尺度测量误差的积分：π_marg(t,x_i)=E_{b_i,z,ε}[1(expit(logit(z)+ε)≥0.04)|t,x_i]；其组内均值π_g(t)=n_g^{-1}Σ_{i∈g}π_marg(t,x_i)即题目所称“Y染色体浓度达标比例”的分组估计，并报告bootstrap置信区间。依据：人类反馈要求π_g为FF分布积分，并对照外部量级；风险：若忽略随机效应或测量误差会高估组达标比例；可验证：提取多组报告π_g(t)及其CI，与Ashoor 2013体重-FF<4%曲线趋势和Wright 2015 FF=4%检出率约62%（该检出率为文献值，非题目常量）做量级对照。
- 假设10【关键假设】：检测误差双通道分离。通道A为浓度测量误差：在logit尺度将ε~N(0,σ²)卷积，σ∈{0,0.5σ_tech,σ_tech,2σ_tech}，进入主达标概率敏感性和错分率κ_σ；σ_tech沿用问题1技术重复估计，定义为同一孕妇同一抽血次数技术重复的logit(y)组内标准差。技术失败类误差（如平台偏差）作独立敏感性，不混入通道A。必须显式声明σ_tech可能低估总误差。依据：palomaki 2018区分低FF与技术失败，人类反馈要求独立处理；风险：只报告技术重复误差会低估真实决策不确定性；可验证：σ_tech与纵向误差/GC技术失败敏感性对比，报告κ_σ随σ变化。
- 假设11【关键假设】：GC含量不采用40.0%–60.0%硬阈值剔除（原文：正常GC含量范围为40% ~ 60%）；因实测GC全距约0.386–0.421，约42%（该比例为数据比例推断值，非题目常量）行低于0.40、无高于0.60，判断为检测平台系统偏差。GC仅作为连续质量协变量或质量权重进入敏感性分析；技术失败类偏差独立敏感性。风险：硬剔除会损失大量记录并引入选择偏差；可验证：纳入GC质量权重前后t_p偏移及技术失败行影响。
- 假设12【关键假设】：推荐窗口端点采用右删失语义。个体t_p,i在[10.0,25.0]内无穿越时记t_p,i=25.0，指标c_i=1表示窗口内未达80%保证；报告n_unsolved和删失比例r_cens。组推荐仍给出具体周数，不降级为方向性结论；若r_cens>0.20，以未删失个体中位数为主口径并显式标注。依据：人类反馈要求保留具体时点并采用右删失；风险：右删失会压缩组内变异；可验证：报告删失/未删失时点差异与上下界25/26周敏感性。
- 假设13【关键假设】：必须报告个体达标时点分布形态，包括直方图、分位数和双峰/多峰检测；若出现双峰，需说明是否由协变量分层或删失造成。依据：q2个体t_p80曾现双峰，人类反馈要求Q3继续报告；风险：双峰分布会使组中位数掩盖重要亚组结构；可验证：直方图与高斯混合模型BIC比较，按BMI/年龄/IVF分层检查峰值。
- 假设14【关键假设】：BMI分组采用数据驱动K∈{2,3,4}，边界和K选择使用一标准误差规则；不得把样本内最优边界表述为精确医学阈值。分组验收必须报告：组间时点差异点估计及bootstrap CI、组内异质性下降、损失改进、边界bootstrap重现频率。依据：q2固定边界30.0 kg/m²重现频率仅24%（该值为历史运行诊断值，非题目常量），人类反馈要求1SE；风险：过细分组会过拟合而不稳定；可验证：1SE与最佳损失选K比较、边界重现频率报告。
- 假设15：组推荐时点t_g采用与个体主交付一致的约束风险最小化语义：t_g=argmin_{t∈[10.0,25.0]}r_delay(t) s.t. π_g(t)≥0.80，在π_g单调下等价于t_g=inf{t∈[10.0,25.0]:π_g(t)≥0.80}；同时报告组内个体t_p,i中位数作为描述性结果。依据：人类反馈要求组级也直接满足达标比例约束；风险：个体中位数与约束最优在右删失或非单调时可能不一致；可验证：比较组中位数与组级首穿一致性，报告差异CI。
- 假设16：延迟风险函数采用连续增函数r_delay(t)=(t−10.0)/17，定义域t∈[10.0,27.0]，分母17由原文窗口10周可检、中期发现13－27周风险高计算；不采用“12周及以下零风险平台”。该函数兼容原文风险等级：早期发现12周以内风险较低、13－27周风险高、28周以后风险极高（原文：早期发现12周以内风险较低；中期发现13－27周风险高；晚期发现28周以后风险极高）。风险：线性风险是最简假设；可验证：指数风险或分段风险与主结果对比。
- 假设17：软损失参考列t*仅作二级风险视角，不作为主交付：L_{ρ,γ}(t;x_i)=r_delay(t)+ρ[1−π_marg(t,x_i)]+γ·κ_σ(t,x_i)，搜索域t∈[10.0,25.0]；ρ∈{0.5,1,2}，γ∈{0,0.5,1,2}扩展敏感性。依据：人类反馈第一条将argmin损失降为二级参考列；风险：权重γ/ρ假设会影响二级参考值；可验证：ρ×γ网格下t*偏移报告。
- 假设18：身高、体重列入模型前必须先画像：检查缺失率、范围、BMI=体重/身高²自洽性；BMI个体代表值取个体内中位数，身高/体重取个体内中位数以减少检测时间点波动。依据：数据智能摘要提示身高体重无缺失但分布非正态，且BMI可直接使用；风险：BMI与体重/身高不一致的记录会污染分组或模型；可验证：输出自洽性异常记录数，剔除前后主结果变化。
- 假设19：AE列“胎儿是否健康”为出生后结果标签，仅作敏感性分层，明确禁止进入问题3实时决策主模型；对其中不健康记录做剔除/保留敏感性。依据：deng 2023显示核型不改变FF，支持保留38条（该数值为标签统计值，非题目常量）做敏感性；风险：不当剔除会改变协变量效应；可验证：剔除AE不健康记录前后t_p及分组结果比较。
- 假设20：Coder必须按parse_hints解析字符串列：孕妇代码用df['孕妇代码'].str.replace('A', '', regex=False).astype(float)，末次月经用pd.to_datetime，染色体的非整倍体用df['染色体的非整倍体'].str.replace('T', '', regex=False).astype(float)，检测孕周自行解析w+d格式。依据：机器数据列解析建议及人类反馈；风险：字符串解析错误会导致样本错误；可验证：输出解析后唯一值与前若干行核对。
- 假设21【关键假设】：结果解释须主动讨论“二元达标事件与连续Y染色体浓度可能受不同因素支配”，以回应mousavi 2022荟萃中BMI/体重对低FF不显著但与FF连续关系可能仍存在的矛盾；方法引用应成组给出Beta回归（Ferrari 2004）、GAMM/Wood 2017、随机效应纵向模型（Laird & Ware 1982）；McKanna FFBR作为按协变量建立FF分布的方法对照；Scheffer 2021作为低FF不良妊娠风险与延迟风险r_delay的临床依据。依据：人类反馈第8条；风险：若仅用一个连续模型或仅用二元模型，会遗漏信号；可验证：主二项与辅助Beta模型中BMI/年龄/体重系数的符号和显著性差异报告。
**符号表**：
- N_rec: 男胎检测记录总数，原始总记录1082行，男胎筛选后待定
- n_male: 男胎孕妇数量，原始总孕妇267位，男胎筛选后待定
- n_i: 第i位孕妇的检测次数
- i: 孕妇索引，i=1,...,n_male
- j: 孕妇内检测索引，j=1,...,n_i
- ga_ij: 数值化孕周，单位：周，ga_ij=w_ij+day_ij/7
- w_ij: 检测孕周文本的整数周部分
- day_ij: 检测孕周文本的天数部分，0≤day_ij<7；纯周格式day_ij=0
- t: 决策用连续孕周变量，单位：周，主决策域[10.0,25.0]
- t_min: 可检测窗口下界，t_min=10.0周（原文：孕期在10周~25周之间可以检测）
- t_max: 可检测窗口上界，t_max=25.0周（原文：孕期在10周~25周之间可以检测）
- t_risk_end: 风险高窗口上界，t_risk_end=27.0周（原文：中期发现13－27周风险高）
- t_risk_extreme: 风险极高起点，t_risk_extreme=28.0周（原文：晚期发现28周以后风险极高）
- t_risk_low_end: 早发现低风险窗口上界，t_risk_low_end=12.0周（原文：早期发现12周以内风险较低）
- y_obs_ij: 第i位孕妇第j次检测的观测Y染色体浓度，无量纲比例，0<y_obs_ij<1
- z_ij: 真实Y染色体浓度随机变量，z_ij∈(0,1)
- D_ij: 二元达标事件指示，D_ij=1(y_obs_ij≥y_thr)
- y_thr: 临床达标阈值，y_thr=0.04，即4.0%
- p: 达标保证概率水平，主值p=0.80
- bmi_ij: 第i位孕妇第j次记录的BMI，单位：kg/m²
- bmi_rep_i: 第i位孕妇个体代表BMI，取个体内中位数，单位：kg/m²
- height_rep_i: 第i位孕妇个体代表身高，取个体内中位数，单位：cm
- weight_rep_i: 第i位孕妇个体代表体重，取个体内中位数，单位：kg
- age_i: 第i位孕妇年龄，单位：岁
- ivf_i: 第i位孕妇IVF妊娠指示，0/1，1表示IVF妊娠
- x_i: 第i位孕妇决策协变量向量，默认包含(age_i, height_rep_i, weight_rep_i, ivf_i)；BMI用于分组轴
- z_mi: 建模协变量分量，第m个固定效应协变量，形态根据四种协变量比较选择
- μ_ij: 辅助Beta分布条件均值，0<μ_ij<1
- φ: 辅助Beta分布精度参数，φ>0
- η_ij: 辅助Beta模型logit线性预测器，η_ij=logit(μ_ij)
- π_ij: 主二项模型中D_ij=1的条件概率
- π_marg(t,x_i): 边缘达标概率，对随机效应和测量误差积分后的P(y_obs≥0.04|t,x_i)
- π_g(t): 第g组组级达标比例估计，π_g(t)=n_g^{-1}Σ_{i∈g}π_marg(t,x_i)
- n_g: 第g组男胎孕妇数
- g: BMI分组索引，g=1,...,K
- K: 分组数量，K∈{2,3,4}
- c_K: K分组所含BMI边界向量，单位：kg/m²
- α0: 主二项模型固定截距
- α_ga: 主二项模型孕周固定效应系数
- α_m: 主二项模型第m个协变量固定效应系数
- β0: 辅助Beta模型固定截距
- β_ga: 辅助Beta模型孕周固定效应系数
- β_m: 辅助Beta模型第m个协变量固定效应系数
- u0i: 主二项模型第i位孕妇随机截距
- u1i: 主二项模型第i位孕妇对中心化孕周(ga_ij−ḡ)的随机斜率
- u_i: 主二项模型随机效应向量，u_i=(u0i,u1i)'
- b0i: 辅助Beta模型第i位孕妇随机截距
- b1i: 辅助Beta模型第i位孕妇对中心化孕周(ga_ij−ḡ)的随机斜率
- b_i: 辅助Beta模型随机效应向量，b_i=(b0i,b1i)'
- Σ_u: 主二项随机效应协方差矩阵
- Σ_b: 辅助Beta随机效应协方差矩阵
- σ_u0²: 主二项随机截距方差
- σ_u1²: 主二项随机斜率方差
- σ_u01: 主二项随机截距与斜率协方差
- σ_b0²: 辅助Beta随机截距方差
- σ_b1²: 辅助Beta随机斜率方差
- σ_b01: 辅助Beta随机截距与斜率协方差
- ḡ: 男胎样本孕周均值，用于随机斜率中心化，单位：周
- ε_ij: 通道A测量误差项，logit尺度，ε_ij~N(0,σ²)
- σ: 通道A测量误差标准差，σ∈{0,0.5σ_tech,σ_tech,2σ_tech}
- σ_tech: 基准技术重复测量误差标准差，logit(y)尺度，沿用问题1技术重复估计
- κ_σ(t,x_i): 在时点t和协变量x_i下，由测量误差σ引起的错分率
- r_delay(t): 延迟风险函数，r_delay(t)=(t−10.0)/17，t∈[10.0,27.0]
- L_{ρ,γ}(t,x_i): 二级软损失，L=r_delay+ρ(1−π_marg)+γ·κ_σ
- ρ: 达标缺失损失权重，ρ∈{0.5,1,2}
- γ: 错分率损失权重，γ∈{0,0.5,1,2}
- t_p,i: 第i位孕妇主交付最佳时点，单位：周
- t_g: 第g组组推荐最佳时点，单位：周
- c_i: 第i位孕妇右删失指示，c_i=1表示[10,25]内无穿越并记25.0，c_i=0表示窗口内有穿越
- n_unsolved: 组内t_p,i右删失人数
- r_cens: 组内删失比例，r_cens=n_unsolved/n_g
- C: BMI边界候选集，单位：kg/m²，C={24,26,28,30,32,34,36}
- Q(c_K): 边界c_K对应的总分组损失
- SE_boot(Q): 总分组损失的bootstrap标准误
- B: 通道B cluster bootstrap次数，B=100
- CI_low_g: 第g组组推荐时点bootstrap区间下界
- CI_high_g: 第g组组推荐时点bootstrap区间上界
- fr(c): 候选BMI边界c的bootstrap重现频率
- AIC: 赤池信息准则
- BIC: 贝叶斯信息准则
- RMSE: 按孕妇分组交叉验证的均方根误差
- VIF: 方差膨胀因子
- edf: 平滑项估计自由度；Q3默认主模型为线性混合，若纳入平滑仅在敏感性中
**公式/方程**：
- 孕周数值解析：ga_ij = w_ij + day_ij/7；纯周格式day_ij=0，0≤day_ij<7。
- logit变换与反变换：logit(u)=ln(u/(1-u))，定义域u∈(0,1)；expit(v)=1/(1+exp(-v))。
- 二元达标事件定义：D_ij=1(y_obs_ij≥y_thr)，y_thr=0.04。
- 主二项模型：D_ij~Bernoulli(π_ij)，ηD_ij=logit(π_ij)=α0+α_ga·ga_ij+Σ_m α_m·z_mi+u0i+u1i·(ga_ij−ḡ)。
- 辅助Beta连续模型：y_obs_ij~Beta(μ_ij,φ)，η_ij=logit(μ_ij)=β0+β_ga·ga_ij+Σ_m β_m·z_mi+b0i+b1i·(ga_ij−ḡ)。
- 随机效应分布：u_i=(u0i,u1i)'~N(0,Σ_u)，b_i=(b0i,b1i)'~N(0,Σ_b)；Σ_u=[[σ_u0²,σ_u01],[σ_u01,σ_u1²]]，Σ_b=[[σ_b0²,σ_b01],[σ_b01,σ_b1²]]。
- 测量误差模型（通道A）：真实浓度z~Beta(μ_ij,φ)，观测值为y_obs=expit(logit(z)+ε)，ε~N(0,σ²)。
- 边缘达标概率：π_marg(t,x_i)=E_{b_i,z,ε}[1(expit(logit(z)+ε)≥0.04)|ga=t,x_i]。
- 数值单调性断言：对t_k∈{10.0,10.1,...,25.0}，要求π_marg(t_k,x_i)−π_marg(t_{k−1},x_i)≥0；若违反，执行假设5回退。
- 个体约束风险最小化等价式：t_p,i=argmin_{t∈[10.0,25.0]} r_delay(t) s.t. π_marg(t,x_i)≥p = inf{t∈[10.0,25.0]:π_marg(t,x_i)≥p}；主值p=0.80。等价成立条件：r_delay(t)严格递增且π_marg(t,x_i)单调不减。
- 个体右删失规则：若{t∈[10.0,25.0]:π_marg(t,x_i)≥0.80}=∅，则t_p,i:=25.0且c_i=1。
- 组级达标比例：π_g(t)=n_g^{-1}Σ_{i∈g}π_marg(t,x_i)。
- 组级约束风险最小化等价式：t_g=argmin_{t∈[10.0,25.0]} r_delay(t) s.t. π_g(t)≥p = inf{t∈[10.0,25.0]:π_g(t)≥p}；主值p=0.80。
- 延迟风险函数：r_delay(t)=(t−10.0)/17，t∈[10.0,27.0]；分母17=27.0−10.0，由原文10周可检与13－27周高风险窗口端点计算。
- 通道A错分率：κ_σ(t,x_i)=E_{b_i,z,ε}[1{(z≥0.04且z_obs<0.04)或(z<0.04且z_obs≥0.04)}]，其中z_obs=expit(logit(z)+ε)。
- 二级软损失参考时点：L_{ρ,γ}(t,x_i)=r_delay(t)+ρ[1−π_marg(t,x_i)]+γ·κ_σ(t,x_i)；t*_i=argmin_{t∈[10.0,25.0]}L_{ρ,γ}(t,x_i)；ρ∈{0.5,1,2}，γ∈{0,0.5,1,2}。
- 分组边界总损失：Q(c_K)=Σ_{g=1}^{K}Σ_{i∈g}L_{ρ,γ}(t_g,x_i)，其中t_g由该边界分组下组级约束最优给出。
- 一标准误差分组选择：在K∈{2,3,4}中选择最小的K，使得Q(c_K)≤min_{K'}Q(c_{K'})+SE_boot(Q(c_{K'})）；边界c_K在相应K下按1SE选择。
- 边界bootstrap重现频率：fr(c)=B^{-1}Σ_{b=1}^{B}1(c*_b=c)，其中c*_b为第b次cluster bootstrap重拟合后的最优边界。
- 通道B组时点区间：CI_low_g=Q_{0.025}({t_g^b}_{b=1}^{B})，CI_high_g=Q_{0.975}({t_g^b}_{b=1}^{B})。
- 误差敏感性时点偏移：delta_t_sigma_tech=t_g,p_σ−t_g,p_0，其中p=0.80，σ∈{0,0.5σ_tech,σ_tech,2σ_tech}。
**解题思路**：步骤1：读取附件.xlsx。Coder必须按parse_hints解析字符串列：孕妇代码用df['孕妇代码'].str.replace('A', '', regex=False).astype(float)，末次月经用pd.to_datetime，染色体的非整倍体用df['染色体的非整倍体'].str.replace('T', '', regex=False).astype(float)，检测孕周自行解析w+d格式。步骤2：筛选男胎样本，以孕妇代码为组单位构建纵向数据；计算ga_ij=w_ij+day_ij/7，并用末次月经与检测日期交叉核对，容差±1周；个体代表值：BMI取个体内中位数，身高/体重取个体内中位数，年龄取个体众数或中位数，IVF取个体众数。步骤3：对身高/体重/BMI做前置画像：缺失率、范围、BMI=体重/身高²自洽性异常数；AE不健康记录仅作敏感性分层，不进入主决策模型。步骤4：拟合主二项模型D_ij~Bernoulli(π_ij)和辅助Beta模型y_obs_ij~Beta(μ_ij,φ)；固定效应包含孕周、协变量（默认体重+身高+年龄+IVF），随机效应为随机截距+中心化孕周随机斜率；报告两种模型系数及显著性。步骤5：显式比较四种协变量形态BMI/体重/身高+体重/三者全入，计算AIC/BIC、孕妇层五折CV RMSE、VIF/concurvity；以1SE规则确定主报告形态；对体重做log变换或分段变换作为对照；残差化身高体重仅作敏感性，且不传播残差不确定性。步骤6：估计边缘达标概率π_marg(t,x_i)，在10.0–25.0周、0.1周网格上做数值单调性断言；若违反，执行随机斜率截断或回退固定效应孕周线性。步骤7：对每个男胎孕妇计算约束风险最小化时点t_p,i=inf{t∈[10,25]:π_marg(t,x_i)≥0.80}；检查与显式argmin r_delay(t) subject to π_marg≥0.80的数值等价性；p敏感性取0.75/0.85/0.90；窗口内无解记t_p,i=25.0、c_i=1，报告n_unsolved和r_cens；绘制个体t_p分布直方图、分位数和双峰检测。步骤8：对BMI进行数据驱动分组K∈{2,3,4}，边界候选按分位数/网格生成，以总分组损失Q(c_K)和一标准误差规则选K与边界；报告组间时点差异点估计及bootstrap CI、组内异质性下降、损失改进、边界bootstrap重现频率；不把最优边界表述为精确医学阈值。步骤9：对每组计算组级达标比例π_g(t)和组推荐时点t_g=inf{t:π_g(t)≥0.80}，输出组级π_g(t)的bootstrap置信区间，并与Ashoor 2013体重-FF<4%曲线趋势、Wright 2015检出率量级进行外部对照。步骤10：通道B用孕妇层cluster bootstrap，B=100，重抽样并重拟合模型后输出每位个体和每组的t_p,t_g的95%分位区间；通道A独立进行σ∈{0,0.5σ_tech,σ_tech,2σ_tech}卷积，输出κ_σ错分率、delta_t_sigma_tech偏移；技术失败/GC平台偏差独立敏感性，明确σ_tech可能低估总误差。步骤11：二级软损失参考列：在t∈[10,25]上对每个个体和组中位数用ρ∈{0.5,1,2}×γ∈{0,0.5,1,2}网格计算t*，作为风险视角参考，不与主交付混淆。步骤12：输出主结果表group,bmi_low,bmi_high,n,median_bmi,t_g,ci_low,ci_high,pi_g_at_tg,median_uncensored,n_unsolved,r_cens,t_star,delta_t_sigma_tech；输出个体时点分布、组达标比例曲线、四种协变量形态比较、误差敏感性和边界bootstrap热图；论文主动讨论二元达标事件与连续Y浓度受不同因素支配。

## 3. 算法与求解
**算法摘要**：
**伪代码/实现步骤**：
1. 步骤1 load_and_parse: 读取附件.xlsx；按parse_hints解析字符串列：df['孕妇代码']=df['孕妇代码'].str.replace('A','',regex=False).astype(float)；末次月经=pd.to_datetime；染色体的非整倍体=df['染色体的非整倍体'].str.replace('T','',regex=False).astype(float)；检测孕周解析w+d格式为w_ij和day_ij，纯周格式day_ij=0，计算ga_ij=w_ij+day_ij/7。
2. 步骤2 male_sample_selection: 以(V列Y染色体浓度有效数值)且(U列Y染色体Z值有效数值)交叉核查男胎，不依据缺失率判定女胎；按孕妇代码分组构建纵向数据；计算个体代表值：BMI、身高、体重取个体内中位数，年龄取个体中位数，IVF取个体众数。
3. 步骤3 covariate_image_QC: 画像身高、体重、BMI缺失率、范围、BMI=体重/身高²自洽性；输出异常记录数；AE列胎儿是否健康仅作敏感性分层，不进入主模型；GC含量不作为硬阈值剔除，按连续质量协变量或权重仅进入敏感性分析。
4. 步骤4 event_and_features: 构建D_ij=1(y_obs_ij≥0.04)；生成建模特征：中心化孕周ga_c=ga_ij−mean(ga)；个体代表身高、体重、BMI、年龄、IVF；构造四种协变量形态：BMI、体重、身高+体重、三项全入。
5. 步骤5 fit_primary_binomial: 对D_ij拟合主二项模型，线性预测器=α0+α_ga·ga_c+协变量项+随机截距+随机斜率；优先使用statsmodels.genmod.BayesMixedGLM.BinomialBayesMixedGLM，失败时退化为GLM+孕妇聚类稳健标准误；报告固定效应系数、显著性、方差组分。
6. 步骤6 fit_auxiliary_beta: 对y_obs_ij拟合辅助Beta混合模型；为满足依赖限制，使用logit(y)的MixedLM（随机截距+随机斜率）作为实现近似，得到固定效应β、随机效应协方差Σ_b和精度φ；报告方差组分、ICC，并比较带/不带随机斜率。
7. 步骤7 covariate_form_comparison: 对四种协变量形态分别计算AIC/BIC、孕妇层五折CV RMSE、VIF；以1SE规则确定主报告形态；默认体重+身高+年龄+IVF但不是强制，需据数据选择；残差化身高体重仅作敏感性且不传播残差不确定性。
8. 步骤8 marginal_probability_integration: 在t∈{10.0,10.1,...,25.0}网格上，对每个个体x_i执行：从N(0,Σ_b)抽取M≤2000个随机效应b；计算μ(t,x_i,b)=expit(β0+β_ga·(t−mean_ga)+协变量项+b0+b1·(t−mean_ga))；抽取z~Beta(μ,φ)，加测量误差ε~N(0,σ²)得到z_obs=expit(logit(z)+ε)；π_marg(t,x_i)=mean(z_obs≥0.04)。
9. 步骤9 monotonicity_assertion_and_fallback: 对每个个体计算网格差分Δπ_marg(t_k)−π_marg(t_{k−1})；若存在负差分则先对孕周随机斜率截断使孕周效应为正；若仍违反，回退为固定效应孕周线性且不引入随机斜率；输出单调性断言统计。
10. 步骤10 individual_t_p_and_equivalence: 对每个个体计算t_p,i=inf{t∈[10.0,25.0]:π_marg(t,x_i)≥0.80}；同时在0.1周网格上显式求解argmin r_delay(t) subject to π_marg(t,x_i)≥0.80；比较两者差异；若窗口内无解则t_p,i=25.0、c_i=1；报告n_unsolved和r_cens。
11. 步骤11 bmi_grouping_1se: 对K∈{2,3,4}按样本BMI分位数或候选集生成边界；对每个分组计算组级π_g(t)和t_g；计算软损失参考L_{ρ,γ}(t,x_i)并聚合形成Q(c_K)；对最佳损失做bootstrap标准误，用1SE规则选择K和边界；输出边界bootstrap重现频率。
12. 步骤12 group_recommendation: 对最终K组计算组达标比例π_g(t)=n_g^{-1}Σπ_marg(t,x_i)，t_g=inf{t∈[10.0,25.0]:π_g(t)≥0.80}；报告组内个体t_p,i中位数、t_g与中位数一致性。
13. 步骤13 channel_A_error_sensitivity: 估计σ_tech=同一孕妇同一抽血次数技术重复logit(y)组内标准差；对σ∈{0,0.5σ_tech,σ_tech,2σ_tech}重算π_marg_σ、t_g,σ和κ_σ错分率；计算δt=t_g,σ−t_g,0。
14. 步骤14 channel_B_bootstrap: 对男胎孕妇做cluster bootstrap B=100次，每次重抽样孕妇单元并重拟合步骤5-6简化版，重算t_p,i、t_g；输出t_g的2.5%–97.5%分位区间和组间差异区间。
15. 步骤15 soft_loss_reference: 对每个个体和组中位数，在t∈[10.0,25.0]上对ρ∈{0.5,1,2}×γ∈{0,0.5,1,2}网格计算L_{ρ,γ}=r_delay(t)+ρ(1−π_marg(t,x_i))+γ·κ_σ(t,x_i)，输出t*作为二级参考列，不与主交付混淆。
16. 步骤16 export_outputs: 写出主结果表results/q3_main.csv；个体时点表results/q3_individual_tp.csv；组达标曲线results/q3_group_prob_curves.csv；误差敏感性results/q3_error_sensitivity.csv；边界热图results/q3_boundary_bootstrap.csv；模型系数与方差组分results/q3_model_coef.csv；所有数值列用np.clip防越界、dropna/fillna防NaN；末尾显式调用df.to_csv(MODELING_OUTPUT_DIR/'results'/'output.csv',index=False)。

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
- 附件为单个Excel表：每行是一次NIPT检测/样本记录，共1082行、31列；列混合孕妇基本信息、测序质量指标、染色体浓度/Z值及异常标签。建模时应按分析问题切分列，不能把全部行当作独立样本。
- 表内是重复测量/纵向结构：同一孕妇代码出现多条记录（孕妇代码唯一数约267，行数约1082），可能存在多次采血或一次采血多次检测。问题2、3的最佳时点/最早达标时间需按孕妇代码分组寻找，不能混行统计。
- 时间列需要预处理：检测孕周为文本格式（如11w+6），必须解析为数值孕周；末次月经和检测日期可用于孕周交叉校验，末次月经存在少量缺失。
- 问题1-3核心变量是Y染色体浓度与检测孕周、孕妇BMI；Y浓度达标阈值为4%。需先筛选男胎样本，并根据Y染色体相关列/性染色体信息区分男女胎。孕妇BMI可直接使用，身高、体重、年龄、怀孕/生产次数、IVF妊娠等可作为协变量。
- 同孕妇多条记录会引入组内相关性；数据概要显示多数数值列非正态且部分列存在较强时序自相关，建模时应考虑按孕妇分组的混合效应、适当变换或稳健方法，不能默认独立正态假设。
- 问题4中目标标签为AB列“染色体的非整倍体”，该列缺失率约88%，且取值为组合文本如T18、T13T18、T21；空白可能表示无异常或未检出，需结合胎儿是否健康AE列谨慎构造标签，不能简单删除或插补。
- 问题4的候选特征应聚焦21/18/13号染色体Z值、X染色体Z值与浓度、总读段数、唯一比对读段数、比对比例、重复比例、过滤比例、GC含量及各染色体GC含量，这些与测序质量和染色体异常信号相关；对问题2、3分析检测误差时，可借助同孕妇多次检测的波动以及测序质量指标度量误差。
- 样本序号、孕妇代码、末次月经、检测日期等主要用于索引、分组和时间计算，不适合作为普通回归特征；胎儿是否健康AE是最终结果，若用于辅助标注需避免与训练过程发生标签泄漏。
**已完成小题的结果（供本题复用）**：
- 小题 1（passed）：C:\Users\lingi\Desktop\Research\test1\outputs\results\q1.csv, C:\Users\lingi\Desktop\Research\test1\outputs\results\output.csv
- 小题 2（passed）：C:\Users\lingi\Desktop\Research\test1\outputs\results\q2.csv

## 5. 预期图表
- fig_roadmap [roadmap] 总体展示问题3从数据到推荐的建模流程与验证链路（数据来源：）
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
- fig_q2_ga_bmi_prob_curves [line] Q2 联合优化结果图（数据来源：）
- fig_q2_bmi_bins_tstar [scatter] Q2 联合优化结果图（数据来源：）
- fig_q2_loss_curves_optimal [line] Q2 联合优化结果图（数据来源：）
- fig_q2_error_shift_sigma [scatter] Q2 联合优化结果图（数据来源：）
- fig_q2_bootstrap_boundary_heatmap [scatter] Q2 联合优化结果图（数据来源：）
- fig_q2_calibration [scatter] Q2 联合优化结果图（数据来源：）
- fig_q2_rho_sensitivity [line] Q2 联合优化结果图（数据来源：）
- fig_q2_fnr_fpr_sigma [scatter] Q2 联合优化结果图（数据来源：）
- fig_q2_monotone_diagnostic [scatter] Q2 联合优化结果图（数据来源：）
- fig_q2_joint_k_selection [scatter] Q2 联合优化结果图（数据来源：）
- fig_q3_data_profile [scatter] 展示数据纵向结构、BMI分布及达标事件随孕周的分布，支持决策窗口与双通道建模（数据来源：results/q3_data_profile.csv）
- fig_q3_covariate_selection [bar] 支持主报告协变量形态选择，回应身高、体重、年龄同时纳入而非仅用BMI的原因（数据来源：results/q3_covariate_forms.csv）
- fig_q3_prob_curves [line] 展示组级达标比例随孕周单调上升及组间差异，为组最佳时点首穿提供直观依据（数据来源：results/q3_group_prob_curves.csv）
- fig_q3_individual_tp_hist [histogram] 报告个体推荐时点分布形态、删失聚集于25周的比例，回答“具体时点”而非方向性结论（数据来源：results/q3_individual_tp.csv）
- fig_q3_group_t_tradeoff [line] 展示每组最佳时点的不确定性和组间差异，验证组级约束风险最小化结果（数据来源：results/q3_main.csv）
- fig_q3_equivalence_check [scatter] 验证单调性条件下主交付可等价为达标保证首穿时点（数据来源：results/q3_equivalence.csv）
- fig_q3_monotone_diagnostic [line] 检验关键假设5（π_marg单调不减）是否满足，若不满足展示回退策略影响（数据来源：results/q3_monotone_diagnostic.csv）
- fig_q3_error_shift_sigma [line] 分析检测误差对分组时点的影响，并量化错分率随误差扩大（数据来源：results/q3_error_sensitivity.csv）
- fig_q3_bmi_boundary_bootstrap [heatmap] 检验数据驱动BMI分组方案的边界稳定性，避免把样本内最优边界表述为精确医学阈值（数据来源：results/q3_boundary_bootstrap.csv）
- fig_q3_p_sensitivity [line] 检验决策保证水平对主推荐时点的影响（数据来源：results/q3_p_sensitivity.csv）
- fig_q3_risk_sensitivity [line] 检验关键假设16（线性延迟风险）对结论的稳健性（数据来源：results/q3_risk_sensitivity.csv）
- fig_q3_model_structure_sensitivity [scatter] 检验关键假设7（二元与连续双通道）是否必要，回应单通道可能遗漏信号的风险（数据来源：results/q3_model_structure_sensitivity.csv）
- fig_q3_model_coef_dualchannel [bar] 支撑假设21：二元达标事件与连续Y浓度可能受不同因素支配，报告双通道系数差异（数据来源：results/q3_model_coef.csv）

## 6. 预期表格
- table_model_comparison：模型比较结果（列：待定；）
- table_smooth_terms：平滑项显著性检验（列：待定；）
- table_random_effects：随机效应方差组分与ICC（列：待定；）
- table_covariate_forms：BMI/体重/身高+体重三种协变量形态比较（列：待定；）
- table_sens_gc：GC处理策略敏感性（列：待定；）
- table_sens_ga_window：孕周外推界限敏感性（列：待定；）
- table_sens_marginal：边缘与条件达标概率差异（列：待定；）
- table_quantile_check：分位数回归与Beta-GAMM达标概率一致性（列：待定；）
- table_q3_main_result：问题3主结果：BMI分组与最佳NIPT时点（列：group, bmi_low, bmi_high, n, median_bmi, t_g, ci_low, ci_high, pi_g_at_tg, median_uncensored, n_unsolved, r_cens, t_star, delta_t_sigma_tech；给出每组BMI区间、最佳NIPT时点、不确定性及删失信息）
- table_q3_covariate_forms：四种协变量形态的模型比较（列：form, AIC, BIC, CV_RMSE, VIF_max, chosen_by_1se；支持主报告协变量形态选择）
- table_q3_model_coef：主二项与辅助Beta模型固定效应及随机效应组分（列：model, term, estimate, ci_low, ci_high, p_value, sigma2_intercept, sigma2_slope, corr_intercept_slope, phi；报告双通道模型系数显著性和随机效应结构）
- table_q3_error_sensitivity：通道A测量误差敏感性（列：group, sigma, t_g_sigma, delta_t_sigma, kappa_sigma, pi_g_at_tg；分析σ∈{0,0.5σ_tech,σ_tech,2σ_tech}对组时点和错分率的影响）
- table_q3_sensitivity_summary：关键假设扰动与对照实验汇总（列：experiment, control_group, perturbation, metric, main_result, sensitivity_result, conclusion；逐条覆盖关键假设的扰动或对照实验，形成‘实验→结论’成对证据）

## 7. 结果契约
```json
{
  "description": "问题3主结果：BMI分组及每组最佳NIPT时点与误差偏移；每组一行，分组数K∈{2,3,4}",
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
      "description": "BMI分组标签，如G1/G2/G3"
    },
    {
      "name": "bmi_low",
      "dtype": "float",
      "min": 15.0,
      "max": 60.0,
      "distinct_required": false,
      "description": "组BMI下界（kg/m²）"
    },
    {
      "name": "bmi_high",
      "dtype": "float",
      "min": 15.0,
      "max": 60.0,
      "distinct_required": false,
      "description": "组BMI上界（kg/m²）"
    },
    {
      "name": "n",
      "dtype": "int",
      "min": 1.0,
      "max": 300.0,
      "distinct_required": false,
      "description": "组内男胎孕妇数"
    },
    {
      "name": "median_bmi",
      "dtype": "float",
      "min": 15.0,
      "max": 60.0,
      "distinct_required": false,
      "description": "组内BMI个体代表值中位数"
    },
    {
      "name": "t_g",
      "dtype": "float",
      "min": 10.0,
      "max": 25.0,
      "distinct_required": true,
      "description": "组最佳NIPT时点（周），必须不同组有差异或给出统计证据支持统一时点"
    },
    {
      "name": "ci_low",
      "dtype": "float",
      "min": 10.0,
      "max": 25.0,
      "distinct_required": false,
      "description": "t_g的bootstrap 2.5%分位"
    },
    {
      "name": "ci_high",
      "dtype": "float",
      "min": 10.0,
      "max": 25.0,
      "distinct_required": false,
      "description": "t_g的bootstrap 97.5%分位"
    },
    {
      "name": "pi_g_at_tg",
      "dtype": "float",
      "min": 0.0,
      "max": 1.0,
      "distinct_required": false,
      "description": "t_g处的组级达标比例，clip到[0,1]"
    },
    {
      "name": "median_uncensored",
      "dtype": "float",
      "min": 10.0,
      "max": 25.0,
      "distinct_required": false,
      "description": "未删失个体t_p中位数"
    },
    {
      "name": "n_unsolved",
      "dtype": "int",
      "min": 0.0,
      "max": 300.0,
      "distinct_required": false,
      "description": "窗口内未达80%保证的人数"
    },
    {
      "name": "r_cens",
      "dtype": "float",
      "min": 0.0,
      "max": 1.0,
      "distinct_required": false,
      "description": "删失比例"
    },
    {
      "name": "t_star",
      "dtype": "float",
      "min": 10.0,
      "max": 25.0,
      "distinct_required": false,
      "description": "二级软损失参考时点"
    },
    {
      "name": "delta_t_sigma_tech",
      "dtype": "float",
      "min": -5.0,
      "max": 5.0,
      "distinct_required": false,
      "description": "σ=σ_tech相对σ=0的t_g偏移（周）"
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
LTM 中出现 '95.0 %'，但题目常量为 [4.0] %（可能 LLM 记错了数值）
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

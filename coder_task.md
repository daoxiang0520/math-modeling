# 建模方案与实现架构说明书

## 0. 当前小题
**当前小题（1/4）**：2025 年高教社杯全国大学生数学建模竞赛题目 
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
风险高；晚期发现（28 周以后）风险极高。 
实践表明，男胎Y 染色体浓度与孕妇孕周数及其身体质量指数（BMI）紧密相关。通常根据孕妇的
BMI 值进行分组（例如：[20,28)，[28,32)，[32,36)，[36,40)，40 以上）分别确定NIPT 的时点（相对孕
期的时间点）。由于每个孕妇的年龄、BMI、孕情等存在个体差异，对所有孕妇采用简单的经验分组和
统一的检测时点进行NIPT，会对其准确性产生较大影响。因此，依据BMI 对孕妇进行合理分组，确定
各不同群组的最佳NIPT 时点，可以减少某些孕妇因胎儿不健康而缩短治疗窗口期所带来的潜在风险。 
为了研究各类孕妇群体合适的NIPT 时点，并对检测的准确性进行分析，附件给出了某地区（大多
为高BMI）孕妇的NIPT 数据。在实际检测中，经常会出现测序失败（比如：检测时点过早和不确定因
素影响等）的情况。同时为了增加检测结果的可靠性，对某些孕妇有多次采血多次检测或一次采血多次
检测的情况。试利用附件提供的数据建立数学模型研究如下问题：
问题1  试分析胎儿Y 染色体浓度与孕妇的孕周数和BMI 等指标的相关特性，给出相应的关系模
型，并检验其显著性。

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
**建模目标**：对男胎样本建立Y染色体浓度与孕周数和BMI等指标之间具有统计显著性的关系模型，以Beta-GAMM为主模型、分位数回归为辅助验证，并输出达标概率曲线、达标孕周反演方法和测量误差估计。

## 2. 建模设定（动态 LTM，编程手不得修改）
**假设**：
- 假设1：分析对象限定为男胎样本，依据孕妇代码、Y染色体浓度（V列）与Y染色体Z值（U列）具有有效数值交叉核查男胎，不能仅凭Y列缺失率判定女胎。依据：数据概要显示Y染色体浓度列缺失率为0，与“女胎Y染色体相关列为空白”存在矛盾；风险：误纳女胎会引入Y≈0噪声；可验证：交叉表检查Y浓度/Y-Z值与孕妇代码、X染色体浓度及胎儿健康列的一致性。
- 假设2：数据记录结构为非独立的纵向重复测量，每行是一次检测记录，同一孕妇可能多次检测；同一孕妇的多条记录必须按组内相关处理，不得把每行当作独立个体。依据：数据智能摘要显示1082条记录/267位孕妇；风险：行级独立假设会低估标准误；可验证：随机效应模型与行级独立模型的参数显著性差异及孕妇层分组交叉验证误差。
- 假设3：孕周变量必须解析为数值周数：ga=w+d/7，兼容纯周格式；用末次月经与检测日期交叉核对，容差±1周写入代码。依据：数据概要显示检测孕周为文本格式（如11w+6，且存在26w、11w等纯周写法）；风险：孕周解析错误会直接污染GA效应；可验证：核对后的异常记录占比及剔除前后模型比较。
- 假设4【关键假设】：Y染色体浓度y在(0,1)内且有界、非正态（Shapiro-Wilk p<0.05），主模型采用Beta分布；实测y范围0.010–0.234，无0/1、无缺失、无零膨胀，因此边界修正仅作注记、不进入主流程。依据：数据概要及正态性检验；风险：Beta分布尾部若失配会影响达标概率估计；可验证：DHARMa残差诊断及Beta模型与高斯/logit-GAM基准比较。
- 假设5【关键假设】：Y染色体浓度的条件均值与孕周、BMI存在非线性和交互关系，设定logit尺度包含s1(ga)、s2(bmi)、ti(ga,bmi)、s3(age)，并用GA分段线性基准检验非线性必要性。依据：文献qiao 2019显示GA分段斜率差异；风险：忽略交互或非线性会扭曲达标孕周反演；可验证：平滑项H0:s_j=0的edf、p值、δAIC及与分段线性基准的比较。
- 假设6【关键假设】：纵向重复测量下随机效应按设计进入，设置随机截距b0i和孕周随机斜率b1i（ga中心化），不依赖ICC显著性检验；报告方差组分与ICC。依据：同一孕妇多次检测的纵向结构；风险：缺失随机斜率会低估个体间达标时间差异；可验证：带/不带随机斜率的AIC/BIC、方差组分及95%置信区间比较。
- 假设7：协变量选择为年龄与IVF指示必入；BMI、体重、身高+体重三种形态显式比较；孕产次候选并在VIF/concurvity检查后决定。依据：qiao 2023显示IVF与胎儿分数显著负相关，hou 2015提示体重可能优于BMI；风险：共线性导致系数不稳定；可验证：VIF、concurvity及三种形态AIC/BIC。
- 假设8【关键假设】：GC含量不采用40%–60%硬阈值剔除（原文：正常GC含量范围为40% ~ 60%）；因实测GC全距0.386–0.421，42%行低于0.4、无高于0.6，判断为检测平台系统偏差；将GC作为连续质量协变量或质量权重并做敏感性分析。依据：原文事实与实测数据；风险：硬剔除会损失451/1082行并引入选择偏差；可验证：GC纳入或加权前后主效应及达标概率变化。
- 假设9：孕周外推界限为10.0–25.0周（原文：通常孕妇的孕期在10 周~25 周之间可以检测胎儿性染色体浓度）；保留超窗样本（最大27.7周）但仅作描述，不在10–25周外做推断。风险：外推可能高估晚期浓度；可验证：排除大于25周样本后的敏感性分析。
- 假设10：技术重复与纵向重复分别处理：同一（孕妇,抽血次数）多行共40组101行（9.3%）用于估计测量误差σ_tech；同一孕妇不同孕周的多次检测作为纵向重复。依据：数据智能摘要；可验证：技术重复组内方差与纵向方差比较。
- 假设11：胎儿不健康记录38条不默认剔除，进行剔除vs保留敏感性分析。依据：核型-胎儿分数关系在文献中有争议（rava 2014显示差异、deng 2023大样本相似）；风险：不当剔除会改变主效应或达标概率；可验证：两种样本集模型效应差异。
- 假设12：序号、年龄的lag-1自相关来自数据按孕妇或检测时间排序，不按时间序列建模；以孕妇随机效应吸收组内相关。依据：数据发现lag-1自相关大于0.5；可验证：按孕妇分组检查残差自相关。
- 假设13【关键假设】：分位数回归作为辅助验证模型，在logit(y)尺度建模，τ∈{0.05,0.10,...,0.95}，施加单调性约束，不放ga×bmi交互，用于与Beta-GAMM后验达标概率交叉验证。依据：Grantz 2023、Rava 2014；风险：尾部分位数稀疏时过拟合；可验证：分位数预测与Beta-GAMM后验达标概率一致性、尾部覆盖。
- 假设14【关键假设】：对新孕妇的达标概率必须采用边缘预测（对随机效应积分），拟合曲线仅作为随机效应=0的条件轨迹；阈值y_thr=4.0%（原文：如果男胎的Y染色体浓度达到或高于4%）与分位数概率τ严格区分。依据：随机效应模型；风险：条件曲线低估个体差异导致达标概率乐观；可验证：边缘与条件达标概率曲线差异。
- 假设15：测序质量相关列（读段数、GC、比对比例等）不作为主模型默认预测变量，仅在质量权重或敏感性分析中使用，以保持主链可解释性。依据：问题1要求刻画Y浓度与孕周、BMI等指标关系；可验证：逐步纳入质量变量后的模型比较与显著性变化。
**符号表**：
- N_rec: 男胎检测记录总数（原始总记录1082行，男胎筛选后待定）
- n_w: 男胎孕妇数量（原始总孕妇267位，男胎筛选后待定）
- n_i: 第i位孕妇的检测次数
- i: 孕妇索引，i=1,...,n_w
- j: 孕妇内检测索引，j=1,...,n_i
- y_ij: 第i位孕妇第j次检测的Y染色体浓度，无量纲比例，0<y<1
- ga_ij: 数值化孕周，单位：周，ga=w+d/7
- w_ij: 检测孕周文本中的整数周部分
- d_ij: 检测孕周文本中的天数部分，0≤d<7
- bmi_ij: 孕妇BMI，单位：kg/m²
- age_i: 孕妇年龄，单位：岁
- height_i: 孕妇身高，单位：cm
- weight_i: 孕妇体重，单位：kg
- ivf_i: IVF妊娠指示，0/1，1表示IVF妊娠
- parity_ac_i: 孕妇怀孕次数，候选协变量
- parity_ad_i: 孕妇生产次数，候选协变量
- gc_ij: GC含量比例，小数；题目原文正常范围为40.0%–60.0%
- y_thr: 临床达标阈值，y_thr=0.04，即4.0%
- μ_ij: Beta分布条件均值，0<μ<1
- φ: Beta分布精度参数，φ>0
- η_ij: logit均值线性预测器，η_ij=logit(μ_ij)
- β0: 固定截距
- β_ivf: IVF固定效应系数
- s1(ga): 孕周平滑函数
- s2(bmi): BMI平滑函数
- s3(age): 年龄平滑函数
- ti(ga,bmi): 孕周与BMI张量积交互平滑项
- b0i: 孕妇随机截距
- b1i: 孕妇随机孕周斜率，对应中心化孕周ga_ij−ḡ
- b_i: 随机效应向量，b_i=(b0i,b1i)'
- Σ_b: 随机效应协方差矩阵
- σ_b0²: 随机截距方差
- σ_b1²: 随机斜率方差
- σ_b01: 随机截距与斜率协方差
- ḡ: 男胎样本孕周均值，用于随机斜率中心化，单位：周
- τ: 分位数概率水平，0<τ<1
- q_logit_y(τ;ga,bmi): logit(y)的条件τ分位数函数
- β0τ: 分位数回归截距项
- g1τ(ga): 分位数回归中孕周平滑项
- g2τ(bmi): 分位数回归中BMI平滑项
- τ_star(ga,bmi): 阈值反演得到的分位数水平，满足q_logit_y(τ_star)=logit(y_thr)
- P_ok(ga,bmi): 达标概率，P(y≥y_thr|ga,bmi)
- P_marg(ga,bmi): 边缘达标概率，对随机效应积分后的预测概率
- f(b;0,Σ_b): 二元正态概率密度，随机效应分布密度
- σ_tech: 技术重复测量误差标准差，logit(y)尺度
- k: 平滑项基维数
- edf: 平滑项估计自由度
- AIC: 赤池信息准则
- BIC: 贝叶斯信息准则
- ΔAIC: 模型与最优模型的AIC差值
- RMSE: 均方根误差
- MAE: 平均绝对误差
- MS_within: 技术重复组内均方，用于估计σ_tech
**公式/方程**：
- 孕周数值解析：ga_ij = w_ij + d_ij/7，其中w_ij为整数周、d_ij为天数；纯周格式d_ij=0。
- logit变换：logit(u)=ln(u/(1-u))，定义域u∈(0,1)。
- Beta-GAMM主模型：y_ij ~ Beta(μ_ij,φ)，η_ij = logit(μ_ij) = β0 + s1(ga_ij) + s2(bmi_ij) + ti(ga_ij,bmi_ij) + s3(age_i) + β_ivf·ivf_i + b0i + b1i·(ga_ij−ḡ)。
- 随机效应分布：b_i=(b0i,b1i)' ~ N(0,Σ_b)，Σ_b=[[σ_b0²,σ_b01],[σ_b01,σ_b1²]]。
- Beta条件方差：Var(y_ij|b_i)=μ_ij(1−μ_ij)/(1+φ)。
- 分位数回归辅助模型：q_logit_y(τ;ga_ij,bmi_ij)=β0τ+g1τ(ga_ij)+g2τ(bmi_ij)，τ∈{0.05,0.10,...,0.95}，对g1τ和g2τ施加单调性约束。
- 达标概率反演：设y_thr=0.04，P_ok(ga,bmi)=1−τ_star(ga,bmi)，其中τ_star满足q_logit_y(τ_star;ga,bmi)=logit(y_thr)=ln(0.04/0.96)。
- 边缘达标概率：P_marg(ga,bmi)=∫ P_ok(ga,bmi|b_i) f(b_i;0,Σ_b) db_i。
- 技术重复测量误差：σ_tech=√MS_within，MS_within为同一(孕妇,抽血次数)技术重复组内logit(y)的REML/单因素方差组内均方。
- 基准对照模型：线性回归y_ij=β0+β_ga·ga_ij+β_bmi·bmi_ij+e_ij；logit(y)高斯GAM；GA分段线性基准（10–12.5/12.5–20/>20周，qiao 2019文献先验）。
- 模型比较：AIC=2k−2lnL，BIC=k·ln(N_rec)−2lnL，δAIC=AIC_model−AIC_min；预测误差采用按孕妇分组的交叉验证RMSE和MAE。
**解题思路**：步骤1：读取附件.xlsx；Coder必须按照parse_hints解析字符串列：孕妇代码用df['孕妇代码'].str.replace('A','',regex=False).astype(float)，末次月经用pd.to_datetime，染色体的非整倍体用df['染色体的非整倍体'].str.replace('T','',regex=False).astype(float)，检测孕周需自行解析w/d格式。步骤2：筛选男胎样本，依据孕妇代码、Y染色体浓度和Y染色体Z值有效数值交叉核查；识别同一孕妇的多行记录，区分纵向重复（同一孕妇不同孕周）与技术重复（同一孕妇同一抽血次数多行）。步骤3：解析孕周ga=w+d/7，兼容纯周格式；用末次月经与检测日期交叉核对，容差±1周。步骤4：构建Beta-GAMM主模型：logit(μ)=β0+s1(ga)+s2(bmi)+ti(ga,bmi)+s3(age)+β_ivf·ivf+b0i+b1i(ga−ḡ)，采用REML选平滑参数，并对基维数k=5/8/10做敏感性检查。步骤5：协变量形态比较BMI、体重、身高+体重三种形式；年龄与IVF必入，孕产次候选，VIF/concurvity检查。步骤6：基准对照：普通线性回归、logit(y)高斯GAM、GA分段线性基准；用AIC/BIC/δAIC和按孕妇分组的交叉验证RMSE/MAE比较，防止行级数据泄漏。步骤7：显著性检验：各平滑项H0:s_j=0，报告edf、p值、δAIC；随机效应报告方差组分、95%置信区间和ICC；诊断包括DHARMa类残差、concurvity、k敏感性，必要时做孕妇层cluster bootstrap。步骤8：分位数回归辅助：在logit(y)尺度拟合单调约束分位数模型，τ∈{0.05,...,0.95}，不放ga×bmi交互；反演达标概率P(y≥0.04|ga,bmi)=1−τ_star。步骤9：输出达标概率曲线、边缘预测P_marg与条件曲线；利用技术重复组估计σ_tech。步骤10：敏感性分析：孕周外推限10.0–25.0周，超窗样本保留但仅描述；GC作为连续质量协变量/权重；胎儿不健康记录做剔除vs保留。衔接说明：本问题只输出达标概率P_ok、达标孕周反演和测量误差σ_tech，不做BMI分组、最佳NIPT时点优化；后续问题2/3将结合题目风险分级（12周内低风险、13–27周高风险、28周以后极高风险）建立决策模型，问题4女胎异常判定不在本问题范围。

## 3. 算法与求解
**算法摘要**：
**伪代码/实现步骤**：
1. 步骤1: load_and_parse_data(filepath): 读取附件.xlsx；孕妇代码列去除非数字字符后astype(float)；末次月经用pd.to_datetime；染色体的非整倍体去除非数字后astype(float)；检测孕周解析为 w=提取整数周，d=提取天数(0<=d<7)，ga=w+d/7；纯周格式d=0；异常格式记入error_log。
2. 步骤2: filter_male_and_replicates(data): 基于孕妇代码、Y染色体浓度列、Y染色体Z值列均有有效数值交叉核查男胎；剔除Y浓度<=0或>=1的异常行（实际应在(0,1)内）；标记技术重复：同一(孕妇代码,抽血次数)多行，技术重复组用于估计σ_tech；同一孕妇不同孕周作为纵向重复。
3. 步骤3: preprocess_for_model(data): 提取核心变量y=Y染色体浓度, ga, bmi, age, ivf, parity_ac, parity_ad, gc；必要时fillna/median或dropna；对y进行clip(y, 1e-6, 1-1e-6)；计算logit_y=log(y/(1-y))；对ga做中心化ga_c=ga-mean(ga)用于随机斜率；构建孕妇分组ID。
4. 步骤4: build_design_matrices(data, k=6): 使用自然三次样条基（或patsy ns函数）构建s1(ga), s2(bmi), s3(age)；用ga和bmi的外积构造ti(ga,bmi)张量积基（低维，如5x5）；设计固定效应矩阵X=[1, s1, s2, ti, s3, ivf]；随机效应矩阵Z=[1, ga_c]；记录样条节点位置。
5. 步骤5: fit_comparison_models(X, Z, logit_y, groups): 拟合a) 主模型：statsmodels.MixedLM(endog=logit_y, exog=X, groups=groups, exog_re=Z)采用REML；若收敛失败或奇异，降级为仅随机截距(statsmodels.MixedLM移除Z第二列)或statsmodels.GEE(exchangeable correlation)；b) 线性回归：statsmodels.OLS(logit_y, [1, ga, bmi])；c) logit高斯GAM：同样用样条但无随机效应；d) GA分段线性基准：设计矩阵含ga分段(10-12.5,12.5-20,>20)和bmi线性。
6. 步骤6: compute_model_metrics(models, data): 对每个模型计算AIC、BIC、log-likelihood、edf(平滑项估计自由度)；按孕妇分组的5折交叉验证计算RMSE和MAE（同一孕妇的所有记录必须在同一折）；记录δAIC=min(AIC)-AIC_model。
7. 步骤7: test_smooth_terms(main_model, reduced_models): 对s1(ga), s2(bmi), ti(ga,bmi), s3(age)分别拟合去掉该项的模型，进行似然比检验或查看anova表；报告各平滑项edf、p值、δAIC；检查concurvity（用模型矩阵相关性近似）。
8. 步骤8: fit_quantile_aux(logit_y, X_qr): 使用statsmodels.QuantReg(logit_y, X_qr)拟合τ∈{0.05,0.10,...,0.95}；X_qr=[1, s1(ga), s2(bmi)]；对g1τ(ga)和g2τ(bmi)施加单调性（通过约束节点系数递增或在数值优化时加罚）；保存分位数预测曲线。
9. 步骤9: invert_threshold_and_marginal(models, ga_grid, bmi_grid): 对阈值y_thr=0.04，logit阈值=ln(0.04/0.96)；从分位数辅助找到满足q_logit_y(τ_star)=logit阈值的τ_star，P_ok=1-τ_star；主模型边缘预测：固定效应预测+随机效应=0条件曲线；再以估计的Σ_b为协方差，Monte Carlo抽样200次随机效应，对每次预测得条件P_ok，平均得边缘P_marg；将P_ok和P_marg输出到网格。
10. 步骤10: estimate_tech_error(data): 对同一(孕妇,抽血次数)的技术重复组，计算组内logit_y的均方MS_within；σ_tech=sqrt(MS_within)；与纵向方差组分σ_b0^2比较，输出到汇总表。
11. 步骤11: sensitivity_analysis(models, data): 执行四类敏感性：a) Beta分布假设：比较主模型与直接GLM(Binomial)或分位数预测的残差分布和达标概率差异；b) 非线性交互：去掉ti(ga,bmi)后模型AIC和达标概率曲线变化；c) GC处理：分别使用GC作为连续线性项、GC质量权重(1/sqrt(GC))、剔除低GC(无硬阈值，仅用分位数阈值如1%或5%极端值)拟合主模型，比较主效应；d) 孕周外推：排除ga>25周的样本重拟合，与全样本模型比较。
12. 步骤12: export_outputs(figs, tables, result_df): 将模型比较、平滑项、随机效应、分位数一致性、敏感性结果按长格式合并为result_df；所有数值进行np.clip或pd.notnull检查，p值clip到[0,1]，概率clip到[0,1]；保存result_df.to_csv(MODELING_OUTPUT_DIR/results/output.csv, index=False)；同时保存图表至指定目录。

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
- 附件.xlsx是单一样本明细表，每行表示一次NIPT检测/采样记录而非一名孕妇；孕妇代码存在重复，说明同一孕妇可能在不同孕周多次检测，建模时应按孕妇代码识别纵向重复测量结构，不能把每行当作独立个体。
- 问题1-3核心变量为Y染色体浓度、检测孕周、孕妇BMI，以及年龄/身高/体重等；其中检测孕周是文本格式（如11w+6），必须解析为数值周数后才能回归或建模。
- 同一孕妇多行、多次采血检测可用于估计Y浓度随孕周变化和最早达标孕周；应优先在孕妇层级聚合或采用混合效应模型，避免把同一孕妇的多条记录视为独立样本。
- 男胎与女胎需在代码中先区分：题目说明女胎Y染色体相关列为空白，但数据概要显示Y染色体浓度缺失率为0，存在矛盾；应检查Y染色体浓度/Y-Z值是否为空、零或非数值，不应仅依据缺失率判断女胎样本。
- AB列（13/18/21染色体非整倍体）缺失率超过50%，但题目定义空白即为无异常，因此该列不应按普通缺失处理，应转换为阴性/阳性标签；该列是问题4的判定标签，且异常类别可能高度不平衡。
- 问题4可用的特征包括X染色体及13/18/21染色体的Z值、X染色体浓度、各染色体GC含量、总读段数、比对比例、重复比例、唯一比对读段数、过滤读段比例、BMI等；读段数和比例之间可能存在函数或强相关关系，建模时需注意共线性并做特征选择。
- L/M/N/O/P/AA及染色体GC含量等测序质量相关列在问题3中可作为检测误差或数据质量的代理变量，需在建模前检查是否要用质量阈值剔除低质量记录，或作为协变量纳入达标时间分析。
- 问题2和3的最佳NIPT时点需把Y染色体浓度≥4%的达标时间与风险等级（12周内低风险、13-27周中风险、28周后高风险）联合建模；从重复检测记录中提取每位男胎孕妇是否达标及最早达标孕周是分组和风险优化的前提。
- 数据概要中某些列显示强时序自相关可能是由于数据按孕妇代码或检测时间排列造成的，不能直接按时间序列处理；应按孕妇分组或考虑重复测量的组内相关性。
- BMI分组不应直接照搬题目示例区间，应基于男胎样本的实际BMI分布和达标时间-风险关系进行数据驱动分组，并注意数据中BMI范围较宽、以高BMI为主，可能影响分组边界。

## 5. 预期图表
- fig_roadmap [roadmap] 展示全文从数据到问题1再到后续问题的总体研究框架，重点突出问题1的主模型与辅助验证关系（数据来源：）
- fig_q1_scatter [scatter] 展示Y浓度与两个核心预测变量的原始关系，支持非线性与交互建模动机（数据来源：results/q1_scatter.csv）
- fig_q1_smooth_ga [line] 展示孕周对logit(Y浓度)的非线性效应及显著性（数据来源：results/q1_smooth_ga.csv）
- fig_q1_smooth_bmi_int [heatmap] 展示交互项空间分布，说明不同BMI水平下孕周效应的变化（数据来源：results/q1_ti_heatmap.csv）
- fig_q1_quantile_curves [line] 展示分位数辅助模型结果，支撑达标概率反演（数据来源：results/q1_quantile_curves.csv）
- fig_q1_prob_curves [line] 展示达标概率随孕周和BMI的变化，并比较边缘与条件预测的差异（数据来源：results/q1_prob_curves.csv）
- fig_diag_resid [scatter] 检验残差正态性和方差齐性，验证模型假设合理性（数据来源：results/q1_resid.csv）
- fig_sens_dist [line] 检验Y浓度分布假设对达标概率估计的影响（数据来源：results/sens_dist.csv）
- fig_sens_interaction [line] 验证交互项对达标概率反演的必要性（数据来源：results/sens_interaction.csv）
- fig_sens_gc [line] 检验GC作为连续协变量或权重与不使用GC的差异，避免硬剔除导致的偏差（数据来源：results/sens_gc.csv）
- fig_sens_marginal [line] 展示忽略个体随机效应导致的条件曲线乐观偏差（数据来源：results/sens_marginal.csv）
- fig_sens_ga_window [line] 检验超窗样本对外推结论的影响（数据来源：results/sens_ga_window.csv）
- fig_anchor_threshold [scatter] 将承重构造y_thr绑定到实际数据分布，说明阈值反演和达标概率的物理意义（数据来源：results/q1_threshold_anchor.csv）

## 6. 预期表格
- table_model_comparison：模型比较结果（列：模型, AIC, BIC, δAIC, RMSE, MAE, 备注；证明Beta/logit-GAMM主模型优于线性回归、logit-GAM和GA分段线性基准）
- table_smooth_terms：平滑项显著性检验（列：平滑项, edf, p值, δAIC(去掉该项), 显著性结论；报告各平滑项的显著性和对模型拟合的贡献）
- table_random_effects：随机效应方差组分与ICC（列：参数, 估计值, 95%置信区间；刻画个体异质性和纵向相关性）
- table_covariate_forms：BMI/体重/身高+体重三种协变量形态比较（列：协变量形态, AIC, BIC, δAIC, 主效应方向；支持假设7中BMI形态的选择）
- table_sens_gc：GC处理策略敏感性（列：策略, AIC, δAIC, s1(ga)最大值变化, 达标概率变化；支持假设8中GC非硬剔除的合理性）
- table_sens_ga_window：孕周外推界限敏感性（列：样本集, N_rec, AIC, δAIC, 10-25周外推P_ok差异；支持假设9中超窗样本保留但仅描述的决策）
- table_sens_marginal：边缘与条件达标概率差异（列：BMI, 孕周, P_marg, P_cond, ΔP；量化条件曲线乐观偏差，支持假设14边缘预测必要性）
- table_quantile_check：分位数回归与Beta-GAMM达标概率一致性（列：BMI, 孕周, 分位数反演P_ok, Beta-GAMM边缘P_ok, 绝对差值；展示辅助模型交叉验证结果）

## 7. 结果契约
```json
{
  "description": "问题1模型汇总与显著性检验结果长表，每行一个模型/指标组合，含模型比较、平滑项、随机效应、分位数一致性和敏感性统计量",
  "allow_single_row": false,
  "min_rows": 5,
  "max_rows": 200,
  "columns": [
    {
      "name": "model",
      "dtype": "category",
      "min": null,
      "max": null,
      "distinct_required": false,
      "description": "模型或分析项名称，如main_model, linear_baseline, s1_ga"
    },
    {
      "name": "metric",
      "dtype": "category",
      "min": null,
      "max": null,
      "distinct_required": false,
      "description": "统计量名称，如AIC, BIC, RMSE, MAE, edf, p_value, sigma_b0_sq, sigma_b1_sq, ICC, sigma_tech, P_ok_diff"
    },
    {
      "name": "value",
      "dtype": "float",
      "min": -10.0,
      "max": 10000.0,
      "distinct_required": false,
      "description": "统计量数值"
    },
    {
      "name": "se",
      "dtype": "float",
      "min": 0.0,
      "max": 100.0,
      "distinct_required": false,
      "description": "标准误或置信区间半宽，如无则填NaN"
    },
    {
      "name": "p_value",
      "dtype": "float",
      "min": 0.0,
      "max": 1.0,
      "distinct_required": false,
      "description": "p值，已clip到[0,1]"
    },
    {
      "name": "n",
      "dtype": "int",
      "min": 0.0,
      "max": 2000.0,
      "distinct_required": false,
      "description": "样本量"
    },
    {
      "name": "bmi",
      "dtype": "float",
      "min": 10.0,
      "max": 60.0,
      "distinct_required": false,
      "description": "对应BMI值（若该指标与BMI相关）"
    },
    {
      "name": "ga",
      "dtype": "float",
      "min": 5.0,
      "max": 30.0,
      "distinct_required": false,
      "description": "对应孕周值（若该指标与孕周相关）"
    },
    {
      "name": "note",
      "dtype": "category",
      "min": null,
      "max": null,
      "distinct_required": false,
      "description": "备注信息"
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

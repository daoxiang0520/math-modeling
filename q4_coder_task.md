# 建模方案与实现架构说明书

## 0. 当前小题
**当前小题（4/4）**：问题4  由于孕妇和女胎都不携带Y 染色体，重要的是如何判定女胎是否异常。试以女胎孕妇的21
号、18 号和13 号染色体非整倍体（AB 列）为判定结果，综合考虑X 染色体及上述染色体的Z 值、GC
含量、读段数及相关比例、BMI 等因素，给出女胎异常的判定方法。 
 
 
 
 
 
 
 

附录1  附件中各列数据的说明 
 
 
附录2  Z 值（Z-score） 
Z 值（Z-score）的计算公式： 
𝑍= 𝑋−𝜇
𝜎
 
其中X 为待检测样本中目标染色体的相对计数比例，𝜇 为正常对照群体中该染色体计数比例的均值，𝜎 
为正常群体中该比例的标准差。在NIPT 中，对于常见染色体非整倍体检测，通常采用Z 值分析方法进
行统计判定。已知染色体非整倍体通常定义为该染色体存在一个或三个拷贝，正常为两个拷贝，且每条
染色体所采集到的读段数量与该染色体长度成正比。 
列 
说明 
列 
说明 
A 
样本序号 
Q 
13 号染色体的Z 值 
B 
孕妇代码 
R 
18 号染色体的Z 值 
C 
孕妇年龄 
S 
21 号染色体的Z 值 
D 
孕妇身高 
T 
X 染色体的Z 值 
E 
孕妇体重 
U 
Y 染色体的Z 值（女胎数据此列为空白） 
F 
末次月经时间 
V 
Y 染色体浓度，即Y 染色体游离DNA 片
段的比例（女胎数据此列为空白） 
G 
IVF 妊娠方式 
W 
X 染色体浓度（其数值是通过生物信息学在
一定假设下通过数据分析估计得出，可能出
现负值） 
H 
检测时间 
X 
13 号染色体的GC 含量 
I 
检测抽血次数 
Y 
18 号染色体的GC 含量 
J 
孕妇本次检测时的孕周（周数+天数） 
Z 
21 号染色体的GC 含量 
K 
孕妇BMI 指标 
AA 
被过滤掉的读段数占总读段数的比例 
L 
原始测序数据的总读段数（个） 
AB 
检测出的13 号，18 号，21 号染色体非整
倍体，即数量异常，空白即为无异常 
M 
总读段数中在参考基因组上比对的比例 
AC 
孕妇的怀孕次数 
N 
总读段数中重复读段的比例 
AD 
孕妇的生产次数 
O 
总读段数中唯一比对的读段数（个） 
AE 
胎儿是否健康（婴儿出生后的结果） 
P 
GC 含量，序列中碱基 G（鸟嘌呤）和 C
（胞嘧啶）所占的比例，是测序数据质量
评估中的一个重要指标，正常 GC 含量范
围为40% ~ 60%，GC 含量过高、过低、
或分布异常可能意味着测序质量存在问题

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
**建模目标**：以女胎孕妇的21号、18号和13号染色体非整倍体AB列结果为判定目标，建立“X染色体浓度可靠性门—逐染色体z阈值—孕妇级保守合并与验证”的三层最小结构判定方法，输出记录级/孕妇级三分类（阳性/阴性/无法判定）及可验证的灵敏度、特异度、PPV、NPV、F1、AUC和覆盖率。

## 2. 建模设定（动态 LTM，编程手不得修改）
**假设**：
- 假设1：分析对象限定为女胎样本；女胎识别直接由工作表特征给出，即Y染色体Z值（U列）与Y染色体浓度（V列）在女胎记录中为空白，不另行设计识别规则。依据：人类架构审核数据画像表明女胎表约605条记录、147位孕妇，且Y相关列U/V全部空白；风险：若存在少量男胎被误标为女胎，会污染可靠性门与z阈值定标；可验证：交叉核查X染色体浓度、X染色体Z值、胎儿是否健康列及孕妇代码的分布一致性。
- 假设2：数据结构为孕妇层纵向重复测量，每行是一次检测记录，同一孕妇可能存在多次检测；女胎异常判定必须按孕妇分块处理，不能把每条记录当作独立样本。依据：人类架构审核数据画像显示女胎表605条记录、147位孕妇，每人1–9次检测（均值约4.1次）；风险：行级独立会导致信息泄漏和标准误低估；可验证：行级独立评估与孕妇分块交叉验证/聚类bootstrap的指标差异。
- 假设3：女胎异常判定标签以AB列“染色体的非整倍体”为判定结果；AB空白即无异常（阴性），AB文本解析为多标签集合，逐染色体构造标签y_ijc=1( c∈AB_ij )。AB列取值如T18、T13T18、T21、T13、T18T21、T13T21。依据：题目问题4明确以AB列为判定结果，附录列说明“空白即为无异常”；风险：空白可能包含未检出或缺失情形而非真阴性；可验证：孕妇内混合标签敏感性分析，比较“空白即阴性”与“剔除孕妇内存在混合标签记录”的指标变化。
- 假设4【关键假设】：女胎表中AE列“胎儿是否健康”全部为“健康”，不能作为金标准或主标签修正依据；主标签采用AB列语义，禁止使用AE列做交叉标签清洗。依据：人类架构审核数据画像显示女胎表605/605条AE均为“健康”，任何“确认阳性=AB异常且AE不健康”的标签方案在女胎表上得到0例；风险：若强制使用AE会抹掉全部阳性信号；可验证：输出AE列取值计数，并展示“AB异常且AE不健康”的确认为0例。
- 假设5【关键假设】：X染色体浓度w_ij不是胎儿浓度代理，而作为NIPT可靠性/平台信号。女胎表w范围约-0.070～0.121，中位数约-0.008，约67%为负值（数据估计值，非题目常量）；负值予以保留，不截断为0。依据：题目原文“女胎的X染色体浓度没有异常，则可认为NIPT的结果是基本准确的，否则难以保证结果准确性要求”以及数据画像显示w呈平台信号特征；风险：将w当作胎儿浓度会锚定错误的拷贝数机制，截断负值会损失平台偏差信息；可验证：保留负值与截断为0两种口径下可靠性门覆盖率-准确率对比，并检查w与孕周/BMI的相关性。
- 假设6【关键假设】：第一层可靠性门以X染色体浓度w_ij为输入，主判据采用w_ij<q10_f（女胎w的经验第10百分位，约-0.0272，该值为数据估计值）；备选界为q5_f、q25_f。门输出为可判定/无法判定（建议重测）。依据：题目原文可靠性条件及人类架构审核数据聚集性：w<q5时24%阳性 vs 3%阴性，w<q10时40% vs 6%，w<q25时61% vs 20%；风险：分位界依赖本数据集且在小样本分组中可能不稳定；可验证：绘制q5/q10/q25下覆盖率-准确率权衡曲线，并用分层bootstrap报告门界稳定性。
- 假设7【关键假设】：门控优先级高于z阈值层；被门判为“无法判定”的记录不进入后续染色体z阈值决策，最终输出为三分类：阳性/阴性/无法判定。依据：题目要求“女胎X染色体浓度没有异常，则可认为NIPT结果基本准确，否则难以保证结果准确性要求”；风险：门界过严会降低覆盖率，过松会纳入不可靠记录；可验证：报告覆盖率（可判定记录占比），并与z阈值层在可判定子集上的混淆结构对照。
- 假设8【关键假设】：第二层逐染色体z阈值仅在“可判定”记录上执行，染色体c∈{13,18,21}，判定式d_ijc=1{z_ijc≥τ_c}；女胎z阈值不能照搬男胎阈值，必须在女胎表上重新定标。依据：人类架构审核数据画像显示女胎z判别力弱：AUC t13约0.42、t18约0.55、t21约0.51，男胎分别为0.74/0.76/0.84（数据估计值）；风险：照搬男胎阈值会高估女胎z信号；可验证：输出女胎vs男胎各染色体z的AUC对比表，并报告女胎可靠子集上z≥3的命中情况。
- 假设9【关键假设】：τ_c选择采用孕妇分层交叉验证内部代价敏感优化；默认漏检:误报=3:1，目标为最小化CV内的代价Cost_c(τ)=3·FN_c(τ)+1·FP_c(τ)；敏感性取1:1、2:1、5:1。若CV选不出稳定操作点，则如实报告临床惯例阈值z≥3在本数据可靠子集上的表现，不伪造性能。依据：临床漏检代价高于误报，且女胎阳性稀疏、z信号弱；风险：代价权重不同会明显改变τ_c；可验证：输出各代价值下τ_c、灵敏度/特异度及bootstrap稳定性区间。
- 假设10【关键假设】：可选校准层不默认启用。校准定义为将z_ijc对X染色体浓度w_ij、染色体特异GC含量gc_ijc、被过滤读段比例filt_rate_ij、读段相关比例及BMI做残差化，得到z_adj_ijc；仅当校准后τ_c及性能有实质改善时才采纳校准版本，否则采用未校正确认版本。依据：人类架构审核要求奥卡姆剃刀，z判别力弱，校准细节应降级；风险：引入残差化可能造成过拟合和不稳定；可验证：比较校准前后CV性能、AIC/BIC和τ_c稳定性。
- 假设11【关键假设】：孕妇级判定与记录级判定必须分开报告。孕妇级在主方案采用保守合并：对染色体c，任一记录阳性即孕妇级阳性D_cons_ic=1(∃j:d_ijc=1)；同时报告多数票与最大风险得分两种合并规则。依据：人类架构审核数据画像显示44位阳性孕妇中43位同时存在AB阴性记录，其中26位仅1条阳性记录；风险：简单合并多数票可能漏检单次阳性；可验证：比较三种孕妇级合并规则（任一阳性、多数票、最大z/最大校准得分）的灵敏度/特异度/阳性一致率。
- 假设12【关键假设】：验证必须采用孕妇分块5折交叉验证（同一孕妇全部记录同折），cluster bootstrap采样次数B≥200，并做置换检验确认信号非随机。依据：重复测量结构与标签稀疏；风险：未分块验证会因同孕妇记录泄漏导致乐观指标；可验证：报告孕妇分块CV指标与行级独立CV指标的差异，并输出置换检验零位分布。
- 假设13【关键假设】：GC含量不采用40.0%–60.0%硬阈值剔除（原文：正常GC含量范围为40% ~ 60%）；女胎表GC范围约0.368–0.443，约36%低于0.40，无高于0.60（数据估计值），判断为平台系统偏差。GC仅作为连续协变量、质量权重或可选校准残差项进入分析，并做独立敏感性；不把低GC行直接剔除。风险：硬剔除会损失约1/3女胎样本并产生选择偏差；可验证：GC连续协变量 vs 质量权重 vs 硬阈值三种处理的覆盖率、混淆指标变化。
- 假设14：数值列多为非正态且部分列存在较强时序自相关（数据发现：多数数值列Shapiro-Wilk p<0.05，序号/年龄lag-1自相关>0.5），因此本问题主模型不假设独立正态；可靠性门采用中位数/MAD等稳健统计，验证采用孕妇分块重采样。依据：数据认知更新与数据智能摘要；风险：若忽视非正态或组内相关会低估不确定性；可验证：残差/得分分布诊断和按孕妇分块bootstrap区间。
- 假设15：临床背景常量仅作解释性边界，不进入三层判定主模型：可检测窗口为10.0–25.0周（原文：通常孕妇的孕期在10周~25周之间可以检测胎儿性染色体浓度）；男胎Y染色体浓度达标阈值4.0%仅适用于男胎，不用于女胎X浓度门（原文：如果男胎的Y染色体浓度达到或高于4%）；早期发现12.0周以内风险较低、中期发现13–27周风险高、晚期发现28.0周以后风险极高（原文：早期发现12周以内风险较低；中期发现13－27周风险高；晚期发现28周以后风险极高）。风险：把时间风险窗误入女胎异常判定会混淆时点优化与异常判定两个决策；可验证：输出无时间变量的主方案与加入孕周敏感性时的指标变化。
- 假设16：Coder必须按parse_hints解析字符串列：孕妇代码用df['孕妇代码'].str.replace('A', '', regex=False).astype(float)，末次月经用pd.to_datetime；AB列“染色体的非整倍体”的原始机器parse_hint为df['染色体的非整倍体'].str.replace('T', '', regex=False).astype(float)，但本问题必须扩展为多标签集合解析，不能将T13T18简单转成数值1318；检测孕周自行解析w+d格式。依据：机器数据列解析建议及AB列组合文本特征；风险：错误解析AB会破坏标签；可验证：输出AB列解析后的唯一值/多标签计数与前若干行核对。
**符号表**：
- N_f_rec: 女胎检测记录总数，数据画像估计约605（该值为数据估计值，非题目常量）
- n_f: 女胎孕妇总数，数据画像估计约147（该值为数据估计值，非题目常量）
- n_i: 第i位女胎孕妇的检测记录数，i=1,...,n_f
- i: 孕妇索引，i=1,...,n_f
- j: 孕妇内检测记录索引，j=1,...,n_i
- c: 目标常染色体索引，c∈{13,18,21}
- w_ij: 第i位孕妇第j次检测的X染色体浓度，无量纲比例，可能为负
- med_w: 女胎样本X染色体浓度的中位数
- MAD_w: 女胎样本X染色体浓度的MAD，MAD_w=median(|w_ij-med_w|)
- w_star_ij: X染色体浓度的稳健参考值，w_star_ij=(w_ij-med_w)/(1.4826·MAD_w)
- q10_f: 女胎样本w_ij的经验第10百分位，主门界，约-0.0272（该值为数据估计值）
- q5_f: 女胎样本w_ij的经验第5百分位，备选门界（数据估计值）
- q25_f: 女胎样本w_ij的经验第25百分位，备选门界（数据估计值）
- gate_ij: 第i位孕妇第j次记录的可靠性门输出，gate_ij=1表示w_ij<门界，判为无法判定；gate_ij=0表示可判定
- z_ijc: 第i位孕妇第j次记录的染色体c Z值，无量纲
- gc_ijc: 第i位孕妇第j次记录的染色体c特异GC含量，小数
- gc_total_ij: 第i位孕妇第j次记录的总GC含量，小数
- filt_rate_ij: 第i位孕妇第j次记录的被过滤读段数占总读段数的比例，小数
- map_rate_ij: 第i位孕妇第j次记录的总读段数中在参考基因组上比对的比例，小数
- unique_reads_ij: 第i位孕妇第j次记录的唯一比对读段数，个
- bmi_ij: 第i位孕妇第j次记录的BMI，单位：kg/m²
- x_ij: 校准协变量向量，默认包含w_ij、gc_ijc、filt_rate_ij、map_rate_ij、bmi_ij；用于可选残差化校准
- z_adj_ijc: 可选校准后的染色体c Z值残差，z_adj_ijc=z_ijc - f_ijc(x_ij)，其中f_ijc为线性回归拟合
- tau_c: 第c号染色体的z阈值，在孕妇分层CV内部按代价敏感选择
- lambda_FNFP: 漏检:误报代价比，主值3:1，敏感性取1:1、2:1、5:1
- d_ijc: 记录级逐染色体判定指示，d_ijc=1{z_ijc≥tau_c}，仅在gate_ij=0时执行
- y_ijc: 记录级AB标签，y_ijc=1(染色体c出现在AB_ij的解析集合中)；AB空白时全部为0
- AB_ij: 第i位孕妇第j次记录的AB列原始文本，如T18、T13T18、T21
- Y_ij: AB列解析出的多标签集合，Y_ij⊆{T13,T18,T21}
- p_hat_ijc: 可选逐记录染色体c阳性风险得分，由z_c或其校准值经回归/单调映射得到；仅用于合并规则比较，不强制进入主判定
- risk_ijc: 逐记录染色体c风险得分，主口径可取z_ijc或p_hat_ijc
- D_cons_ic: 孕妇级保守合并判定：任一记录阳性则阳性，D_cons_ic=1(Σ_j d_ijc≥1)
- D_maj_ic: 孕妇级多数票合并判定：D_maj_ic=1(Σ_j d_ijc > n_i/2)
- D_max_ic: 孕妇级最大风险得分合并判定：D_max_ic=1(max_j risk_ijc ≥ tau_c)
- coverage: 覆盖率，可判定记录数占女胎总记录数比例，coverage=(Σ_ij 1(gate_ij=0))/N_f_rec
- TP_c: 染色体c真阳性数
- FN_c: 染色体c假阴性数
- FP_c: 染色体c假阳性数
- TN_c: 染色体c真阴性数
- sens_c: 染色体c灵敏度，sens_c=TP_c/(TP_c+FN_c)
- spec_c: 染色体c特异度，spec_c=TN_c/(TN_c+FP_c)
- ppv_c: 染色体c阳性预测值，ppv_c=TP_c/(TP_c+FP_c)
- npv_c: 染色体c阴性预测值，npv_c=TN_c/(TN_c+FN_c)
- f1_c: 染色体c的F1分数，f1_c=2·TP_c/(2·TP_c+FP_c+FN_c)
- auc_c: 染色体c在可靠子集上的AUC，数据估计值见表
- Cost_c(tau): 染色体c阈值tau对应的CV代价，Cost_c=λ_FNFP·FN_c+FP_c，λ_FNFP默认3
- B: cluster bootstrap次数，B≥200
- CI_low_c: 染色体c评估指标的bootstrap 2.5%分位
- CI_high_c: 染色体c评估指标的bootstrap 97.5%分位
- perm_p_c: 染色体c置换检验p值，衡量指标是否高于随机零位
- K_fold: 孕妇分块交叉验证折数，K_fold=5
**公式/方程**：
- AB多标签解析：Y_ij=parse_set(AB_ij)⊆{T13,T18,T21}；AB_ij空白时Y_ij=∅。标签y_ijc=1(c∈Y_ij)，c∈{13,18,21}。
- X浓度稳健参考：w_star_ij=(w_ij-med_w)/(1.4826·MAD_w)，med_w=median_ij(w_ij)，MAD_w=median_ij(|w_ij-med_w|)。
- 可靠性门主判据：gate_ij=1{w_ij<q10_f}；q10_f为女胎w的经验第10百分位，主界约-0.0272（数据估计值）。gate=1输出无法判定，不进入z阈值层。
- 可靠性门敏感性：门界q∈{q5_f,q10_f,q25_f}，gate_ij(q)=1{w_ij<q}；绘制coverage-准确率权衡曲线。
- 逐染色体记录级判定（仅在gate_ij=0）：d_ijc=1{z_ijc≥tau_c}，c∈{13,18,21}。
- 可选校准层：拟合z_ijc=α_c+x_ij'γ_c+e_ijc，在校正集上得到z_adj_ijc=z_ijc-x_ij'γ_hat_c；仅当校准后稳定性能和CV代价有实质改善时，tau_c用z_adj_ijc定标，否则采用未校准z_ijc。
- 阈值代价函数：Cost_c(tau)=λ_FNFP·FN_c(tau)+FP_c(tau)，λ_FNFP默认3（漏检:误报=3:1）。在孕妇分块CV内对tau进行网格搜索，tau_c=argmin_tau Cost_c(tau)。
- 敏感性代价比：λ_FNFP∈{1,2,3,5}对应1:1、2:1、3:1、5:1，输出不同λ下tau_c及指标。
- 记录级混淆指标：sens_c=TP_c/(TP_c+FN_c)，spec_c=TN_c/(TN_c+FP_c)，ppv_c=TP_c/(TP_c+FP_c)，npv_c=TN_c/(TN_c+FN_c)，f1_c=2·TP_c/(2·TP_c+FP_c+FN_c)。
- 覆盖率：coverage=(Σ_ij 1(gate_ij=0))/N_f_rec；三分类输出分布：n_invalid=Σ_ij gate_ij，n_positive_c=Σ_ij 1(gate_ij=0,d_ijc=1)，n_negative_c=Σ_ij 1(gate_ij=0,d_ijc=0)。
- 孕妇级保守合并：D_cons_ic=1(Σ_j d_ijc≥1)，主口径用此规则；同时报告多数票D_maj_ic=1(Σ_j d_ijc > n_i/2)和最大风险得分D_max_ic=1(max_j risk_ijc ≥ tau_c)。
- 风险得分定义：risk_ijc=z_ijc（主口径）或risk_ijc=p_hat_ijc（可选校准映射）。最大风险得分合并仅在记录级tau_c下比较，不单独重选阈值。
- 孕妇分块5折CV：将n_f位孕妇随机均分为5折，每折包含该孕妇所有记录；训练折内选tau_c，验证折内计算标签级指标，合并为CV estimate。
- Cluster bootstrap区间：对n_f位孕妇有放回重采样B≥200次，每次在重采样样本内执行门界+CV阈值+合并，得到指标集，CI_low_c=Q_{0.025}({metric_c^b})，CI_high_c=Q_{0.975}({metric_c^b})。
- 置换检验：随机打乱AB标签（保持孕妇结构）重复计算指标，构造零位分布，perm_p_c=零位中指标≥观测指标的比例；若perm_p_c≥0.05，信号不显著。
- GC连续协变量/权重敏感性：主模型按连续gc_total_ij/gc_ijc纳入校准或加权；硬阈值对照仅报告，不采用；硬阈值判据依据原文正常GC含量范围为40.0% ~ 60.0%，但女胎表GC范围约0.368–0.443，因此硬阈值会损失大量记录。
- X浓度处理敏感性：主口径保留负值；对照口径w_trunc_ij=max(w_ij,0)，比较两种口径的coverage和指标。
- 女胎vs男胎对照：对同一z阈值法在女胎可靠子集与男胎表数据上分别估计AUC并报告差异，作为女胎必须单独定方法的证据。
**解题思路**：步骤1：读取附件.xlsx。Coder必须按照parse_hints解析字符串列：孕妇代码用df['孕妇代码'].str.replace('A', '', regex=False).astype(float)；末次月经用pd.to_datetime；AB列“染色体的非整倍体”的机器parse_hint为df['染色体的非整倍体'].str.replace('T', '', regex=False).astype(float)，但本问题必须扩展为多标签集合解析，不能将T13T18简单转为数值1318；检测孕周自行解析w+d格式。步骤2：筛选女胎子集：U列Y染色体Z值与V列Y染色体浓度为空白/缺失的记录构成女胎表，交叉核查孕妇代码、X染色体浓度W列及AE列；输出女胎记录数约605、孕妇数约147的数据画像。步骤3：数据画像：报告年龄、身高、体重、BMI、GC、各染色体GC、读段比例、被过滤比例、z值、X浓度的缺失率/范围/非正态性；检查BMI=体重/身高²自洽性；统计AB组合计数（T13/T18/T21/多标签）和孕妇内标签稳定性；确认AE列全部“健康”。步骤4：构造逐染色体标签：解析AB→y_ijc，空白全0；明确不采用AE交叉清洗。步骤5：第一层可靠性门：对女胎全样本w_ij计算med_w、MAD_w、w_star_ij及q5_f/q10_f/q25_f；主界采用w_ij<q10_f（约-0.0272）判无法判定，q5_f/q25_f为敏感性；输出gate_ij与覆盖率。步骤6：第二层逐染色体z阈值：仅在gate_ij=0的记录上，对c∈{13,18,21}在孕妇分块5折CV内部按代价敏感目标选择τ_c，默认漏检:误报=3:1，敏感性1:1、2:1、5:1；可选校准：将z_ijc对w_ij、染色体GC、被过滤比例、比对比例、BMI做线性残差化，比较校准前后CV代价和τ_c稳定性，若无实质改善则采用未校准z。输出d_ijc及可靠子集混淆矩阵、sens/spec/ppv/npv/f1/auc。诚实边界：若CV选不出稳定操作点，报告临床惯例z≥3在可靠子集上的表现。步骤7：孕妇级合并：对每染色体分别计算任一记录阳性（主）、多数票、最大风险得分三种合并规则；报告孕妇级混淆指标。步骤8：验证：孕妇分块5折CV防止同孕妇记录泄漏；cluster bootstrap≥200给出指标区间；置换检验检验信号非随机。步骤9：敏感性分析：门界q5/q10/q25的coverage-准确率曲线；AB标签空白=阴性 vs 剔除孕妇内混合标签；GC连续协变量/质量权重/硬阈值；X浓度保留负值 vs 截断0；主方案 vs 全特征惩罚逻辑回归 vs 男胎同一方法对照。步骤10：交付图表：三层判定流程图、可靠性门覆盖率-准确率曲线、各染色体阈值表和可靠子集混淆矩阵、记录级与孕妇级对比表、女胎/男胎z判别力对比表、孕妇级合并规则比较表；论文主动讨论“门价值＞阈值层”以及女胎必须单独定标的证据。

## 3. 算法与求解
**算法摘要**：
**伪代码/实现步骤**：
1. 步骤1: 读取附件.xlsx；按parse_hints解析：孕妇代码=str.replace('A','',regex=False).astype(float)，末次月经=pd.to_datetime，检测孕周解析w+d，AB列不做数值转换而保留原文本。
2. 步骤2: 筛选女胎子集 df_female = df[df['Y染色体Z值'].isna() & df['Y染色体浓度'].isna()]；输出记录数N_f_rec=len(df_female)、孕妇数n_f=df_female['孕妇代码'].nunique()；交叉核查孕妇代码、X染色体浓度W列、AE列，确认AE列全部为'健康'。
3. 步骤3: 构造多标签：定义 parse_ab(s) -> set，若s为空返回空集合；对c∈{13,18,21}设置 y_ijc = 1(c in parse_ab(AB_ij))；输出AB列唯一组合计数与前若干行核对。
4. 步骤4: 数据画像：计算w_ij范围、中位数med_w、MAD_w、q5_f/q10_f/q25_f；统计GC总/染色体特异范围及低于0.40比例；统计AB组合计数和孕妇内标签稳定性；AE列取值计数。
5. 步骤5: 第一层可靠性门：gate_ij = 1(w_ij < q10_f)；coverage = 1 - mean(gate_ij)；输出q5/q10/q25下覆盖率和可靠子集阳性/阴性标签比例。
6. 步骤6: 准备第二层：仅使用gate_ij==0的记录；定义候选tau网格 = np.arange(-3.0, 6.0, 0.5)，保证步长≥0.5；染色体列表c_list=[13,18,21]。
7. 步骤7: 孕妇分块5折CV阈值选择：使用sklearn.model_selection.GroupKFold(n_splits=5)，groups=df_female['孕妇代码']；在每个训练折内对每个染色体c，遍历tau网格，计算标签级FN_c和FP_c（记录级），代价Cost=3*FN_c+FP_c；选择最小代价tau_c_train；在验证折上应用d_ijc=1(z_ijc>=tau_c_train)并与y_ijc比较，累计混淆计数。
8. 步骤8: 若某折某染色体训练集中无阳性导致所有tau代价相同，则tau_c_train设为临床惯例3.0，并在结果中标注；若验证折AUC无法计算（单类），使用try-except将auc设为np.nan后填充为0.5。
9. 步骤9: 可选校准层仅作敏感性：对可判定记录拟合线性回归 z_ijc ~ w_ij + gc_ijc + gc_total_ij + filt_rate_ij + map_rate_ij + bmi_ij，得到残差z_adj_ijc；在训练折内重选tau_c_adj；比较未校准与校准CV代价及稳定性；主方案采用未校准。
10. 步骤10: 记录级判定整合：对全部女胎记录，若gate_ij=1则输出无法判定；若gate_ij=0且d_ijc=1输出阳性；否则阴性；生成记录级三分类表df_pred_record。
11. 步骤11: 孕妇级合并：对每染色体c，计算任一记录阳性D_cons_ic=1(Σ_j d_ijc>=1)；多数票D_maj_ic=1(Σ_j d_ijc>n_i/2)；最大风险得分D_max_ic=1(max_j z_ijc>=tau_c)；分别计算孕妇级混淆矩阵和指标。
12. 步骤12: 指标计算：定义函数metrics(y_true, y_pred)返回TP,FN,FP,TN,sens,spec,ppv,npv,f1；所有概率/比例指标用np.clip限制在[0,1]；AUC使用roc_auc_score，单类时返回0.5。
13. 步骤13: Cluster bootstrap：B=200；对n_f个孕妇有放回重采样（保持孕妇所有记录），每次重采样在重采样集上重新计算q10_f、CV内tau_c、记录级d及孕妇级合并，记录coverage、sens、spec、ppv、f1、auc等指标；得到B个bootstrap样本。
14. 步骤14: 置换检验：保持孕妇结构，随机打乱AB标签（在孕妇代码内重排标签或跨孕妇保持块结构），重复计算记录级AUC或F1，构造零位分布；perm_p = (1+Σ(metric_perm>=metric_obs))/(1+B_perm)；B_perm=100；若perm_p>=0.05则不宣称信号非随机。
15. 步骤15: 敏感性分析循环：依次改变门界q∈{q5_f,q10_f,q25_f}，代价比λ∈{1,2,3,5}，X浓度截断vs保留，GC硬阈值剔除vs连续协变量vs质量权重，AB混合标签剔除，校准层开关，加入孕周变量，男胎同一方法对照；每个场景输出关键指标到sensitivity_results列表。
16. 步骤16: 导出结果：df_summary（含层级、染色体、规则、tau、coverage、sens、spec、ppv、npv、f1、auc及bootstrap上下限）保存到 MODELING_OUTPUT_DIR/results/q4.csv；同时生成 output.csv 为兼容输出，包含主要指标；确保无NaN/Inf，数值列clip。

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
- 小题 3（passed）：outputs\results\q3.csv

## 5. 预期图表
- fig_roadmap [roadmap] 总体技术路线图，展示数据流与三层判定结构（数据来源：）
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
- fig_q4_data_profile [boxplot] 展示女胎X浓度分布特征、负值比例、AB标签组合及AE列计数，支撑数据构造物证（数据来源：results/q4_data_profile.csv）
- fig_q4_gate_curve [line] 展示不同门界下覆盖率和可靠子集准确率的变化，验证门界选择及其不确定性（数据来源：results/q4_gate_sensitivity.csv）
- fig_q4_z_roc [line] 展示女胎z值判别力弱，说明必须单独定标，且与男胎对照（数据来源：results/q4_z_roc.csv）
- fig_q4_cost_sens [line] 展示代价权重λ对tau_c和代价曲线的影响，验证关键假设9（数据来源：results/q4_cost_sensitivity.csv）
- fig_q4_merge_compare [bar] 比较保守合并与其他合并规则，支持任一阳性作为主规则（数据来源：results/q4_merge_metrics.csv）
- fig_q4_calibration [scatter] 验证校准层是否需要开启，展示未校准与校准的差异（数据来源：results/q4_calibration_sens.csv）
- fig_q4_gc_sens [line] 验证GC不硬剔除假设，展示硬剔除会损失样本并改变指标（数据来源：results/q4_gc_sens.csv）
- fig_q4_w_trunc [line] 验证关键假设5和17，展示负值信息不可丢弃（数据来源：results/q4_w_trunc_sens.csv）
- fig_q4_cv_bootstrap [boxplot] 展示不确定性区间和信号显著性，验证孕妇分块假设（数据来源：results/q4_bootstrap_metrics.csv）
- fig_male_female_auc [bar] 证明女胎必须单独定标，不照搬男胎阈值（数据来源：results/q4_male_female_auc.csv）
- fig_q4_time_sens [line] 验证临床时间窗变量不进入主模型假设（数据来源：results/q4_time_sens.csv）

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
- tab_q4_data_profile：表1 女胎数据画像与标签构造核查（列：指标, 值, 备注；呈现女胎记录数、孕妇数、w中位数/范围、负值比例、AB组合计数、AE列计数、GC范围，支撑数据画像和构造验证）
- tab_q4_gate：表2 X浓度可靠性门覆盖率与可靠子集标签分布（列：门界, q值, 覆盖率, 可靠子集阳性比例, 可靠子集阴性比例, 无法判定记录数；展示q5/q10/q25三种门界下的覆盖率和标签结构，验证门控优先级与门界敏感性）
- tab_q4_thresholds：表3 各染色体z阈值及记录级指标（列：染色体, tau_c, 代价比λ, TP, FN, FP, TN, sens, spec, ppv, npv, f1, auc, bootstrap_CI；报告第二层阈值选择结果和记录级性能，体现CV代价敏感选择）
- tab_q4_merge：表4 孕妇级合并规则指标对比（列：染色体, 合并规则, TP, FN, FP, TN, sens, spec, ppv, f1；比较任一阳性、多数票、最大风险得分三种合并规则，支持保守合并主规则）
- tab_q4_sensitivity：表5 关键假设敏感性实验汇总（列：关键假设, 扰动/对照, 比较指标, 期望结论；逐条覆盖关键假设4,5,6,7,8,9,11,13，提供扰动实验证据）
- tab_q4_bootstrap：表6 cluster bootstrap和置换检验结果（列：染色体, 指标, 观测值, CI_low, CI_high, perm_p；报告不确定性和信号显著性，验证孕妇分块与置换检验）

## 7. 结果契约
```json
{
  "description": "女胎异常判定方法的关键指标与阈值表，包含记录级和孕妇级结果，按染色体分层输出三分类指标和不确定性",
  "allow_single_row": false,
  "min_rows": 12,
  "max_rows": 200,
  "columns": [
    {
      "name": "level",
      "dtype": "category",
      "min": null,
      "max": null,
      "distinct_required": false,
      "description": "记录级record或孕妇级pregnant"
    },
    {
      "name": "chrom",
      "dtype": "category",
      "min": null,
      "max": null,
      "distinct_required": false,
      "description": "目标染色体c，取13/18/21"
    },
    {
      "name": "rule",
      "dtype": "category",
      "min": null,
      "max": null,
      "distinct_required": false,
      "description": "孕妇级合并规则，conservative/majority/max_risk；记录级填na"
    },
    {
      "name": "tau",
      "dtype": "float",
      "min": -5.0,
      "max": 10.0,
      "distinct_required": false,
      "description": "该染色体z阈值，临床惯例3.0为备选"
    },
    {
      "name": "coverage",
      "dtype": "float",
      "min": 0.0,
      "max": 1.0,
      "distinct_required": false,
      "description": "可靠性门覆盖率，可判定记录占比"
    },
    {
      "name": "sens",
      "dtype": "float",
      "min": 0.0,
      "max": 1.0,
      "distinct_required": false,
      "description": "灵敏度/召回率"
    },
    {
      "name": "spec",
      "dtype": "float",
      "min": 0.0,
      "max": 1.0,
      "distinct_required": false,
      "description": "特异度"
    },
    {
      "name": "ppv",
      "dtype": "float",
      "min": 0.0,
      "max": 1.0,
      "distinct_required": false,
      "description": "阳性预测值"
    },
    {
      "name": "npv",
      "dtype": "float",
      "min": 0.0,
      "max": 1.0,
      "distinct_required": false,
      "description": "阴性预测值"
    },
    {
      "name": "f1",
      "dtype": "float",
      "min": 0.0,
      "max": 1.0,
      "distinct_required": false,
      "description": "F1分数"
    },
    {
      "name": "auc",
      "dtype": "float",
      "min": 0.0,
      "max": 1.0,
      "distinct_required": false,
      "description": "ROC曲线下面积，单类时为0.5"
    },
    {
      "name": "ci_low",
      "dtype": "float",
      "min": 0.0,
      "max": 1.0,
      "distinct_required": false,
      "description": "cluster bootstrap 2.5%分位"
    },
    {
      "name": "ci_high",
      "dtype": "float",
      "min": 0.0,
      "max": 1.0,
      "distinct_required": false,
      "description": "cluster bootstrap 97.5%分位"
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

# 本批趋势分析交接

## 数据与业务口径
本项目只分析一批单变量观测，未知参数为截距和斜率。源观测满足 y_i=截距+斜率*t_i+误差_i；不同id的误差相互独立、均值为零、服从正态分布，variance是已知测量方差。相同id无论被引用多少次，始终是同一次测量。

数据包为JSON对象，observations中每行有唯一字符串id及time、value、variance；windows中每项有唯一id及weights对象，键为源id，值为非负权重且和为1。数字用整数或[分子,正分母]表达，不接受近似小数。本题只覆盖2～16条源观测、2～12个窗口、整体窗口时间不全相同的合法输入；无需扩展解析器。窗口时间、读数均按同一组权重加权。plans把窗口id分为非空块，每个窗口恰好出现一次；分块是处理安排，不代表重新采样。

完整批次按已知测量误差做广义最小二乘，不能改用残差估计噪声或t分布。窗口可线性相关，重复窗口不增加独立测量信息。对奇异协方差，可选取窗口权重行的任一极大线性无关子集进行等价拟合，不要求特定伪逆实现；不能丢弃仍独立的信息。information_rank定义为完整窗口误差协方差的秩，不是参数个数或窗口个数。95%区间固定为斜率±1.96×标准误。源误差独立和线性传播确定全部窗口协方差，分块不能改变同一批次的规范结果。

## 已有入口
python3 pipeline.py data/incident.json cut
从指定包读取plans中的计划，标准输出一个JSON对象。result有intercept、slope、slope_variance、standard_error、ci95、information_rank、window_count；trace保留窗口权重和加权值、每块采用的basis_windows、协方差、information、rhs、trace_loadings以及合并摘要。trace_loadings是该块rhs两分量相对各源value的线性系数，供血缘追溯。精确数输出为[分子,分母]，标准误及区间为浮点数。Python入口pipeline.run(package,plan)返回同构字典，精确数使用Fraction。

python3 checks.py normal
正常包是4个互不共享来源的窗口，三种分块均应返回intercept=1、slope=2、slope_variance=1/5、information_rank=4。

python3 checks.py local
检查两个包全部公开计划的块内协方差对称性、可逆性及信息矩阵的非负性。这些局部检查在当前实现下通过，不构成完整批次正确性的保证。

python3 checks.py incident
事故包有6个来源、8个窗口及5种计划。入口在计划间比较slope_variance与information_rank，不一致时以TREND_PARTITION_DRIFT失败。logs/是实际执行所得业务JSON和异常末行摘录，不包含诊断结论。

## 本轮交付
保留matrix.py、resample.py、summary.py、trend.py、pipeline.py、checks.py及原始data/、logs/，不实施修复。新增REPORT.md、analysis.py和evidence.json，分析只需覆盖本包和下述有限对照，不需要通用统计框架或全组合搜索。

REPORT.md引用源码行号与日志计划名，定位首次偏离业务口径的中间状态，说明根因、必要触发条件及本包影响；分别核定各计划的截距、斜率、信息秩、报告标准误及区间。不能只说“共享观测导致相关”，也不能把斜率接近视为结果正确。

analysis.py从原始数据独立重建窗口协方差和规范拟合；可复用matrix.py的基础算术，但独立期望不得调用resample.materialize、summary、trend或把all计划当真值。对每种计划，还要从原始观测推导现有斜率估计量的线性系数和实际方差，与其报告方差分别比较；现有估计量实际方差不能直接拿规范拟合方差代替。用精确有理数核对权重、协方差、秩、线性系数和方差，标准误及区间的绝对误差不超过1e-10。

对照只改变一个因素：原包仅换分块；为完整批次追加一个与w3共用来源及权重的重复窗口；把每个窗口引用的源替换为该引用独有的新id，保留原time、value、variance和权重，代表数值相同但误差独立的新测量。明确这些操作的区别及其估计、秩和区间影响；派生数据在分析脚本中构造，不覆盖原件。追加重复窗口比较all与全部单窗口分块，独立id对照运行原有五种计划。无需讨论未提供的历史批次、在线更新或性能问题。

evidence.json保存可重算数值、对照变换和两次原事故失败记录，REPORT.md给出与证据相符的结论；无需固定JSON层级。python3 analysis.py verify应重算并核对证据及报告中的数值锚点，运行原正常与局部回归，并再次连续两次运行原事故入口，确认均非零退出且含TREND_PARTITION_DRIFT；任一不符则失败。python3 analysis.py demo展示规范结果、现有结果及一个真实单因素差异。两种模式均不修改原始材料。

## 运行交接
.pointwise-runtime.json是工具读取的运行清单：version为1，kind为cli，setup是argv数组列表，verify与demo各为一个argv数组，不使用shell命令串。当前清单展示既有业务；交付时将verify接至analysis.py verify、demo接至analysis.py demo，使清单直接重放报告证据。

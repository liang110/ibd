# 下一轮论文重构方案（基于 V24 两份审稿意见）

本方案只处理审稿意见指向的科学断点，不扩充模型目录，不引入新的深度学习架构，也不改变冻结的主要估计对象。所有修改必须回写唯一构建源，再整包生成主文、Supplementary、图表和审计文件。

## 审稿意见的共同结论与取舍

两份报告虽然对投稿成熟度的评分不同，但对真正影响可信度的问题高度一致：

- 直接优先级最高的是主 benchmark 冻结预测与 calibration 分析的估计口径一致性；
- Nielsen 的 country × phenotype 构成必须先做可识别性审计，不能用合成 batch control 代替真实混杂诊断；
- 126 物种契约是标签无关的 transductive 坐标契约，不是 strict source-only deployment，也不证明跨平台生物学等价；
- participant aggregation、source perturbation ensemble 和模型比较必须拆成不同 estimand；
- `nested CV`、模型纳入时序和点估计排名的统计措辞必须收窄；
- 代码、输入哈希、冻结预测、环境和重建入口需要形成第三方可运行的最小闭环。

以下建议是条件性 P1，而非阻塞投稿的理由：来源端 participant-equal weighting、半合成效应梯度、matched-source-size 分析和更多图形类型。若原始输入或样本量不能支持这些分析，应将其标为未估计，不通过补充叙述替代证据。两份报告均不支持引入 Transformer、图神经网络或更多分类器作为“创新”替代方案。

## 1. 先锁定估计对象与分析分层

主分析继续限定为固定有序 126 物种 HMP2–PRISM transductive contract、来源端拟合、6 个有向 source-to-target edge、10 次 source-side perturbation 和参与者聚簇 bootstrap。把以下对象分开：

- sample-level split-averaged metric；
- participant-aggregated metric；
- source perturbation SD；
- target participant-cluster CI；
- source-model perturbation ensemble；
- participant-level aggregation。

不得把这些层级混合后再比较模型。`nested cross-validation` 改称为 `repeated source-side model-selection perturbations with grouped inner CV`，因为 20% reserve 不承担传统 outer-test 估计。

## 2. P0：校准与主 benchmark 必须使用同一冻结预测

对 42 个主方向–模型组合，直接从完整 source refit 产生的冻结 target predictions 重新计算 Brier、calibration intercept、calibration slope、ECE 和 reliability bins。每个 bootstrap replicate 使用同一套目标参与者重采样权重，同时应用于所有模型和 source perturbation。

prevalence baseline 为常数概率时 calibration slope 标为 `not estimable`，不得报告不可识别的有限拟合值。Table 2、Supplementary calibration tables、Figure 4 和正文只能引用同一估计口径。Brier skill 可以作为派生描述量，但不得替代原始 Brier。

验收：每个方向–模型组合有唯一 prediction hash、唯一 bootstrap weight hash、统一行数和统一指标定义。

## 3. P0：明确 participant aggregation 与 ensemble aggregation

分别报告：

1. 先对每个 source perturbation 形成样本级预测，再计算指标；
2. 先对同一 target participant 的重复样本聚合概率，再计算 participant-level 指标；
3. 多个 source perturbation 的概率集成；
4. 对不同模型的配对比较。

这四者不得在一个“ensemble AUROC”字段中混用。正文以 observed cohort-edge/model-conditional ordering 描述点估计；只有配对区间不跨 0 的比较才可称为有区间支持的差异。

## 4. P0：Nielsen 地域与疾病构成审计

先输出 country × phenotype 的样本数和参与者数列联表，再判断子组是否可识别。只有同时包含病例和对照且数量足够的分层，才计算分层 AUROC、PR-AUC、Brier 和 participant-cluster CI；缺失类别的分层标为 `not estimable`。不得把合成 batch-only control 当作真实 Nielsen 混杂已被控制的证据。

结果只能被表述为地域/疾病构成敏感性诊断，不能延伸为临床混杂控制、因果解释或广泛人群泛化。

## 5. P0：特征契约与历史不等价运行

补充每队列 126 个特征的非零率、保留丰度质量占比、canonicalization 碰撞数、排除数和原始 header→canonical ID crosswalk。明确“计算坐标一致”不等于“生物学测量等价”。

严格 source-specific feature-universe deployment 只有在恢复原始 HMP2/PRISM 矩阵后才能重跑；恢复前不得猜测，也不得把 Table 11 归档运行称为定量敏感性分析。若历史 `cap126_prev0p10` 仍与主流程不等价，必须说明其输入/执行状态差异并从定量稳健性结论中排除。

## 6. P1：模型纳入时序和权重敏感性

记录正式模型家族和超参数网格在任何目标评估前已锁定。探索性图模型或未纳入模型应单列为 exploratory，不得写成事后依据外部结果排除。若计算资源允许，增加 participant-equal weighting 或每人单次抽样的来源端敏感性；否则把训练样本权重限制写成未估计。

## 7. 结果和图表重构顺序

- Figure 1：信息边界和固定契约；
- Figure 2：已知真值、身份泄漏和 clinical-proxy rejection；
- Figure 3：六方向 discrimination（AUROC、PR-AUC、配对差值）；
- Figure 4：同一冻结预测下的 Brier skill、校准曲线和 calibration slope；
- Figure 5：participant aggregation、ensemble aggregation 与来源/地域敏感性；SIAMCAT 移到 Supplementary protocol reference。

图注必须标明统计单位、中心值、不确定性、higher/lower-is-better 和不可估计项。图表、正文和机器可读 CSV 必须从同一 registry 生成。

## 8. 最终门禁

- 126 物种唯一性、crosswalk、契约哈希：pass；
- 42 组主冻结预测与校准指标可追溯：pass；
- 共同 bootstrap 权重：pass；
- prevalence baseline calibration slope 正确标记：pass；
- Nielsen country × phenotype 可识别性：pass；
- Table 11 不进入定量稳健性结论：pass；
- citation、公式、图注、表格字段、摘要和 Supplementary 一致性：pass；
- GitHub 最小包可在干净环境中运行 unit tests 和 smoke test：pass。

另设一个不改变科学估计对象的投稿合规门：根据作者实际使用过程核对 AI 辅助写作、代码和图表工具的披露是否符合 BIB/ISCB 当期政策。不能通过删除披露句来规避政策，也不把该合规问题包装成统计修复。

只有 P0 门禁全部通过后，才允许进行正文重构。若 P0 证据不能支持更强的主线，论文应保持“固定契约下的可审计三队列 protocol benchmark”定位，而不通过措辞包装成广泛迁移可靠性发现。

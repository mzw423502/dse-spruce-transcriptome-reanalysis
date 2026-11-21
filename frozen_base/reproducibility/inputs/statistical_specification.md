# 宿主四组因子检验预冻结规格

状态：主统计问题与阈值已在正式表达结果产生前写定；R 环境尚受系统执行许可阻塞，以下代码尚未通过真实 R 包执行 QA。系统 NumPy 设计代数 QA 与 R 包、负二项模型、长度 offset 验证是不同状态，不能互相替代。冻结时间、SHA256 以本目录的冻结回执为准。后续只可记录有原因的实现修正，不可按显著数修改检验方案。

## 1. 分析对象和可识别量

仅分析 PRJNA734471 的宿主 RNA。四组 Pa、PaPs、PaHp、PaPsHp 各有三个独立的生物学混样，统计 n=12；每个混样包含的幼苗数不增加 n。DSE=Phialocephala sphaeroides，病原=Heterobasidion parviporum。样品是单终点整体幼苗。每次接种的因果解释还取决于原实验随机化信息，表达结果本身不识别免疫机制、定殖量或抗病因果中介。

主参考为 64,238 个 exact CDS sequence equivalence units；每个 target_id 与 feature_id 恒等。它们是序列计数单元，不保证分别对应一个独立基因。含多个原 MA 标识的序列单元不拆分或复制计数。主注释采用共同条目，union 仅作标明的注释敏感性分析；KEGG map 编号不是 KO。

目标量是经宿主文库组成和有效长度调整的相对宿主 RNA 丰度变化。没有外部 spike-in，不能解释成每株绝对转录量；宿主归一化也不能估计真菌生物量。原文 5,269 个差异条目属于经筛选的结果，不能用作检验总体。

## 2. 主模型、五个对比和检验总体

对序列单元 j、样品 s 拟合负二项模型：

`log2(mu_js) = log2(o_js) + beta0_j + betaD_j D_s + betaP_j P_s + betaI_j D_s P_s`。

其中 o_js 为标准归一化因子，包含文库组成及 tximport 提供的样品有效长度信息。计算软件的内部链接函数为自然对数，输出效应转换为 log2 fold change。设计式固定为 `~DSE*pathogen`，0 为参考。不得把 betaD 或 betaP 误写成跨另一个因子平均的主效应。

| 名称 | 系数 | 组间表达比较 |
|---|---|---|
| 主交互 | betaI | (PaPsHp−PaPs)−(PaHp−Pa)，各差值位于 log2 丰度尺度 |
| 无病原时 DSE 效应 | betaD | PaPs−Pa |
| 有病原时 DSE 效应 | betaD+betaI | PaPsHp−PaHp |
| 无 DSE 时病原效应 | betaP | PaHp−Pa |
| 有 DSE 时病原效应 | betaP+betaI | PaPsHp−PaPs |

主交互等价于倍数变化的比值的 log2：log2[(PaPsHp/PaPs)/(PaHp/Pa)]。它判断在 log 尺度上偏离加和，而不是原始计数的加法偏差。正交互只表示有 DSE 时病原效应向正方向移动，不能单凭正号称“协同促生”或“免疫增强”。

预过滤为至少三个样品的 kallisto estimated fragments ≥10。阈值按全体样品定义，独立于分组的效应大小和显著性；保留原始计数是小数这一事实，DESeq2 导入器按其标准流程处理。过滤后所有单元组成固定检验总体，输出全部通过和不通过单元及原因。

主分析使用 DESeq2、betaPrior=FALSE、Wald 检验、parametric dispersion trend（软件若自动改用 local，记录警告和实际拟合类型）、parallel=FALSE。不自动加入批次、潜变量或 SVA。已知批次若在分析前获得，须另立有时间戳的可识别性修订，不能按 PCA 效果选择。

## 3. 多重检验、区间和异常值

- 主检验族：一个交互对比 × 全部固定合格序列单元，BH FDR<0.05。
- 次级简单效应检验族：四个简单效应 × 全部固定合格单元，共同 BH FDR<0.05。另输出各对比 BH 和所有五个对比合并 BH，字段名明确，不能挑较小值宣称通过。
- NA 检验仍占计划检验总数；不可通过删除 Cook 异常行缩小 BH 分母。DESeq2 默认 padj 列保留为软件输出，主要判断使用明确命名的 fixed-universe 列。
- 效应采用未收缩 log2 fold change、标准误和点态 95% Wald 区间。区间不是同时区间，也未经显著后选择校正；不把未显著当等效或无效。|log2FC| 不另设主检验门槛；实际大小、基线丰度和区间共同解释。
- 维持默认 Cook 判定，但禁止自动替换计数，`minReplicatesForReplace=Inf`。不因结果“不好看”删样。保存 Cook 矩阵、迭代收敛字段、归一化因子和 NA 原因。任何技术排除必须来自独立的文库完整性/质量证据并单独报告。

## 4. 导入和成熟方法敏感性

主导入：`tximport(type="kallisto", countsFromAbundance="no")`，随后 `DESeqDataSetFromTximport`，由其处理平均有效长度的 normalization factors。不得直接对 TPM/FPKM 拟合 DESeq2，也不得把长度调整重复添加两遍。

预设敏感性：同一固定单元总体，重新导入 `countsFromAbundance="lengthScaledTPM"`，使用 edgeR TMM + robust dispersion/QL 拟合、同一因子设计与五个对比，不再加长度 offset。比较交互效应方向、量级及检验一致性；不是两个独立数据集，不能称外部复验。方法不一致时保留分歧，不选显著者作为主结果。

完整参考、去歧义/低质量参考、文库质控和链方向的敏感性由参考与流水线协议另行冻结。统计函数不自行改变其规则。逐样本留一如后续实施，只评价估计受单库影响的程度；不得要求 n=11 重算后仍显著作为稳定性门槛，也不把它称未见样本预测验证。

## 5. 预注释模块和响应方向

采用 limma voom + mroast（9,999 次旋转、固定种子、mean statistic、midp=FALSE），在注释而非表达结果定义的集合上检验交互，并平行检验两种 DSE 背景下的病原效应。只纳入过滤后含 15–500 个序列单元的集合；去重单元与成员，不按显著性选择集合。输入主集合须在正式结果前冻结，记录来源、集合数、成员、共享条目规则与集合重叠。旋转检验保留基因间相关结构，但序列碎片和同源成员仍可能导致功能证据重复，不能把集合显著解释为一条独立机制被证实。

三个对比 × 合格集合组成一个 directional 检验族，BH；mixed 检验另列为探索性检验族，明确区分。mroast 是 self-contained 检验，回答集合内有无表达响应，不等于相对全基因组过度富集。仅原理性功能注释不支持特定植物通路“激活”；动物命名 KEGG 条目须按共享同源注释解释或排除于植物机制叙述。

响应方向图的横轴为 `PaHp−Pa`，纵轴为 `PaPsHp−PaPs`，都为 log2 丰度尺度。两轴来自不重叠的组对，避免直接共享均值产生的主要算术耦合；它们仍共享归一化与离散度估计，因此不能把跨单元相关系数当独立证据。呈现效应和不确定性、零线与同效应线。

每个单元同时对 P0、P1、I 输出 Bonferroni 校正的三个边际 Wald 区间（各 98.3333%，形成近似至少 95% 的单元内联合覆盖；不同时覆盖所有单元，也未校正结果筛选）。方向类别仅作以下保守、可回溯的描述：先要求 I 的主族 FDR<0.05、P0/P1 均通过四简单效应共同 FDR<0.05、三个联合区间均排除零；P0/P1 异号为“有支持的方向相反”，同号且 I 与 P0 同号为“同向、幅度较大”，同号且 I 与 P0 异号为“同向、幅度较小”。这些检验的合取不自动具有类别层面的 FDR 保证，必须保留此限制。其余标记“不确定性或效应改变尚未解决”，不能仅凭一个显著、另一个不显著称新生响应、消失响应或反转。没有进行预设容许界的等效检验，因此不得把未显著 I 称同等效应或方向稳定。该响应平面是解释和展示方案，不是重新命名交互系数的新算法。

禁止使用 betaD 与 betaI 的朴素负相关证明“基线预激越强、病原反应越弱”。在四组估计误差独立等方差时，Cov(betaD,betaI)=−Var(PaPs)−Var(Pa)，即使真实效应皆为零，相关期望也为 −1/sqrt(2)。本项目合成设计 QA 专门检验这一点。

## 6. 实现状态、验证和来源

`stat_primary.R` 已准备输入断言、固定总体、五对比、CI、两套成熟方法、模块、全量诊断导出；正式执行需要 R 包加载和合成数据 QA 先通过。`stat_synthetic_R_QA.R` 生成明示的模拟文库，检验正负交互和仅长度改变的零效应；其数据绝不可纳入论文结果。当前无法运行时必须保留 `R_ENVIRONMENT_UNVALIDATED` 状态，不生成 `environment_ready.json`。

官方来源访问日：2026-10-17（Asia/Shanghai；远程 UTC 日期仍为 2026-10-16）。当期在线文档可能高于拟安装版本，因此最终还需以安装版本帮助页和实际执行核验 API。

- DESeq2 官方说明（交互设计、计数导入、Cook、区间）：https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
- tximport 官方说明（original counts + offset 与 length-scaled counts 的互斥用法）：https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
- limma 官方手册 roast/mroast/voom：https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf
- micromamba 官方安装：https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html

完整方法报告须引用实际安装版本；软件输出完成不代表获得支持投稿的生物学结论。

补充的结果前注释核查：主注释 shared 与 union 在四类条目上的成员差异实际为零，因此不重复完全相同的集合分析充当敏感性。候选主题集合另存 `modules_precommit.json` 和 `modules_candidate_audit.tsv`，缺失条目不根据表达结果替换。

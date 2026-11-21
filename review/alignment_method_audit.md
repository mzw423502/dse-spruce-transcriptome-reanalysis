# 固定 MA—来源序列对的比对方法追溯审计

## 结论

冻结记录支持把本研究的序列指标明确写为：**对固定历史 MA—来源 accession.version 映射的氨基酸序列局部比对描述值**。比对不是核酸比对，也不是重新检索数据库所得的最佳命中；固定映射不构成直系同源、同一基因或功能等价的证据。

实际比对输入为当前 `Paab.1_0b` 肽序列和由精确 accession.version 锁定的 GBSeq 来源蛋白序列。宿主肽序列在此前的 CDS—肽序列身份审计中，曾按标准核遗传密码表以正向 frame 0 翻译 CDS（模糊密码子记为 X，残余不完整密码子另行标记，至多去除一个末端终止符，不搜索替代阅读框）进行核对；但**该翻译不是本次局部比对时重新执行的步骤**，本次比对直接读取冻结肽序列。

## 已核实的方法

- 引擎：Python 3.12.3、Biopython 1.88 `Bio.Align.PairwiseAligner`、NumPy 1.24.2。
- 报告值：Gotoh 局部比对；BLOSUM62；缺口开放罚分 −11、延伸罚分 −1，即长度为 *n* 的缺口罚分为 −11 − (*n*−1)。另存有全局 Gotoh 结果，但正文六候选使用的是局部结果。
- 序列处理：肽序列大写；宿主可至多去除一个末端 `*`，内部 `*` 和 `X` 不被静默替换；来源序列为 GBSeq 蛋白序列，accession.version、原始 XML 和序列哈希均锁定。
- 一致率：`local_identity_pct_paired = 100 × 相同的非缺口残基对 / paired_residues`。因此分母不含任一侧为 gap 的列；另有仅限 20 个标准氨基酸的身份计数，但不是此处报告的一致率。
- 双侧覆盖度：宿主或来源覆盖度均为 `100 × paired_residues / 各自完整肽序列长度`。配对残基数不含 gap 列；不可用跨度覆盖度 `(end-start)/length` 代替。
- 并列最优：至多读取固定库顺序中的前 32 个最优解，并在该集合中按扁平化坐标字典序选取一个；若枚举被截断，不代表生物学对应唯一。
- 固定链：保留既有所有 MA/unit/source accession.version 行及其行哈希，以 `(当前宿主肽 SHA-256, source accession.version, 来源肽 SHA-256)` 构成去重后的序列对键；不按本轮分数删行，不产生 best hit，不进行直系同源或功能转移判定。

## 六个展示候选的冻结局部指标复核

| 来源描述简称 | MA | source accession.version | 局部一致率（%） | 宿主覆盖度（%） | 来源覆盖度（%） |
|---|---|---:|---:|---:|---:|
| ERF017-like | MA_81979g0010 | XP_016682358.1 | 44.83 | 86.57 | 87.88 |
| NHL13-like | MA_71687g0010 | XP_014500081.1 | 43.07 | 91.33 | 54.80 |
| CSLD3-like | MA_595526g0010 | XP_022718576.1 | 79.86 | 98.95 | 98.51 |
| GST-like | MA_10426104g0010 | XP_010923805.1 | 55.83 | 92.61 | 73.76 |
| XTH23-like | MA_10432086g0010 | XP_012574185.1 | 70.94 | 97.07 | 92.33 |
| CHS4-like | MA_5735g0010 | XP_027350008.1 | 86.79 | 100.00 | 85.17 |

六行均在冻结 `all27745_current_unit_source_alignment_metrics.tsv` 中为 `local_status=ALIGNED`，此表的数值未被本次审计重算。

## 可直接用于正文（约 1.5 段）

来源蛋白描述仅用于为候选建立可追溯的序列来源标签。我们保留既有的 MA—来源 accession.version 历史映射，并以当前 Paab.1_0b 肽序列与相应 GBSeq 来源蛋白进行固定对局部氨基酸比对；本步骤未进行全库检索或重新选择最佳命中。比对采用 Biopython 1.88 `PairwiseAligner` 的 Gotoh 局部算法、BLOSUM62 矩阵、缺口开放和延伸罚分分别为 −11 和 −1。报告的一致率为非缺口残基对中字符完全相同的比例，宿主和来源覆盖度均以非缺口配对残基数除以各自完整蛋白长度计算。因此，这些数值描述的是固定序列对的局部相似性和双侧覆盖范围，而非直系同源、同一基因或功能等价的判据。

为保留来源链的历史可追溯性，所有既有 MA/unit/source accession.version 行及其来源记录均被保留；不因比对得分或覆盖度删除来源，也不把任一来源称为本研究重新确定的 best hit。对于六个展示候选，CSLD3-like 的局部一致率及宿主/来源覆盖度为 79.86%/98.95%/98.51%，XTH23-like 为 70.94%/97.07%/92.33%，而 ERF017-like 和 NHL13-like 的一致率分别为 44.83% 和 43.07%（NHL13-like 来源覆盖度为 54.80%）。这些差异支持将标签作为候选优先级的来源线索，而不能替代云杉中独立的序列身份或功能验证。

## 可直接用于补充材料 S4（详细方法段）

**固定 MA—来源蛋白序列对及局部比对。** 为避免把来源蛋白描述误作重新检索得到的最佳命中，本研究不进行全数据库相似性搜索。我们从既有来源链中保留全部 MA/unit/source accession.version 行及其行哈希；每个可用对由当前宿主肽序列 SHA-256、来源 accession.version 和来源蛋白 SHA-256 固定，随后仅对这些预先存在的序列对计算指标。宿主输入为冻结的 Paab.1_0b 肽序列（大写、至多去除一个末端终止符，保留内部终止符和 X）；来源输入为对应 accession.version 的 GBSeq 蛋白序列，并保留原始 XML 与序列哈希。此前的 CDS—肽序列身份审计以标准核遗传密码表的正向 frame 0 翻译核对肽序列，其中模糊密码子转为 X、末端不完整密码子单独标记、不搜索替代阅读框；该审计用于来源追溯，不是在本次比对时重新翻译 CDS。

局部指标用 Python 3.12.3 中 Biopython 1.88 的 `Bio.Align.PairwiseAligner` 完成，采用 Gotoh 局部算法、BLOSUM62、缺口开放罚分 −11、延伸罚分 −1（长度 *n* 的缺口罚分为 −11 − (*n*−1)）。对每个 `ALIGNED` 记录，程序重建两条比对字符串并复核 BLOSUM62/仿射缺口总分。`local_identity_pct_paired` 的分子为字符完全相同且两侧均非缺口的残基对，分母为全部非缺口配对残基；`local_host_paired_residue_coverage_pct` 和 `local_source_paired_residue_coverage_pct` 分别为该配对残基数除以完整宿主或来源蛋白长度。故 gap 列不进入一致率分母或覆盖度分子；跨度覆盖度为另一独立指标，未用于替代双侧残基覆盖度。对于并列最优局部解，程序至多读取固定库顺序中的前 32 个，并在该集合内选择扁平化坐标字典序最小者；被截断的并列集合不被解释为唯一生物学对应。全局 Gotoh 指标被存档作结构性比较，但本文候选表和正文报告局部结果。整个流程不设以一致率或覆盖度赋予功能资格的阈值，也不据此产生 best hit、直系同源、宿主 GO 归属或功能等价结论。

## 冻结证据定位

| 证据 | 冻结位置 |
|---|---|
| 比对引擎、参数、指标定义 | `/home/11/dse_audit_20261016T1226Z/12_formal_RNA_20261016T1640Z/03_reference/annotation_extension/current_sequence_evidence_20261028_1212/align_current_pairs.py` |
| 完整局部/全局协议与软件版本 | `/home/11/dse_audit_20261016T1226Z/12_formal_RNA_20261016T1640Z/03_reference/annotation_extension/current_sequence_evidence_20261028_1212/all_fixed_pairs/protocol.json` |
| 完整执行回执 | `/home/11/dse_audit_20261016T1226Z/12_formal_RNA_20261016T1640Z/03_reference/annotation_extension/current_sequence_evidence_20261028_1212/all_fixed_pairs/execution_receipt.json` |
| 固定历史 MA/source 链与序列对构造规则 | `/home/11/dse_audit_20261016T1226Z/12_formal_RNA_20261016T1640Z/03_reference/annotation_extension/current_sequence_evidence_20261028_1212/freeze_assets.py`、`all_legacy_chains.tsv`、`unique_current_pairs.tsv` |
| 27745 个 unit/source 对的冻结局部和全局指标 | `/home/11/dse_audit_20261016T1226Z/12_formal_RNA_20261016T1640Z/03_reference/annotation_extension/current_sequence_evidence_20261028_1212/combined_sequence_source_evidence/all27745_current_unit_source_alignment_metrics.tsv` |
| 肽序列直接输入和 CDS→肽身份核对边界 | `/home/11/dse_audit_20261016T1226Z/12_formal_RNA_20261016T1640Z/03_reference/annotation_extension/hourly_20261017_0944/transfer_policy_review/peptide_flag_addendum.json`、`.../hourly_20261017_0643/cds_peptide_identity/` |


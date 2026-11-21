本analysis目录包含两轮补分析。本轮GO功能富集及图3复现入口见 enrichment/README_REPRODUCE.txt；下文沿用上一轮候选蛋白定向结构域分析说明。

候选蛋白 Pfam 域证据包：复现说明
====================================

一、范围和结论边界

本包是对 6 个候选宿主蛋白及其 6 条冻结来源蛋白的定向 Pfam 域证据检查，
不是全 Pfam-A 搜索、不是新的转录组重分析，也不验证具体直系同源关系或蛋白功能。
正式模型集为 11 个经官方 InterPro/Pfam 元数据核对的预选 HMM；正式目标为 12 条
固定氨基酸序列（6 条宿主当前肽和 6 条冻结来源蛋白）。模型清单、版本、哈希和
下载来源见 models/model_manifest.tsv；序列的固定来源、长度和哈希见
inputs/scanned_sequence_manifest.tsv。

二、正式分析

正式分析在银河端项目隔离目录执行，使用 Python 3.12.3、PyHMMER 0.12.3，以 1 CPU
逐个 HMM 对固定的 12 条蛋白进行 hmmsearch。正式判定使用 Pfam gathering threshold
（GA）：全序列分数和相应域分数均须达到该模型的 GA 门槛。完整域坐标、E 值、
全序列分数、两类 GA 门槛、HMM/蛋白覆盖度和 GA 通过状态见
results/domain_hit_summary_with_ga.tsv。

正式输出为 17 个域命中；其中宿主当前肽有 7 个 GA 域命中，6 个候选中有 5 个至少
带一个宿主 GA 命中。E 值仅对应这一次固定 12 条目标蛋白的 hmmsearch 搜索空间，
不能替代全 Pfam-A 搜索的 E 值。零显示的 E 值是浮点数值下溢的显示结果，不应解释为
期望随机命中数严格为零。未命中不构成“没有该域”的证据；域命中亦不构成具体直系同源或功能
已验证的证据。

PF01616（官方元数据为 Orbi_NS3）和 PF01436（NCL-1/HT2A/LIN-41 NHL repeat 命名
体系）在正式扫描前已排除，不能用作植物 XTH 或 NHL13-like 的阴性证据。其排除理由
记录在 receipts/complete_targeted_scan_receipt.json。

三、门槛外诊断

results/ga_outside_diagnostic.tsv 是为解释 ERF 宿主 AP2 及 GST 宿主 N 端近阈值情况
保留的探索性诊断。它采用宽松报告阈值（E/domE/incE/incdomE=1000）读取最高报告分数
及坐标，但不改变正式 GA 分类，不能把其中结果提升为正式阳性。原始诊断脚本
scripts_candidate_domain_diagnostic_20261110.py 保留其历史路径和运行证据；可移植的复现
入口为 rerun_diagnostic_from_frozen_inputs.py。它同样依赖项目内的 PyHMMER 环境和同包模型/
序列，若环境不可用，不应在 Windows 上尝试替代计算。

四、离线复现正式结果

本交付包提供 HMM、固定 FASTA、模型清单和离线复现脚本，但刻意不包含 PyHMMER 的
安装目录或其他大型工具缓存。必须在银河端或等价的受控 Linux 项目环境中准备项目
独立的 PyHMMER 0.12.3，且不在 Windows 上进行科学计算。

在包根目录下运行：

  python3 rerun_from_frozen_inputs.py --root <本包根目录> --output <新的输出目录>

脚本仅读取本包 inputs/candidate_host_and_frozen_source_peptides.faa、
models/model_manifest.tsv 所列的 11 个 HMM 及项目独立 PyHMMER，使用 1 CPU 与
gathering 阈值，不下载模型、不访问早期冻结参考目录。输出 rerun_domain_hits.tsv。

已在银河端实际复跑；复跑得到 17 行，且与正式增强审核表的共同核心字段
（模型 accession、序列 header、HMM/目标坐标）核对一致，QA 结论为 PASS。复跑表
SHA-256 为 c64689ce13b88292ab47589849da23c24cdf91e56aafa70b26e92ba6c7a82633。

门槛外诊断可在同一受控银河环境中运行：

  python3 rerun_diagnostic_from_frozen_inputs.py --root <本包根目录> --output <新的诊断输出目录>

该脚本仅复现 AP2、GST_N 和 GST_N_3 的 6 条固定诊断记录；其输出不得替代正式 GA 表。
已在银河端实际复跑，复跑记录与原诊断在模型、序列、GA 门槛、最高报告分数、E 值及
坐标等共同字段完全一致，QA 结论为 PASS。

五、应优先读取的状态记录

receipts/complete_targeted_scan_receipt.json 是当前正式完成记录，优先于早期
BLOCKED_NO_LOCAL_SCANNER_AND_SECURE_MODEL_RETRIEVAL_FAILED 历史收据。后者保留是为了
完整记录最初工具/旧 Pfam URL 受阻的环境事实，不表示最终扫描失败。正式完成记录
包含实际运行方法、退出状态、版本、CPU、资源摘要及模型排除理由。

本包只用于可复现的域证据审计。引用或稿件表述时，应使用“域/家族支持”“历史来源
别名”，并保留上述分析范围限制。

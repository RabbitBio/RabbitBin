# 在 RabbitBin 现有架构内继续定位问题

基线为最新 `origin/main`：`bf9e1013941bf7864501e1ac771c77db73525d6e`，
本轮结束前再次核对远端未变。独立实验工作树，默认算法未改；阶段结束时
未提交，随后按用户要求归档至[独立实验分支](../README.md)。

**结论：这轮没有找到可直接启用的通用改进。** 数值饱和不是当前三套结果
差距的主要解释；候选容量确实浪费，但修正阶段先后次序只改善连通性，
同时引入大量跨 genome 边。接下来更值得检验的是：如何用现有信息学到
可靠、各阶段一致的“同源归属度量”，而非继续叠加硬过滤或只增加候选。
这仍是工作假说，不能据此认定 joint representation 已被证明是唯一或主要病因。

## 统一设置与防泄漏

- 现有 BAM 输入，不组装、不重比对；LPA 消融复用同一 BAM 生成的 cache。
- contig >=1000 bp；其中 >=2500 bp 建图，1000–2499 bp 后续招募。
- bin >=200000 bp，RabbitBin `amber`；seed=42，64 线程，PMH m=500，
  mutual top-200，dual-depth=5，既有 cutoff=0.7153318629591614。
- 不新增生物学阈值、dataset 分支或调参搜索；三套数据执行完全相同规则。
  Gold 只给 AMBER 和本目录离线统计，不进入候选筛选或 LPA。
- HQ/MQ 沿用当前 RabbitBin AMBER 定义：HQ completeness>90%、purity>95%；
  MQ completeness>70%、purity>90%，MQ 包含满足条件的 HQ，不是互斥计数。

预先固定的实验设计见 [DESIGN.md](DESIGN.md)，完整指标见
[summary.tsv](results/summary.tsv)，边统计见 [graphs.tsv](results/graphs.tsv)。

## 1. LPA 数值问题：存在，但本次解释力很小

保持原始图及边权，比较原始 Fisher CDF、已有 upper-tail、稳定 log-tail 和
已有 log-sum。log-tail 实现的是 `-log Q(m, -sum log(1-w))`，在最大 Poisson
项附近缩放求和，避免先算 exp(-u) 的下溢和 1-Q 的抵消。
其平局规则采用 log 域一 ULP；因此保留数学排序，但不声称与旧 `1e-12`
CDF 绝对容差在决策上等价。它不使原始边权自动成为校准 p-value，也不
修复 Fisher 对相邻边独立性的假设。

| 原始图决策诊断 | Marine | Plant-associated | Strain-madness |
|---|---:|---:|---:|
| 活跃节点决策次数 | 305854 | 62608 | 35704 |
| 最佳 CDF 已等于 1 | 14520 | 5771 | 2428 |
| 两个以上候选标签同时等于 1 | 67 | 0 | 0 |
| 同一轨迹上 log-tail 会改变的决策次数 | 36 | 0 | 0 |
| 涉及独立 contig 数 | 15 | 0 | 0 |
| 重复标签 guard 冻结的节点数 | 112 | 42 | 9 |

影子诊断使用**同一运行轨迹上的相同邻居标签状态**，不改变生产决策；
三个影子基线的 members 与原始基线逐字节一致。冻结计数不是冻结影响的
因果消融，不能仅凭数量少就断言其完全无害。

| 原图聚合 | Marine HQ/MQ/加权纯度 | Plant HQ/MQ/加权纯度 | Strain HQ/MQ/加权纯度 |
|---|---|---|---|
| Fisher 基线 | 286 / 353 / .9264 | 84 / 94 / .9797 | 33 / 38 / .4038 |
| upper-tail | 286 / 353 / .9264 | 84 / 94 / .9797 | 33 / 38 / .4038 |
| stable log-tail | 286 / 353 / .9264 | 84 / 94 / .9797 | 33 / 38 / .4038 |
| log-sum | 291 / 358 / .9237 | 81 / 92 / .9571 | 35 / 41 / .3840 |

Stable log-tail 的 Marine members 有变化，但上述指标未变；Plant、Strain
members 逐字节不变。不能把“CDF 会饱和”直接推成当前主要精度瓶颈。

## 2. 将既有完整 coverage 准入移到 PMH top-200 截断之前

原来先以 Spearman 做保守初筛，再按 PMH 排 top-200、取 mutual，最后检查
`min(Spearman, Jcov) >= cutoff`。本实验要求最后这个条件先通过，才允许
候选占据 PMH heap 名额。PMH 仍负责排序，mutual 和最终边权计算均不变。
实现上仍利用 PMH 的 heap 上界剪枝，再对有竞争力的 pair 做完整 coverage
检查；这里“前移”指先于有界 heap 的占位，不要求先于每次 PMH 计算。

这**改变了邻居集合**，不是一个与旧图等价的性能优化。开关为
`RABBIT_CANDIDATE_COVERAGE=1`，仅允许 fresh、>=3 samples、>25000 长 contig 的默认
mutual PMH + fused-coverage 路径，拒绝旧 cache/外部图/LSH/GFA/SNV 和
candidate-stage gold 审计组合，避免无效或误标的实验。

| 数据集 | HQ 基线→前移 | MQ 基线→前移 | 加权纯度基线→前移 |
|---|---|---|---|
| Marine | 286 → 272 | 353 → 340 | .9264 → .8823 |
| Plant-associated | 84 → 84 | 94 → 94 | .9797 → .9795 |
| Strain-madness | 33 → 32 | 38 → 47 | .4038 → .3924 |

离线集合检查确认：**三套新图都包含全部旧有效边，且这些边的权重逐项相同**。
因此这是一个很干净的“保留旧边并补充候选”的干预，但补边不等于提升 bin 质量。

| 数据集 | 有同-genome 邻居的长 contig 碱基比例 | 新增同-genome 边 | 新增跨 genome 边 |
|---|---|---:|---:|
| Marine | 92.37% → 97.12% | 2084141 | 325419 |
| Plant-associated | 80.50% → 81.32% | 211508 | 8008 |
| Strain-madness | 50.09% → 93.19% | 105995 | 360410 |

该比例按碱基计，不是边召回率；分母限 gold 可标注且该 genome 至少有两条
建图长 contig 的碱基。Strain 新增边中 77.3% 跨 genome；补充正确关系的
同时也放大了错误关系。其平均 genome completeness 从 .3690 升到 .7759，
但 HQ/纯度没有相应提升，也说明仅报 completeness 或 MQ 容易掩盖混合 bin。

### 新增边是否还能靠现有单项分数简单排序？

在 Strain **新增的已选边**中：

| 分数 | 同 genome 中位数 | 跨 genome 中位数 | 高分预测同 genome 的 AUC |
|---|---:|---:|---:|
| PMH | .6875 | .7708 | .3493 |
| coverage weight | .7873 | .7952 | .4759 |

这是受候选选择影响的条件统计，不是全体 pair 的 AUC，也不证明任何联合
函数都不可能区分它们。但它说明：单纯“高 PMH 加高 coverage”不是充分的
正确归属证据；把边扩大后全部交给 coverage 标量，风险很明确。

## 3. 扩大图 × LPA 聚合的交叉消融

复用上述新图 cache，先确认 Fisher replay members 与从 BAM 完整运行逐字节
一致，再仅替换聚合规则，排除输入差异。

| 新图聚合 | Marine HQ/MQ/加权纯度 | Plant HQ/MQ/加权纯度 | Strain HQ/MQ/加权纯度 |
|---|---|---|---|
| Fisher | 272 / 340 / .8823 | 84 / 94 / .9795 | 32 / 47 / .3924 |
| stable log-tail | 272 / 340 / .8823 | 85 / 94 / .9796 | 32 / 47 / .3923 |
| log-sum | 273 / 343 / .8674 | 80 / 90 / .9533 | 43 / 66 / .3702 |
| raw sum | 269 / 339 / .8609 | 80 / 91 / .9524 | 42 / 64 / .3703 |

可以把 Strain HQ 拉高，但代价是纯度和另外两套数据。不能挑出 Strain 的
43 HQ 宣称成功，也不能按数据集选择不同聚合规则。

## 性能与验证

这轮 full-BAM 单次墙钟（包含保存 cache 和导出图，AMBER 另计）：

| 数据集 | 基线/前移，秒 | 基线/前移峰值，KiB | pair-pass 基线/前移，秒 |
|---|---|---|---|
| Marine | 14.58 / 10.67 | 1499944 / 1478004 | 2.472 / 2.778 |
| Plant | 13.44 / 12.24 | 2767124 / 2694376 | .212 / .213 |
| Strain | 75.53 / 76.84 | 591100 / 637088 | 2.000 / 3.924 |

**不是加速证明**：BAM 页缓存/读盘状态不同，Marine/Plant 总墙钟下降主要
伴随 depth 耗时变化；Strain 候选 pair-pass 接近翻倍。质量已不满足通用
改进要求，因此未将它作为性能优化推广，也未对失败配置做参数调优。

通过的检查：

- 3400 个 log-tail 测例，与 long-double `boost::math::gamma_q` 对比，
  最大绝对误差约 9.1e-13；包含 CDF 饱和和 exp(-u) 下溢区间。
- 合成候选测试：后排正确幅度邻居被找回，旧有效边及权重完整保留，
  1/4 线程图与 members 一致。
- 三套 full-BAM baseline、影子诊断 baseline 均与先前基线 members 相同；
  三套新图 cache replay 与 fresh 输出相同。
- `test/test_coverage_edges.py`、`test/test_abundance_filter_graph.py` 通过。
- `git diff --check` 通过。没有将未注册测试的 ctest 当作验证依据。

## 下一步：仍基于这套架构，而不是再换一套 binner

本轮支持把下一步具体化为两个可证伪方向，并不提前保证有效：

1. **学习现有 PMH/composition 与 coverage 的联合可靠性。** 保留 BAM 输入、
   紧凑特征和 bounded graph；用同一长 contig 的独立片段构造自监督正对，
   从真实分段 BAM 计数获得片段 coverage，不直接复制父 contig 的 coverage
   冒充独立样本。用随机跨 contig 对作受污染的负对时必须显式处理“其实
   同 genome”的可能。学习结果同时用于候选排序与边评分，不能前一阶段
   学联合关系、后一阶段又丢掉 composition。长度用于表达估计可靠性，
   不手工调 dataset 权重；gold 只作最后盲评。
2. **在统一度量上做 bin 级一致性检验。** 现在 connected-but-fragmented
   的 contig 常已获得当前 fragment 的局部支持；仅增加邻居或换个 LPA
   加和式不够。下一步可检查多个 fragment 的联合分布是否相容，再决定
   是否值得实现合并/再分配，不能直接把互相连接当成可合并证据。

与上一轮线性 [joint representation 实验](../joint_representation/README.md)
结合，优先级应从“特征拼接/硬门槛堆叠”转到“可靠性学习 + 全流程使用同一
度量”。本轮默认路线和用户其他工作树均未修改。

## 复现

```bash
# 当前图 LPA；诊断只输出计数，不影响选择
bash experiments/architecture_local/run_lpa.sh
SCORES=logtail bash experiments/architecture_local/run_lpa.sh
OUTROOT=experiments/architecture_local/results_shadow SCORES=fisher RB_LPA_AUDIT=1 \
  bash experiments/architecture_local/run_lpa.sh

# BAM 起点对比
bash experiments/architecture_local/run_candidate.sh

# 扩大图上复查聚合
CACHE_ROOT=experiments/architecture_local/results CACHE_NAME=coverage_first.r1 \
LABEL_PREFIX=coverage_lpa SCORES='fisher logtail logsum sum' \
  bash experiments/architecture_local/run_lpa.sh

# 离线统计；此脚本才读取 gold
experiments/candidate_stage/venv/bin/python experiments/architecture_local/summarize.py

c++ -std=c++17 -O2 experiments/architecture_local/test_fisher_log.cpp -o build/test_fisher_log
build/test_fisher_log
python3 experiments/architecture_local/test_candidate.py build/src/rabbitbin build/test-candidate-coverage
```

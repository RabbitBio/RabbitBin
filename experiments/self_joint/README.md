# 自监督联合度量：阶段性总结（已暂停测试）

按用户要求，已完成已启动的基础度量这一轮，停止后续测试。
**没有启动 signal 系列的 CAMI2 实验，没有修改默认分箱路线。**
阶段结束时尚未提交；随后按用户要求归档到独立实验分支，见
[实验入口](../README.md)。代码归档不代表恢复测试或启用默认设置。
基于 origin/main `bf9e1013941bf7864501e1ac771c77db73525d6e` 的独立工作树。

## 当前结论

1. 从 BAM 同一次扫描提取两半实际 coverage 的基础设施已完成并验证。
2. 自监督线性度量可以明显改善 Strain，但 Marine/Plant 下降，不能通用替换。
3. 完整 joint 与去掉跨特征块项的 block 结果接近；Strain 上 coverage-only
   更好。**目前没有证据把增益归因于跨 composition/coverage 的联合学习。**
4. 这不是对所有 joint representation 的否定：本阶段只验证了正对驱动的
   线性 Mahalanobis 度量，没有验证非线性表征、可靠性分层或负向约束。
5. 精确全局邻居检索仍有明显开销，尚不满足直接替代 RabbitBin 的质量/速度要求。

全部结果在 [all_results.tsv](results/all_results.tsv)：3 种表示 × 8 个阶段/
聚合配置 × 3 个数据集 = 72 个实验配置；加上 6 个基线 replay 对照、
3 条原始基线，共 81 行。没有任何非恒等配置能同时保持三套 HQ 和
加权纯度均不低于原始基线。

## 做了什么，没有做什么

保持 BAM 输入，未组装、未重比对、未添加 taxonomy/marker 数据库。
contig >=1000 bp；>=2500 bp 建图，其余后续招募；bin >=200 kb 用 RabbitBin
AMBER，seed42、64 threads、mutual top-200、dual-depth5。
所有配置在三套数据上统一执行，gold 只用于 AMBER 和结束后的边统计。

借鉴 [SemiBin 的同-contig 切分正对思路](https://www.nature.com/articles/s41467-022-29843-y)，
但本轮没有复制其神经网络或负对/marker 方法：

- >=5000 bp 的 contig 等分为两段；5000=2×现有建图长度下限，不另调阈值。
- 两半分别统计 canonical 4-mer（现有 PMH 所用的 composition 信息），
  使用已有 per-base GC 背景校正，不是添加一类外部特征。
- 两半 coverage 从真实 BAM CIGAR 对各自区间累计，保留 all-read 与 MAPQ>=5
  两通道、现有 read identity 过滤、各自边缘裁剪；绝不复制父 contig coverage。
- 对两半特征差学习二阶矩 W，再用 W^-1 定义 Mahalanobis 距离。
  [OAS](https://scikit-learn.org/1.8/modules/covariance.html) 的收缩量由数据估计，
  没有手调 shrinkage、embedding 维数或数据集权重。
- joint 使用完整 W；block 使用同一 W 的对角块，去掉跨 composition/coverage
  项；coverage 只保留其 coverage 子块。不是分别调三套数据的最优参数。
- 先在旧 PMH 图上只换边权，再检查同一 learned embedding 同时选择邻居和
  赋边权的结果。后者是质量验证原型，不是已集成的生产 PMH 检索替代品。

训练正对：Marine **51240**、Plant **29001**、Strain **22238**。
按父 contig 做二折检查，未把同一父 contig 的两半拆到训练/验证两侧。
随机其他父片段只用于 proxy 排序诊断，不作为训练的负标签；它可能来自
同 genome，因此 proxy 不是 genome 分类准确率。前处理使用全数据的无标签
列均值和尺度，是 transductive 设置；交叉检查隔离的是 W 拟合。

两个不重叠区域不代表 reads/mates 完全统计独立，assembly chimera 也可能
污染正对；现有 W 是共享噪声模型，尚未解决长度/深度异方差。这些是限制，
不是已验证的后续改进收益。

## 主结果：同一表示用于邻居与边权，保留 Fisher-LPA

每格为 **HQ / MQ / 加权纯度**。MQ 使用当前 RabbitBin 定义，包含 HQ，
不是与 HQ 互斥的计数。

| 模式 | Marine | Plant-associated | Strain-madness |
|---|---|---|---|
| 原始 RabbitBin | 286 / 353 / .9264 | 84 / 94 / .9797 | 33 / 38 / .4038 |
| 自监督 joint | 160 / 191 / .7883 | 74 / 81 / .9533 | 67 / 106 / .5796 |
| 自监督 block | 159 / 191 / .7983 | 74 / 84 / .9546 | 66 / 108 / .5801 |
| 自监督 coverage-only | 157 / 178 / .8107 | 67 / 77 / .9190 | 75 / 112 / .6208 |

Strain 中 joint 的 HQ 从 33 到 67，纯度从 .4038 到 .5796，这是真实改善；
但是 block 几乎相同，coverage-only 还更好，不能据此说“跨特征联合表示
已解决主问题”。Marine 的代价很大，不能按数据集挑模式隐藏退化。

Strain 的同-genome 邻居碱基比例从 50.09% 到 joint 的 94.16%，同时跨
genome 边从 132270 到 477578；Marine 的该比例只从 92.37% 到 94.21%，
跨 genome 边却从 115880 到 411452。连通率增加仍不等于分箱精度增加。

## 保留现有架构的阶段对照

不仅测试了整个 graph 替换，也保留了原 PMH 图/原 coverage 打分的单阶段
对照，并检查 Fisher、sum、原 coverage 准入及既有 neutral-point 映射。
这些复用上一轮同一组控制规则，不按结果扫参数。

例如，旧 PMH 图上 learned joint weight + sum-LPA 在 Marine 达到
292 HQ（原286），但纯度略降；Plant 只有81 HQ（原84）。所有配置合起来
仍没有一个能在三套同时维持 HQ 与加权纯度，所以目前全部保持 opt-in。

## 性能边界

以下是研究原型的分段计时，不是一次集成生产 pipeline 的端到端 benchmark：

| 数据集 | BAM分箱并导出两半 coverage | 两半 composition 导出 | 特征装载/变换 | W拟合+二折检查 | joint精确邻居搜索 |
|---|---:|---:|---:|---:|---:|
| Marine | 15.94 s | 2.90 s | 11.47 s | 2.55 s | 25.33 s |
| Plant | 30.37 s | 2.94 s | 4.35 s | 1.52 s | 9.37 s |
| Strain | 62.50 s | 2.05 s | 8.07 s | 6.22 s | 3.32 s |

学习矩阵不是最大瓶颈；Marine 的 exact top-200 搜索本身已高于原始 RabbitBin
此前整轮约14–15秒的观测量级。特征装载也有 Python/CSV 原型开销。
不能声称当前版本“精度提升且性能接近”。BAM时间受缓存和机器其他I/O影响，
未在暂停前继续做成对重复性能测试。训练过程的峰值（依次运行三个模式，
不含 BAM）：约 1.93 / 1.28 / 0.96 GiB。

## 已完成的验证与文件完整性

- 合成 BAM 独立区域计数 oracle：M/= /X、I/D/N、MAPQ 双通道、偶数/奇数
  长度、过滤长度、1/4 线程、不同 shard 数均逐 float32 相同。
- 强制 B2 重同步失败后回退，确认丢弃部分 fragment sums，不会重复计数。
- 启用新导出时，合成 BAM 的普通 depth/members 与关闭导出一致。
- 三套真实 CAMI2 导出后的 baseline members 与既有基线逐字节一致。
- learned metric 的白化恒等式、Mahalanobis 距离等价、block/coverage 对照
  和 held-parent proxy 的代数单元检查通过。
- 原图导入 identity 的三套 members 与基线一致；每套 26 份 AMBER 完整，
  三套共78份，均已完成图文件解析及离线金标准统计。
- `bash -n` 和 `git diff --check` 通过。

执行记录说明：为加入后续 signal 版本，曾在当前 driver 尚未完全退出时
更新脚本。三个数据集及所有输出完成后，该 driver 尾部出现 shell EOF
解析报错；不是某个 binning/AMBER 失败。已逐项核对78份结果，当前磁盘
脚本语法检查通过，未因该尾部错误重新启动任何测试。

## 暂停点

- 本轮的训练、分箱、AMBER 已结束，没有继续运行测试。
- 信噪比/likelihood-ratio 度量仅完成公式、代码和合成代数检查，**未运行
  signal_joint / signal_block / signal_coverage 的 CAMI2 实验**。
- 没有开始 bin-level 合并试验，也没有修改默认算法。
- [DESIGN.md](DESIGN.md) 保留下一阶段定义；等待用户明确继续后再运行。

## 文件与复现入口（暂停期间不会自动执行）

- [extract.sh](extract.sh)：BAM 同扫导出，两半 composition；检查原始成员不变。
- [learn.py](learn.py)：无 gold 的度量学习、父-contig 二折检查及 graph 输出。
- [run.sh](run.sh)：统一三个模型、三个数据集及阶段/聚合对照。
- [report.py](report.py)：唯一读取 gold 的离线汇总入口之一，另一个是 AMBER。
- [test_fragment_depth.py](test_fragment_depth.py)、[test_metric.py](test_metric.py)：
  独立计数与代数检查。依赖复用上一轮
  [requirements.txt](../joint_representation/requirements.txt)。

# RabbitBin 联合表示实验

基于 GitHub `origin/main` 的 `bf9e1013941bf7864501e1ac771c77db73525d6e`，
本轮结束时再次核对远端，提交未变。实验在独立工作树进行；阶段结束时
未提交，随后按用户要求归档至[独立实验分支](../README.md)，默认设置不变。

这轮结果支持继续检查阶段之间的信息损失，但**尚不支持把缺少 joint
representation 定为三套数据共同的主要原因**。轻量联合表示大幅改善了
Strain，而该数据上 coverage-only 更强；Marine、Plant 的结果没有稳定改善。
全部配置中，没有一个同时保持三套数据的 HQ 和加权纯度不下降。

## 设置与范围

- 复用现有 BAM 计算得到的同一份 RabbitBin v4 cache，避免反复读取 BAM
  和改变输入。没有重新比对或组装，也没有引入 marker、分类数据库或金标准特征。
- 保持 contig >=1,000 bp、初始图 contig >=2,500 bp、bin >=200,000 bp，
  seed 42、64 线程、mutual top-200，使用 RabbitBin `amber`。
- 训练、归一化、邻居搜索和 binning 不读取 CAMI gold。Gold 只用于 AMBER
  及完成后的边正确性统计。
- 共 9 种表示、每种 8 个阶段/聚合配置、3 套数据，即 216 个实验配置；
  加上原始基线、原图重导入和基线 sum-LPA 对照，汇总表共 225 行。
- 这些是线性表示与图流程的机制实验，没有训练 SemiBin/COMEBin 式神经网络。

完整协议见 [DESIGN.md](DESIGN.md)，全部指标见
[all_results.tsv](results/all_results.tsv)。失败配置也保留。

## 联合表示

composition 使用 136 维 canonical 4-mer；coverage 使用现有 all-read 和
MAPQ>=5 两个通道。coverage 按列归一化后取 log1p，composition 使用平方根
频率；两个块均中心化并按总方差归一化，避免维数直接决定相对权重。

1. `concat`：归一化后拼接，使用共同空间中的欧氏距离。
2. `whiten`：对拼接向量估计完整的 OAS 收缩协方差，使用其逆平方根变换，
   因此跨模态协方差进入距离。
3. `cca`：同一 contig 的 composition/coverage 是配对视图，通过正则 CCA
   学习对齐坐标，保留所有数值上非零方向，按 canonical correlation 连续加权。

另对这三种方法加入 RabbitBin 原有的碱基背景校正，得到 `gc_*`；并加入
原始 composition、校正后 composition、coverage 三个单模态对照。
没有按数据集选择权重、维数或邻居数。

联合距离用于精确 top-200 邻居搜索。边权为局部尺度 Gaussian affinity：
`exp(-d² / sqrt(s_i*s_j))`，其中 `s_i` 为第 200 邻居的平方距离。
这是局部密度归一化，没有另外搜索带宽。

OAS 实现使用 [scikit-learn](https://scikit-learn.org/stable/modules/generated/sklearn.covariance.OAS.html)，
精确邻居搜索使用 [FAISS IndexFlatL2](https://github.com/facebookresearch/faiss/wiki/Faiss-indexes)。
依赖版本在 [requirements.txt](requirements.txt) 中固定。

## 主结果：同一距离控制候选和边权

下表统一保留 RabbitBin 默认 Fisher-LPA、split 和 recruit，仅替换图表示。
每格为 HQ / MQ，不是按数据集选择各方法的最佳聚合配置。

| 表示 | Marine | Plant-associated | Strain-madness |
| --- | ---: | ---: | ---: |
| RabbitBin 原始基线 | 286 / 353 | 84 / 94 | 33 / 38 |
| concat | 159 / 187 | 71 / 80 | 52 / 76 |
| whiten | 157 / 190 | 76 / 85 | 45 / 73 |
| cca | 169 / 202 | 70 / 82 | 34 / 54 |
| gc_concat | 163 / 198 | 76 / 83 | 51 / 72 |
| gc_whiten | 166 / 197 | 72 / 84 | 47 / 72 |
| gc_cca | 159 / 190 | 77 / 85 | 33 / 55 |
| composition-only | 137 / 172 | 34 / 60 | 6 / 9 |
| GC composition-only | 145 / 177 | 38 / 57 | 11 / 12 |
| coverage-only | 164 / 181 | 71 / 79 | 62 / 95 |

Strain 的加权纯度：基线 **0.4038**，concat **0.5786**，coverage-only
**0.6488**。这说明已有 coverage 信息确实能支持更好的分箱，但等权联合表示
不是这一提升的必要条件。换成已有 GC 校正也没有消除 Marine 的退步。

## 阶段、过滤和聚合对照

每种表示分别测试新候选+原 coverage 评分、原有存活边+联合权重、完整联合图。
又测试标准 sum-LPA、保留原 coverage 拒绝条件，以及将 affinity 映射到
Fisher 原有中性点以上。后者只用现有解析常数，不按 gold 拟合。

| concat 对照 | Marine HQ/MQ | Plant HQ/MQ | Strain HQ/MQ |
| --- | ---: | ---: | ---: |
| 只换候选，使用原 coverage 权重/过滤 | 275 / 346 | 85 / 95 | 35 / 53 |
| 保持原存活边，换联合权重，Fisher-LPA | 243 / 272 | 74 / 80 | 33 / 41 |
| 保持原存活边，换联合权重，sum-LPA | 291 / 358 | 81 / 91 | 35 / 42 |
| 同时使用联合候选与权重，Fisher-LPA | 159 / 187 | 71 / 80 | 52 / 76 |

因此，不能只把 Gaussian affinity 填入 Fisher 聚合就宣称整个流程已校准；
但改成 sum-LPA 也没有得到三套共同改善。相同旧图和相同旧权重使用 sum-LPA，
本身就得到 293/81/35 HQ，说明部分小幅收益来自聚合变化。

Strain 上 Fisher 中性点映射可让 concat 达到 80 HQ、coverage-only 达到
89 HQ，但加权纯度分别为 0.4382、0.5307，低于其未映射版本。
这类 HQ/纯度取舍必须同时报告，不能只展示最高 HQ。

唯一在三套上 HQ 都不低于基线的非恒等配置，是原图+GC composition 权重+
sum-LPA：287/85/35 HQ，但纯度 0.9229/0.9781/0.3924 全部低于基线，
且它并不是联合表示的证据，因此没有提升为默认。

## 图层面的解释

以下比例为有至少一条存活同-genome 边的长 contig 碱基比例；分母只包含至少
有两条长 contig 的已标注 genome。它衡量邻居覆盖，不等于最终 bin 召回率。

| 数据/图 | 同 genome 边 | 跨 genome 边 | 有正确邻居的长 contig bp |
| --- | ---: | ---: | ---: |
| Marine 基线 | 2,729,504 | 115,880 | 92.37% |
| Marine concat | 5,877,764 | 553,909 | 97.55% |
| Strain 基线 | 56,206 | 132,270 | 50.09% |
| Strain concat | 279,871 | 522,525 | 96.31% |
| Strain coverage-only | 418,842 | 501,471 | 97.29% |

Strain 的候选关系覆盖得到很大改善；Marine 同时引入了约 4.8 倍的错误边。
这是统一空间在两套数据表现不同的直接图层证据，但仍不是对最终 HQ 变化的
单因素因果分解。表示、过滤和后续聚合存在相互作用。

## 性能与验证

单个表示变换通常低于 0.5 秒。主要成本是精确全量邻居搜索，原始三种联合
表示的搜索耗时为：Marine 18.76–26.91 秒，Plant 8.69–9.83 秒，
Strain 1.00–2.22 秒。此外还有读取/归一化、4-mer 提取、图输出和 RabbitBin
尾部处理。不能把这些缓存实验的尾部耗时当作从 BAM 开始的完整耗时。
本轮没有候选同时通过三数据集质量条件，故没有进一步集成为默认 BAM 路线，
也没有宣称已经满足端到端性能要求。

- C++ 编译通过；Python 和 shell 脚本语法检查通过。
- 三套原图重导入后的最终 members TSV 与原基线逐字一致。
- 现有 `test_coverage_edges.py` 和 `test_abundance_filter_graph.py` 通过。
- 新增碱基计数输出前后的三套原始 4-mer 文件 SHA256 全部一致。
- 默认 binning 路径没有更改；新增 `--external-graph` 和
  `--external-graph-coverage` 仅供显式实验。

导入器要求无重复的无向边按 `(i,j)` 排序、`i<j`，节点名称/长度/顺序必须匹配
cache；本目录的图生成脚本已满足这些约束。原图恒等对照仅重排边记录后导入。

## 对下一步的判断

本轮支持“已有信息在候选/评分流程中未被充分利用”，尤其是 Strain 中的
coverage。它不支持“把两类向量拼接或统一为一个线性距离就能普遍解决问题”。
CCA 强调两视图共有结构，可能丢失 coverage 独有的 strain 判别方向；等权
拼接则无法表达两个模态在不同邻域中的可靠性。这两点是待验证的机制解释，
不是已完成的因果证明。

后续更值得测试的是：用真实同-contig 片段约束学习联合度量及其不确定性，
保留各模态独有的判别方向，并让候选、边权和招募共享该表示。这样的训练应
有独立片段验证和明确的成本预算；不应按 CAMI 数据集名称选择算法或权重。

## 复现

需要先编译本工作树的 RabbitBin，并在脚本所用的 Python 环境中安装
`requirements.txt` 中的依赖。

```bash
mkdir -p experiments/joint_representation/build
g++ -O3 -fopenmp -std=c++17 experiments/joint_representation/k4_features.cpp \
  -o experiments/joint_representation/build/k4_features

bash experiments/joint_representation/run.sh
bash experiments/joint_representation/run_controls.sh

METHODS='comp coverage gc_concat gc_whiten gc_cca gc_comp' \
  RUN_LABEL=representation_controls bash experiments/joint_representation/run.sh
METHODS='comp coverage gc_concat gc_whiten gc_cca gc_comp' \
  bash experiments/joint_representation/run_controls.sh

OPENBLAS_NUM_THREADS=1 experiments/candidate_stage/venv/bin/python \
  experiments/joint_representation/summarize.py \
  --output experiments/joint_representation/results/all_results.tsv
```

默认脚本路径指向本工作区已有 BAM cache 和 CAMI2 assembly；可通过脚本中的
PY/RB/BASE_ROOT/ASSEMBLY_ROOT/GOLD_ROOT/OUTROOT 等环境变量指定相同格式的输入输出。
生成的图、日志、AMBER 报告及计时文件保存在被 git 忽略的 `results/` 中。

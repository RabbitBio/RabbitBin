# 迁移后续跑说明（当前暂停）

代码分支：`experiments/joint-representation-20260925`。
本文件只提供命令，不会自动启动实验。先读 `self_joint/README.md` 和
`self_joint/DESIGN.md`；用户明确要求继续后才运行实验命令。

## 1. 构建与依赖

在这个分支的仓库根目录执行。不要复制旧机器的 build 或 venv。

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target rabbitbin -j 16
g++ -O3 -fopenmp -std=c++17 experiments/joint_representation/k4_features.cpp \
  -o build/k4_features
mkdir -p experiments/joint_representation/build
cp build/k4_features experiments/joint_representation/build/k4_features

python3 -m venv experiments/candidate_stage/venv
experiments/candidate_stage/venv/bin/python -m pip install \
  -r experiments/joint_representation/requirements.txt
```

原始环境为 Python 3.13；固定版本见 requirements。若系统 Python 太老，
使用已有 conda 的兼容 Python 环境或另建环境，不要静默改依赖版本。
C++17/OpenMP、CMake >=3.16、Boost >=1.66 的 program_options/filesystem/system/
graph/serialization/iostreams/regex 是原工程要求。HTSlib、libdeflate、zlib
可用现有系统开发包；原 CMake 也有部分依赖的下载回退。迁移没有运行构建
或实验；目标机的实际依赖可用性需要在恢复时确认。

## 2. 恢复中间状态

传输包中的 `state/experiments/` 保持相对于仓库的原目录布局。Git clone 后
把它复制到仓库 `experiments/` 下（迁移时已执行则不必重复）。其中包含：

- candidate_stage：三套 baseline.cache、原 retained 图及节点表、原 members；
- joint_representation：整条 contig 的 features.f32/base.f32/nodes.tsv；
- self_joint：真实 BAM 两半 coverage、两半 composition、训练矩阵/元数据；
- 历轮汇总表、AMBER、计时、日志和成员表（不含重复的实验大图）。

缓存是当前固定输入的派生状态，不是任意其他 assembly 的通用模型。cache
版本为 v4，float32/64-bit size_t，小端；不得与另一套 contig 名称、长度、
样本顺序或 BAM 过滤条件混用。不要将它作为“从 BAM 起算”的性能耗时。

没有传 BAM/FASTQ/assembly/gold 数据集。沿用默认路径：

- assembly: `/home/bigssd/zt/cami2/DATASET/CAMI2_DATASET_GoldStandardAssembly.fasta`
- BAM list: `/home/bigssd/zt/CAMI2_remap_bams/DATASET/bam.list`
- AMBER gold: `/home/bigssd/zt/runs/cami2_benchmark/DATASET/prep/gold_len.binning`

DATASET 为 `marine`、`plant_associated`、`strain_madness`。恢复前核对 BAM list
内部路径及样本次序。使用缓存做信噪比下一轮时，不需要重新读取 BAM/assembly；
但 AMBER 仍需要匹配原输入的 gold。路径不同可设置 `ASSEMBLY_ROOT`、`BAM_ROOT`、
`GOLD_ROOT`，不必修改算法。另支持 `BASE_ROOT`、`FEATURE_ROOT`、`OUTROOT`、
`PY`、`RB`、`K4`、`THREADS`。

## 3. 用户批准后才运行的下一轮

建议先复制本阶段目录到一个新的输出目录，以保留暂停时的全部记录。若目标
目录已存在，先确认是否是另一轮结果，不要覆盖。

```bash
test ! -e experiments/self_joint/results_signal
cp -a experiments/self_joint/results experiments/self_joint/results_signal

OUTROOT="$PWD/experiments/self_joint/results_signal" \
METHODS='signal_joint signal_block signal_coverage' \
RUN_LABEL=signal_learning TRAINING_LABEL=signal \
  bash experiments/self_joint/run.sh marine plant_associated strain_madness
```

所有三套、三种结构和相同阶段控制一起报告，不能根据 gold 选择一个数据集
专用的获胜模式。`run.sh` 会执行学习、分箱及 AMBER，绝不是仅准备文件。

若要重新提取两半特征（已有缓存完整时无需执行）：

```bash
bash experiments/self_joint/extract.sh marine plant_associated strain_madness
```

若原始 baseline.cache 缺失，先按 candidate_stage/run_clustering.sh 中
baseline 的 `--save-cache` / `--export-retained-graph` 命令重建；该脚本的
后半还会运行 Infomap/Leiden，若只需缓存，不要直接运行整个脚本。

## 4. 汇总与验证边界

传输时未携带数 GB 的重复 learned/scored graph，已有 all_results.tsv
已经保留其离线边统计。因此旧轮的 report.py 不能在缺图时原样重算。
完整重跑某轮后、该轮图齐全时可调用 report.py。新的输出目录若复制了旧
AMBER 但没复制旧图，应只对新方法重新汇总，或先恢复旧图。例如 signal
轮完成后，用迁移时补充的纯汇总过滤参数（不影响训练或选择模型）：

```bash
OUTROOT="$PWD/experiments/self_joint/results_signal" \
  experiments/candidate_stage/venv/bin/python experiments/self_joint/report.py \
  --method-prefix sj_signal_ \
  --output experiments/self_joint/results_signal/all_results_signal.tsv
```

不要把缺图解释成算法失败，也不要覆盖冻结的旧汇总表。

可在恢复工作后执行的验证入口（本次迁移不会执行）：

```bash
python3 experiments/self_joint/test_fragment_depth.py build/src/rabbitbin build/test-fragment-depth
OPENBLAS_NUM_THREADS=1 experiments/candidate_stage/venv/bin/python experiments/self_joint/test_metric.py
```

当前结果并未证明 joint representation 是主要病因。joint 与 block 接近、
Strain coverage-only 更好；Marine/Plant 退化，exact kNN 也有明显成本。
signal 系列尚无 CAMI 结果；bin-level 合并未开始。下一轮也必须按这几个
反证对照解释结果，并完整报告精度、纯度、端到端耗时和内存。

# 自监督联合度量：预先固定的试验

基线 origin/main bf9e101；三套 CAMI2 同一协议，1000 bp 输入，2500 bp
长 contig 图，200 kb bin，seed42，64 threads，top-200，dual-depth5。
不读 gold 训练/选择模型，不添加参考数据库，不按数据集选择参数。

## 数据与学习

每条长度 >=2×现有 graph minContig 的 contig 等分为两个不重叠片段。
composition 是现有 PMH 所用 canonical 4-mer 的完整频率（136 维），
使用已有 per-base GC 背景校正；coverage 使用同一次 BAM 扫描中对两半
分别做 CIGAR 区域计数的 all-read 与 MAPQ>=5 通道，不复制父 contig 覆盖度。
两半分别执行现有 edge trim。区域不重叠不意味着 reads/mates 在统计上
完全独立，且 assembly chimera 会污染正对；不隐瞒这些假设。

沿用前轮变换：composition 平方根频率，coverage 列均值归一化后 log1p；
两个特征块按全体长 contig 的 RMS 尺度归一化。无需随机负标签：
用两半特征差 `d = z_left-z_right` 学习 within-contig 二阶矩 `W`。
差分符号任意，因此固定零均值，用 OAS 数据驱动收缩保证矩阵可逆，
不手调 shrinkage。距离为 `(z_i-z_j)^T W^-1 (z_i-z_j)`。

固定三个模式：

1. `joint`：完整 W，包含 composition/coverage 跨块项。
2. `block`：使用同一 W 的对角块，跨块置零；控制是否真有 joint coupling 收益。
3. `coverage`：仅用 coverage 块；控制增益是否只是 coverage 噪声校准。

这是一阶线性/二次度量试验，不是复现 SemiBin/COMEBin 神经网络。
不构造随机“负对”，因此不存在把同 genome 随机对强行监督为负的问题；
代价是没有直接学习区分近缘 genome 的负向约束。

采用按父 contig 名称固定哈希的二折交叉检查，报告未见父 contig 两半的
距离压缩/相对于随机其他片段的排序，不用 gold 选择模式。最终度量使用
所有正对拟合；不将该 proxy 当成 genome-level accuracy。

## 图与评估

先在现有 PMH retained graph 上仅换边权，隔离 learned weight 的作用；
再用同一 embedding 选 exact mutual top-200，并用于边权，检查统一度量。
保留同一候选图+原 coverage 分数以及保留原 coverage 准入的对照。
权重仍用前轮 self-tuning kNN kernel，带宽由已有第200邻居决定，不调参。
Fisher 与 raw-sum LPA 分别报告，防止 kernel 尺度与 Fisher 偏好混淆。
所有表示/消融在三套数据上完整报告，不 cherry-pick。

FAISS exact search 是质量验证原型，不预先宣称满足生产性能目标；BAM
计数、特征、学习、检索、binning 和内存分别计时。若没有跨数据集的
可靠质量改善，不用金标准挑配置或把原型并入默认流程。

## 训练诊断后补充、查看本轮 AMBER 前固定的信噪比对照

纯 within whitening 会给每个方向单位噪声，未必保留判别信号。增加由
同一加性噪声模型推导的 signal 版本，同样 joint/block/coverage 三种结构：

`W = E[(left-right)(left-right)^T]`

`C = Cov((left+right)/2)`

独立父片段差的协方差为 `U = 2*C + W/2`。同源与独立源高斯密度的
对数比，其距离部分为 `M = W^-1 - U^-1`。在 W 白化域估计
`V = 2*W^-1/2 C W^-1/2 - I/2`，将负信号特征值截为零（PSD 约束），
每个方向的度量权重自动为 `lambda/(1+lambda)`。没有人工选保留维数，
没有标签负对；OAS 仍由数据估计。额外的 Gaussian/同方差/独立噪声假设
必须在结果中说明，不把 log likelihood 当成已经校准的同 genome 后验。

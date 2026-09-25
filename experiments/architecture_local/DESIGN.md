# RabbitBin 架构内的局部实验

基线：origin/main bf9e1013941bf7864501e1ac771c77db73525d6e。
不替换 PMH，不引入外部训练/标记基因，不读取 gold 作模型输入。
沿用 >=1000 bp contig 输入，>=2500 bp 长 contig 建图与后续短片段招募，
>=200 kb bin 参与 RabbitBin AMBER；三个 CAMI2 数据集统一 seed=42、64 线程、
PMH m=500、top-200、dual-depth=5、现有 Fisher neutral-point coverage cutoff。

## 预先固定的问题和干预

1. 同一图上检查 Fisher CDF 饱和、绝对平局容差和重复标签冻结。
   先用已有 fisher / tail / logsum 开关；再测试数学上保持 Fisher 排序的
   log-survival 实现及明确的浮点平局规则。后者是数值/决策规则消融，
   不声称与包含旧绝对容差的生产输出等价。
2. 让已有 coverage 完整准入条件在 PMH top-200 截断前生效。
   原有 Spearman 初筛保留；只把本来稍后会拒绝的边提前排除，给其他候选
   留出容量。PMH 排序、mutual-kNN、coverage 分数及 cutoff 均不变。
   这会改变候选集合，不是与旧图等价的加速优化。

每项分别在三套数据上报告 HQ / MQ / weighted purity、时间、内存。
不按数据集选模式，不根据 gold 调阈值，不将 HQ 单项增加当作整体改善。
默认行为保持不变，所有试验 opt-in。

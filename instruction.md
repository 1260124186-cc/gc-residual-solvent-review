# 完成残留溶剂乙酸乙酯的序列审核

车间刚刚完成一批药物中间体的顶空气相色谱检查。分析员把五个校准点、两个方法空白、两份平行样品和对应的加标样品放在同一序列里，现在需要把原始色谱和序列表整理成可以直接交给QA复核的JSON报告。`/workspace` 下已经有一套早期脚本，但它只适用于无基线漂移、固定稀释倍数的演示数据，不能直接用于这次序列。请保留现有命令行入口和包名，补齐并修正实现。

命令行格式：

```bash
python -m gc_review.cli SEQUENCE.csv CHROMATOGRAMS.csv OUTPUT.json
```

`sequence.csv` 固定包含 `run_id,sample_id,role,replicate,dilution_factor,nominal_mg_l,spike_added_mg_l`。`role` 只允许 `calibration`、`blank`、`sample`、`spike`。校准行必须填写正数 `nominal_mg_l`，加标行必须填写正数 `spike_added_mg_l`；其余不适用字段留空。`chromatograms.csv` 固定包含 `run_id,time_s,signal_mv`，同一进样的时间必须严格递增。进样号在序列中必须唯一，且色谱与序列的进样号集合必须完全一致。

乙酸乙酯的积分窗口为39.2至44.8秒，积分前锚点区间为36.4至38.0秒，积分后锚点区间为46.0至47.6秒。内标峰使用54.4至59.6秒的积分窗口、52.0至53.6秒的前锚点和60.4至62.0秒的后锚点。每个锚点区间的基线值取该区间内信号算术平均值，时间取区间中心；把前后两个基线点连成直线，作为整个积分窗口的基线。面积用相邻采样点梯形法累加，纵坐标是扣除该点线性基线后的信号，不做截断。

乙酸乙酯面积除以内标面积得到响应比 `r`。校准曲线采用带截距模型 `r = a + b*x`，每个校准点的权重为 `1/x^2`。请按这些权重求加权最小二乘解，并用同一组权重计算加权决定系数 `R2` 和加权残差标准误 `s_yx = sqrt(sum(w_i*(r_i-a-b*x_i)^2)/(n-2))`。`LOD = 3.3*s_yx/b`，`LOQ = 10*s_yx/b`。校准至少要有五个不同浓度且 `b > 0`、`R2 >= 0.995`。每个校准点的回算回收率为 `100*((r-a)/b)/x`，最低浓度允许80%至120%，其余浓度允许85%至115%；任一点都不合格时校准不通过。

样品和加标的仪器浓度先按 `(r-a)/b` 计算。方法空白也按同样方式计算，并乘以各自的稀释倍数；空白校正值取空白结果均值和零中的较大者。样品最终浓度等于 `max(0, 仪器浓度*稀释倍数 - 空白校正值)`。同一 `sample_id` 的两条 `sample` 行是平行样，报告其平均浓度，并计算 `RPD=100*abs(c1-c2)/((c1+c2)/2)`；两个浓度均为零时RPD记为0。加标回收率为 `100*(加标样最终浓度-平行样平均浓度)/spike_added_mg_l`。

空白校正值不得超过0.10 mg/L，平行样RPD不得超过10%，加标回收率必须在80%至120%之间。若校准或空白不合格，所有样品状态为 `qc_failed`；否则浓度低于LOQ的样品状态为 `below_loq`，平行样或加标不合格的样品状态为 `qc_failed`，其余为 `pass`。整批状态在有任一 `qc_failed` 样品时为 `fail`，否则为 `pass`；`below_loq` 本身不导致整批失败。

输出的根对象固定包含 `calibration`、`blank_mean_mg_l`、`blank_correction_mg_l`、`qc`、`samples`、`overall_status`。`calibration` 包含 `slope`、`intercept`、`weighted_r2`、`sdyx`、`lod_mg_l`、`loq_mg_l`、`accepted`、`standards`；`standards` 是按浓度升序排列的对象数组，每个对象包含 `nominal_mg_l`、`response_ratio`、`recovery_percent`。`qc` 包含 `calibration`、`blank`、`duplicate`、`spike_recovery`，值为 `pass` 或 `fail`，其中后两项表示该批所有相关样品是否合格。`samples` 以样品名为键，每个值包含 `concentration_mg_l`、`replicates_mg_l`、`duplicate_rpd_percent`、`spike_recovery_percent`、`status`。所有浓度单位为mg/L，面积使用色谱文件中的信号单位和秒。数值不要提前取整。

输入缺列、空字段不符合角色要求、数值非有限、稀释倍数非正、未知角色、同一进样内时间重复或倒序、序列与色谱进样号不一致，以及样品缺少两条平行样或缺少加标样时，程序应以非零状态退出，不能生成一个看似完整的报告。`/workspace/data` 中有一次完整的模拟序列，可用于本地运行。该目录是教学和交接用的模拟记录，不是原始实验记录。

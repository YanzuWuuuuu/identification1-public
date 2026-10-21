# 非参数 FRF 对照图

`plot_frf.py --z all --panels all` 读取已保存的 FRF 与参数模型，不重新估计、不插值、不丢弃频点。`--z 0` 只生成中心位置；`--reverse` 可额外生成后半模型对前半 FRF 的反向检查。

每个位置包括：

- `full_vs_full`：8 周期模型与同记录 8 周期 FRF，仅为拟合检查。
- `front_vs_back`：前 4 周期模型与后 4 周期 FRF。原固定时延仍共享先前估计，因此不是全流程盲测。
- `main`：M1、M3、M5n 与灰色 LPM 参照。
- `all_methods`：M1、M2、M3、M3w、M5、M5n 全部方法。未收敛端点图例明确标为 failed endpoint。
- `magnitude`：幅值；`phase`：各曲线独立展开相位；`complex_difference`：相对于 LPM 的逐频点复相对差。

每图均有 PNG、PDF、SVG。`source.csv` 保留画图用的复数值和随机标准误；`figure_manifest.json` 保留源文件摘要、频点数、收敛状态以及相位整体整周移动；`figure_QA.json` 为显示恒等式与文字可编辑性检查。

灰色带是 ±2 倍随机标准误，不含 LPM 的平滑偏差及反演偏差，不当作 95% 置信带。误差图中的灰色点线是 2×随机复标准误除以参照幅值。相位展开仅用于展示，所有误差都由原始复数 FRF 计算。相位局部整周不确定和低幅值附近的宽误差带均保留。

原数据的幅值单位尚未重新校准，所以幅值标签采用 recorded units。FRF 是估计参照，不是真值；曲线更接近参照不等于该方法更准确。

## LRM 诊断补图

`plot_frf.py --z all --reference lrm` 另存带 `lrm_` 前缀的主图。采用已有一阶、31 谱线的 LRM 结果，配置通过前半数据的纯数据交叉验证固定。随机标准误来自逐周期 jackknife，仍不含平滑和反演偏差。这是 LPM 分辨率问题的诊断补充，不能覆盖、删除 LPM 结果，也不能将更接近 LRM 的方法判定为更准确。

复现说明（Python 3, numpy, scipy）
SISO：
  python3 fd_cl.py 20000 0.5 200        # 9.1(b) 表
  python3 pred_siso.py 0.15 100         # 9.1(a) 一阶偏差公式 vs 蒙特卡洛
MIMO 支承模态（9.2 节，需设置环境变量）：
  export NOFLEX=1 FC=60 AL=30 AH=100 FHI=100 NN=300000 SW=0.03 NPER=10
  python3 predict2.py 0.2 3e-7          # 定理1一阶偏差 + 定理2协方差预测
  python3 exact_limit.py                # 极限方程(4.1)的精确根
  python3 run_mc4.py 0.2 3e-7 25 101 mc.npz dir_band_W1,dir_band_Wg,ref_band_W1,ref_band_Wg
  python3 design3.py 0.2 3e-7           # 各逻辑坐标的反馈噪声占比与输出信噪比
MIMO 刚体+柔性对照（9.3 节）：python3 predict.py 1e-7   （使用 mimo_cl.py）
结果：results/mcA..D.npz（15+15+35+35=100 次蒙特卡洛；列 = log频率误差(3), log阻尼误差(3), 迭代次数）

v2 新增：
  chk_C_siso.py                          # SISO：M^-1 b, (M-C)^-1 b, M0^-1 b 对照
  HH=1e-4 python3 chk_C.py               # MIMO：同上（需 9.2 的环境变量；外层差分步长须 >=1e-4，1e-6 数值不稳定）
  数据长度扫描（9.3 节）：
  export NOFLEX=1 FC=60 AL=30 AH=100 FHI=100 SW=0.03 FSPR=0.01 ZSPR=0.4
  NN=90000 NPER=3 python3 run_mc5.py 0.2 3e-7 40 1001 sw_090k.npz dir_band_W1,ref_band_W1,dir_band_Wm,ref_band_Wm
  （300k: NN=300000 NPER=10；900k: NN=900000 NPER=30）
  python3 analyze_sweep.py               # 偏差/标准误/RMSE/零空间角

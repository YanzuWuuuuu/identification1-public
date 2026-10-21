实测对比（第 10.1 节）复现
1) MATLAB 导出各位置数据（已执行，输出在 experiment_data/srivc_no_lrm_20261015/closed_loop_consistency_20261019/cl_data_z*.mat）：
   读 auto_identify_motor1234_z*.mat 的 u_Motor1..4、pre_y_0..3__0_、u_z/u_Rx/u_Ry（第 2 列），D = U - V*Q_alloc.'，存 -v7。
2) Python：
   ZPOS=0 python3 run_real.py 1 9 dir_all,dir_lines,ref_K,ref_loo     # 全 8 周期
   ZPOS=0 python3 run_real.py 1 5 ...  / 5 9 ...                       # 两半
   ZPOS=0 python3 xval.py                                              # 跨半段验证
   （load.py 中 DIR 指向 cl_data_z*.mat 所在目录，需按本机路径修改）
方法名：dir_all=直接法全部频点；dir_lines=直接法只用激励谱线；ref_K=参考IV(控制器由 v_fb<-y ARX 精确辨识)；ref_loo=留一法两阶段工具变量(不稳定，仅参考)

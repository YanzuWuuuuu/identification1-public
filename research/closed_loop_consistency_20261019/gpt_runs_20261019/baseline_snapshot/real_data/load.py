import scipy.io as sio, numpy as np, os
F='/mnt/user-data/uploads/identifiation1/experiment_data/srivc_no_lrm_20261015/full_coordinate_calibration_20261016/recovered_reference_z0_candidate.mat'
DIR='/mnt/user-data/uploads/identifiation1/experiment_data/srivc_no_lrm_20261015/closed_loop_consistency_20261019/'
def load():
    z=os.environ.get('ZPOS')
    if z is None:
        R=sio.loadmat(F)
        return dict(U=R['U'],Y=R['Y'],D=R['D_candidate'],V=R['V_logic'],Q=R['Q_allocation'],P=int(R['period_samples'][0,0]),fs=float(R['fs'][0,0]))
    R=sio.loadmat(DIR+f'cl_data_z{z}.mat')
    return dict(U=R['U'],Y=R['Y'],D=R['D'],V=R['V'],Q=R['Q'],P=int(R['P'][0,0]),fs=float(R['fs'][0,0]))

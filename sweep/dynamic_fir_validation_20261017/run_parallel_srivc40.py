"""Independent source/configuration jobs; writes disjoint candidate files.
The main run subsequently loads these frozen training-only results.
"""
from scipy.io import loadmat
from run_validation import OUT, read_coords, project, run_candidate

if __name__=='__main__':
    raw=loadmat(OUT/'input_z0.mat',squeeze_me=True,variable_names=['f','Htrain','Hmodel'])
    c=read_coords();f=raw['f'];Pt=project(raw['Htrain'],c);Pm=project(raw['Hmodel'],c)
    for delay in [0,2,5]:run_candidate('srivc',40,delay,f,Pm,Pt,c)

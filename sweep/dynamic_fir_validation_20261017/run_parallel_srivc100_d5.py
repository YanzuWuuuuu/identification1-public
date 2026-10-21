"""Disjoint delay-five job; main reuses its completed candidate."""
from scipy.io import loadmat
from run_validation import OUT,read_coords,project,run_candidate
if __name__=='__main__':
    raw=loadmat(OUT/'input_z0.mat',squeeze_me=True,variable_names=['f','Htrain','Hmodel'])
    c=read_coords();f=raw['f']
    run_candidate('srivc',100,5,f,project(raw['Hmodel'],c),project(raw['Htrain'],c),c)

"""Audit saved sensitivity results and bound their engineering interpretation."""
from pathlib import Path
import json,csv,hashlib
import numpy as np

HERE=Path(__file__).resolve().parent
# Reuse only model setup, not campaign execution. MATLAB independently rebuilt P.
setup=(HERE/'perturbation_sweep.py').read_text(encoding='utf-8').split('rows=[]')[0]
scope={'__file__':str(HERE/'perturbation_sweep.py')}
exec(compile(setup,'model_setup','exec'),scope)
B,C,B2,C2,M,ny,nu,dc,db=[scope[k] for k in ['B','C','B2','C2','M','ny','nu','dc','db']]
freq0,zeta=scope['freq0'],scope['zeta']; OUT=scope['OUT']
rows=list(csv.DictReader((OUT/'cases.csv').open(encoding='utf-8-sig')))
mat=json.loads((OUT/'audit_matlab.json').read_text())

def model(f):
    s=2j*np.pi*np.asarray(f)
    F=sum(np.outer(C2[:,i],B2[i])[None]/(s*s+2*zeta[i]*2*np.pi*freq0[i]*s+(2*np.pi*freq0[i])**2)[:,None,None] for i in range(5))
    return M[None]/s[:,None,None]**2+F,F

expected=model(mat['hz'])[0]
check=np.moveaxis(np.array(mat['frf_real'])+1j*np.array(mat['frf_imag']),2,0)
frf_err=float(np.linalg.norm(check-expected)/np.linalg.norm(expected))
assert frf_err<1e-8,frf_err

def geom(e):
    cp=C+e*dc;bp=B+e*db
    y=np.linalg.svd(cp,full_matrices=True)[0][:,:,-1]
    u=np.linalg.svd(bp,full_matrices=True)[2][:,-1,:]
    y*=np.sign(y@ny)[:,None];u*=np.sign(u@nu)[:,None]
    return cp,bp,y,u

# Independently use analytic integral of f^-4 for every nonzero static case.
errors=[]; geometry_errors=[]
for level in [0.0001,0.0003,.001,.003,.01,.03,.1]:
    cp,bp,y,u=geom(level)
    geometry_errors.extend(abs(np.linalg.norm(cp-C,axis=(1,2))/np.linalg.norm(C)-level))
    geometry_errors.extend(abs(np.linalg.norm(bp-B,axis=(1,2))/np.linalg.norm(B)-level))
    for lo in [2,20]:
        val=1e6*abs(np.einsum('na,ab,nb->n',y,M,u))*np.sqrt((lo**-3-800**-3)/(3*(800-lo)))/(2*np.pi)**2
        selected=sorted([r for r in rows if r['kind']=='geometry_relative' and r['method']=='both_nominal' and float(r['level'])==level and int(r['band_low'])==lo],key=lambda r:int(r['seed']))
        reported=np.array([float(r['rigid_nm']) for r in selected])
        errors.extend(abs(reported-val)/val)

# Fine quadrature check on ideal flexible output, both evaluation bands.
fine=np.linspace(2,800,100001)
_,F=model(fine); z=np.einsum('a,fab,b->f',ny,F,nu)
quad=[]
for lo in [2,20]:
    ff=np.r_[lo,fine[fine>lo]]
    zz=np.r_[np.einsum('a,fab,b->f',ny,model([lo])[1],nu),z[fine>lo]]
    val=1e6*np.sqrt(np.trapezoid(abs(zz)**2,ff)/(800-lo))
    old=float(next(r for r in rows if r['kind']=='geometry_relative' and r['method']=='both_nominal' and float(r['level'])==0 and int(r['band_low'])==lo)['flex_nm'])
    quad.append({'band_low':lo,'relative_error':abs(old-val)/val,'old_nm':old,'fine_nm':val})
assert max(errors)<1e-3
assert max(q['relative_error'] for q in quad)<1e-3

# Equal actual mechanical rigid motion, only compare output extraction.
# Each of 4 sensors has 1um spatial-RMS mechanical displacement at the test tone.
# Synthetic single-tone test, not a closed-loop trajectory or identified model.
_,_,y0,_=geom(.01); _,_,yr,_=geom(.001)
rigid_tests=[]
for hz in [20,100]:
    s=2j*np.pi*hz
    for delay in [0,10,100,1000]:
        dy=100e-6+delay*1e-6*np.array([0,.3,.7,1])
        hy=np.exp(-s*dy)/(1+s/(2*np.pi*2000))
        delta=.9*(dy-dy[0]);L=max(delta)+20e-6
        cy=np.exp(-s*(L-delta))
        for axis in range(3):
            physical=C[:,axis]/np.sqrt(np.mean(C[:,axis]**2))*1000
            for method,y,cc in [('nominal_1pct',y0,np.ones(4)),('static_0p1pct',yr,np.ones(4)),('dynamic_prescribed',yr,cy)]:
                leak=abs(np.einsum('na,a->n',y,cc*hy*physical))
                rigid_tests.append(dict(hz=hz,delay_spread_us=delay,axis=['rx','ry','z'][axis],method=method,leak_nm_p10=float(np.percentile(leak,10)),leak_nm_p50=float(np.median(leak)),leak_nm_p90=float(np.percentile(leak,90))))

# Retention per modal coefficient, align to ideal output/input directions;
# at each pole frequency, divide out estimated common delays/reference filters.
# This is a spatial-direction comparison: common attenuation removed explicitly.
_,_,yr,ur=geom(.001)
retention=[]
for delay in [10,100,1000]:
    dy=100e-6+delay*1e-6*np.array([0,.3,.7,1]);du=100e-6+delay*1e-6*np.array([1,.6,.2,0])
    for method,fraction in [('static',0.),('dynamic_prescribed',.9)]:
        per=[]
        for i,hz in enumerate(freq0):
            s=2j*np.pi*hz
            ey=np.exp(-s*(1-fraction)*(dy-dy[0]));eu=np.exp(-s*(1-fraction)*(du-du[0]))
            actual=(yr@(C2[:,i]*ey))*(ur@(B2[i]*eu))
            target=(ny@C2[:,i])*(B2[i]@nu)
            rel=abs(actual-target)/abs(target)
            per.append(dict(mode=i+1,complex_relative_error_p50=float(np.median(rel)),amplitude_ratio_p50=float(np.median(abs(actual/target)))))
        retention.append(dict(delay_spread_us=delay,method=method,modes=per))

result=dict(matlab_frf_relative_error=frf_err,geometry_relative_norm_max_error=float(max(geometry_errors)),static_analytic_integral_max_relative_error=float(max(errors)),flex_quadrature=quad,
saved_input_data={k:v for k,v in mat.items() if not k.startswith('frf_')},
rigid_motion_tests=rigid_tests,modal_projection_checks=retention,
noise_line=dict(N=20000,fs=5000,sample_noise_nm=100,complex_peak_phasor_noise_rms_nm=2*100/np.sqrt(20000),equivalent_sinusoidal_rms_noise_nm=100*np.sqrt(2/20000),assumptions='Exact excited frequency, rectangular 4s record, independent white output noise, unit output vector, magnitude-one delay compensator. Not valid automatically for closed-loop colored noise.'),
scope='Audits existing exact frequency sensitivity. Added rigid motion and modal projection checks are synthetic and do not rerun identification.')
(OUT/'audit_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k not in ['rigid_motion_tests','modal_projection_checks']},indent=2))
print('Z-axis synthetic 1um-per-sensor RMS motion:')
for r in rigid_tests:
    if r['hz']==20 and r['axis']=='z':print(r)
print('Modal checks:')
for r in retention:
    if r['delay_spread_us']==100:print(r)

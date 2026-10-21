"""Controlled stage-0 sensitivity, not an identification performance experiment."""
from pathlib import Path
import re,json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUT=HERE/'perturbation_outputs'; OUT.mkdir(exist_ok=True)
txt=(ROOT/'matlab/0/Multiple_main_Control.m').read_text(encoding='utf-8',errors='replace')
def matrix(name):
    s=re.search(r'^\s*'+name+r'\s*=\s*\[([^\]]+)\]',txt,re.M).group(1)
    return np.array([[float(x) for x in row.replace(',',' ').split()] for row in s.split(';')])
B=matrix('B1'); C=matrix('C1'); B2=matrix('B2'); C2=matrix('C2')
freq0=matrix('f').ravel(); zeta=matrix('kesai').ravel()
mass=[float(re.search(r'^\s*'+n+r'\s*=\s*([\d.]+)',txt,re.M).group(1)) for n in ['Jrx','Jry','mz']]
M=C@np.diag(1/np.array(mass))@B
ny=np.linalg.svd(C,full_matrices=True)[0][:,-1]
nu=np.linalg.svd(B,full_matrices=True)[2][-1]
f=np.unique(np.r_[np.geomspace(2,200,700),np.linspace(200,800,2401),20])
s=2j*np.pi*f
A=np.einsum('ai,ib->iab',C2,B2)
F=np.einsum('fi,iab->fab',1/(s[:,None]**2+2*zeta*(2*np.pi*freq0)*s[:,None]+(2*np.pi*freq0)**2),A)
R=M[None,:,:]/s[:,None,None]**2
seeds=80
rng=np.random.default_rng(20261014)
dc=rng.normal(size=(seeds,4,3)); db=rng.normal(size=(seeds,3,4))
# Relative column/row error preserves units of rotational vs translational ports.
dc=dc/np.linalg.norm(dc,axis=1)[:,None,:]*np.linalg.norm(C,axis=0)[None,None,:]
db=db/np.linalg.norm(db,axis=2)[:,:,None]*np.linalg.norm(B,axis=1)[None,:,None]
ey=rng.normal(size=(seeds,4)); eu=rng.normal(size=(seeds,4))
ey-=np.outer(ey@ny,ny); ey/=np.linalg.norm(ey,axis=1)[:,None]
eu-=np.outer(eu@nu,nu); eu/=np.linalg.norm(eu,axis=1)[:,None]

def geom(e):
    cp=C+e*dc;bp=B+e*db
    y=np.linalg.svd(cp,full_matrices=True)[0][:,:,-1]
    u=np.linalg.svd(bp,full_matrices=True)[2][:,-1,:]
    y*=np.sign(y@ny)[:,None];u*=np.sign(u@nu)[:,None]
    return y,u

rows=[]
def evaluate(kind,x,method,ys,us,hy=None,hu=None,cy=None,cu=None):
    hy=np.ones((len(f),4),complex) if hy is None else hy
    hu=np.ones_like(hy) if hu is None else hu
    cy=np.ones_like(hy) if cy is None else cy
    cu=np.ones_like(hu) if cu is None else cu
    # Match actual actuator command-vector norm at each frequency between methods.
    vv=us[:,None,:]*cu[None,:,:]
    vv/=np.linalg.norm(vv,axis=2)[:,:,None]
    left=ys[:,None,:]*cy[None,:,:]*hy[None,:,:]
    right=vv*hu[None,:,:]
    zr=np.einsum('nfa,fab,nfb->nf',left,R,right,optimize=True)
    zf=np.einsum('nfa,fab,nfb->nf',left,F,right,optimize=True)
    gain=np.sum(abs(ys[:,None,:]*cy[None,:,:])**2,axis=2)
    for lo in (2,20):
        ix=f>=lo;ff=f[ix];bw=800-lo
        # Total norm of input RMS=1 mN; flat one-sided PSD over each band.
        rigid=1e6*np.sqrt(np.trapezoid(abs(zr[:,ix])**2,ff,axis=1)/bw)
        flex=1e6*np.sqrt(np.trapezoid(abs(zf[:,ix])**2,ff,axis=1)/bw)
        # 100nm/sample at fs=5kHz, independent output-added white sensor noise.
        noise=np.sqrt(np.trapezoid(gain[:,ix]*100**2/2500,ff,axis=1))
        for n in range(len(ys)):
            rows.append(dict(kind=kind,level=float(x),method=method,band_low=lo,seed=n,
                rigid_nm=rigid[n],flex_nm=flex[n],noise_nm=noise[n],
                rigid_flex=rigid[n]/flex[n],rigid_noise=rigid[n]/noise[n],
                sinr_db=20*np.log10(flex[n]/np.hypot(rigid[n],noise[n])),
                input_mN_for_leak_equals_noise=noise[n]/max(rigid[n],1e-25),
                angle_y_deg=np.degrees(np.arccos(np.clip(abs(ys[n]@ny),0,1))),
                angle_u_deg=np.degrees(np.arccos(np.clip(abs(us[n]@nu),0,1)))))

levels=[0,0.0001,0.0003,0.001,0.003,0.01,0.03,0.1]
for e in levels:
    y,u=geom(e)
    evaluate('geometry_relative',e,'both_nominal',y,u)
    evaluate('geometry_relative',e,'output_only',y,np.tile(np.array([1.,0,0,0]),(seeds,1)))
    evaluate('geometry_relative',e,'one_side_exact',np.tile(ny,(seeds,1)),u)
for deg in [0,.001,.003,.01,.03,.1,.3,1,3,10]:
    a=np.deg2rad(deg)
    evaluate('null_angle_deg',deg,'both_nominal',np.cos(a)*ny+np.sin(a)*ey,np.cos(a)*nu+np.sin(a)*eu)

pattern_y=np.array([0,.3,.7,1]);pattern_u=np.array([1,.6,.2,0])
# Fixed directions across all channel mismatch levels; 1% nominal vs .1% residual.
y0,u0=geom(.01);yr,ur=geom(.001)
for kind,vals in [('delay_spread_us',[0,1,3,10,30,100,300,1000]),('timeconstant_spread_fraction',[0,.01,.03,.1,.3,1])]:
    for x in vals:
        ty=np.full(4,1/(2*np.pi*2000));tu=ty.copy()
        dy=np.full(4,100e-6);du=dy.copy()
        if kind=='delay_spread_us':dy+=x*1e-6*pattern_y;du+=x*1e-6*pattern_u
        else:ty*=1+x*pattern_y;tu*=1+x*pattern_u
        hy=np.exp(-s[:,None]*dy)/(1+s[:,None]*ty)
        hu=np.exp(-s[:,None]*du)/(1+s[:,None]*tu)
        for method,y,u,fraction in [('B0_nominal',y0,u0,None),('B1_static_exact',np.tile(ny,(seeds,1)),np.tile(nu,(seeds,1)),None),('B1_static_residual',yr,ur,None),('B4_oracle',np.tile(ny,(seeds,1)),np.tile(nu,(seeds,1)),1.),('B4_residual',yr,ur,.9)]:
            cy=cu=None
            if fraction is not None:
                # Estimate 90% of each departure from reference; 10% remains.
                dhy=dy[0]+fraction*(dy-dy[0]);dhu=du[0]+fraction*(du-du[0])
                thy=ty[0]+fraction*(ty-ty[0]);thu=tu[0]+fraction*(tu-tu[0])
                def comp(d,t):
                    delta=d-d[0];L=max(0,max(delta))+20e-6
                    return np.exp(-s[:,None]*(L-delta))*(1+s[:,None]*t)/(1+s[:,None]*t[0])
                cy=comp(dhy,thy);cu=comp(dhu,thu)
            evaluate(kind,x,method,y,u,hy,hu,cy,cu)

with (OUT/'cases.csv').open('w',newline='',encoding='utf-8-sig') as h:
    w=csv.DictWriter(h,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
summary=[]
keys=sorted(set((r['kind'],r['level'],r['method'],r['band_low']) for r in rows))
for k in keys:
    rr=[r for r in rows if (r['kind'],r['level'],r['method'],r['band_low'])==k]
    q=dict(zip(['kind','level','method','band_low'],k))
    for metric in ['rigid_nm','flex_nm','noise_nm','rigid_flex','rigid_noise','sinr_db','input_mN_for_leak_equals_noise','angle_y_deg','angle_u_deg']:
        v=[r[metric] for r in rr]
        for pct in [10,50,90]:q[f'{metric}_p{pct}']=float(np.percentile(v,pct))
    summary.append(q)
(OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
contract=dict(source='matlab/0/Multiple_main_Control.m',modes=freq0.tolist(),D=0,seeds=seeds,
    seed=20261014,input_rms_mN=1,noise_sample_rms_nm=100,fs_hz=5000,
    noise_location='independent additive measured-output noise before output compensation',
    channel_base_bandwidth_hz=2000,channel_base_delay_us=100,
    note='Synthetic exact frequency-response propagation; no estimator rerun, no closed-loop trajectory claim. B4 residual is prescribed 0.1% geometry and 10% channel-departure error, not measured estimation accuracy.',
    geometry='Every C column and B row perturbed by specified relative Euclidean norm; true plant unchanged.',
    normalization='Equal total actual command RMS 1mN and flat PSD per evaluation band for every method; continuous frequency quadrature.')
(OUT/'contract.json').write_text(json.dumps(contract,indent=2),encoding='utf-8')
# Verification: exact one-side null cancels rigid term; equal channels exact static cancels.
assert max(r['rigid_nm'] for r in rows if r['method']=='one_side_exact')<1e-7
assert max(r['rigid_nm'] for r in rows if r['method']=='B4_oracle')<1e-7
assert max(r['rigid_nm'] for r in rows if r['method']=='B1_static_exact' and r['level']==0)<1e-7
assert np.max(abs(ny@M))<1e-12 and np.max(abs(M@nu))<1e-12
fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for lo,col in [(2,'tab:blue'),(20,'tab:orange')]:
    for kind,method,a in [('geometry_relative','both_nominal',ax[0,0]),('null_angle_deg','both_nominal',ax[0,1])]:
        rr=[q for q in summary if q['kind']==kind and q['method']==method and q['band_low']==lo and q['level']>0]
        xx=[q['level']*(100 if kind=='geometry_relative' else 1) for q in rr]
        a.loglog(xx,[q['rigid_flex_p50'] for q in rr],color=col,label=f'{lo}-800 Hz')
        a.fill_between(xx,[q['rigid_flex_p10'] for q in rr],[q['rigid_flex_p90'] for q in rr],color=col,alpha=.15)
for kind,a in [('delay_spread_us',ax[1,0]),('timeconstant_spread_fraction',ax[1,1])]:
    for method in ['B0_nominal','B1_static_residual','B4_residual']:
        rr=[q for q in summary if q['kind']==kind and q['method']==method and q['band_low']==20]
        xx=[q['level']*(100 if kind=='timeconstant_spread_fraction' else 1) for q in rr]
        a.semilogy(xx,[max(q['rigid_nm_p50'],1e-12) for q in rr],marker='o',label=method)
    a.axhline(np.sqrt((800-20)/2500)*100,color='k',ls='--',label='Static output noise')
    rr=[q for q in summary if q['kind']==kind and q['method']=='B4_residual' and q['band_low']==20]
    a.semilogy([q['level']*(100 if kind=='timeconstant_spread_fraction' else 1) for q in rr],[q['noise_nm_p50'] for q in rr],ls=':',color='red',label='B4 output noise')
ax[0,0].set(xlabel='C-column / B-row error (%)',ylabel='Rigid / target flexible RMS')
ax[0,1].set(xlabel='Null vector angular error (deg)',ylabel='Rigid / target flexible RMS')
ax[1,0].set(xlabel='Per-side delay spread (microseconds)',ylabel='RMS displacement (nm), input 1 mN')
ax[1,1].set(xlabel='Per-side time-constant spread (%)',ylabel='RMS displacement (nm), input 1 mN')
for a in ax.flat:a.grid(True,which='both',alpha=.3);a.legend(fontsize=8)
fig.savefig(OUT/'perturbation_summary.png',dpi=180)
print(json.dumps({'cases':len(rows),'groups':len(summary),'checks':'passed','output':str(OUT)}))

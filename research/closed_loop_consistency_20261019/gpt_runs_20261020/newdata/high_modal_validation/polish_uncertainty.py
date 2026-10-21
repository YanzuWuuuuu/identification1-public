from run import *
from analyze import loadfit,uncertainty,csvout,match
import itertools

def main():
    dest=HERE/'diagnostics';pairs=[]
    for w in ('H100','H300','H493'):
        s,m,_,_=setup(w);anchor=np.asarray(loadfit('A_nominal',w,'M3')['beta'])
        for method in ('M3','M3w'):
            estimates={}
            for ctrl in ('lowgain','nominal','highgain'):
                rid='A_'+ctrl;p=dest/f'{rid}_{w}_{method}.json'
                if not p.exists():continue
                polish=json.loads(p.read_text());r=loadfit(rid,w,method);r.update(beta=polish['beta'],strict_score_pass=polish['polished_strict_pass'])
                u=uncertainty(rid,w,method,r,dest)
                mt=match(m,anchor,np.asarray(r['beta']),s);perm=[a['fit_index'] for a in mt]
                ff,zz=m.modal(np.asarray(r['beta']));estimates[ctrl]=(ff[perm],zz[perm],np.asarray(u['frequency_u_hz'])[perm],np.asarray(u['damping_u_pp'])[perm],u,mt)
            for ca,cb in itertools.combinations(estimates,2):
                aa=estimates[ca];bb=estimates[cb]
                for j in range(m.nm):
                    for k,par in enumerate(('f','zeta')):
                        du=float(np.hypot(aa[k+2][j],bb[k+2][j]));delta=float(bb[k][j]-aa[k][j])
                        ok=aa[4]['local_numerical_conditions'] and bb[4]['local_numerical_conditions'] and aa[5][j]['candidate_match'] and bb[5][j]['candidate_match']
                        pairs.append(dict(window=w,method=method,component=j+1,background=s['background'][j],parameter=par,a=ca,b=cb,difference=delta,combined_u=du,abs_z=abs(delta)/du,eligible=bool(ok),below_2u=abs(delta)<2*du))
    csvout('polished_controller_pairs.csv',pairs)
    print('POLISHED_UNCERTAINTY_COMPLETE',len(pairs))

if __name__=='__main__':main()

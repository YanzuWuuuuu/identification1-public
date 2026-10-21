from run import *
sys.path.insert(0,str(PARENT/'baseline'))
import realfit_extended as old

def validate():
    rng=np.random.default_rng(20261020);f=np.geomspace(3,35,401);c=old.Ctx(f,TAU,'zoh');m=ModalModel(3,TAU,25);cc=m.context(f)
    b=rng.normal(size=old.NB)*1e-6;b[old.IW:old.IZ]=np.log(2*np.pi*np.array([10.6,11.7,16.1]));b[old.IZ:old.ITR]=np.log([.03,.04,.05]);b[old.ITR:]*=.01
    order=[]
    for o in range(4):
        for j in range(3):order.extend(list(range(32*j+4*o,32*j+4*o+4))+list(range(32*j+16+4*o,32*j+16+4*o+4))+[old.ITR+8*j+2*o,old.ITR+8*j+2*o+1])
        order.extend(list(range(96+4*o,96+4*o+4))+list(range(112+4*o,112+4*o+4)))
    order+=list(range(old.IW,old.ITR));bn=b[order]
    U=rng.normal(size=(len(f),4))+1j*rng.normal(size=(len(f),4));Uh=U+.1*rng.normal(size=U.shape)
    yy,jj=m.output(cc,bn,U);y0,j0=old.out_jac(c,b,U)
    rep=dict(output_vs_old_relative=float(np.linalg.norm(yy-y0)/np.linalg.norm(y0)),
             jac_vs_old_forward_difference_relative=float(np.linalg.norm(jj-j0[:,:,order])/np.linalg.norm(j0)),
             transfer_vs_old_relative=float(np.linalg.norm(m.resp(cc,bn)-old.resp(c,b))/np.linalg.norm(old.resp(c,b))))
    L=rng.normal(size=(len(f),4,4))+1j*rng.normal(size=(len(f),4,4));W=L@L.conj().transpose(0,2,1)+np.eye(4)
    Y=yy+(rng.normal(size=yy.shape)+1j*rng.normal(size=yy.shape))*1e-6
    zz=m.output(cc,bn,Uh)[1];jw=np.einsum('koi,kip->kop',W,jj);ew=np.einsum('koi,ki->ko',W,Y-yy)
    A=np.einsum('kop,koq->pq',zz.conj(),jw,optimize=True).real;g=np.einsum('kop,ko->p',zz.conj(),ew).real
    Ac,gc,_,_=equations(m,cc,U,Y,Uh,bn,W,chunk=200)
    sc=1/np.sqrt(abs(np.diag(A)));rep['scaled_moment_matrix_relative']=float(np.linalg.norm((Ac-A)*sc[:,None]*sc[None,:])/np.linalg.norm(A*sc[:,None]*sc[None,:]))
    rep['scaled_moment_gradient_relative']=float(np.linalg.norm((gc-g)*sc)/np.linalg.norm(g*sc))
    dv=rng.normal(size=m.nb)*np.maximum(abs(bn),1e-8);dv[m.iw:]=rng.normal(size=6)*.2;h=2e-5
    fd=(m.output(cc,bn+h*dv,U,False)[0]-m.output(cc,bn-h*dv,U,False)[0])/(2*h)
    analytic=np.einsum('kop,p->ko',jj,dv)
    rep['directional_derivative_relative']=float(np.linalg.norm(fd-analytic)/np.linalg.norm(fd))
    G,Gj=m.gmodel(cc,bn);fdg=(m.resp(cc,bn+h*dv)-m.resp(cc,bn-h*dv))/(2*h)
    rep['transfer_derivative_relative']=float(np.linalg.norm(fdg-np.einsum('koip,p->koi',Gj,dv))/np.linalg.norm(fdg))
    lb,ub=m.bounds([8,9,13],[13,15,20]);start=bn.copy();start[m.iw:]+=rng.normal(size=6)*1e-4
    # Noiseless recovery, independent instrument and measured-input designs.
    recovered,info=iv_solve(m,cc,U,yy,Uh,None,start,lb,ub)
    rep['noiseless_max_log_modal_error']=float(np.max(abs(recovered[m.iw:]-bn[m.iw:])));rep['noiseless_converged']=info['success']
    js(HERE/'validation.json',rep);print(json.dumps(rep,indent=2))
    assert rep['output_vs_old_relative']<1e-7 and rep['transfer_vs_old_relative']<1e-7
    assert rep['jac_vs_old_forward_difference_relative']<1e-4
    assert rep['scaled_moment_matrix_relative']<1e-10 and rep['scaled_moment_gradient_relative']<1e-10
    assert rep['directional_derivative_relative']<1e-5 and rep['transfer_derivative_relative']<1e-5
    assert rep['noiseless_max_log_modal_error']<1e-9 and rep['noiseless_converged']

if __name__=='__main__':validate()

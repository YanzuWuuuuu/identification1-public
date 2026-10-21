"""Damped IV-Gauss-Newton for the SAME reference-IV moment equations.

This is not output-error least squares. The moment Jacobian includes the
derivative of the changing instrument and model-output Hessian contraction.
A candidate is accepted only if the actual frozen-scale IV moment norm decreases.
LM damping replaces pole-step truncation. A 0.25 log-pole trust radius is enforced
by increasing damping, not by clipping an accepted step. No parameter penalty.
"""
import numpy as np

def exact_moment_jacobian(beta,U,Y,model,instrument,W,iw,itr,K=None,response=None):
    """Derivative of g including the changing reference instrument.

    Only six modal coordinates are nonlinear at fixed input. Their Hessian
    columns use centered differences with outer step 1e-4; mixed linear-modal
    entries follow symmetry. duhat/db=-S K (dG/db) uhat is analytic.
    Returns minus dg/db, for compatibility with A*step=g.
    """
    A,g,score=moments(beta,U,Y,model,instrument,W)
    yh,_=model(beta,U);uh=instrument(beta);e=Y-yh
    we=e if W is None else np.einsum('kij,kj->ki',W,e)
    n=len(beta); H=np.zeros((n,n)); ids=np.arange(iw,itr)
    for j in ids:
        bp=beta.copy();bm=beta.copy();bp[j]+=1e-4;bm[j]-=1e-4
        dz=(model(bp,uh)[1]-model(bm,uh)[1])/2e-4
        H[:,j]=np.real(np.einsum('kop,ko->p',dz.conj(),we))
    linear=np.r_[0:iw,itr:n]
    H[np.ix_(ids,linear)]=H[np.ix_(linear,ids)].T
    if K is not None and response is not None:
        zero=np.zeros_like(U);_,jzero=model(beta,zero)
        Gj=np.empty((len(U),4,4,n),complex)
        for q in range(4):
            unit=zero.copy();unit[:,q]=1
            Gj[:,:,q,:]=model(beta,unit)[1]-jzero
        SK=np.linalg.solve(np.eye(4)[None]+K@response(beta),K)
        du=-np.einsum('kab,kbcp,kc->kap',SK,Gj,uh,optimize=True)
        B=np.einsum('koip,ko->kip',Gj.conj(),we,optimize=True)
        H+=np.real(np.einsum('kip,kiq->pq',B,du.conj(),optimize=True))
    else:
        raise ValueError('Reference-IV exact moment derivative needs K and response')
    return A-H,g,score

def moments(beta,U,Y,model,instrument,W):
    yh,J=model(beta,U); uh=instrument(beta); _,Z=model(beta,uh)
    e=Y-yh
    if W is None:
        ZW=Z.conj()
    else:
        ZW=np.einsum('koi,kop->kip',W,Z.conj(),optimize=True)
    A=np.real(ZW.reshape(-1,len(beta)).T@J.reshape(-1,len(beta)))
    g=np.real(np.einsum('kop,ko->p',ZW,e))
    if W is None:
        zn=np.sum(np.abs(Z)**2,axis=(0,1));en=float(np.sum(np.abs(e)**2))
    else:
        zn=np.real(np.einsum('kop,kop->p',ZW,Z));en=float(np.real(np.einsum('ki,kij,kj->',e.conj(),W,e)))
    score=float(np.max(np.abs(g)/np.maximum(np.sqrt(np.maximum(zn,0)*en),1e-300)))
    return A,g,score

def solve(beta0,U,Y,model,instrument,W,iw,itr,maxiter=80,tol=1e-5,K=None,response=None):
    beta=beta0.copy(); damping=1e-3; history=[]; accepted=0; evaluations=0
    for iteration in range(maxiter):
        A,g,score=exact_moment_jacobian(beta,U,Y,model,instrument,W,iw,itr,K,response);evaluations+=1
        scale=1/np.sqrt(np.maximum(np.abs(np.diag(A)),1e-300))
        As=A*scale[:,None]*scale[None,:];gs=g*scale
        # Diagonal two-sided scaling uses the moment Jacobian here.
        merit=float(gs@gs);normal=As.T@As;rhs=As.T@gs
        diag=np.maximum(np.diag(normal),1e-24)
        success=False
        for attempt in range(20):
            try:
                dx=np.linalg.solve(normal+damping*np.diag(diag),rhs)
                step=scale*dx
                if not np.all(np.isfinite(step)) or np.max(np.abs(step[iw:itr]))>.25:
                    damping*=10;continue
                candidate=beta+step
                Ac,gc,sc=moments(candidate,U,Y,model,instrument,W);evaluations+=1
                candidate_merit=float((gc*scale)@(gc*scale))
                if np.isfinite(candidate_merit) and candidate_merit<merit:
                    success=True;break
                damping*=10
            except (np.linalg.LinAlgError,ValueError,FloatingPointError,OverflowError):
                damping*=10
        history.append({'iteration':iteration+1,'score':score,'merit':merit,'damping':damping,
                        'accepted':success,'trials':attempt+1})
        if not success:
            # A stationary IV root can have no representable merit decrease.
            if score<1e-8:
                return beta,dict(success=iteration<79,status='pole_step_converged',iterations=iteration+1,score_max=score,accepted=accepted,evaluations=evaluations,history=history)
            return beta,dict(success=False,status='no_moment_descent',iterations=iteration+1,score_max=score,accepted=accepted,evaluations=evaluations,history=history)
        beta=candidate;accepted+=1;damping=max(damping/3,1e-14)
        # Keep A's numerical stopping criterion, and report the full score too.
        if np.max(np.abs(step[iw:itr]))<tol:
            return beta,dict(success=iteration<79,status='pole_step_converged' if iteration<79 else 'iteration_limit',iterations=iteration+1,score_max=sc,accepted=accepted,evaluations=evaluations,history=history)
    return beta,dict(success=False,status='iteration_limit',iterations=maxiter,score_max=sc,accepted=accepted,evaluations=evaluations,history=history)

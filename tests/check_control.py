#!/usr/bin/env python3
"""Recompute finite-window and bounded-control photon-receiver checks.

The no-go proof is analytical; numerical schedules are tests, not an exhaustive
control search. Reported rates are explicit coherent-information lower bounds,
not optimized quantum capacities or rates per second. All data are regenerated.
"""
from __future__ import annotations
import argparse,hashlib,json,platform,time
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import scipy
import mpmath as mp
from scipy.integrate import quad,solve_ivp
from scipy.linalg import eigh
from capture.finite_window import (FinitePulse,finite_eigenvalue,kernel,exact_cell_matrix,
    minimum_window,tail_coefficient,truncated_exponential)
from capture.capacity import boundary_metrics,threshold_r,coherent_information,coherent_information_mp
from capture.receiver import Schedule,matrix_fock,kernel_metrics,split_quad
ROOT=Path(__file__).resolve().parents[1]


def spectral_checks():
    normerr=[];eigenerr=[];odeerr=[];rows=[];gaps=[]
    for r,S in [(0.,1.),(.1,.2),(.3,.7),(.5,2.),(.735,5.),(1.,3.),(3.,.2),(3.,8.)]:
        f=FinitePulse(r,S);met=f.metrics();lam=met['largest_kernel_eigenvalue']
        normerr.append(abs(quad(lambda u:float(f(u))**2,0,S,epsabs=2e-12)[0]-1))
        for u in [0.,S/5,S/2,S]:
            val=quad(lambda v:kernel(u,v,r)*f(v),0,u,epsabs=2e-12)[0]+quad(lambda v:kernel(u,v,r)*f(v),u,S,epsabs=2e-12)[0]
            eigenerr.append(abs(val-lam*f(u)))
        # Independent backward shooting starts at the Robin endpoint. Forward
        # shooting is exponentially ill-conditioned at large r*S.
        sol=solve_ivp(lambda u,y:[y[1],(r*r-2*r*np.exp(-u)/lam)*y[0]],(S,0.),[1.,-r],method='DOP853',rtol=3e-12,atol=2e-14)
        assert sol.success
        odeerr.append(abs(sol.y[1,-1]/sol.y[0,-1]))
        last=-1
        for n in [32,64,128,256,512]:
            K,h=exact_cell_matrix(r,S,n);val=eigh(K+np.outer(h,h),subset_by_index=[n-1,n-1],eigvals_only=True)[0]
            assert val>=last-2e-10 and val<=lam+2e-10,(r,S,n,val,lam)
            gaps.append({'r':r,'S':S,'cells':n,'continuum_minus_Galerkin':float(lam-val)})
            last=val
        rows.append(met)
    assert max(normerr+eigenerr)<1e-10
    assert max(odeerr)<5e-9
    # At long windows the principal Bessel root approaches its half-line
    # limit below floating precision. Starting a bracket at that rounded
    # root used to skip the principal mode and falsely report Q=0, e.g.
    # r=.5,S=20 gave Lambda=.1407 instead of 1.17996.
    long_window=[]
    for r in [1e-6,1e-4,.01,.02,.1,.3,.5,20.]:
        limit=boundary_metrics(r)['kernel_largest_eigenvalue'];previous=0.
        for S in [2.,10.,20.,30.,40.,50.]:
            value=finite_eigenvalue(r,S)
            # Roundoff is allowed; these floats do not certify which side of
            # an unresolved strict continuum boundary a parameter lies on.
            assert previous<=value+2e-9 and value<=limit+2e-9,(r,S,value,previous,limit)
            if S>=30:assert abs(value-limit)<2e-9,(r,S,value,limit)
            long_window.append({'r':r,'S':S,'largest_kernel_eigenvalue':value,
                                'half_line_minus_finite':float(limit-value)})
            previous=value
    assert finite_eigenvalue(.5,20.)>1.
    # The Y admixture can underflow while Y itself overflows in the far
    # tail. The normalized principal pulse and its channel stay finite.
    far_tail=FinitePulse(20.,50.);far_metrics=far_tail.metrics()
    assert np.isfinite([far_tail(0.),far_tail(50.),far_metrics['population'],
                        far_metrics['coherence_amplitude']]).all()
    assert abs(quad(lambda u:float(far_tail(u))**2,0,50.,epsabs=2e-12)[0]-1)<1e-10
    return {'max_pulse_norm_error':float(max(normerr)),'max_integral_eigenfunction_error':float(max(eigenerr)),
            'max_independent_ODE_normalized_Neumann_error':float(max(odeerr)),'finite_pulses':rows,
            'nested_exact_cell_lower_bounds':gaps,'long_window_principal_mode_regression':long_window}


def controlled_dynamics_checks():
    schedules=[Schedule([0,1,2,3],[1,.25,.8],[.3,-.4,.7],.2,.4,1.),
               Schedule([0,1,2,3],[1,0,.5],[.3,-.4,.7],0.,.8,1.),
               Schedule([0,.3,1.3,2.8,4.],[.15,.95,.4,1.],[-.6,.7,.15,-.3],.25,.3,1.),
               Schedule([0,1,2],[1,1],[0,0],0.,.5,1.)]
    rows=[];errs=[];point=[]
    for j,sch in enumerate(schedules):
        a=1.3;N=np.sqrt(a/(-np.expm1(-a*sch.T)))
        pulse=lambda x:complex(N*np.exp(-a*x/2+1j*.4*x*x))
        dynamics=matrix_fock(sch,pulse);integ=kernel_metrics(sch,pulse)
        ep=abs(dynamics['p']-integ['p']);em=abs(dynamics['m']-integ['m']);errs.extend([ep,em])
        bound=sch.beta_max*finite_eigenvalue(sch.r_min,sch.Gamma_max*sch.T)
        coarse=sch.top_margin_objective(24);fine=sch.top_margin_objective(48)
        # Nyström diagonal-cusp integration can approach the exact value from
        # above, so only the distinctly nonsaturating schedules get this check.
        if j<3:assert fine<bound-1e-3
        else:assert abs(fine-bound)<1e-4
        x,w,K,h=sch.discrete(16);z=sch.exposure(x);g,G,_=sch.rates(x);mask=G>0
        xm=x[mask];zm=z[mask];bet=g[mask]/G[mask]
        Lclock=np.exp(-(zm[:,None]+zm[None,:])/2)*np.sqrt(bet[:,None]*bet[None,:])*(
                np.exp(-sch.gamma_phi*abs(xm[:,None]-xm[None,:]))+np.exp(-sch.gamma_phi*(xm[:,None]+xm[None,:])))
        ref=sch.beta_max*kernel(zm[:,None],zm[None,:],sch.r_min)
        violation=float(np.max(Lclock-ref));point.append(violation);assert violation<1e-12
        rows.append({'schedule':j,'lookback_edges':sch.edges.tolist(),'incident_rates':sch.incident.tolist(),
            'detuning_rates':sch.detuning.tolist(),'kappa':sch.kappa,'gamma_phi':sch.gamma_phi,
            'p':dynamics['p'],'m_real':dynamics['m'].real,'m_imag':dynamics['m'].imag,
            'p_kernel_vs_matrix_Fock_error':ep,'m_kernel_vs_matrix_Fock_error':em,
            'top_objective_Nystrom_24_per_segment':coarse,'top_objective_Nystrom_48_per_segment':fine,
            'exact_maximum_rate_bound':bound,'clock_pointwise_max_violation':violation})
    assert max(errs)<2e-9
    return {'max_matrix_Fock_vs_kernel_error':float(max(errs)),'max_clock_inequality_violation':float(max(point)),
            'schedules':rows,'scope':'Control optimality follows from the time-reparameterization proof, not the sampled schedules.'}


def positive_sqrt(A):
    e,V=eigh((A+A.conj().T)/2)
    if e[0]<-1e-10:raise ValueError('Matrix is not positive semidefinite.')
    return (V*np.sqrt(np.maximum(e,0.)))@V.conj().T


def multimode_kraus(A,b,broadcast=False):
    n=len(b);F=np.zeros((4 if broadcast else 2,n+1),complex);F[0,0]=1
    F[1,1:]=b.conj()
    if broadcast:F[2,1:]=b.conj()
    out=[F]
    if broadcast:components=[(3,A-np.outer(b,b.conj())),(0,np.eye(n)-A-np.outer(b,b.conj()))]
    else:components=[(1,A-np.outer(b,b.conj())),(0,np.eye(n)-A)]
    for target,C in components:
        for row in positive_sqrt(C):
            f=np.zeros_like(F);f[target,1:]=row;out.append(f)
    return out


def broadcast_checks():
    sch=Schedule([0,.7,2.],[1.,.35],[.2,-.7],.2,1.,1.)
    _,_,A,b=sch.discrete(4);E=multimode_kraus(A,b);F=multimode_kraus(A,b,True);n=len(b)+1
    tp=np.linalg.norm(sum(f.conj().T@f for f in F)-np.eye(n));errs=[]
    for i in range(n):
        for j in range(n):
            unit=np.zeros((n,n),complex);unit[i,j]=1
            ref=sum(k@unit@k.conj().T for k in E)
            joint=sum(k@unit@k.conj().T for k in F);x=joint.reshape(2,2,2,2)
            left=np.einsum('abcb->ac',x);right=np.einsum('abac->bc',x)
            errs.extend([np.linalg.norm(ref-left),np.linalg.norm(ref-right)])
    assert max(errs+[tp])<2e-11
    return {'input_dimension':n,'trace_preservation_error':float(tp),'max_marginal_error':float(max(errs)),
            'scope':'Complex-phase finite-mode test of the generic broadcast construction; the all-mode result is analytical.'}


def deadline_checks():
    rows=[];errs=[]
    for beta,rs in [(1.,[0.,.1,.3,.5,.6,.7,.72,.73,.735,.7357]),(.9,[0.,.3,.5]),(.75,[0.,.1,.3]),(.6,[0.,.1])]:
        for r in rs:
            S=minimum_window(r,beta);lam=finite_eigenvalue(r,S);errs.append(abs(beta*lam-1))
            assert beta*finite_eigenvalue(r,S+1e-4)>1
            assert beta*finite_eigenvalue(r,S-1e-4)<1
            rows.append({'beta':beta,'r':r,'infimum_dimensionless_window':S,'boundary_residual':float(beta*lam-1)})
    assert np.isinf(minimum_window(0.,.5)) and np.isinf(minimum_window(1.,1.))
    assert max(errs)<1e-10
    rc=threshold_r();h=1e-5
    deriv=(boundary_metrics(rc+h)['kernel_largest_eigenvalue']-boundary_metrics(rc-h)['kernel_largest_eigenvalue'])/(2*h)
    constant=float(np.log(tail_coefficient(rc)/(-deriv)));asym=[]
    for delta in [1e-3,1e-4,1e-5,1e-6,1e-7,1e-8]:
        S=minimum_window(rc-delta);val=(1+2*rc)*S+np.log(delta)
        asym.append({'rc_minus_r':delta,'window':S,'scaled_window_plus_log_delta':float(val)})
    assert abs(asym[-1]['scaled_window_plus_log_delta']-constant)<6e-4
    tail=[]
    for r in [.1,.5,.7]:
        li=boundary_metrics(r)['kernel_largest_eigenvalue'];theory=tail_coefficient(r)
        for S in [4.,6.,8.]:
            coefficient=(li-finite_eigenvalue(r,S))*np.exp((1+2*r)*S)
            tail.append({'r':r,'S':S,'scaled_eigenvalue_tail':float(coefficient),'asymptotic_coefficient':theory})
        assert abs(coefficient/theory-1)<1e-3
    return {'windows':rows,'max_threshold_error':float(max(errs)),'beta_one_critical_r':rc,
            'critical_T2_over_T1':float(1/(.5+rc)), 'logarithmic_window_coefficient':float(1/(1+2*rc)),
            'asymptotic_log_constant':constant,'near_boundary_windows':asym,'tail_checks':tail,
            'scope':'Window at the listed infimum still has Q=0; every strictly longer window has Q>0 under the stated model.'}


def finite_rate_checks():
    rows=[];errs=[]
    for S in [2.,3.,4.,6.]:
        a=np.sqrt(3);r=.5;met=truncated_exponential(a,r,S);f=lambda x:complex(met['normalization']*np.exp(-a*x/2))
        sch=Schedule([0.,S],[1.],[0.],0.,r,1.)
        exact=matrix_fock(sch,f);errs.extend([abs(exact['p']-met['population']),abs(exact['m']-met['coherence_amplitude'])])
        ci=coherent_information(1/8,met['population'],met['coherence_amplitude']);assert ci>0
        rows.append({**met,'input_photon_probability':1/8,'achievable_CI_bits_per_use':float(ci)})
    precision=[]
    for digits in [70,110]:
        with mp.workdps(digits):
            a=mp.sqrt(3);r=mp.mpf('.5');S=mp.mpf(4);q=mp.mpf(1)/8;N2=a/(-mp.expm1(-a*S));c=(a+1)/2
            p=2*N2/(c-r)*(-mp.expm1(-(c+r)*S)/(c+r)+mp.expm1(-2*c*S)/(2*c))
            m2=N2*(-mp.expm1(-(c+r)*S)/(c+r))**2
            ci=coherent_information_mp(q,p,m2)
            precision.append({'decimal_precision':digits,'p':mp.nstr(p,45),'m_squared':mp.nstr(m2,45),'coherent_information':mp.nstr(ci,45)})
    x=precision[0].copy();y=precision[1].copy();x.pop('decimal_precision');y.pop('decimal_precision');assert x==y
    assert max(errs)<1e-9
    return {'max_closed_form_vs_matrix_Fock_error':float(max(errs)),'finite_pulse_witnesses':rows,'precision_checks':precision,
            'scope':'Explicit diagonal-input CI lower bounds, not a globally optimized channel capacity, code, memory or rate per second.'}


def main(output):
    if not __debug__:raise RuntimeError('Keep Python assertions enabled.')
    start=time.time()
    result={'status':'PASS','scope':'Analytical finite-window and scalar-control theorem supported by independent integrations and finite approximations. No experimental data or device-noise fitting.',
        'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mp.__version__},
        'finite_spectrum':spectral_checks(),'time_dependent_controls':controlled_dynamics_checks(),
        'multimode_broadcast':broadcast_checks(),'minimum_duration':deadline_checks(),'finite_rate_witness':finite_rate_checks()}
    result['elapsed_seconds']=time.time()-start
    result['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(p for d in ('capture','tests','tools') for p in (ROOT/d).glob('*.py'))}
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':'PASS','groups':5,'elapsed_seconds':result['elapsed_seconds'],'output':str(output)},indent=2))
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=ROOT/'artifacts/control.json');a=parser.parse_args();main(a.output)

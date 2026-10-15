"""Checks for exact population/coherence optimization of a dephasing receiver.

All data are recomputed. The finite-Galerkin checks are lower bounds on the
continuum variational optimum, not replacements for its analytical proof.
No experimental data, fitted parameters, or external repository is used.
Run with assertions enabled. Approximate computational errors are reported.
"""
from __future__ import annotations
import argparse, hashlib, json, platform, time
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.linalg import eigh
from scipy.optimize import brentq, minimize_scalar
from capture.objectives import (Pulse, population_pulse, coherence_pulse, weighted_pulse,
                     kernel, integrate_metrics)
ROOT = Path(__file__).resolve().parents[1]


def kernel_apply(pulse, u):
    r=pulse.r
    fn=lambda v:np.exp(-(u+v)/2-r*abs(u-v))*pulse(v)
    return quad(fn,0,u,epsabs=2e-12,epsrel=2e-11)[0]+quad(fn,u,60,epsabs=2e-12,epsrel=2e-11)[0]


def analytic_checks():
    eigen_errors=[];stationary_errors=[];norm_errors=[];ode_errors=[];metric_rows=[]
    for r in [0.,.05,.1,.5,1.,3.,10.,50.,100.]:
        for beta in ([1.,.6] if r in [.5,3.] else [1.]):
            for pulse in [population_pulse(r,beta),coherence_pulse(r,beta),weighted_pulse(r,beta)]:
                exact=pulse.metrics();calc=integrate_metrics(pulse)
                norm_errors.append(abs(calc['pulse_norm_on_truncated_domain']-1))
                ode_errors.extend([abs(calc['population']-exact['population']),abs(calc['coherence_amplitude']-exact['coherence_amplitude'])])
                assert -1e-10<=exact['cp_slack_p_minus_m2'] and 0<=exact['population']<=beta+1e-10
                if beta==1. and pulse.kind!='coherence':metric_rows.append(exact)
    for r in [.1,.5,3.]:
        pp=population_pulse(r);mu=pp.metrics()['population']
        for u in [0.,.07,.5,2.,6.]:eigen_errors.append(abs(kernel_apply(pp,u)-mu*pp(u)))
        for beta,w in [(1.,1.),(.6,.2),(.6,5.)]:
            p=weighted_pulse(r,beta,w);z=8*r*beta/p.x**2
            for u in [0.,.07,.5,2.,6.]:
                stationary_errors.append(abs(z*p(u)-beta*kernel_apply(p,u)-w*np.sqrt(beta)*np.exp(-(r+.5)*u)))
    gram_errors=[]
    for r in [.1,1.,3.]:
        for u,v in [(0.,1.),(.2,.8),(2.,.3),(1.,1.)]:
            lhs=kernel(u,v,r)-np.exp(-(r+.5)*(u+v))
            rhs=2*r*np.exp(-(u+v)/2)*quad(lambda s:np.exp(-r*(u+v-2*s)),0,min(u,v),epsabs=1e-13)[0]
            gram_errors.append(abs(lhs-rhs))
    assert max(norm_errors)<2e-9
    assert max(ode_errors)<2e-9
    assert max(eigen_errors+stationary_errors+gram_errors)<2e-9
    return {'max_norm_error':float(max(norm_errors)),'max_scalar_Fock_ODE_metric_error':float(max(ode_errors)),
            'max_Bessel_eigenfunction_error':float(max(eigen_errors)),'max_variational_equation_error':float(max(stationary_errors)),
            'max_CP_Gram_identity_error':float(max(gram_errors)),'metrics':metric_rows}


def matrix_fock(pulse):
    """Independently arranged four-block Fock master equation (2x2 matrices).

    The input convention gives rho10=-m|e><g|. A known final phase flip is
    applied to convert it to the desired optical-to-atomic logical basis.
    """
    sm=np.array([[0.,1.],[0.,0.]],complex);pe=np.diag([0.,1.]).astype(complex)
    lin=np.sqrt(pulse.beta)*sm
    def dissip(l,x):
        n=l.conj().T@l
        return l@x@l.conj().T-(n@x+x@n)/2
    def L(x):return dissip(sm,x)+2*pulse.r*dissip(pe,x)
    def comm(a,b):return a@b-b@a
    def rhs(t,y):
        q00,q10,q01,q11=y.reshape(4,2,2);xi=pulse(-t)
        return np.stack([L(q00),L(q10)+xi*comm(q00,lin.conj().T),
                         L(q01)+xi*comm(lin,q00),
                         L(q11)+xi*comm(q01,lin.conj().T)+xi*comm(lin,q10)]).ravel()
    gg=np.diag([1.,0.]).astype(complex);y0=np.stack([gg,np.zeros((2,2)),np.zeros((2,2)),gg]).ravel()
    sol=solve_ivp(rhs,[-40.,0.],y0,method='DOP853',rtol=2e-11,atol=1e-13)
    if not sol.success:raise RuntimeError(sol.message)
    q00,q10,q01,q11=sol.y[:,-1].reshape(4,2,2)
    phase=np.diag([1.,-1.]);maps=[phase@x@phase for x in [q00,q10,q01,q11]]
    q00,q10,q01,q11=maps
    # Reference-system order: |0,g>,|0,e>,|1,g>,|1,e>.
    choi=np.block([[q00,q01],[q10,q11]])/2
    bell=np.array([1.,0.,0.,1.])/np.sqrt(2)
    fe=float(np.vdot(bell,choi@bell).real)
    states=[np.array([1.,0.]),np.array([0.,1.]),np.array([1.,1.])/np.sqrt(2),np.array([1.,-1.])/np.sqrt(2),np.array([1.,1j])/np.sqrt(2),np.array([1.,-1j])/np.sqrt(2)]
    vals=[]
    for psi in states:
        d=np.outer(psi,psi.conj());out=d[0,0]*q00+d[1,0]*q10+d[0,1]*q01+d[1,1]*q11
        vals.append(float(np.vdot(psi,out@psi).real))
    p=float(q11[1,1].real);m=float(q10[1,0].real)
    return {'population':p,'coherence_amplitude':m,'entanglement_fidelity':fe,
            'six_state_average_fidelity':float(np.mean(vals)), 'choi_min_eigenvalue':float(np.linalg.eigvalsh(choi)[0]),
            'trace_error':float(max(abs(np.trace(q00)-1),abs(np.trace(q11)-1))),
            'off_diagonal_trace_error':float(max(abs(np.trace(q10)),abs(np.trace(q01))))}


def fock_checks():
    rows=[]
    for r,beta in [(0.,1.),(.5,1.),(3.,1.),(3.,.6)]:
        for pulse in [population_pulse(r,beta),weighted_pulse(r,beta)]:
            calc=matrix_fock(pulse);a=pulse.metrics()
            err=max(abs(calc['population']-a['population']),abs(calc['coherence_amplitude']-a['coherence_amplitude']),abs(calc['entanglement_fidelity']-a['entanglement_fidelity']),abs(calc['six_state_average_fidelity']-a['average_qubit_fidelity']))
            assert err<2e-9 and calc['choi_min_eigenvalue']>-1e-10
            rows.append({'r':r,'beta':beta,'kind':pulse.kind,'max_analytical_difference':float(err),**calc})
    return rows


def power_integral(a,b,exponent):
    if abs(exponent+1)<1e-12:
        if a==0:raise ValueError('Divergent integral')
        return np.log(b/a)
    return (b**(exponent+1)-a**(exponent+1))/(exponent+1)


def galerkin(r,n):
    """Exact cell integrals on q=exp(-u) in [0,1], not sampled diagonal kernels.

    F(q)=f(-log(q))/sqrt(q); K(q,q')=(min/max)^r, h(q)=q^r.
    The n equal-bin normalized indicator functions give nested subspaces.
    """
    widths=1/n;grid=np.linspace(0.,1.,n+1)
    aa=grid[:-1];bb=grid[1:]
    pos=(bb**(r+1)-aa**(r+1))/(r+1)
    neg=np.zeros(n)
    for j in range(1,n):neg[j]=power_integral(aa[j],bb[j],-r)
    M=np.zeros((n,n))
    for i in range(n):
        M[i,i]=2/(widths*(r+1))*((bb[i]**2-aa[i]**2)/2-(aa[i]**(r+1)*neg[i] if i else 0.))
        if i+1<n:
            M[i,i+1:]=pos[i]*neg[i+1:]/widths
            M[i+1:,i]=M[i,i+1:]
    return M,pos/np.sqrt(widths)


def galerkin_checks():
    rows=[];certificate_errors=[];rng=np.random.default_rng(20261015)
    for r,beta in [(.5,1.),(3.,1.),(3.,.6)]:
        prev_p=prev_f=-1.
        exact_p=population_pulse(r,beta).metrics()['population']
        exact_f=weighted_pulse(r,beta).metrics()['average_qubit_fidelity']
        for n in [64,128,256,512]:
            M,h=galerkin(r,n);ev,V=eigh(beta*M);b=np.sqrt(beta)*h;bb=V.T@b
            lo=ev[-1]+max(1e-13,abs(bb[-1])*1e-7);hi=ev[-1]+np.linalg.norm(bb)+1
            z=brentq(lambda z:np.sum((bb/(z-ev))**2)-1,lo,hi,xtol=2e-14)
            f=V@(bb/(z-ev));p=float(f@(beta*M)@f);m=float(np.sqrt(beta)*h@f);fav=(3+p+2*m)/6
            assert ev[-1]>=prev_p-1e-11 and fav>=prev_f-1e-11
            assert ev[-1]<=exact_p+1e-10 and fav<=exact_f+1e-10
            assert np.linalg.norm((z*np.eye(n)-beta*M)@f-b)<2e-11
            for k in range(3):
                test=rng.normal(size=n);test/=np.linalg.norm(test)
                gap=p+2*m-(test@(beta*M)@test+2*b@test)
                d=f-test;cert=d@(z*np.eye(n)-beta*M)@d
                certificate_errors.append(abs(gap-cert));assert gap>=-1e-12
            rows.append({'r':r,'beta':beta,'bins':n,'population_lower_bound':float(ev[-1]),'population_continuum_gap':float(exact_p-ev[-1]),'fidelity_lower_bound':float(fav),'fidelity_continuum_gap':float(exact_f-fav)})
            prev_p,prev_f=ev[-1],fav
    assert max(certificate_errors)<2e-10
    return {'nested_finite_subspace_bounds':rows,'max_global_certificate_identity_error':float(max(certificate_errors))}


def exponential_metrics(a,r,beta=1.):
    p=4*beta*a/((1+a)*(1+a+2*r));m=2*np.sqrt(beta*a)/(1+a+2*r)
    return {'exponential_rate_over_Gamma':float(a),'population':float(p),'coherence_amplitude':float(m),'average_qubit_fidelity':float((3+p+2*m)/6)}


def comparisons():
    rows=[]
    for r in [.1,.5,1.,3.,10.,100.]:
        result=minimize_scalar(lambda z:-exponential_metrics(np.exp(z),r)['average_qubit_fidelity'],bounds=(-8.,np.log(20+20*r)),method='bounded',options={'xatol':1e-12})
        a=float(np.exp(result.x));m=exponential_metrics(a,r);exact=weighted_pulse(r).metrics()
        assert exact['average_qubit_fidelity']>=m['average_qubit_fidelity']-1e-11
        rows.append({'r':r,**m,'gap_to_unrestricted_optimum':float(exact['average_qubit_fidelity']-m['average_qubit_fidelity'])})
    pop_threshold=brentq(lambda r:population_pulse(r).metrics()['average_qubit_fidelity']-2/3,2.,4.,xtol=1e-12)
    f_threshold=brentq(lambda r:weighted_pulse(r).metrics()['average_qubit_fidelity']-2/3,2.,4.,xtol=1e-12)
    # Coherence-only pulse supplies a lower bound on the optimal fidelity;
    # combining two independent maxima supplies an upper bound.
    asym=[]
    for r in [10.,50.,100.]:
        pmax=population_pulse(r).metrics()['population'];mmax=coherence_pulse(r).metrics()['coherence_amplitude']
        low=coherence_pulse(r).metrics()['average_qubit_fidelity'];opt=weighted_pulse(r).metrics()['average_qubit_fidelity'];up=(3+pmax+2*mmax)/6
        assert low<=opt+1e-11 and opt<=up+1e-11
        asym.append({'r':r,'lower_bound_fidelity':low,'exact_optimum_fidelity':opt,'upper_bound_fidelity':up,'leading_strong_noise_coefficient':1/(3*np.sqrt(2*r)),'exact_advantage_above_reset_baseline':opt-.5})
    return {'best_exponential_comparator':rows,'two_thirds_crossing_r_population_objective':float(pop_threshold),'two_thirds_crossing_r_fidelity_objective':float(f_threshold),'asymptotic_sandwich_checks':asym,
            'benchmark_scope':'Uniform unknown single-rail qubit, deterministic channel; below 2/3 is not a proof of entanglement breaking or zero quantum capacity.'}


def main(output):
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    start=time.time()
    data={'status':'PASS','date':'2026-10-15','scope':'Fixed Markov two-level receiver, single-rail qubit, unconditioned state at one capture time. Analytical identities, independent ODEs, and nested finite-space optimization. No device or experimental model fit.',
          'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
          'continuum_and_scalar_ODE':analytic_checks(),'independent_matrix_Fock_equations':fock_checks(),
          'independent_global_optimization':galerkin_checks(),'strong_comparators':comparisons()}
    data['elapsed_seconds']=time.time()-start
    data['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(p for d in ('capture','tests','tools') for p in (ROOT/d).glob('*.py'))}
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'status':data['status'],'groups':4,'elapsed_seconds':data['elapsed_seconds'],'output':str(output)},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=ROOT/'artifacts/objectives.json')
    main(parser.parse_args().output)

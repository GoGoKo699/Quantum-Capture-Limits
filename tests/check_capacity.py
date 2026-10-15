#!/usr/bin/env python3
"""Recompute checks for the photon receiver's zero/positive capacity boundary.

Exact symbolic channel identities, independent continuum/Fock checks, nested
finite-subspace optima, and high-precision achievable-rate witnesses are kept
separate. A positive coherent information is a lower bound, not the capacity.
No experiments or device error budgets are represented. Keep assertions on.
"""
from __future__ import annotations
import argparse, hashlib, json, platform, time
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import scipy
import sympy as sp
import mpmath as mp
from scipy.integrate import quad, solve_ivp
from scipy.linalg import eigh
from scipy.optimize import brentq, minimize_scalar
from capture.objectives import population_pulse, weighted_pulse, kernel, integrate_metrics
from capture.capacity import (boundary_pulse, boundary_metrics, threshold_r,
    exponential_metrics, best_exponential, kraus, binary_broadcast_kraus,
    coherent_information, coherent_information_mp)
ROOT=Path(__file__).resolve().parents[1]


def channel(K,rho):return sum((x@rho@x.conj().T for x in K))

def partials(rho):
    x=rho.reshape(2,2,2,2)
    return np.einsum('abcb->ac',x),np.einsum('abac->bc',x)


def entropy(vals):
    a=np.real(np.asarray(vals));a=a[a>1e-15]
    return float(-a@np.log2(a))


def algebra_checks():
    p,s=sp.symbols('p s',real=True)
    J=sp.Matrix([[1,0,0,sp.sqrt(s)],[0,0,0,0],[0,0,1-p,0],[sp.sqrt(s),0,0,p]])/2
    B=sp.diag(1-p/2,p/2)
    gap=sp.simplify(sp.trace(B*B)-sp.trace(J*J))
    assert gap==(1-p-s)/2 and sp.det(J)==0
    a,c=sp.symbols('a c',positive=True)
    n=a*(2*a+c+1);d=(a+1)*(a+c)**2
    assert sp.expand(sp.diff(n,a)*d-n*sp.diff(d,a)+(a+c)*(2*a**3+2*a**2+(1-3*c)*a-c*c-c))==0
    errors=[];tp=[];purity_errors=[];ci_errors=[]
    for p0,m in [(.1,0.),(.5,np.sqrt(.5)),(.6,.6),(.8,.3),(1.,0.),(.02,.1)]:
        kk=kraus(p0,m);bb=binary_broadcast_kraus(p0,m)
        tp.append(np.linalg.norm(sum(x.T@x for x in bb)-np.eye(2)))
        for i in range(2):
            for j in range(2):
                unit=np.zeros((2,2),complex);unit[i,j]=1
                aa,ab=partials(channel(bb,unit));ref=channel(kk,unit)
                errors.extend([np.linalg.norm(aa-ref),np.linalg.norm(ab-ref)])
        choi=np.zeros((4,4),complex)
        for i in range(2):
            for j in range(2):
                unit=np.zeros((2,2));unit[i,j]=1
                choi[2*i:2*i+2,2*j:2*j+2]=channel(kk,unit)/2
        out=channel(kk,np.eye(2)/2)
        purity_errors.append(abs(np.trace(out@out)-np.trace(choi@choi)-(1-p0-m*m)/2))
    for p0,m in [(.6,.6),(.8,.85),(.9,.8),(.99,.99),(.2,.4)]:
        kk=kraus(p0,m)
        for q in [.01,.1,.3,.8]:
            rho=np.diag([1-q,q]);out=channel(kk,rho)
            env=np.array([[np.trace(x@rho@y.T) for y in kk] for x in kk])
            direct=entropy(np.linalg.eigvalsh(out))-entropy(np.linalg.eigvalsh(env))
            ci_errors.append(abs(direct-coherent_information(q,p0,m)))
    assert max(errors+tp+ci_errors+purity_errors)<1e-12
    return {'symbolic_symmetric_extension_gap':str(gap),
            'symbolic_exponential_stationary_polynomial':'2a^3+2a^2+(1-3c)a-c^2-c',
            'max_binary_broadcast_marginal_error':float(max(errors)),
            'max_broadcast_trace_preserving_error':float(max(tp)),
            'max_choi_purity_check_error':float(max(purity_errors)),
            'max_independent_complement_entropy_error':float(max(ci_errors))}


def matrix_fock(pulse):
    sm=np.array([[0.,1.],[0.,0.]],complex);pe=np.diag([0.,1.]).astype(complex)
    lin=np.sqrt(pulse.beta)*sm
    def dis(l,x):
        n=l.conj().T@l
        return l@x@l.conj().T-(n@x+x@n)/2
    def L(x):return dis(sm,x)+2*pulse.r*dis(pe,x)
    def comm(a,b):return a@b-b@a
    def rhs(t,y):
        q00,q10,q01,q11=y.reshape(4,2,2);xi=pulse(-t)
        return np.stack([L(q00),L(q10)+xi*comm(q00,lin.conj().T),
            L(q01)+xi*comm(lin,q00),L(q11)+xi*comm(q01,lin.conj().T)+xi*comm(lin,q10)]).ravel()
    gg=np.diag([1.,0.]).astype(complex)
    y0=np.stack([gg,np.zeros((2,2)),np.zeros((2,2)),gg]).ravel()
    sol=solve_ivp(rhs,[-45.,0.],y0,method='DOP853',rtol=1e-11,atol=1e-13)
    assert sol.success,sol.message
    q00,q10,q01,q11=sol.y[:,-1].reshape(4,2,2)
    return float(q11[1,1].real),float(-q10[1,0].real)


def continuum_checks():
    eigen=[];norm=[];metrics=[];odes=[];matrix=[];rows=[]
    for r in [0.,.02,.1,.5,.72,1.,3.,10.]:
        for beta in ([1.,.6] if r in [.5,.72] else [1.]):
            f=boundary_pulse(r,beta);b=boundary_metrics(r,beta);old=f.metrics()
            norm.append(abs(quad(lambda u:f(u)**2,0,60,epsabs=1e-12)[0]-1))
            m=np.sqrt(beta)*quad(lambda u:np.exp(-(r+.5)*u)*f(u),0,60,epsabs=1e-12)[0]
            metrics.extend([abs(m-b['coherence_amplitude']),abs(old['population']-b['population'])])
            ode=integrate_metrics(f,45.)
            odes.extend([abs(ode['population']-b['population']),abs(ode['coherence_amplitude']-b['coherence_amplitude'])])
            if r in [.1,.72,3.]:
                pp,mm=matrix_fock(f)
                matrix.extend([abs(pp-b['population']),abs(mm-b['coherence_amplitude'])])
                for u in [0.,.2,1.,3.,6.]:
                    z=quad(lambda v:kernel(u,v,r)*f(v),0,u,epsabs=1e-12)[0]+quad(lambda v:kernel(u,v,r)*f(v),u,60,epsabs=1e-12)[0]
                    z+=np.exp(-(r+.5)*u)*m/np.sqrt(beta)
                    eigen.append(abs(z-b['kernel_largest_eigenvalue']*f(u)))
            rows.append(b)
    assert max(eigen+norm+metrics+odes+matrix)<3e-9
    return {'max_norm_error':float(max(norm)),'max_eigenfunction_integral_error':float(max(eigen)),
            'max_independent_overlap_error':float(max(metrics)), 'max_scalar_Fock_error':float(max(odes)),
            'max_matrix_Fock_error':float(max(matrix)), 'boundary_pulse_metrics':rows}


def power_int(a,b,s):
    return np.log(b/a) if abs(s+1)<1e-12 else (b**(s+1)-a**(s+1))/(s+1)


def galerkin(r,n):
    """Exact double cell integrals in y=e^-u, orthonormal indicators in L²[0,1]."""
    grid=np.linspace(0,1,n+1);a=grid[:-1];b=grid[1:];w=1/n
    v=(b**(r+1)-a**(r+1))/(r+1);neg=np.zeros(n)
    for j in range(1,n):neg[j]=power_int(a[j],b[j],-r)
    M=np.zeros((n,n))
    for i in range(n):
        M[i,i]=2/(w*(r+1))*((b[i]**2-a[i]**2)/2-(a[i]**(r+1)*neg[i] if i else 0.))
        M[i,i+1:]=v[i]*neg[i+1:]/w;M[i+1:,i]=M[i,i+1:]
    return M,v/np.sqrt(w)


def psqrt(A):
    ev,V=eigh(A)
    if ev[0]<-2e-11:raise ValueError('Not positive semidefinite.')
    return (V*np.sqrt(np.maximum(ev,0)))@V.T


def multimode_kraus(A,b,broadcast=False):
    n=len(b);bb=np.outer(b,b);basis=np.eye(n)
    if not broadcast:
        F=np.zeros((2,n+1));F[0,0]=1;F[1,1:]=b
        out=[F]
        for target,C in [(0,np.eye(n)-A),(1,A-bb)]:
            for row in psqrt(C):
                f=np.zeros((2,n+1));f[target,1:]=row;out.append(f)
    else:
        F=np.zeros((4,n+1));F[0,0]=1;F[1,1:]=b;F[2,1:]=b
        out=[F]
        for target,C in [(3,A-bb),(0,np.eye(n)-A-bb)]:
            for row in psqrt(C):
                f=np.zeros((4,n+1));f[target,1:]=row;out.append(f)
    return out


def optimization_and_multimode_checks():
    rows=[];certs=[];rng=np.random.default_rng(20261015)
    for r in [.1,.5,.72,1.,3.]:
        exact=boundary_metrics(r)['kernel_largest_eigenvalue'];last=-1
        for n in [32,64,128,256,512]:
            K,h=galerkin(r,n);T=K+np.outer(h,h)
            ev,V=eigh(T);lam=ev[-1];f=V[:,-1]
            assert last<=lam+1e-12 and lam<=exact+1e-11
            for _ in range(2):
                x=rng.normal(size=n);x/=np.linalg.norm(x);d=x-f
                certs.append(abs((lam-x@T@x)-d@(lam*np.eye(n)-T)@d))
            rows.append({'r':r,'basis_size':n,'largest_eigenvalue_lower_bound':float(lam),
                         'continuum_gap':float(exact-lam)})
            last=lam
    errors=[];tp=[]
    for r,beta in [(0.,.5),(.5,.6),(1.,1.),(3.,1.)]:
        n=7;K,h=galerkin(r,n);A=beta*K;b=np.sqrt(beta)*h
        assert eigh(A+np.outer(b,b),eigvals_only=True)[-1]<1+1e-12
        KK=multimode_kraus(A,b);FF=multimode_kraus(A,b,True)
        tp.append(np.linalg.norm(sum(f.T@f for f in FF)-np.eye(n+1)))
        for i in range(n+1):
            for j in range(n+1):
                rho=np.zeros((n+1,n+1),complex);rho[i,j]=1
                ref=channel(KK,rho);x,y=partials(channel(FF,rho))
                errors.extend([np.linalg.norm(x-ref),np.linalg.norm(y-ref)])
    assert max(errors+tp+certs)<2e-10
    return {'nested_eigenvalue_lower_bounds':rows,'max_Rayleigh_gap_identity_error':float(max(certs)),
            'multimode_input_dimension':8,'max_multimode_broadcast_marginal_error':float(max(errors)),
            'max_multimode_trace_preserving_error':float(max(tp)),
            'scope':'Finite mode checks support the analytical all-mode operator construction; finite truncation is not the proof.'}


def thresholds_and_competitors():
    popcrit=brentq(lambda r:(lambda x:x['population']+x['coherence_amplitude']**2-1)(population_pulse(r).metrics()),.1,1.,xtol=1e-12)
    expcrit=brentq(lambda r:best_exponential(r)['p_plus_m2']-1,.1,1.,xtol=1e-12)
    qcrit=threshold_r(1.)
    fidelitycrit=brentq(lambda r:(lambda x:x['population']+x['coherence_amplitude']**2-1)(weighted_pulse(r).metrics()),.1,1.,xtol=1e-12)
    bp=[]
    for beta in [.5,.55,.6,.75,.9,1.]:bp.append({'beta':beta,'maximum_r_for_positive_capacity':threshold_r(beta)})
    rows=[];errors=[]
    for r in [.1,.3,.5,.72,1.,3.,10.]:
        exp=best_exponential(r);ref=minimize_scalar(lambda z:-exponential_metrics(np.exp(z),r)['p_plus_m2'],bounds=(-6,np.log(30+20*r)),method='bounded',options={'xatol':1e-12})
        errors.append(abs(-ref.fun-exp['p_plus_m2']))
        assert exp['p_plus_m2']<=boundary_metrics(r)['p_plus_m2']+1e-12
        rows.append({'r':r,'best_exponential':exp,'population_optimal':population_pulse(r).metrics(),
                     'fidelity_optimal':weighted_pulse(r).metrics(),'boundary_optimal':boundary_metrics(r)})
    assert max(errors)<1e-10 and popcrit<expcrit<qcrit
    return {'beta1_threshold_r':{'population_objective':popcrit,'fidelity_objective':fidelitycrit,
            'best_exponential_boundary_objective':expcrit,'unrestricted_boundary_objective':qcrit},
            'incident_fraction_thresholds':bp,'max_independent_exponential_objective_error':float(max(errors)),
            'comparisons':rows,
            'scope':'Zero/positive capacity is exact in the model. No quantum-capacity value is optimized here.'}


def high_precision_witnesses():
    outputs=[]
    for digits in [90,140]:
        with mp.workdps(digits):
            r=mp.mpf('0.72');p=mp.mpf(200)/333;s=mp.mpf(5000)/12321;q=mp.mpf('1e-26')
            ciexp=coherent_information_mp(q,p,s)
            assert ciexp>0 and abs(p+s-1-mp.mpf(79)/12321)<mp.mpf(10)**(-digits+5)
            nu=2*r;seed=boundary_pulse(.72).x
            x=mp.findroot(lambda x:(mp.besselj(nu-1,x)-mp.besselj(nu+1,x))/2,(seed*.99,seed*1.01))
            lam=8*r/x**2;s2=lam*lam/(4-2*r*lam);p2=lam-s2;delta=p2+s2-1
            C=-p2*mp.log(p2)+(1-p2)*mp.log(1-p2)+(p2-s2)*mp.log(p2-s2)
            q2=mp.exp(C/delta);ci=coherent_information_mp(q2,p2,s2)
            assert ci>0
            # Eigenvalue and threshold root at high precision, independent of SciPy root refinement.
            rr=mp.findroot(lambda rr:(mp.besselj(2*rr-1,mp.sqrt(8*rr))-mp.besselj(2*rr+1,mp.sqrt(8*rr)))/2,(mp.mpf('.73'),mp.mpf('.74')))
            outputs.append({'decimal_precision':digits,'exact_exponential_example':{
                'r':'18/25','a':'2','population':'200/333','coherence_squared':'5000/12321','margin':'79/12321',
                'input_one_photon_probability':mp.nstr(q,32),'coherent_information_lower_bound_bits_per_use':mp.nstr(ciexp,32)},
                'boundary_pulse_example':{'r':'0.72','p':mp.nstr(p2,32),'m_squared':mp.nstr(s2,32),
                'margin':mp.nstr(delta,32),'input_one_photon_probability':mp.nstr(q2,32),
                'coherent_information_lower_bound_bits_per_use':mp.nstr(ci,32)},
                'unrestricted_threshold_r':mp.nstr(rr,32)})
    a=outputs[0].copy();b=outputs[1].copy();a.pop('decimal_precision');b.pop('decimal_precision');assert a==b
    # Small-q expansion coefficient tested on channels on both sides of the boundary.
    asym=[]
    with mp.workdps(90):
        for pp,ss in [('0.8','0.7'),('0.6','0.4'),('0.6','0.3')]:
            p=mp.mpf(pp);s=mp.mpf(ss);q=mp.mpf('1e-25');D=p+s-1
            C=-p*mp.log(p)+(1-p)*mp.log(1-p)+(p-s)*mp.log(p-s)
            approx=q*(D*mp.log(1/q)+D+C)/mp.log(2)
            exact=coherent_information_mp(q,p,s)
            error=abs(exact-approx)/q
            assert error<mp.mpf('1e-20')
            asym.append({'p':pp,'m_squared':ss,'error_divided_by_q':mp.nstr(error,16)})
    return {'repeated_precision':outputs,'small_input_expansion':asym,
            'scope':'Explicit achievable-rate lower bounds only. These small numbers are not upper bounds on multi-use capacity or finite-block performance.'}


def main(output):
    if not __debug__:raise RuntimeError('Keep assertions enabled.')
    start=time.time();result={'status':'PASS','date':'2026-10-15',
        'scope':'Fixed Markov receiver, arbitrary temporal states in vacuum plus one-photon sector, independent channel uses. Exact boundary, not a capacity formula or physical memory-rate prediction.',
        'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__,'mpmath':mp.__version__},
        'channel_algebra':algebra_checks(),'continuum_and_dynamics':continuum_checks(),
        'finite_space_and_multimode':optimization_and_multimode_checks(),
        'thresholds_and_comparators':thresholds_and_competitors(),'high_precision_witnesses':high_precision_witnesses()}
    result['elapsed_seconds']=time.time()-start
    result['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(p for d in ('capture','tests','tools') for p in (ROOT/d).glob('*.py'))}
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','groups':5,'elapsed_seconds':result['elapsed_seconds'],'output':str(output)},indent=2))
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=ROOT/'artifacts/capacity.json');args=ap.parse_args();main(args.output)

"""Finite-duration noisy photon capture and its rate-limited control boundary.

r=gamma_phi/Gamma_max and S=Gamma_max*T. gamma_phi is EXTRA optical
coherence decay, not population decay. beta is the maximum incident fraction.
Only number-preserving coupling/phase controls and at most one photon per use
are covered. The optimized eigenvalue decides zero vs positive unassisted
capacity; it does not optimize the capacity value or a rate per second.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import jv, yv, jvp, yvp, gamma
from scipy.optimize import brentq
from scipy.integrate import quad
from scipy.linalg import eigh
from capture.capacity import boundary_metrics, first_derivative_zero, threshold_r


def validate(r:float,S:float,beta:float=1.) -> None:
    if not np.isfinite([r,S,beta]).all() or r<0 or S<=0 or not 0<beta<=1:
        raise ValueError('Require finite r>=0, S>0 and 0<beta<=1.')


def kernel(u,v,r:float):
    u=np.asarray(u);v=np.asarray(v)
    return np.exp(-(u+v)/2)*(np.exp(-r*abs(u-v))+np.exp(-r*(u+v)))


def boundary_determinant(k:float,r:float,S:float) -> float:
    nu=2*r;y=k*np.exp(-S/2)
    return float(jvp(nu,k)*yv(nu+1,y)-yvp(nu,k)*jv(nu+1,y))


def finite_root(r:float,S:float) -> float:
    validate(r,S)
    if r==0:raise ValueError('The r=0 case is rank one; no root is needed.')
    if r>20 or S>50 or S<1e-4 or r<1e-6:
        raise ValueError('Floating-point Bessel solver tested only for 1e-6<=r<=20, 1e-4<=S<=50.')
    lo=first_derivative_zero(2*r)
    f0=boundary_determinant(lo,r,S)
    if abs(f0)<1e-300:raise ArithmeticError('Underflow in boundary determinant.')
    step=.05;hi=lo+step
    limit=lo+100+10*np.sqrt(r/S)
    while np.sign(boundary_determinant(hi,r,S))==np.sign(f0):
        hi+=step
        if hi>limit:raise ArithmeticError('First finite-window root not bracketed.')
    return float(brentq(lambda x:boundary_determinant(x,r,S),hi-step,hi,xtol=2e-13,rtol=1e-14))


def finite_eigenvalue(r:float,S:float) -> float:
    validate(r,S)
    if r==0:return float(-2*np.expm1(-S))
    k=finite_root(r,S)
    return float(8*r/(k*k))


@dataclass
class FinitePulse:
    r:float
    S:float
    beta:float=1.
    def __post_init__(self):
        validate(self.r,self.S,self.beta)
        self.k=None if self.r==0 else finite_root(self.r,self.S)
        if self.r==0:
            self.coefficient=0.;self.norm=np.sqrt(-np.expm1(-self.S))
        else:
            y=self.k*np.exp(-self.S/2);nu=2*self.r
            self.coefficient=-jv(nu+1,y)/yv(nu+1,y)
            norm2=quad(lambda u:float(self.raw(u))**2,0,self.S,epsabs=3e-13,epsrel=3e-13)[0]
            self.norm=np.sqrt(norm2)
        if self.raw(0)<0:self.norm=-self.norm
    def raw(self,u):
        u=np.asarray(u,dtype=float)
        if self.r==0:return np.exp(-u/2)
        x=self.k*np.exp(-u/2)
        return np.exp(-u/2)*(jv(2*self.r,x)+self.coefficient*yv(2*self.r,x))
    def __call__(self,u):
        x=np.asarray(u,dtype=float)
        return np.where((x>=0)&(x<=self.S),self.raw(x)/self.norm,0.)
    def metrics(self):
        lam=finite_eigenvalue(self.r,self.S)
        m=np.sqrt(self.beta)*quad(lambda u:np.exp(-(self.r+.5)*u)*float(self(u)),0,self.S,epsabs=2e-12,epsrel=2e-12)[0]
        p=self.beta*lam-m*m
        return {'r':self.r,'S':self.S,'beta':self.beta,'largest_kernel_eigenvalue':lam,'population':float(p),
                'coherence_amplitude':float(m),'boundary_margin':float(self.beta*lam-1),
                'root':self.k,'Y_coefficient':float(self.coefficient)}


def minimum_window(r:float,beta:float=1.) -> float:
    if not np.isfinite([r,beta]).all() or r<0 or not 0<beta<=1:
        raise ValueError('Require r>=0, 0<beta<=1.')
    if beta<=.5:return float('inf')
    if r==0:return float(np.log(2*beta/(2*beta-1)))
    if beta*boundary_metrics(r)['kernel_largest_eigenvalue']<=1:return float('inf')
    nu=2*r;k=np.sqrt(8*beta*r)
    # This branch stays before zeros of Y_(nu+1), for the below-threshold
    # receiver parameters. Dividing avoids large Y values at long windows.
    def residual(S):
        y=k*np.exp(-S/2)
        return jvp(nu,k)-yvp(nu,k)*jv(nu+1,y)/yv(nu+1,y)
    hi=1.
    while residual(hi)>0:
        hi*=2
        if hi>100:raise ArithmeticError('Threshold window not bracketed.')
    return float(brentq(residual,0.,hi,xtol=2e-12,rtol=1e-13))


def tail_coefficient(r:float) -> float:
    """A(r) in Lambda_infty-Lambda_S ~ A exp(-(1+2r)S)."""
    if r==0:return 2.
    b=boundary_metrics(r);nu=2*r;j=b['derivative_zero'];lam=b['kernel_largest_eigenvalue']
    D=jv(nu,j)**2-jv(nu-1,j)*jv(nu+1,j)
    c2=(j/2)**(2*nu)/(gamma(nu+1)**2*D)
    return float(lam*c2/(nu+1))


def exact_cell_matrix(r:float,S:float,n:int):
    """Exact cell integrals in y=exp(-u), a nested piecewise-constant subspace.

    Eigenvalues are Galerkin lower bounds, not exact continuum eigenvalues.
    """
    validate(r,S)
    if n<2:raise ValueError('n>=2 is required.')
    edges=np.linspace(np.exp(-S),1.,n+1);a=edges[:-1];b=edges[1:];w=b-a
    plus=(b**(r+1)-a**(r+1))/(r+1)
    minus=np.log(b/a) if abs(r-1)<1e-10 else (b**(1-r)-a**(1-r))/(1-r)
    mat=np.zeros((n,n))
    for i in range(n):
        mat[i,i]=2/(r+1)*((b[i]**2-a[i]**2)/2-a[i]**(r+1)*minus[i])/w[i]
        mat[i,i+1:]=plus[i]*minus[i+1:]/np.sqrt(w[i]*w[i+1:]);mat[i+1:,i]=mat[i,i+1:]
    h=plus/np.sqrt(w)
    return mat,h


def truncated_exponential(a:float,r:float,S:float,beta:float=1.):
    validate(r,S,beta)
    if not np.isfinite(a) or a<=0:raise ValueError('a>0 required.')
    N2=a/(-np.expm1(-a*S));c=(a+1)/2
    if abs(c-r)<1e-9:
        integral=(1-np.exp(-2*r*S)*(1+2*r*S))/(2*r)**2
        p=2*beta*N2*integral
    else:
        p=2*beta*N2/(c-r)*(-np.expm1(-(c+r)*S)/(c+r)+np.expm1(-2*c*S)/(2*c))
    m=np.sqrt(beta*N2)*(-np.expm1(-(c+r)*S))/(c+r)
    return {'a':a,'r':r,'S':S,'beta':beta,'population':float(p),'coherence_amplitude':float(m),
            'boundary_margin':float(p+m*m-1),'normalization':float(np.sqrt(N2))}

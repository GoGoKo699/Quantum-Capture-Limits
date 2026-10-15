"""Number-preserving, time-dependent single-receiver controls.

Lookback time x>=0, with capture at x=0. Rates and detunings are in any common
inverse-time unit. g(x) is the INCIDENT population-decay coupling, kappa the
fixed inaccessible population-decay rate, and gamma_phi extra coherence decay.
No coherent transverse drive, bath record, or shelving level is included.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import quad,solve_ivp
from scipy.linalg import eigh


@dataclass
class Schedule:
    edges:np.ndarray
    incident:np.ndarray
    detuning:np.ndarray
    kappa:float
    gamma_phi:float
    gmax:float
    def __post_init__(self):
        self.edges=np.asarray(self.edges,float);self.incident=np.asarray(self.incident,float);self.detuning=np.asarray(self.detuning,float)
        if len(self.edges)!=len(self.incident)+1 or self.detuning.shape!=self.incident.shape:raise ValueError('Length mismatch.')
        if self.edges[0]!=0 or np.any(np.diff(self.edges)<=0):raise ValueError('Edges must begin at zero and increase.')
        if self.kappa<0 or self.gamma_phi<0 or self.gmax<=0 or np.any(self.incident<0) or np.any(self.incident>self.gmax+1e-14):raise ValueError('Invalid rates.')
        if not np.isfinite(np.r_[self.edges,self.incident,self.detuning,self.kappa,self.gamma_phi,self.gmax]).all():raise ValueError('Nonfinite parameters.')
        self.total=self.incident+self.kappa;self.R_edges=np.r_[0.,np.cumsum(self.total*np.diff(self.edges))]
        self.D_edges=np.r_[0.,np.cumsum(self.detuning*np.diff(self.edges))]
    @property
    def T(self):return float(self.edges[-1])
    @property
    def Gamma_max(self):return self.gmax+self.kappa
    @property
    def beta_max(self):return self.gmax/self.Gamma_max
    @property
    def r_min(self):return self.gamma_phi/self.Gamma_max
    def index(self,x):return np.clip(np.searchsorted(self.edges,x,side='right')-1,0,len(self.incident)-1)
    def rates(self,x):
        j=self.index(x);return self.incident[j],self.total[j],self.detuning[j]
    def exposure(self,x):
        j=self.index(x);return self.R_edges[j]+self.total[j]*(np.asarray(x)-self.edges[j])
    def phase(self,x):
        j=self.index(x);return self.D_edges[j]+self.detuning[j]*(np.asarray(x)-self.edges[j])
    def h(self,x):
        g,_,_=self.rates(x)
        return np.sqrt(g)*np.exp(-self.exposure(x)/2-self.gamma_phi*np.asarray(x)+1j*self.phase(x))
    def kernel(self,x,y):
        gx,_,_=self.rates(x);gy,_,_=self.rates(y)
        return np.sqrt(gx*gy)*np.exp(-(self.exposure(x)+self.exposure(y))/2-self.gamma_phi*abs(np.asarray(x)-np.asarray(y))+1j*(self.phase(x)-self.phase(y)))
    def nodes(self,n_per_interval=24):
        z,w=roots_legendre(n_per_interval);xx=[];ww=[]
        for a,b in zip(self.edges[:-1],self.edges[1:]):
            xx.extend((a+b)/2+(b-a)*z/2);ww.extend(w*(b-a)/2)
        return np.array(xx),np.array(ww)
    def discrete(self,n_per_interval=24):
        x,w=self.nodes(n_per_interval);sw=np.sqrt(w);h=sw*self.h(x)
        K=sw[:,None]*self.kernel(x[:,None],x[None,:])*sw[None,:]
        return x,w,K,h
    def top_margin_objective(self,n_per_interval=24):
        x,w,K,h=self.discrete(n_per_interval);n=len(x)
        return float(eigh(K+np.outer(h,h.conj()),subset_by_index=[n-1,n-1],eigvals_only=True)[0])


def quad_complex(fn,lo,hi,**kw):
    return quad(lambda x:np.real(fn(x)),lo,hi,**kw)[0]+1j*quad(lambda x:np.imag(fn(x)),lo,hi,**kw)[0]


def split_quad(fn,edges):
    return sum(quad_complex(fn,a,b,epsabs=3e-12,epsrel=3e-12,limit=150) for a,b in zip(edges[:-1],edges[1:]))


def matrix_fock(schedule:Schedule,pulse):
    """Independent 2x2 four-block Fock-wavepacket evolution, including phase."""
    sm=np.array([[0.,1.],[0.,0.]],complex);pe=np.diag([0.,1.]).astype(complex)
    gg=np.diag([1.,0.]).astype(complex)
    def dissipator(L,X):
        n=L.conj().T@L
        return L@X@L.conj().T-(n@X+X@n)/2
    def comm(A,B):return A@B-B@A
    state=np.stack([gg,np.zeros((2,2)),np.zeros((2,2)),gg]).astype(complex).ravel()
    # Physical propagation is in increasing t=-lookback, so reverse the segments.
    for j in range(len(schedule.incident)-1,-1,-1):
        g=schedule.incident[j];G=schedule.total[j];det=schedule.detuning[j];lin=np.sqrt(g)*sm
        def L(X):return -1j*det*comm(pe,X)+G*dissipator(sm,X)+2*schedule.gamma_phi*dissipator(pe,X)
        def rhs(t,y):
            q00,q10,q01,q11=y.reshape(4,2,2);xi=pulse(-t)
            return np.stack([L(q00),L(q10)+xi*comm(q00,lin.conj().T),
                L(q01)+xi.conjugate()*comm(lin,q00),
                L(q11)+xi*comm(q01,lin.conj().T)+xi.conjugate()*comm(lin,q10)]).ravel()
        sol=solve_ivp(rhs,(-schedule.edges[j+1],-schedule.edges[j]),state,method='DOP853',rtol=3e-12,atol=2e-14)
        if not sol.success:raise RuntimeError(sol.message)
        state=sol.y[:,-1]
    q00,q10,q01,q11=state.reshape(4,2,2)
    return {'p':float(q11[1,1].real),'m':complex(-q10[1,0]),'trace_error':float(abs(np.trace(q11)-1)),
            'adjoint_error':float(np.linalg.norm(q01-q10.conj().T))}


def kernel_metrics(schedule:Schedule,pulse):
    m=split_quad(lambda x:np.conj(schedule.h(x))*pulse(x),schedule.edges)
    def outer(x):
        ed=np.r_[schedule.edges[schedule.edges<x],x]
        inner=split_quad(lambda y:schedule.kernel(x,y)*pulse(y),ed) if len(ed)>1 else 0j
        return np.conj(pulse(x))*inner
    p=2*split_quad(outer,schedule.edges).real
    return {'p':float(p),'m':complex(m)}

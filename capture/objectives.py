"""Noisy one-photon capture in a fixed two-level Markov receiver.

Time u is look-back time in units of inverse population decay rate Gamma.
r = gamma_phi / Gamma, where gamma_phi is the extra optical-coherence decay.
beta = Gamma_in / Gamma. Pulse functions are normalized on u in [0, infinity).
No monitoring, postselection, time-dependent coupling, or active correction.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import jv
from scipy.optimize import brentq
from scipy.integrate import quad


def _validate(r: float, beta: float = 1.0) -> None:
    if not np.isfinite([r, beta]).all() or r < 0 or not (0 < beta <= 1):
        raise ValueError('Require finite r >= 0 and 0 < beta <= 1.')


def first_bessel_zero(order: float) -> float:
    """First strictly positive zero for order > -1 (including -1 < order < 0)."""
    if not np.isfinite(order) or order <= -1:
        raise ValueError('Order must be finite and greater than -1.')
    lo = max(1e-7, order)
    old = float(jv(order, lo))
    if old <= 0:
        raise ArithmeticError('Unexpected sign before the first Bessel zero.')
    hi = lo + 0.1
    while jv(order, hi) > 0:
        lo = hi
        hi += min(0.25, 0.1 + 0.001 * max(order, 0))
        if hi > 3 * max(order, 1) + 100:
            raise ArithmeticError('Bessel root not bracketed.')
    return float(brentq(lambda z: jv(order, z), lo, hi, xtol=5e-14))


def population_eigenvalue(r: float) -> float:
    _validate(r)
    if r == 0:
        return 1.0
    j = first_bessel_zero(2*r - 1)
    return float(8*r/(j*j))


@dataclass(frozen=True)
class Pulse:
    r: float
    beta: float
    x: float | None
    kind: str

    def __call__(self, u):
        uu = np.asarray(u, dtype=float)
        if np.any(uu < 0):
            raise ValueError('Look-back time u must be nonnegative.')
        if self.r == 0:
            out = np.exp(-uu/2)
        elif self.kind == 'coherence':
            out = np.sqrt(1+2*self.r)*np.exp(-(self.r+.5)*uu)
        else:
            nu, x = 2*self.r, float(self.x)
            a, b, c = jv(nu-1,x), jv(nu,x), jv(nu+1,x)
            norm2 = b*b-a*c
            if norm2 <= 0:
                raise ArithmeticError('Bessel mode normalization is not positive.')
            out = np.exp(-uu/2)*jv(nu, x*np.exp(-uu/2))/np.sqrt(norm2)
        return float(out) if uu.ndim == 0 else out

    def metrics(self) -> dict:
        r, beta = self.r, self.beta
        if r == 0:
            p, m = beta, np.sqrt(beta)
        elif self.kind == 'coherence':
            p, m = beta/(1+r), np.sqrt(beta/(1+2*r))
        elif self.kind == 'population':
            mu = population_eigenvalue(r)
            p, m = beta*mu, np.sqrt(beta)*mu
        else:
            nu, x = 2*r, float(self.x)
            a, b, c = jv(nu-1,x), jv(nu,x), jv(nu+1,x)
            norm2 = b*b-a*c
            hd = 2*c/(x*np.sqrt(norm2))
            source = 2*a/(x*np.sqrt(norm2))
            p = beta*(8*r/x**2-source*hd)
            m = np.sqrt(beta)*hd
        fe = (1+p+2*m)/4
        return {'r':r, 'beta':beta, 'kind':self.kind, 'x':self.x,
                'population':float(p), 'coherence_amplitude':float(m),
                'entanglement_fidelity':float(fe), 'average_qubit_fidelity':float((2*fe+1)/3),
                'cp_slack_p_minus_m2':float(p-m*m)}


def population_pulse(r: float, beta: float=1.) -> Pulse:
    _validate(r,beta)
    return Pulse(r,beta,None if r==0 else first_bessel_zero(2*r-1),'population')


def coherence_pulse(r: float,beta: float=1.) -> Pulse:
    _validate(r,beta)
    return Pulse(r,beta,None,'coherence')


def weighted_pulse(r:float,beta:float=1.,weight:float=1.) -> Pulse:
    """Maximize p+2*weight*m at fixed pulse norm. weight=1 maximizes qubit fidelity.

    The closed Bessel expression is evaluated over a moderate parameter range.
    A failure at extreme parameters raises rather than silently changing precision.
    """
    _validate(r,beta)
    if not np.isfinite(weight) or weight<=0:
        raise ValueError('weight must be finite and positive.')
    if r==0:
        return Pulse(r,beta,None,'fidelity' if weight==1 else f'weight_{weight:g}')
    nu=2*r
    j=first_bessel_zero(nu-1)
    def equation(x):
        a,b,c=jv(nu-1,x),jv(nu,x),jv(nu+1,x)
        if a==0 or not np.isfinite([a,b,c]).all():
            raise ArithmeticError('Bessel quotient underflow or singularity.')
        return weight*weight*x*x*((b/a)**2-c/a)/(4*beta)-1
    lo=j*.2
    while equation(lo)>0:
        lo*=.5
    hi=j*(1-1e-9)
    x=float(brentq(equation,lo,hi,xtol=4e-13))
    return Pulse(r,beta,x,'fidelity' if weight==1 else f'weight_{weight:g}')


def kernel(u,v,r:float):
    return np.exp(-(np.asarray(u)+np.asarray(v))/2-r*np.abs(np.asarray(u)-np.asarray(v)))


def integrate_metrics(pulse:Pulse,cutoff:float=50.) -> dict:
    """Independent scalar ODE of one-photon Fock-state equations, forward in time.

    c'=-(1/2+r)c+sqrt(beta)xi, p'=-p+2sqrt(beta)xi*c.
    The pulse is truncated at time -cutoff without renormalization; its omitted
    tail probability is reported separately. All pulse functions here are real.
    """
    from scipy.integrate import solve_ivp
    r,beta=pulse.r,pulse.beta
    def rhs(t,y):
        xi=pulse(-t)
        return [-(.5+r)*y[0]+np.sqrt(beta)*xi,
                -y[1]+2*np.sqrt(beta)*xi*y[0]]
    sol=solve_ivp(rhs,[-cutoff,0],[0.,0.],method='DOP853',rtol=3e-11,atol=2e-13)
    if not sol.success:
        raise RuntimeError(sol.message)
    norm=quad(lambda u:pulse(u)**2,0,cutoff,epsabs=2e-11,epsrel=2e-11)[0]
    return {'population':float(sol.y[1,-1]),'coherence_amplitude':float(sol.y[0,-1]),
            'pulse_norm_on_truncated_domain':float(norm)}

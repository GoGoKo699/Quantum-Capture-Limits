"""Exact positive-capacity boundary of the fixed noisy single-photon receiver.

r = extra optical-coherence decay / total radiative population decay.
beta = incident radiative fraction. Only vacuum and <=1 photon per independent
channel use are admitted. Rates are per use, not per second. The boundary pulse
maximizes p+|m|^2; it is NOT asserted to maximize the quantum-capacity value.
"""
from __future__ import annotations
import numpy as np
from scipy.special import jvp
from scipy.optimize import brentq
from capture.objectives import Pulse, _validate


def first_derivative_zero(order: float) -> float:
    """First strictly positive zero of J_order' for order>0.

    The root tends to zero as order->0+. This is not the first positive zero
    of J_0', so order=0 must be handled by the rank-one noise-free limit.
    Floating-point implementation is checked at r=.02 through 10.
    """
    if not np.isfinite(order) or order <= 0:
        raise ValueError('Require positive finite order.')
    lo = max(1e-7, order)
    if jvp(order, lo) <= 0:
        raise ArithmeticError('Unexpected derivative sign below the first root.')
    hi = lo + .1
    while jvp(order, hi) > 0:
        lo, hi = hi, hi + .1
        if hi > 3*max(order, 1) + 100:
            raise ArithmeticError('Derivative root not bracketed.')
    return float(brentq(lambda x: jvp(order, x), lo, hi, xtol=2e-14))


def boundary_pulse(r: float, beta: float = 1.) -> Pulse:
    _validate(r, beta)
    return Pulse(r, beta, None if r == 0 else first_derivative_zero(2*r), 'capacity_boundary')


def boundary_metrics(r: float, beta: float = 1.) -> dict:
    pulse = boundary_pulse(r, beta)
    lam = 2. if r == 0 else 8*r/pulse.x**2
    m2 = beta*lam*lam/(4-2*r*lam)
    p = beta*lam-m2
    return {'r': r, 'beta': beta, 'derivative_zero': pulse.x,
            'kernel_largest_eigenvalue': float(lam), 'population': float(p),
            'coherence_amplitude': float(np.sqrt(m2)), 'p_plus_m2': float(beta*lam),
            'positive_capacity_margin': float(beta*lam-1),
            'average_fidelity': float((3+p+2*np.sqrt(m2))/6),
            'scope': 'Exact boundary objective, not the quantum capacity or its rate-optimal pulse.'}


def threshold_r(beta: float = 1.) -> float | None:
    _validate(0., beta)
    if beta <= .5:
        return None  # No positive-capacity region even at r=0.
    hi = 1.
    while boundary_metrics(hi, beta)['positive_capacity_margin'] > 0:
        hi *= 2
    return float(brentq(lambda r: boundary_metrics(r, beta)['positive_capacity_margin'], 1e-8, hi, xtol=5e-13))


def exponential_metrics(a: float, r: float, beta: float = 1.) -> dict:
    _validate(r, beta)
    if not np.isfinite(a) or a <= 0:
        raise ValueError('Require a>0.')
    p = 4*beta*a/((1+a)*(1+a+2*r))
    m = 2*np.sqrt(beta*a)/(1+a+2*r)
    return {'a': float(a), 'r': r, 'beta': beta, 'population': float(p),
            'coherence_amplitude': float(m), 'p_plus_m2': float(p+m*m),
            'average_fidelity': float((3+p+2*m)/6)}


def best_exponential(r: float, beta: float = 1.) -> dict:
    """Global p+m² maximum within normalized rising exponentials.

    The stationary cubic has exactly one positive root by its sign pattern.
    Both endpoints of the objective tend to zero.
    """
    _validate(r, beta)
    c = 1+2*r
    poly = lambda a: 2*a**3+2*a*a+(1-3*c)*a-c*c-c
    hi = c+2
    a = brentq(poly, 0., hi, xtol=2e-14)
    return exponential_metrics(a, r, beta)


def kraus(p: float, m: float) -> list[np.ndarray]:
    if not (0 <= p <= 1 and 0 <= m*m <= p+1e-13):
        raise ValueError('Require 0<=m²<=p<=1.')
    return [np.diag([1., m]), np.array([[0., np.sqrt(1-p)], [0., 0.]]),
            np.diag([0., np.sqrt(max(0., p-m*m))])]


def binary_broadcast_kraus(p: float, m: float) -> list[np.ndarray]:
    """CPTP extension with two identical receiver marginals when p+m²<=1."""
    if not (0 <= p <= 1 and m*m <= p+1e-13 and p+m*m <= 1+1e-13):
        raise ValueError('Outside the broadcastable/antidegradable region.')
    f0 = np.array([[1., 0.], [0., m], [0., m], [0., 0.]])
    f1 = np.zeros((4, 2)); f1[3, 1] = np.sqrt(max(0., p-m*m))
    f2 = np.zeros((4, 2)); f2[0, 1] = np.sqrt(max(0., 1-p-m*m))
    return [f0, f1, f2]


def coherent_information(q: float, p: float, m: float) -> float:
    """CI for diag(1-q,q); not a global single-letter or regularized optimum.

    Stable for moderate q. For very small positive witnesses use the arbitrary-
    precision function below; cancellation can defeat ordinary floating point.
    """
    if not 0 <= q <= 1:
        raise ValueError('Require q in [0,1].')
    if q in (0., 1.):
        return 0.
    a = (1-p)*q; t = 1-a; det = (p-m*m)*q*(1-q)
    lp = (t+np.sqrt(max(0., t*t-4*det)))/2
    lm = det/lp if lp else 0.
    xlog = lambda x: 0. if x <= 0 else -x*np.log2(x)
    missing = a+lm
    henv = xlog(a)+xlog(lm)-(1-missing)*np.log1p(-missing)/np.log(2)
    y = p*q
    hout = xlog(y)-(1-y)*np.log1p(-y)/np.log(2)
    return float(hout-henv)


def coherent_information_mp(q, p, m2):
    """Arbitrary precision evaluation. Caller sets mpmath precision."""
    import mpmath as mp
    q,p,m2 = mp.mpf(q),mp.mpf(p),mp.mpf(m2)
    if not (0 < q < 1 and 0 <= m2 <= p <= 1):
        raise ValueError('Invalid channel/input parameters.')
    a=(1-p)*q; t=1-a; det=(p-m2)*q*(1-q)
    lp=(t+mp.sqrt(t*t-4*det))/2; lm=det/lp
    xlog=lambda x: -x*mp.log(x,2) if x else mp.mpf('0')
    missing=a+lm; y=p*q
    henv=xlog(a)+xlog(lm)-(1-missing)*mp.log1p(-missing)/mp.log(2)
    hout=xlog(y)-(1-y)*mp.log1p(-y)/mp.log(2)
    return hout-henv

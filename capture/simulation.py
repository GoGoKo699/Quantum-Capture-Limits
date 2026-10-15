"""Exact finite-mode tests of the continuum capture preprocessing identity.

This module does not replace the continuum proof by a discretization. The
finite matrices are Nyström models, checked to define CPTP channels in each
reported instance. `capture.receiver` contains the independently organized
Fock-dynamics implementation. All rates use one common inverse-time unit.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.linalg import eigh
from scipy.integrate import quad
from capture.receiver import Schedule, quad_complex


def psd_sqrt(a: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    a = np.asarray(a, complex)
    eig, v = eigh((a + a.conj().T) / 2)
    if eig[0] < -tol:
        raise ValueError(f'Matrix not PSD: minimum eigenvalue {eig[0]}')
    return (v * np.sqrt(np.maximum(eig, 0))) @ v.conj().T


def receiver(a: np.ndarray, b: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Vacuum-plus-n-mode -> receiver qubit map E_(a,b). Linear on all inputs."""
    n = len(b)
    if a.shape != (n, n) or rho.shape != (n + 1, n + 1):
        raise ValueError('Receiver dimension mismatch')
    p = np.trace(a @ rho[1:, 1:])
    return np.array([[rho[0, 0] + np.trace(rho[1:, 1:]) - p, rho[0, 1:] @ b],
                     [b.conj() @ rho[1:, 0], p]], complex)


def receiver_kraus(a: np.ndarray, b: np.ndarray) -> list[np.ndarray]:
    """Separate Kraus construction of E_(a,b); tests its CP domain."""
    n = len(b)
    f = np.zeros((2, n + 1), complex)
    f[0, 0] = 1
    f[1, 1:] = b.conj()
    out = [f]
    for target, residue in [(1, a - np.outer(b, b.conj())), (0, np.eye(n) - a)]:
        for row in psd_sqrt(residue):
            f = np.zeros((2, n + 1), complex)
            f[target, 1:] = row
            out.append(f)
    return out


def apply_kraus(kraus: list[np.ndarray], rho: np.ndarray) -> np.ndarray:
    return sum(k @ rho @ k.conj().T for k in kraus)


def preprocess(rho: np.ndarray, c: np.ndarray, d: np.ndarray,
               phase: np.ndarray, gamma: float) -> np.ndarray:
    """Trace-preserving attenuation + extra-time Wiener phase channel.

    Deterministic input phase is exp(-i phase). Lost photons go to vacuum;
    no conditional normalization and no success postselection are performed.
    """
    n = len(c)
    if rho.shape != (n + 1, n + 1) or len(d) != n or len(phase) != n:
        raise ValueError('Preprocessor dimension mismatch')
    if gamma < 0 or np.min(c) < -1e-12 or np.max(c) > 1 + 1e-12 or np.min(d) < -1e-12:
        raise ValueError('Invalid attenuation, delay, or noise rate')
    amp = c * np.exp(-1j * phase)
    out = np.empty_like(rho, complex)
    out[0, 0] = rho[0, 0] + np.sum((1 - c*c) * np.diag(rho)[1:])
    hnoise = np.exp(-gamma * np.maximum(d, 0))
    out[0, 1:] = rho[0, 1:] * amp.conj() * hnoise
    out[1:, 0] = rho[1:, 0] * amp * hnoise
    correlation = np.exp(-gamma * np.abs(d[:, None] - d[None, :]))
    out[1:, 1:] = rho[1:, 1:] * amp[:, None] * amp.conj()[None, :] * correlation
    return out


def preprocessing_kraus(c: np.ndarray, d: np.ndarray,
                         phase: np.ndarray, gamma: float) -> list[np.ndarray]:
    """Independent loss Kraus + Schur-product Kraus construction."""
    n = len(c)
    loss = [np.diag(np.r_[1., c * np.exp(-1j * phase)])]
    for i in range(n):
        k = np.zeros((n + 1, n + 1), complex)
        k[0, i + 1] = np.sqrt(max(0., 1 - c[i]*c[i]))
        loss.append(k)
    labels = np.r_[0., d]
    correlation = np.exp(-gamma * np.abs(labels[:, None] - labels[None, :]))
    eig, v = eigh(correlation)
    if eig[0] < -1e-10:
        raise ValueError('Noise correlation matrix is not PSD')
    dephase = [np.diag(v[:, i] * np.sqrt(eig[i])) for i in range(n + 1) if eig[i] > 1e-14]
    return [k @ l for k in dephase for l in loss]


def discrete_factorization(schedule: Schedule, order: int = 3) -> dict:
    """Two independently expressed kernels on paired physical/exposure grids.

    At Gamma=0 a grid cell never couples. Its input is deterministically lost to
    vacuum by the preprocessor; it is represented by an unused dark dimension
    in the finite reference model. The continuum implementation simply discards
    that subspace before the time-coordinate isometry.
    """
    x, w, ag, bg = schedule.discrete(order)
    g, rate, _ = schedule.rates(x)
    u = schedule.exposure(x)
    d = x - u / schedule.Gamma_max
    assert np.min(d) >= -1e-12
    beta = np.divide(g, rate, out=np.zeros_like(g), where=rate > 0)
    c = np.sqrt(beta / schedule.beta_max)
    # du=Gamma dx. Zero-width entries are unused auxiliary dimensions only.
    wu = rate * w
    ref_a = schedule.beta_max * np.sqrt(wu[:, None] * wu[None, :]) * np.exp(
        -(u[:, None] + u[None, :]) / 2 - schedule.r_min * np.abs(u[:, None] - u[None, :]))
    ref_b = np.sqrt(schedule.beta_max * wu) * np.exp(-u/2 - schedule.r_min*u)
    return {'x': x, 'w': w, 'u': u, 'd': np.maximum(d, 0), 'c': c,
            'phase': schedule.phase(x), 'A_control': ag, 'b_control': bg,
            'A_reference': ref_a, 'b_reference': ref_b,
            'unused_zero_rate_dimensions': int(np.count_nonzero(rate == 0))}


def apply_with_reference(map_fn, rho: np.ndarray, reference_dim: int, input_dim: int) -> np.ndarray:
    """I_reference tensor map, evaluating all reference blocks explicitly."""
    if rho.shape != (reference_dim * input_dim, reference_dim * input_dim):
        raise ValueError('Reference input dimension mismatch')
    rows = []
    for i in range(reference_dim):
        row = []
        for j in range(reference_dim):
            row.append(map_fn(rho[i*input_dim:(i+1)*input_dim, j*input_dim:(j+1)*input_dim]))
        rows.append(row)
    return np.block(rows)


def _split_integral(fn, intervals):
    return sum(quad_complex(fn, lo, hi, epsabs=3e-12, epsrel=3e-12, limit=150)
               for lo, hi in intervals if hi > lo)


def continuum_simulated_metrics(schedule: Schedule, pulse) -> dict:
    """Compute E_max(P(rho)) by independent integration in exposure time u.

    The physical-time master equation does not enter this implementation.
    The one-photon input may have amplitude in uncoupled intervals. Its lost
    weight is kept as vacuum probability in the trace-preserving preparation.
    """
    segments = []
    for j, rate in enumerate(schedule.total):
        if rate > 0:
            segments.append((float(schedule.R_edges[j]), float(schedule.R_edges[j+1]), j))
    if not segments:
        return {'p': 0., 'm': 0j, 'prepared_photon_probability': 0.}
    intervals = [(a, b) for a, b, _ in segments]

    def parameters(u):
        # At measure-zero shared endpoints any adjacent continuous integral
        # convention is valid. Zero-rate intervals become jumps in x(u).
        index = len(segments)-1
        for idx, (_, hi, _) in enumerate(segments):
            if u < hi:
                index = idx
                break
        lo, _, j = segments[index]
        rate = schedule.total[j]
        x = schedule.edges[j] + (u-lo)/rate
        beta = schedule.incident[j]/rate
        c = np.sqrt(beta/schedule.beta_max)
        d = max(0., float(x-u/schedule.Gamma_max))
        psi = pulse(x)*np.exp(-1j*schedule.phase(x))/np.sqrt(rate)
        return c, d, psi

    def coherence(u):
        c, d, psi = parameters(u)
        h0 = np.sqrt(schedule.beta_max)*np.exp(-u/2-schedule.r_min*u)
        return h0*c*np.exp(-schedule.gamma_phi*d)*psi

    m = _split_integral(coherence, intervals)
    photon = _split_integral(lambda u: parameters(u)[0]**2*abs(parameters(u)[2])**2, intervals).real

    def outer(u):
        cu, du, pu = parameters(u)
        clipped = [(a, min(b, u)) for a, b in intervals if a < u]
        def inner(v):
            cv, dv, pv = parameters(v)
            k0 = schedule.beta_max*np.exp(-(u+v)/2-schedule.r_min*abs(u-v))
            return pu.conjugate()*cu*k0*cv*pv*np.exp(-schedule.gamma_phi*abs(du-dv))
        return _split_integral(inner, clipped)

    p = 2*_split_integral(outer, intervals).real
    return {'p': float(p), 'm': complex(m), 'prepared_photon_probability': float(photon)}

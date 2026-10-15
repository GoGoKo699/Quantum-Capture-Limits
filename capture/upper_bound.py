"""Regularized quantum-capacity upper bounds from an antidegradable comparison.

The analytical bound applies to the full temporal input channel, not only a
selected pulse or a single-letter optimization. Numerical receiver spectra use
floating-point Bessel roots; their evaluated bounds are not interval-certified.
See docs/capacity-upper-bound.md for the proof and continuity attribution.
"""
from __future__ import annotations

import numpy as np

from capture.capacity import boundary_metrics
from capture.finite_window import finite_eigenvalue, validate


def continuity_penalty(epsilon: float) -> float:
    """2 epsilon + g(epsilon), in bits, for a qubit output.

    epsilon is HALF the diamond distance. g(epsilon) is evaluated with log1p
    to retain small positive values; its continuous value at zero is zero.
    """
    if not np.isfinite(epsilon) or not 0 <= epsilon <= 1:
        raise ValueError('Require finite 0<=epsilon<=1.')
    if epsilon == 0:
        return 0.
    g = ((1 + epsilon) * np.log1p(epsilon)
         - epsilon * np.log(epsilon)) / np.log(2)
    return float(2 * epsilon + g)


def spectral_capacity_bound(mu: float) -> dict:
    """Bound Q for E_(A,b) given mu >= ||A+|b><b|||.

    The channel must satisfy 0<=|b><b|<=A<=I. A certified spectral UPPER
    estimate may replace the exact norm; a Galerkin lower bound may not.
    No claim that the chosen antidegradable comparison is nearest is made.
    """
    if not np.isfinite(mu) or mu < 0:
        raise ValueError('Require a finite nonnegative spectral objective.')
    tau = 1. if mu <= 1 else 1 / mu
    epsilon = 0. if mu <= 1 else (mu - 1) / mu
    return {
        'spectral_objective': float(mu),
        'output_transmission': float(tau),
        'diamond_distance_upper': float(epsilon),
        'capacity_upper_bound': min(1., continuity_penalty(epsilon)),
        'scope': ('Upper bound on regularized unassisted capacity in qubits per '
                  'independent use; diamond_distance_upper is HALF the diamond '
                  'norm. The exact capacity and nearest comparison are not asserted.'),
    }


def receiver_capacity_bound(r: float, beta: float = 1., S: float | None = None) -> dict:
    """Evaluate the analytical upper bound using the continuum receiver spectrum.

    S=None denotes the unlimited window. Finite S is dimensionless Gamma_M*T.
    Values near the sign boundary inherit floating-point root uncertainty.
    """
    if S is None:
        lam = boundary_metrics(r, beta)['kernel_largest_eigenvalue']
    else:
        validate(r, S, beta)
        lam = finite_eigenvalue(r, S)
    result = spectral_capacity_bound(beta * lam)
    result.update(r=float(r), beta=float(beta), S=None if S is None else float(S),
                  numerical_status='Floating-point evaluation, not an interval certificate.')
    return result

#!/usr/bin/env python3
"""Independent finite checks for the all-mode quantum-capacity upper bound.

The continuum comparison and regularized-capacity continuity theorem are
analytical inputs, not proved by these finite tests. Galerkin eigenvalues are
not used as upper bounds. All distances below use HALF the diamond norm.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import mpmath as mp
import numpy as np
import scipy
from scipy.optimize import minimize_scalar
import sympy as sp

from capture.upper_bound import (continuity_penalty, receiver_capacity_bound,
                                 spectral_capacity_bound)

ROOT = Path(__file__).resolve().parents[1]


def apply_channel(operators, rho):
    return sum(k @ rho @ k.conj().T for k in operators)


def damping_operators(tau):
    return [np.diag([1., np.sqrt(tau)]),
            np.array([[0., np.sqrt(1-tau)], [0., 0.]])]


def choi(operators):
    """Unnormalized Choi matrix, input index first, output index second."""
    d = operators[0].shape[1]
    out = operators[0].shape[0]
    result = np.zeros((d*out, d*out), complex)
    for i in range(d):
        for j in range(d):
            unit = np.zeros((d, d), complex)
            unit[i, j] = 1.
            result[out*i:out*(i+1), out*j:out*(j+1)] = apply_channel(operators, unit)
    return result


def diamond_certificate_checks():
    # Let s=sqrt(tau), d=1-s^2, a=1-s. The dual feasible Z obeys
    # Z >= +/-J(Id-AD_tau) and ||Tr_output Z||_infty = 2d.
    # Factorizations below establish positivity for every 0<=s<=1,
    # including s=1 without the removable 0/0 in 2a^2/d.
    s = sp.symbols('s', real=True, nonnegative=True)
    d, a = 1-s*s, 1-s
    J = sp.Matrix([[0, 0, 0, a], [0, 0, 0, 0],
                   [0, 0, -d, 0], [a, 0, 0, d]])
    z00 = 2*(1-s)/(1+s)
    Z = sp.Matrix([[z00, 0, 0, a], [0, 0, 0, 0],
                   [0, 0, d, 0], [a, 0, 0, d]])
    minus = sp.diag(z00, 0, 2*d, 0)
    v = sp.Matrix([1, 0, 0, 1+s])
    plus = z00 * v * v.T
    assert all(sp.simplify(x) == 0 for x in Z-J-minus)
    assert all(sp.simplify(x) == 0 for x in Z+J-plus)
    assert sp.simplify(2*d-z00-2*s*(1-s)*(2+s)/(1+s)) == 0

    # A matching primal feasible point uses rho=|1><1| and
    # W=rho tensor (|1><1|-|0><0|), so Tr(JW)=2d.
    rho = np.diag([0., 1.])
    W = np.kron(rho, np.diag([-1., 1.]))
    primal_box = np.kron(rho, np.eye(2))
    assert np.linalg.eigvalsh(primal_box-W)[0] >= 0
    assert np.linalg.eigvalsh(primal_box+W)[0] >= 0
    assert sp.simplify(sp.trace(J*sp.Matrix(W))-2*d) == 0

    rows, errors, lowest = [], [], []
    for tau in [0., 1e-8, .01, .25, .5, .9, .999999, 1.]:
        root = np.sqrt(tau)
        delta = 1-tau
        offdiag = 1-root
        actual_J = choi([np.eye(2)])-choi(damping_operators(tau))
        witness = np.zeros((4, 4))
        witness[0, 0] = 2*(1-root)/(1+root)
        witness[0, 3] = witness[3, 0] = offdiag
        witness[2, 2] = witness[3, 3] = delta
        partial = np.trace(witness.reshape(2, 2, 2, 2), axis1=1, axis2=3)
        upper = float(np.linalg.eigvalsh(partial)[-1])
        lower = float(np.trace(actual_J@W).real)
        achieved = float(np.linalg.norm(rho-apply_channel(damping_operators(tau), rho), ord='nuc'))
        errors.extend([abs(upper-2*delta), abs(lower-2*delta), abs(achieved-2*delta)])
        lowest.extend([float(np.linalg.eigvalsh(witness-actual_J)[0]),
                       float(np.linalg.eigvalsh(witness+actual_J)[0])])
        rows.append({'output_transmission': tau, 'half_diamond_distance': delta,
                     'dual_objective_full_norm': upper,
                     'primal_objective_full_norm': lower})
    assert min(lowest) > -2e-14
    assert max(errors) < 2e-14
    return {'symbolic_dual_factorizations': 'Z-J=diag(2(1-s)/(1+s),0,2(1-s^2),0); Z+J=2(1-s)/(1+s) vv^T, v=(1,0,0,1+s)',
            'max_primal_dual_gap_or_distance_error': max(errors),
            'minimum_dual_slack_eigenvalue': min(lowest), 'examples': rows,
            'scope': 'Explicit primal and dual feasible certificates; no SDP solver or state sampling is used to infer the exact diamond norm.'}


def psqrt(matrix):
    values, vectors = np.linalg.eigh(matrix)
    assert values[0] > -3e-13
    return (vectors*np.sqrt(np.maximum(values, 0))) @ vectors.conj().T


def independent_capture_operators(A, b, broadcast=False):
    """Build Kraus matrices directly from the positive operator remainders."""
    n = len(b)
    output = 4 if broadcast else 2
    F = np.zeros((output, n+1), complex)
    F[0, 0] = 1.
    F[1, 1:] = b.conj()
    bb = np.outer(b, b.conj())
    if broadcast:
        F[2, 1:] = b.conj()
        remainders = [(3, A-bb), (0, np.eye(n)-A-bb)]
    else:
        remainders = [(1, A-bb), (0, np.eye(n)-A)]
    operators = [F]
    for target, remainder in remainders:
        for row in psqrt(remainder):
            operator = np.zeros((output, n+1), complex)
            operator[target, 1:] = row
            operators.append(operator)
    return operators


def capture_blocks(A, b, rho):
    population = np.trace(A@rho[1:, 1:])
    return np.array([[np.trace(rho)-population, rho[0, 1:]@b],
                     [b.conj()@rho[1:, 0], population]])


def multimode_comparison_checks():
    rng = np.random.default_rng(2026101505)
    errors, trace_errors, minimums, rows = [], [], [], []
    specifications = [([.3], [.4]), ([1.], [1.]),
                      ([.91, .7, .4], [.88, .18, .12]),
                      ([.95, .85, .7, .6, .4, .2, .05], [.7, .3, .2, .15, .1, .08, .02])]
    for eigenvalues, weights in specifications:
        n = len(eigenvalues)
        unitary, _ = np.linalg.qr(rng.normal(size=(n, n))+1j*rng.normal(size=(n, n)))
        A = (unitary*np.asarray(eigenvalues))@unitary.conj().T
        b = unitary@(np.sqrt(eigenvalues)*np.asarray(weights))
        mu = float(np.linalg.eigvalsh(A+np.outer(b, b.conj()))[-1])
        result = spectral_capacity_bound(mu)
        tau, epsilon = result['output_transmission'], result['diamond_distance_upper']
        original = independent_capture_operators(A, b)
        compared = independent_capture_operators(tau*A, np.sqrt(tau)*b)
        broadcast = independent_capture_operators(tau*A, np.sqrt(tau)*b, True)
        composed = [damp@operator for damp in damping_operators(tau) for operator in original]
        for operators in [original, compared, broadcast, composed]:
            trace_errors.append(float(np.linalg.norm(sum(k.conj().T@k for k in operators)-np.eye(n+1))))
        minimums.append(float(np.linalg.eigvalsh(np.eye(n)-tau*(A+np.outer(b, b.conj())))[0]))
        for i in range(n+1):
            for j in range(n+1):
                unit = np.zeros((n+1, n+1), complex)
                unit[i, j] = 1.
                target = capture_blocks(tau*A, np.sqrt(tau)*b, unit)
                full = apply_channel(broadcast, unit).reshape(2, 2, 2, 2)
                for actual in [apply_channel(compared, unit), apply_channel(composed, unit),
                               np.trace(full, axis1=1, axis2=3),
                               np.trace(full, axis1=0, axis2=2)]:
                    errors.append(float(np.linalg.norm(actual-target)))
                errors.append(float(np.linalg.norm(apply_channel(original, unit)-capture_blocks(A, b, unit))))
        max_reference_distance = 0.
        for rank in [1, 3]:
            reference_dimension = 3
            X = rng.normal(size=(reference_dimension*(n+1), rank))+1j*rng.normal(size=(reference_dimension*(n+1), rank))
            rho = X@X.conj().T
            rho /= np.trace(rho)
            def with_reference(operators):
                return apply_channel([np.kron(np.eye(reference_dimension), k) for k in operators], rho)
            target = with_reference(composed)
            extension = with_reference(broadcast).reshape(reference_dimension, 2, 2, reference_dimension, 2, 2)
            for actual in [with_reference(compared),
                           np.trace(extension, axis1=2, axis2=5).reshape(2*reference_dimension, 2*reference_dimension),
                           np.trace(extension, axis1=1, axis2=4).reshape(2*reference_dimension, 2*reference_dimension)]:
                errors.append(float(np.linalg.norm(actual-target)))
            distance = float(np.linalg.norm(with_reference(original)-target, ord='nuc')/2)
            assert distance <= epsilon+2e-13
            max_reference_distance = max(max_reference_distance, distance)
        rows.append({'photon_modes': n, 'spectral_objective': mu, 'output_transmission': tau,
                     'half_diamond_upper_bound': epsilon,
                     'largest_sampled_reference_trace_distance': max_reference_distance})
    assert max(errors+trace_errors) < 3e-12
    assert min(minimums) > -3e-13
    return {'max_channel_identity_or_marginal_error': max(errors),
            'max_trace_preservation_error': max(trace_errors),
            'minimum_broadcast_remainder_eigenvalue': min(minimums), 'examples': rows,
            'scope': 'Complex finite-mode Kraus and all-matrix-unit checks with entangled references. Sampling does not prove antidegradability or diamond contraction in the continuum.'}


def mp_binary_entropy(value):
    if value in (0, 1):
        return mp.mpf('0')
    return -value*mp.log(value, 2)-(1-value)*mp.log(1-value, 2)


def independent_penalty(epsilon):
    return 2*epsilon+(1+epsilon)*mp_binary_entropy(epsilon/(1+epsilon))


def capacity_function_checks():
    penalty_errors, bound_errors, rows = [], [], []
    with mp.workdps(90):
        for epsilon in [0., 1e-16, 1e-12, 1e-8, .001, .1, .5, 1.]:
            expected = float(independent_penalty(mp.mpf(epsilon)))
            actual = continuity_penalty(epsilon)
            penalty_errors.append(abs(actual-expected)/max(expected, 1e-300))
        previous = -1.
        for mu in [0., .25, .9, 1., 1+1e-12, 1.001, 1.05, 1.2, 1.5, 2., 4.]:
            value = spectral_capacity_bound(mu)
            tau = min(mp.mpf(1), 1/mp.mpf(mu)) if mu else mp.mpf(1)
            epsilon = 1-tau
            expected = min(mp.mpf(1), independent_penalty(epsilon))
            bound_errors.extend([abs(value['capacity_upper_bound']-float(expected)),
                                 abs(value['diamond_distance_upper']-float(epsilon)),
                                 abs(value['output_transmission']-float(tau))])
            assert previous <= value['capacity_upper_bound'] <= 1.
            previous = value['capacity_upper_bound']
    assert max(penalty_errors) < 5e-15
    assert max(bound_errors) < 5e-15

    # In the zero-extra-dephasing case the receiver reduces, after selecting
    # its single occupied temporal mode, to amplitude damping with eta.
    # Compare against its established exact capacity: max_q[h2(eta q)
    # -h2((1-eta)q)] for eta>=1/2, and zero below. This is only a special case.
    def entropy(x):
        return float(-x*np.log2(x)-(1-x)*np.log2(1-x)) if 0 < x < 1 else 0.
    for eta in [0., .2, .5, .500001, .51, .55, .6, .7, .9, 1.]:
        if eta <= .5:
            exact = 0.
        else:
            optimum = minimize_scalar(lambda q: -(entropy(eta*q)-entropy((1-eta)*q)),
                                      bounds=(0., 1.), method='bounded',
                                      options={'xatol': 1e-13})
            assert optimum.success
            exact = max(0., -float(optimum.fun))
        upper = spectral_capacity_bound(2*eta)['capacity_upper_bound']
        assert exact <= upper+2e-13
        rows.append({'noiseless_population_transmission': eta,
                     'amplitude_damping_exact_capacity_numerically_optimized': exact,
                     'general_capacity_upper_bound': upper})
    for function, bad_values in [(continuity_penalty, [-1., 1.01, np.inf, np.nan]),
                                 (spectral_capacity_bound, [-1., np.inf, np.nan])]:
        for value in bad_values:
            try:
                function(value)
            except ValueError:
                pass
            else:
                raise AssertionError(f'Invalid value accepted: {function.__name__}({value})')
    return {'max_relative_penalty_error_vs_90_digit_entropy_expression': max(penalty_errors),
            'max_spectral_bound_error_vs_90_digit_expression': max(bound_errors),
            'amplitude_damping_comparisons': rows,
            'scope': 'Only amplitude damping has its exact capacity used here. The general bound does not optimize a noisy-receiver capacity.'}


def mp_principal_root(r, S=None):
    """Independent high-precision sign-bracketing of continuum Bessel roots."""
    order = 2*r
    derivative = lambda z: (mp.besselj(order-1, z)-mp.besselj(order+1, z))/2
    if S is None:
        residual = derivative
        left = order
    else:
        def residual(z):
            endpoint = z*mp.exp(-S/2)
            yd = (mp.bessely(order-1, z)-mp.bessely(order+1, z))/2
            return derivative(z)*mp.bessely(order+1, endpoint)-yd*mp.besselj(order+1, endpoint)
        left = mp_principal_root(r)-mp.mpf('1e-40')
    left_value = residual(left)
    right = left+mp.mpf('.1')
    for _ in range(1000):
        if residual(right)*left_value < 0:
            break
        left, right = right, right+mp.mpf('.1')
    else:
        raise AssertionError('Independent principal Bessel root not bracketed.')
    root = mp.findroot(residual, (left, right), solver='anderson', tol=mp.mpf('1e-60'))
    assert left < root < right
    return root


def receiver_and_boundary_checks():
    rows, errors, scaling = [], [], []
    with mp.workdps(70):
        for r, beta, S in [(0., 1., None), (0., .8, .7), (0., .4, 4.),
                           (.1, 1., None), (.1, 1., .3), (.5, 1., None),
                           (.5, 1., 4.), (.72, 1., None), (.72, 1., 6.),
                           (1., .8, 3.), (3., 1., None)]:
            rr, bb = mp.mpf(r), mp.mpf(beta)
            SS = None if S is None else mp.mpf(S)
            if r == 0:
                eigenvalue = mp.mpf(2) if S is None else 2*(1-mp.exp(-SS))
            else:
                root = mp_principal_root(rr, SS)
                eigenvalue = 8*rr/(root*root)
            mu = bb*eigenvalue
            expected = min(mp.mpf(1), independent_penalty(max(mp.mpf(0), 1-1/mu)))
            result = receiver_capacity_bound(r, beta, S)
            errors.extend([abs(result['spectral_objective']-float(mu)),
                           abs(result['capacity_upper_bound']-float(expected))])
            rows.append({'r': r, 'beta': beta, 'S': S,
                         'continuum_spectral_objective_70_digits': mp.nstr(mu, 32),
                         'capacity_upper_bound': result['capacity_upper_bound']})

        # Follow an actual continuum receiver to its threshold from the
        # positive side. No finite-space eigenvalue supplies this bound.
        def threshold_equation(rr):
            z = mp.sqrt(8*rr)
            return (mp.besselj(2*rr-1, z)-mp.besselj(2*rr+1, z))/2
        critical = mp.findroot(threshold_equation, (mp.mpf('.73'), mp.mpf('.74')))
        for displacement in ['0.01', '0.001', '0.0001', '0.00001']:
            rr = float(critical-mp.mpf(displacement))
            result = receiver_capacity_bound(rr)
            mu = 8*mp.mpf(rr)/mp_principal_root(mp.mpf(rr))**2
            epsilon = 1-1/mu
            expected = independent_penalty(epsilon)
            errors.append(abs(result['capacity_upper_bound']-float(expected)))
            scaling.append({'r': rr, 'spectral_excess': float(mu-1),
                            'capacity_upper_bound': result['capacity_upper_bound']})
        assert all(scaling[i+1]['capacity_upper_bound'] < scaling[i]['capacity_upper_bound']
                   for i in range(len(scaling)-1))

    assert max(errors) < 2e-11
    residual_errors = []
    for delta in [1e-3, 1e-5, 1e-7, 1e-9]:
        result = spectral_capacity_bound(1+delta)
        epsilon = result['diamond_distance_upper']
        remainder_coefficient = result['capacity_upper_bound']/epsilon+np.log2(epsilon)
        residual_errors.append(abs(remainder_coefficient-(2+1/np.log(2))))
    assert all(residual_errors[i+1] < residual_errors[i] for i in range(len(residual_errors)-1))
    assert residual_errors[-1] < 2e-8
    return {'max_error_vs_independent_continuum_roots': max(errors), 'receiver_examples': rows,
            'continuum_approach_to_boundary': scaling,
            'small_excess_expansion': 'U=epsilon log2(1/epsilon)+(2+log2(e))epsilon+O(epsilon^2), epsilon=(mu-1)/mu',
            'last_linear_remainder_coefficient_error': residual_errors[-1],
            'scope': 'High-precision continuum root comparisons are numerical checks, not certified intervals; strict decisions at an unresolved floating-point boundary are not claimed.'}


def main(output):
    if not __debug__:
        raise RuntimeError('Keep assertions enabled.')
    start = time.time()
    result = {'status': 'PASS', 'date': '2026-10-15',
              'scope': 'Upper bounds on regularized unassisted capacity per independently reset one-photon capture use; no exact noisy capacity or experimental validation.',
              'versions': {'python': platform.python_version(), 'numpy': np.__version__,
                           'scipy': scipy.__version__, 'sympy': sp.__version__, 'mpmath': mp.__version__},
              'diamond_primal_dual_certificate': diamond_certificate_checks(),
              'multimode_antidegradable_comparison': multimode_comparison_checks(),
              'capacity_function': capacity_function_checks(),
              'continuum_receivers_and_boundary': receiver_and_boundary_checks()}
    result['elapsed_seconds'] = time.time()-start
    result['source_sha256'] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in sorted(p for directory in ('capture', 'tests', 'tools')
                                               for p in (ROOT/directory).glob('*.py'))}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', 'groups': 4, 'elapsed_seconds': result['elapsed_seconds'],
                      'output': str(output)}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'artifacts/upper_bound.json')
    main(parser.parse_args().output)

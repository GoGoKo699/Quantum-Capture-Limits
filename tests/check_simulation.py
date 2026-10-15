#!/usr/bin/env python3
"""Recompute bounded-coupling capture-channel simulation checks.

The continuum identity and capacity ordering are analytical. Finite-mode CPTP
maps, arbitrary matrix units, reference-entangled states, independent physical-
time Fock ODEs, and a fixed-source counterexample test the construction. No
experimental data, hardware fit, or quantum-capacity numerical estimate is used.
"""
from __future__ import annotations
import argparse, hashlib, json, platform, time
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import eigvalsh
import sympy as sp
from capture.receiver import Schedule, matrix_fock
from capture.simulation import (discrete_factorization, preprocess, preprocessing_kraus,
    receiver, receiver_kraus, apply_kraus, apply_with_reference,
    continuum_simulated_metrics)

ROOT = Path(__file__).resolve().parents[1]


def schedules():
    return [
        Schedule([0, .7, 2.2, 3.], [1., .25, .8], [.3, -.4, .7], .2, .4, 1.),
        Schedule([0, 1., 2., 3.], [1., 0., .5], [.3, -.4, .7], 0., .8, 1.),
        Schedule([0, .5, 1.4, 2.], [.2, .9, 0.], [0., .7, -.2], .25, 0., 1.),
        Schedule([0, 2.], [1.], [.4], 0., .5, 1.),
        Schedule([0, 2.], [0.], [-.6], 0., .4, 1.),
    ]


def factorization_checks():
    rng = np.random.default_rng(154404)
    rows = []
    for k, sch in enumerate(schedules()):
        d = discrete_factorization(sch, order=3)
        ac, bc = d['A_control'], d['b_control']
        ar, br = d['A_reference'], d['b_reference']
        c, delay, phase = d['c'], d['d'], d['phase']
        dim = len(c) + 1
        # These must be real CPTP finite models before any channel comparison.
        kc, kr = receiver_kraus(ac, bc), receiver_kraus(ar, br)
        kp = preprocessing_kraus(c, delay, phase, sch.gamma_phi)
        tp = np.linalg.norm(sum(v.conj().T @ v for v in kp) - np.eye(dim))
        rc_tp = np.linalg.norm(sum(v.conj().T @ v for v in kc) - np.eye(dim))
        rr_tp = np.linalg.norm(sum(v.conj().T @ v for v in kr) - np.eye(dim))
        direct = lambda x: receiver(ac, bc, x)
        simulated = lambda x: receiver(ar, br, preprocess(x, c, delay, phase, sch.gamma_phi))
        e_units = e_prep = e_kraus = 0.
        for i in range(dim):
            for j in range(dim):
                unit = np.zeros((dim, dim), complex)
                unit[i, j] = 1
                e_units = max(e_units, float(np.linalg.norm(direct(unit) - simulated(unit))))
                e_prep = max(e_prep, float(np.linalg.norm(preprocess(unit, c, delay, phase, sch.gamma_phi) - apply_kraus(kp, unit))))
                e_kraus = max(e_kraus, float(np.linalg.norm(direct(unit) - apply_kraus(kc, unit))))
        # Include coherences with an untouched reference and a mixed entangled input.
        refdim = 3
        raw = rng.normal(size=(refdim*dim, 4)) + 1j*rng.normal(size=(refdim*dim, 4))
        rho = raw @ raw.conj().T
        rho /= np.trace(rho)
        out = apply_with_reference(direct, rho, refdim, dim)
        sim = apply_with_reference(simulated, rho, refdim, dim)
        ent_error = np.linalg.norm(out - sim)
        # The construction must not increase the expected number of input photons.
        reduced = rho.reshape(refdim, dim, refdim, dim).trace(axis1=0, axis2=2)
        pre = preprocess(reduced, c, delay, phase, sch.gamma_phi)
        n_before = float(np.trace(reduced[1:, 1:]).real)
        n_after = float(np.trace(pre[1:, 1:]).real)
        assert n_after <= n_before + 1e-12
        assert max(tp, rc_tp, rr_tp, e_units, e_prep, e_kraus, ent_error) < 3e-11
        rows.append({'schedule': k, 'input_dimension': dim,
            'edges': sch.edges.tolist(), 'incident_rates': sch.incident.tolist(),
            'detuning_rates': sch.detuning.tolist(), 'kappa': sch.kappa,
            'gamma_phi': sch.gamma_phi,
            'unused_zero_rate_dimensions': d['unused_zero_rate_dimensions'],
            'preprocessor_TP_error': float(tp), 'receiver_TP_error': float(max(rc_tp, rr_tp)),
            'all_matrix_units_simulation_error': e_units,
            'preprocessor_Kraus_vs_formula_error': e_prep,
            'receiver_Kraus_vs_formula_error': e_kraus,
            'reference_entangled_state_error': float(ent_error),
            'minimum_prepared_density_eigenvalue': float(eigvalsh(pre)[0]),
            'input_photon_probability': n_before, 'prepared_photon_probability': n_after,
            'reference_kernel_max_eigenvalue': float(eigvalsh(ar)[-1])})
    return {'finite_modes_only': True, 'cases': rows,
            'scope': 'Finite Nyström channel instances are checked CPTP. The proof, not the samples, establishes the continuum/all-code statement.'}


def continuum_checks():
    rows = []
    for i, sch in enumerate(schedules()):
        a = 1.3
        norm = np.sqrt(a/(-np.expm1(-a*sch.T)))
        pulse = lambda x: complex(norm*np.exp(-a*x/2 + .4j*x*x))
        # Four full 2x2 Fock-state master-equation blocks, in physical time.
        direct = matrix_fock(sch, pulse)
        # The simulation is integrated in compressed exposure time, including
        # loss-to-vacuum and Wiener phase correlations. It does not call the ODE.
        sim = continuum_simulated_metrics(sch, pulse)
        pe = abs(direct['p'] - sim['p'])
        me = abs(direct['m'] - sim['m'])
        assert max(pe, me, direct['trace_error'], direct['adjoint_error']) < 2e-9
        assert -1e-12 <= sim['prepared_photon_probability'] <= 1 + 1e-10
        rows.append({'schedule': i, 'direct_population': direct['p'],
            'direct_coherence_real': direct['m'].real,
            'direct_coherence_imag': direct['m'].imag,
            'independent_population_error': float(pe), 'independent_coherence_error': float(me),
            'prepared_photon_probability': sim['prepared_photon_probability'],
            'trace_error': direct['trace_error'], 'adjoint_error': direct['adjoint_error']})
    return {'cases': rows, 'scope': 'Physical-time Fock ODE versus independently integrated source preprocessing and maximum-rate receiver. No temporal-bin approximation is used in these integrals.'}


def fixed_source_check():
    """Prescribed square pulse: a counterexample to misreading source freedom."""
    T = 4.
    switch = 2*np.log(2)
    def g(t):
        return 1. if t <= switch else 1/(1+t-switch)
    def rhs(t, y):
        rate = g(t)
        return [-rate*y[0]/2 + np.sqrt(rate/T)]
    sol = solve_ivp(rhs, (0, switch), [0.], method='DOP853', rtol=1e-12, atol=1e-14)
    assert sol.success
    sol = solve_ivp(rhs, (switch, T), sol.y[:, -1], method='DOP853', rtol=1e-12, atol=1e-14)
    assert sol.success
    p_modulated = 1-(2*np.log(2)-1)/T
    p_constant = 4*(-np.expm1(-T/2))**2/T
    p_free_shape = -np.expm1(-T)
    err = abs(sol.y[0, -1]**2 - p_modulated)
    assert err < 1e-10 and p_constant < p_modulated < p_free_shape
    # Noiseless amplitude damping channel fidelities, no postselection.
    fidelity = lambda p: (3+p+2*np.sqrt(p))/6
    return {'dimensionless_duration': T, 'kappa': 0., 'gamma_phi': 0., 'gmax': 1.,
        'switch_time': float(switch), 'constant_coupling_fixed_square_probability': float(p_constant),
        'bounded_modulation_fixed_square_probability': float(p_modulated),
        'maximum_coupling_free_source_probability': float(p_free_shape),
        'direct_ODE_vs_modulation_formula_error': float(err),
        'fixed_square_constant_average_qubit_fidelity': float(fidelity(p_constant)),
        'fixed_square_modulated_average_qubit_fidelity': float(fidelity(p_modulated)),
        'scope': 'An explicit noiseless fixed-source example, not a new pulse protocol and not a violation of the joint-source optimization theorem.'}


def memory_boundary_check():
    """A colored-noise test of the *construction*, not a capacity counterexample."""
    # With max rate 1, reference times (0,1,2), actual times (0,1,3) are
    # consistent with a rate-1 segment followed by rate 1/2. Label zero is vacuum.
    u = np.array([0., 1., 2.])
    x = np.array([0., 1., 3.])
    delay = x-u
    markov = np.exp(-np.abs(delay[:, None]-delay[None, :]))
    # Quasistatic Gaussian detuning with variance sigma^2=1: C(t)=exp(-t^2/2).
    ratio = np.exp(-((x[:, None]-x[None, :])**2-(u[:, None]-u[None, :])**2)/2)
    a, b = sp.symbols('a b', real=True)
    symbolic = sp.Matrix([[1,1,a],[1,1,b],[a,b,1]]).det()
    assert sp.expand(symbolic+(a-b)**2) == 0
    assert eigvalsh(markov)[0] > -1e-12 and eigvalsh(ratio)[0] < -1e-3
    return {'markov_ratio_min_eigenvalue': float(eigvalsh(markov)[0]),
        'quasistatic_ratio_matrix': ratio.tolist(),
        'quasistatic_ratio_min_eigenvalue': float(eigvalsh(ratio)[0]),
        'quasistatic_ratio_determinant': float(np.linalg.det(ratio)),
        'exact_determinant_formula': '-(a-b)^2',
        'scope': 'The same time-warp/phase-noise construction is not automatically CP for quasistatic noise. This does NOT rule out other simulators or prove a capacity advantage for modulation in that changed model.'}


def main(output: Path):
    if not __debug__:
        raise RuntimeError('Keep assertions enabled; do not run python -O')
    start = time.time()
    result = {'status': 'PASS',
        'scope': 'Continuum channel ordering is proved analytically. Finite CPTP/reference tests and independently organized ODE/quadrature calculations are consistency checks, not experimental validation or a numerical capacity computation.',
        'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__, 'sympy': sp.__version__},
        'full_channel_simulation': factorization_checks(),
        'independent_continuum_dynamics': continuum_checks(),
        'source_freedom_boundary': fixed_source_check(),
        'Markov_assumption_boundary': memory_boundary_check()}
    result['elapsed_seconds'] = time.time()-start
    result['source_sha256'] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(p for d in ('capture','tests','tools') for p in (ROOT/d).glob('*.py'))}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status':'PASS', 'groups':4, 'elapsed_seconds':result['elapsed_seconds'], 'output':str(output)}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'artifacts/simulation.json')
    args = parser.parse_args()
    main(args.output)

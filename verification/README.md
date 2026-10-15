# Verification and evidence

Run `python tools/verify.py` from the repository root. Dependencies are pinned in requirements.txt. The local record uses Python 3.12; GitHub Actions uses Python 3.13. Assertions must remain enabled. No external dataset or account is used.

| Suite | Groups | Independently checked content |
|---|---:|---|
| `tests/check_simulation.py` | 4 | CPTP preparation and receiver Kraus maps; every finite matrix unit; reference-entangled inputs; physical-time Fock equations versus exposure-time simulation; fixed-source example; failure of the proposed colored-noise residual |
| `tests/check_control.py` | 5 | Finite Bessel spectrum versus integral/ODE calculations, including long-window principal roots; arbitrary schedules; finite-mode extensions; minimum time and tail asymptotics; finite-pulse coherent information |
| `tests/check_capacity.py` | 5 | Channel extension algebra; continuum eigenfunctions; full matrix Fock equations; nested Galerkin bounds; all-mode extensions; high-precision small-input witnesses |
| `tests/check_objectives.py` | 4 | Pulse eigenfunctions; scalar and full matrix dynamics; independent constrained optimization; exponential comparators and asymptotics |
| `tests/check_upper_bound.py` | 4 | Exact amplitude-damping distance certificates; multimode antidegradable comparisons; continuity-bound evaluation and amplitude-damping capacity checks; finite/infinite-window examples and boundary scaling |

Raw JSON reports are written to artifacts/ and include detailed tolerances and source hashes. `validated.json` records a successful run; `python tools/verify.py --record` refreshes it after all checks succeed. The verifier also checks Python syntax, local Markdown links and documentation formatting.

## What the evidence means

The continuum factorization, countable channel extensions and variational optima rest on analytical proofs. The checks compare independent finite formulations of those results: Nyström channels are checked to be CPTP before comparison, Galerkin eigenvalues supply spectral lower bounds, and high-precision coherent information supplies achievable rates for specified inputs. Capacity upper bounds use continuum spectral formulas or justified upper estimates.

## Reproducible record

The scientific implementations are in capture/. The five suites use independent ODE, operator, quadrature, algebraic and finite-space formulations. `validated.json` records the dependency versions and source hashes. Numerical evaluations of analytical capacity upper bounds use floating-point continuum roots; they are numerical evaluations rather than interval certificates.

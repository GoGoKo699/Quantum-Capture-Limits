# Verification and evidence

Run `python tools/verify.py` from the repository root. Dependencies are pinned in requirements.txt; the reference environment is Python 3.13. Assertions must remain enabled. No external dataset or account is used.

| Suite | Groups | Independently checked content |
|---|---:|---|
| `tests/check_simulation.py` | 4 | CPTP preparation and receiver Kraus maps; every finite matrix unit; reference-entangled inputs; physical-time Fock equations versus exposure-time simulation; fixed-source example; failure of the proposed colored-noise residual |
| `tests/check_control.py` | 5 | Finite Bessel spectrum versus integral/ODE calculations; arbitrary schedules; finite-mode extensions; minimum time and tail asymptotics; finite-pulse coherent information |
| `tests/check_capacity.py` | 5 | Channel extension algebra; continuum eigenfunctions; full matrix Fock equations; nested Galerkin bounds; all-mode extensions; high-precision small-input witnesses |
| `tests/check_objectives.py` | 4 | Pulse eigenfunctions; scalar and full matrix dynamics; independent constrained optimization; exponential comparators and asymptotics |

Raw JSON reports are written to artifacts/ and include detailed tolerances and source hashes. They are regenerated, not substituted from a cache. `validated.json` is a compact record of an actual run; `python tools/verify.py --record` refreshes it only after all checks succeed. The verifier also checks syntax, local Markdown links, and stale exploration labels.

## What the evidence means

The continuum factorization, countable channel extensions and variational optima rest on analytical proofs. Finite Nyström instances are checked to be CPTP before channel comparisons. Galerkin eigenvalues are lower bounds, not exact continuum eigenvalues. A numerical root residual is not exhaustive originality clearance. High-precision coherent information is an achievable lower bound for the specified input, not an optimized quantum capacity. None of these tests validates an experimental noise model or a bandwidth-constrained implementation.

The code distinguishes radiative population decay from extra optical-coherence decay; trace loss from postselection; channel simulation from excitation matching; and internal numerical evidence from independent scientific review.

## Reproducible record

The scientific engines are deduplicated in capture/. The four suites retain independently organized ODE, operator, quadrature, algebraic, and finite-space tests. `validated.json` records the versions and source hashes used. Changes in elapsed runtime do not signify a change in the scientific outputs. The full original-result comparison was performed during consolidation, excluding only metadata such as filenames, paths, hashes, and runtimes.

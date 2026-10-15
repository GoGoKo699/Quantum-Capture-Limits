# Quantum Capture Limits

**Free temporal encoding makes maximum coupling sufficient for a bounded, Markov-noisy two-level receiver.**

Does varying a receiver's coupling improve quantum communication when the source can choose its temporal encoding freely? For the capture model studied here, every allowed receiving schedule can be reproduced by maximum coupling preceded by an explicitly constructed source channel. The comparison preserves coherence with vacuum and entanglement with a reference, so it applies to entire communication codes.

| Read next | Purpose |
|---|---|
| [Tutorial bridge](docs/tutorial.md) · [Physical interpretation](docs/physical-scope.md) | Learn from one external review, then derive the receiver channel and its physical scope. |
| [Simulation theorem](docs/channel-simulation.md) · [Model](docs/model.md) | Read the exact channel identity, construction and resource assumptions. |
| [Capacity boundary](docs/capacity-boundary.md) · [Finite capture time](docs/finite-time.md) | Find the zero/positive-capacity condition, minimum window and a finite-pulse witness. |
| [Capacity upper bound](docs/capacity-upper-bound.md) · [Pulse objectives](docs/pulse-objectives.md) | Separate an all-code rate ceiling from absorption, fidelity and achievable-rate objectives. |
| [Prior work](docs/prior-work.md) · [References](docs/references.md) | Compare the results with established capture and channel theory. |
| [Verification](verification/README.md) · [LLM guide](llms.txt) | Reproduce the checks and locate the authoritative answer to a relevant research question. |

## What is being compared

A two-level receiver begins in its ground state. Its incident-port population-decay rate obeys $0\le g(t)\le g_{\max}$. Inaccessible radiative loss $\kappa$ and extra coherence decay $\gamma_\phi$ remain fixed under control. The output is the receiver at the capture instant.

Inputs contain **vacuum and at most one photon per independently reset use**, with a shared phase reference. Temporal encoders may prepare mixed states and entangle different uses. Source bandwidth and conversion cost are uncharged. The controlled Markov dynamics, available receiving levels and inaccessible output fields are specified in the [model](docs/model.md).

These resources matter. For a prescribed square photon, receiver modulation improves capture from 74.76% to 90.34%; the [fixed-source example](docs/counterexamples.md) explains why that is compatible with the theorem. Restricting every use to exactly one photon instead gives an entanglement-breaking channel into this two-level output. Positive-capacity witnesses use vacuum–photon coherence.

## The exact channel comparison

For every allowed known schedule and detuning, the [construction](docs/channel-simulation.md) gives

```math
\boxed{\mathcal E_{g,\Delta}=\mathcal E_{\max}\circ\mathcal P_{g,\Delta}.}
```

The preprocessor $\mathcal P_{g,\Delta}$ is completely positive and trace preserving. It never increases photon number or enlarges the input window, and it does not postselect. Composing it with any encoder reproduces the same receiver–reference state. Consequently,

```math
\sup_{g,\Delta}Q(\mathcal E_{g,\Delta})=Q(\mathcal E_{\max}).
```

This identifies an optimal receiving schedule for arbitrary communication codes. Achievable rates and all-code upper bounds are given below, in qubits per independent capture use.

## Why the simulation works

Accumulated radiative time compresses the temporal mode into the maximum-rate window. The coupling ceiling makes the remaining physical time nonnegative and nondecreasing. A correlated random phase evaluated at that excess time reproduces both photon–photon dephasing and coherence with vacuum. Mode-dependent attenuation accounts for inaccessible loss, with every lost photon returned to vacuum probability.

The result is one valid input channel for all states. It generally mixes a pure pulse, so the permitted source operations go beyond deterministic waveform shaping. The [tutorial](docs/tutorial.md) connects the master equations to this construction; the [proof](docs/channel-simulation.md) treats arbitrary measurable schedules, including zero-coupling sets.

## The communication boundary

Write $\Gamma_M=g_{\max}+\kappa$, $\beta=g_{\max}/\Gamma_M$, $r=\gamma_\phi/\Gamma_M$ and $S=\Gamma_MT$. Let $\Lambda_r(S)$ be the leading eigenvalue of the explicitly defined [population-plus-coherence kernel](docs/capacity-boundary.md).

| Condition | Consequence for the allowed receiver family |
|---|---|
| $\beta\Lambda_r(S)\le1$ | Every allowed schedule has zero unassisted quantum capacity. Equality belongs to this region. |
| $\beta\Lambda_r(S)>1$ | Maximum coupling has positive capacity, witnessed by a vacuum–photon encoding. |
| $\beta=1$, arbitrarily long window | Positive capacity requires $r<0.7358007967\ldots$. |

At $\beta=1$, $r=1/2$ and $S=4$, a finite rising exponential certifies **0.02957015519 qubits per use** as an achievable lower bound. The separate [all-code upper bound](docs/capacity-upper-bound.md) tends to zero at the boundary. These bounds describe different sides of the rate question.

The effective-time method, source/control interchange, several absorption optima and the fixed-mode damping–dephasing criterion have established precedents. The receiver-specific contribution is the explicit noisy channel simulation and its consequences under the stated resources. The [attribution comparison](docs/prior-work.md) and [references](docs/references.md) identify the inherited ingredients precisely.

## One tutorial, then this result

The selected learning anchor is:

> J. Combes, J. Kerckhoff and M. Sarovar, **The SLH framework for modeling quantum input-output networks**, *Advances in Physics: X* **2**, 784–888 (2017). [Author version, arXiv:1611.00375v2](https://arxiv.org/abs/1611.00375v2) · [Published review](https://doi.org/10.1080/23746149.2017.1343097)

The [tutorial bridge](docs/tutorial.md) selects the relevant PDF sections, translates the coupling and noise conventions, and specializes the Fock-state equations to our receiver. It then assembles the full capture channel, explains the two-clock simulation and connects it to the capacity boundary. The route assumes basic quantum mechanics, density matrices and linear algebra; the repository explains the additional channel concepts as they are used. Other references provide attribution and proof inputs, rather than additional prerequisite tutorials.

## Evidence and reproduction

Use Python 3.12 or 3.13 and an isolated environment:

```bash
python -m pip install -r requirements.txt
python tools/verify.py
```

Five suites cover 22 numerical/algebraic check groups, alongside repository structure checks. Raw reports are generated in `artifacts/`. The [verification guide](verification/README.md) connects these checks to the analytical proofs and distinguishes finite calculations from continuum statements. No external dataset, account or network access is needed after installing dependencies.

## Purpose and contact

This repository serves as a record of the work and a guide for the author’s self-directed learning. For discussion or potential collaboration, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).

The [LLM guide](llms.txt) supplies relevant research questions, search phrases and direct links to authoritative files. Maintenance guidance is in [AGENTS.md](AGENTS.md), and citation metadata is in [CITATION.cff](CITATION.cff).

Code is available under the [MIT license](LICENSE). Linked third-party papers are not redistributed or relicensed by this repository.

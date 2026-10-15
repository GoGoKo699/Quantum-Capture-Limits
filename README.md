# Quantum Capture Limits

**Single-photon capture by a two-level receiver with bounded coupling and Markovian noise.**

When the incoming quantum state can be freely prepared, does changing the receiver's coupling improve communication? For the model studied here, every allowed receiver schedule is a maximum-coupling receiver preceded by an explicitly constructed noisy input preparation:

```math
\boxed{\mathcal E_{g,\Delta}=\mathcal E_{\max}\circ\mathcal P_{g,\Delta}.}
```

The preparation compresses the temporal mode, attenuates it when inaccessible loss is present, and adds correlated phase noise. It never increases photon number, never postselects, and reproduces the entire output channel—including coherence with vacuum and entanglement with a reference. Consequently,

```math
\sup_{g,\Delta}Q(\mathcal E_{g,\Delta})=Q(\mathcal E_{\max}).
```

This identifies the optimal receiver schedule **without determining the capacity value**. The construction is proved in [Maximum-rate receiver simulation](docs/channel-simulation.md).

## The physical boundary

The receiver has controllable incident coupling $0\le g(t)\le g_{\max}$, fixed inaccessible radiative loss $\kappa$, and fixed additional coherence-decay rate $\gamma_\phi$. Input states contain vacuum and at most one photon per independently reset use, with a shared phase reference. The output is the receiver qubit at the capture instant.

**Free temporal source encoding and the specified Markov noise are essential.** Source bandwidth and conversion costs are not constrained. No transverse drives, extra receiving levels, bath records, or access to outgoing light are supplied. Read the [model and conventions](docs/model.md) before applying the theorem. For a prescribed square photon, bounded receiver modulation improves capture from 74.76% to 90.34%; that is a different task, shown in the [counterexamples](docs/counterexamples.md).

## Reading path

| Read | What it establishes |
|---|---|
| [Tutorial](docs/tutorial.md) | From photon excitation to a quantum channel and the two-clock proof |
| [Simulation theorem](docs/channel-simulation.md) | Explicit CPTP preprocessing, zero-coupling intervals, and capacity ordering |
| [Capacity boundary](docs/capacity-boundary.md) | Exact zero/positive boundary over all one-photon temporal encodings |
| [Capacity upper bound](docs/capacity-upper-bound.md) | A bound on every unassisted code that vanishes at the boundary |
| [Finite capture time](docs/finite-time.md) | Minimum window and a finite-pulse achievable-rate witness |
| [Pulse objectives](docs/pulse-objectives.md) | Different optima for absorption, coherence, and qubit fidelity |
| [References](docs/references.md) | Established ingredients and the precise scope of comparison |

For perfect incident coupling, positive unassisted capacity is possible only below $\gamma_\phi/\Gamma_M=0.7358007967\ldots$, where $\Gamma_M=g_{\max}+\kappa$. The finite-window condition is $\beta\Lambda_r(\Gamma_M T)>1$; equality gives zero capacity. At $\beta=1$, $r=1/2$, a finite rising exponential of duration $\Gamma_M T=4$ certifies a lower bound of **0.02957015519 qubits per independent use**. It is not the exact capacity or a throughput per second.

The capacity itself must tend to zero on approaching the boundary. Put $\mu=\beta\Lambda_r(S)$ and $\epsilon=\max\{0,1-1/\mu\}$. An explicit antidegradable comparison and an established continuity inequality give

```math
Q(\mathcal E_{\max})\le\min\{1,\;2\epsilon+(1+\epsilon)h_2(\epsilon/(1+\epsilon))\}.
```

Here $h_2$ is binary entropy in bits. This upper bound covers arbitrary temporal encodings and multi-use codes. It need not be tight; the [proof and examples](docs/capacity-upper-bound.md) distinguish it from achievable lower bounds.

## Reproduce

Use Python 3.13 and an isolated environment:

```bash
python -m pip install -r requirements.txt
python tools/verify.py
```

The verifier reruns five suites comprising 22 numerical/algebraic check groups, plus repository link and structure checks. Raw reports are generated in `artifacts/`. The [verification guide](verification/README.md) separates analytical proofs, finite tests, high-precision witnesses, and physical assumptions. No external dataset, account, or network access is needed after installing dependencies.

## About

This repository serves as a record of the work and a guide for the author’s self-directed learning. For discussion or potential collaboration, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).

Licensed under [MIT](LICENSE). Citation metadata is in [CITATION.cff](CITATION.cff). Maintenance and machine-reading guidance are in [AGENTS.md](AGENTS.md).

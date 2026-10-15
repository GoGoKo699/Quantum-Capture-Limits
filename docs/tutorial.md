# A route from photon capture to channel simulation

The central question is not whether a photon can excite an atom. It is whether an unknown quantum state can be transferred, and which control resources can improve that transfer. This tutorial connects those questions to the [main proof](channel-simulation.md).

## 1. Specify a channel use

At the start, a two-level receiver is in its ground state. An input during $[-T,0]$ contains vacuum and at most one photon. The output is the receiver at time zero, not the whole outgoing field and not a permanently stored qubit. Single-rail coherence means a superposition $\alpha|\mathrm{vac}\rangle+\zeta|1_f\rangle$, so sender and receiver need a phase reference.

The basic distinction is between **radiative population decay** and **additional coherence decay**. If total radiative decay is $\Gamma$ and the extra dephasing rate is $\gamma_\phi$, optical coherence decays at $\Gamma/2+\gamma_\phi$. Confusing these conventions changes every threshold.

For the dynamical foundation, use one primary anchor: Baragiola, Cook, Brańczyk and Combes, [N-Photon wave packets interacting with an arbitrary quantum system](references.md#b12). Focus on the four auxiliary density operators for photon indices zero and one. `capture.receiver.matrix_fock` implements those equations directly with two-by-two matrices.

## 2. Keep both population and coherence

A single chosen input mode gives

```math
\mathcal E_f(\rho)=\begin{pmatrix}
\rho_{00}+(1-p)\rho_{11}&m^*\rho_{01}\\
m\rho_{10}&p\rho_{11}
\end{pmatrix}.
```

$p$ is excitation probability; $m$ is coherent transfer amplitude. They are not independent: complete positivity requires $|m|^2\le p\le1$. A receiver with high $p$ but small $m$ can be a poor qubit receiver.

The distinction between particle survival and coherence with vacuum is established [AO07](references.md#ao07). Indeed, if every use is restricted to exactly one photon, this receiver produces only diagonal energy-basis outputs and the restricted channel is entanglement breaking. The vacuum sector supplies the coherent qubit degree of freedom in the positive-rate witnesses; see the [model](model.md#why-the-vacuum-sector-matters).

After removing a known phase, the uniformly averaged qubit fidelity is $(3+p+2|m|)/6$. This is a single-use performance measure. The condition $p+|m|^2>1$ instead determines positive unassisted quantum capacity for this channel family. These criteria answer different questions; the latter uses coding over independent uses and is attributed in [Capacity boundary](capacity-boundary.md).

**Check:** at $r=3$, a fidelity slightly above $2/3$ coexists with zero unassisted capacity for the entire allowed temporal receiver. This is not a contradiction and does not imply an entanglement-breaking channel.

## 3. Replace one pulse with the full temporal space

A pulse is a vector in $L^2([0,T])$. An arbitrary input is a density operator on $\mathbb C|\mathrm{vac}\rangle\oplus L^2([0,T])$. Its one-photon block can mix temporal modes; its off-diagonal block carries coherence with vacuum.

Two objects characterize capture: a positive kernel $K$ for population and a vector $h$ for coherence. Their inequalities $0\preceq|h\rangle\langle h|\preceq K\preceq I$ make the entire receiver a valid quantum channel. Read the block map in [Model and conventions](model.md), then compare its explicit finite Kraus representation in `capture.simulation.receiver_kraus`.

This step is what prevents an optimum over one pulse family from being advertised as a theorem about every encoding.

## 4. Use two clocks

Physical time measures noise exposure. Accumulated radiative time measures how much coupling has acted:

```math
u=\int_0^x[g(s)+\kappa]ds,\qquad
\Gamma_M=g_{\max}+\kappa,\qquad
d(u)=x(u)-u/\Gamma_M.
```

The coupling ceiling makes $d$ nonnegative and nondecreasing. Therefore

```math
|x(u)-x(v)|=|u-v|/\Gamma_M+|d(u)-d(v)|.
```

That equality permits the residual noise to be represented by an actual quantum channel. A mere inequality between kernels would be insufficient to establish full channel simulation.

Integrated-control clocks already appear in optimized quantum memories [G07C](references.md#g07c). The additional step here is constructing the residual noise on vacuum and all temporal modes simultaneously, as explained in the [prior-work comparison](prior-work.md).

Temporal reshaping preserves norm through the Jacobian $f(x(u))/\sqrt{g(x(u))+\kappa}$. Completely uncoupled components become vacuum; they are not conditionally discarded and renormalized. The generalized inverse retains time spent in gaps, so the proof does not forget dephasing while coupling is off.

## 5. Construct the preparation, not only its desired output

After temporal compression, attenuate amplitudes by the useful-port ratio and add Wiener phases evaluated at $d(u)$. The phase correlations reproduce both photon–photon coherence and vacuum–photon coherence. Each realization is unitary; averaging it is completely positive. Attenuation sends the missing photon weight to vacuum, making the preparation trace preserving.

The resulting identity is

```math
\mathcal E_{g,\Delta}=\mathcal E_{\max}\circ\mathcal P_{g,\Delta}.
```

To verify it, substitute the prepared photon and coherence blocks into the maximum-rate map. To understand its consequence, take any encoder and absorb $\mathcal P$ into it. The decoder receives exactly the same joint state, even with an external reference. No capacity formula or additivity assumption is required.

**Check:** this is input preprocessing. It does not assert that a qubit-output degradation map exists in the opposite position.

## 6. Test the resource assumptions

For a fixed square photon, tunable coupling improves capture. Source-optimized maximum coupling does better still. The three numbers compare different resource sets, as shown in [Counterexamples](counterexamples.md).

For quasistatic Gaussian detuning, the proposed residual phase-correlation matrix has a negative determinant. That invalidates this specific extension of the proof, not every conceivable simulator and not all channel ordering under colored noise.

## 7. Separate a rate witness from a rate ceiling

A positive coherent-information example proves that some communication rate is achievable. A very small example does not prove that every code is slow. To obtain a ceiling, let $\mu=\|A+|b\rangle\langle b|\|>1$ and attenuate the receiver's output with photon transmission $1/\mu$. The resulting comparison channel lies exactly on the antidegradable boundary. Its distance from the original receiver is at most $1-1/\mu$ in half diamond norm, including any reference system.

An established continuity inequality then bounds the original receiver's full quantum capacity by a function that tends to zero with this distance. Read [Capacity upper bound](capacity-upper-bound.md) for the finite-code argument: it allows mixed encoders and infinitely many temporal modes without an additivity assumption. The comparison is a mathematical tool; it does not claim that adding attenuation improves capture.

The reading order is [model](model.md) → [proof](channel-simulation.md) → [counterexamples](counterexamples.md) → [capacity boundary](capacity-boundary.md) → [capacity upper bound](capacity-upper-bound.md) → [finite-time examples](finite-time.md). The [verification guide](../verification/README.md) connects each analytical claim to its independent checks.

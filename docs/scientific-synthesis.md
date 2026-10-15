# A channel comparison for bounded noisy capture

Photon capture can be improved by changing either the incoming waveform or the receiver control. Earlier memory and absorption results already establish this interchange in important regimes. The question addressed here is whether receiver modulation adds a communication capability when temporal source encoding is unrestricted and the receiver remains subject to fixed Markov noise. The answer follows from an explicit comparison of the entire noisy channel, including its action on an unknown state entangled with a reference.

## Central statement and resources

Consider a two-level receiver initially in its ground state, with incident population-decay rate $0\le g(t)\le g_{\max}$, fixed inaccessible vacuum loss $\kappa$, and fixed extra coherence decay $\gamma_\phi$. Its output is evaluated at the end of a window of length $T$. The input is vacuum plus at most one photon per independently reset use, with a shared phase reference and unrestricted temporal encoders. The controlled Markov dilation, coupling ceiling and noise rates are held fixed across the comparison.

For every allowed known schedule and longitudinal detuning, the [simulation theorem](channel-simulation.md) constructs a normal CPTP map satisfying

```math
\mathcal E_{g,\Delta}=\mathcal E_{\max}\circ\mathcal P_{g,\Delta}.
```

The preprocessor neither creates photons nor enlarges the window. It can be absorbed into any source encoder, including one used in an entangled multi-use code. Hence maximum coupling attains the supremum of unassisted capacity over the allowed schedules. The theorem identifies an optimal receiver, not a capacity-achieving pulse or the positive capacity value.

## Why the construction works

Accumulated radiative time $u=\int_0^x[g(s)+\kappa]ds$ compresses each schedule into the maximum-rate window. Put $\Gamma_M=g_{\max}+\kappa$ and $d(u)=x(u)-u/\Gamma_M$. The ceiling makes $d$ nonnegative and nondecreasing. Wiener phases at this excess time reproduce both residual photon–photon dephasing and residual coherence with vacuum. A mode-dependent attenuation channel accounts for inaccessible loss, returning all lost weight to vacuum. The generalized inverse retains elapsed physical time across zero-coupling sets.

Clock reshaping is an established technique. The substantive construction is one completely positive residual map that matches both capture kernels simultaneously. It generally mixes pure pulses, so arbitrary source channels are a stronger resource than deterministic waveform shaping. The [prior-work comparison](prior-work.md) identifies the precise antecedents and overlapping formulas.

## Consequences for communication

Write the maximum-rate channel's excitation effect as $A$ and vacuum-coherence vector as $b$. An explicit symmetric two-output extension gives the [all-mode boundary](capacity-boundary.md):

```math
Q>0\quad\Longleftrightarrow\quad
\mu=\|A+|b\rangle\langle b|\|>1.
```

The positive direction uses the established fixed-mode damping–dephasing criterion. For this receiver, $\mu=\beta\Lambda_r(S)$, with $\beta=g_{\max}/\Gamma_M$, $r=\gamma_\phi/\Gamma_M$, and $S=\Gamma_MT$. The spectral solution gives an exact finite-window condition and the half-line threshold $r=0.7358007967\ldots$ at $\beta=1$. Equality gives zero capacity. A pulse maximizing absorption can have zero capacity even when another pulse has positive capacity.

Explicit finite pulses give achievable lower bounds. Separately, output attenuation by $1/\mu$ when $\mu>1$ produces an antidegradable comparison channel. Established finite-output continuity yields the [all-code upper bound](capacity-upper-bound.md)

```math
Q\le\min\{1,2\epsilon+g_2(\epsilon)\},\qquad
\epsilon=\max\{0,1-1/\mu\},\qquad
g_2(\epsilon)=(1+\epsilon)\log_2(1+\epsilon)-\epsilon\log_2\epsilon.
```

Thus full capacity tends to zero at the boundary, even though its positive value and exact asymptotic scaling remain undetermined. All rates are qubits per reset use.

## Physical interpretation

Within this source-flexible model, modulation cannot restore positive capacity in the region $\mu\le1$. Improvements there require changing a stated resource or assumption. Fixed sources, finite conversion bandwidth, control-dependent noise, extra receiving levels, outgoing-field access and storage after capture describe different comparisons. Exactly-one-photon inputs in every use give a binary measure–prepare channel into this two-level output; the positive-rate witnesses use vacuum–photon coherence. The [physical scope](physical-scope.md) relates the idealization to direct waveguide coupling without treating energy-loading measurements as capacity validation. The [verification record](../verification/README.md) supports reproducibility; it is neither an experimental test nor a certificate of priority.

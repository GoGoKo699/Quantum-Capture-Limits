# Maximum-rate receiver simulation

Under the [specified control and source resources](model.md), every controlled capture channel factors as

```math
\boxed{\mathcal E_{g,\Delta}=\mathcal E_{\max}\circ\mathcal P_{g,\Delta}.}
```

The input map is normal, completely positive and trace preserving, never increases photon number, and preserves the available input time window. It simulates vacuum coherence and reference entanglement, not just capture probability. Therefore maximum coupling attains the supremum of the unassisted capacity over allowed receiver schedules when temporal source encoding is free. This orders capacity values without computing the value.

Assume measurable bounded incident coupling, integrable known detuning, fixed vacuum radiative loss and fixed Markov dephasing. Completely uncoupled measurable sets are discarded to vacuum; the accumulated-rate generalized inverse retains elapsed time across gaps. The finite-interval implementation below uses piecewise-constant schedules as tests, not as a restriction of the theorem.

## Constructing the simulator

### Compress to accumulated radiative time

Change coordinate to

```math
u=R(x),\qquad 0\le u\le U=R(T)\le S.
```

Where $\Gamma>0$, use the inverse $x=x(u)$ and the isometry

```math
\psi(u)=\frac{f(x(u))}{\sqrt{\Gamma(x(u))}}.
```

Extend its output by zero on $(U,S]$. This reshapes the incoming photon without adding an excitation. A deterministic phase mask multiplies $\psi(u)$ by $e^{-iD(x(u))}$.

If $g=\kappa=0$ on an interval, incoming components there never interact with the receiver. Discard that one-photon subspace to vacuum, including its coherence with the retained subspace, before making the isometry. The lost weight is retained in the vacuum population; it is not conditionally normalized away. The generalized inverse of R has jumps across those intervals. Those jumps are retained in the extra time below, so dephasing during an idle period is not omitted.

### The two remaining penalties

Define

```math
\beta(u)=\frac{g(x(u))}{g(x(u))+\kappa},\qquad
c(u)=\sqrt{\frac{\beta(u)}{\beta_M}}\le1,
```

and the excess physical time

```math
\boxed{d(u)=x(u)-\frac{u}{\Gamma_M}.}
```

The coupling cap implies that d is nonnegative and nondecreasing. On smooth portions,

```math
d'(u)=\frac1{\Gamma(x(u))}-\frac1{\Gamma_M}\ge0.
```

Thus for any ordered u and v,

```math
|x(u)-x(v)|=\frac{|u-v|}{\Gamma_M}+|d(u)-d(v)|.
```

This additive identity is why the residual dephasing is an actual channel rather than merely a positive numerical penalty in one objective.

The factor c is attenuation of the photon amplitude. Its trace-preserving channel on vacuum plus one photon is

```math
\rho_{11}\mapsto C\rho_{11}C,\qquad
\rho_{01}\mapsto\rho_{01}C,\qquad
\rho_{00}\mapsto\rho_{00}+\operatorname{Tr}[(I-C^2)\rho_{11}],
```

where C multiplies wavefunctions by c(u). This is ordinary loss into an unobserved ancillary port, not postselection. If $\kappa=0$, c=1 wherever g>0: in that case the only loss added by this step is from completely uncoupled temporal components.

### Extra waiting is a valid phase-noise channel, including vacuum coherence

Let $B_s$ be a standard real Wiener process with $B_0=0$, and let

```math
V_B|\mathrm{vac}\rangle=|\mathrm{vac}\rangle,\qquad
(V_B\psi)(u)=e^{i\sqrt{2\gamma_\phi}B_{d(u)}}\psi(u).
```

Each $V_B$ is unitary on vacuum plus the one-photon space. Averaging $V_B\rho V_B^\dagger$ therefore gives a normal CPTP map. Its matrix elements satisfy

```math
\rho_{11}(u,v)\mapsto e^{-\gamma_\phi|d(u)-d(v)|}\rho_{11}(u,v),
```

```math
\rho_{01}(u)\mapsto e^{-\gamma_\phi d(u)}\rho_{01}(u).
```

The second line is essential. Preserving photon-sector populations alone would not simulate the single-rail quantum channel. The covariance with the vacuum is fixed by taking its phase label to be d=0.

This random-unitary representation is a mathematical implementation of the required dephasing channel; it does not assert that a particular laboratory bath is classical or reveal a noise record to the decoder. Attenuation remains a separate process and is not reclassified as classical phase noise.

The simulator $\mathcal P_g$ is the composition of this phase channel, attenuation, known phase mask, and time-coordinate map. Every step is CPTP and photon-number nonincreasing.

### Direct proof of the channel identity

After time reshaping and removal of the known phase, the controlled kernels are

```math
\widetilde h_g(u)=c(u)e^{-\gamma_\phi d(u)}h_M(u),
```

```math
\widetilde K_g(u,v)=c(u)c(v)
 e^{-\gamma_\phi|d(u)-d(v)|}K_M(u,v).
```

Applying the adjoints of the attenuation and phase-noise channels to $h_M,K_M$ gives precisely these two expressions. Transforming back restores $D(x)$ and the $\sqrt\Gamma$ Jacobians. Therefore every block in the full capture map agrees:

```math
\mathcal E_g(\rho)=\mathcal E_{\max}(\mathcal P_g(\rho)).
```

Both sides are normal channels on the same finite-window vacuum-plus-one-photon space. Their equality extends to arbitrary trace-class input states and, by tensoring with identity, to inputs entangled with an arbitrary reference. The proof is independent of any finite temporal discretization.

For a constant submaximum rate $g_0$, the interpretation is especially direct. Put $\Gamma_0=g_0+\kappa$, $\beta_0=g_0/\Gamma_0$. The input is compressed by $\Gamma_0/\Gamma_M$, its photon transmission probability is $\beta_0/\beta_M$, and the extra dephasing in exposure time has coefficient

```math
\gamma_\phi\left(\frac1{\Gamma_0}-\frac1{\Gamma_M}\right).
```

No numerical optimization is needed to identify the comparison channel.

## Capacity values and finite-code consequences

The information-theoretic inference is standard channel simulation/data processing [J21](references.md#j21). It must be distinguished from the receiver-specific construction above.

Suppose an encoder $\mathcal A_n$, n controlled receiver uses, and decoder $\mathcal D_n$ achieve a specified entanglement fidelity. Replace each controlled use by maximum coupling and incorporate $\mathcal P_g^{\otimes n}$ into the encoder. The final reference–decoder state is unchanged, exactly. If different deterministic schedules are used in different slots, use the corresponding product of preprocessors.

The preprocessing does not increase the input photon number or support outside the available window, so it respects both the at-most-one-photon-per-use restriction and an additional upper bound on mean photon number. The encoder may discard systems, so the fact that preprocessing can convert a pure pulse to a mixed temporal state is allowed. No classical side information about the added phases is delivered to the receiver.

Consequently every unassisted quantum code for the controlled receiver can be reproduced using the maximum-rate receiver. Hence

```math
Q(\mathcal E_g)\le Q(\mathcal E_{\max}).
```

Maximum coupling is itself an allowed schedule, proving equality after the control supremum. This is stronger than optimizing $p+|m|^2$, stronger than a single-letter coherent-information comparison, and independent of possible nonadditivity or an unknown optimal code.

It does **not** imply that one fixed deterministic wavepacket at maximum coupling reproduces every controlled one-mode code. The simulator can mix temporal modes. The capacity comparison concerns the full input space and standard encoding freedom specified in the [model](model.md). Nor does it compare the same fixed incoming state at two schedules.

The same simulation lets the sender reproduce classical messages and other receiver-output communication tests when their encoding resources include this preprocessing. No general claim is made about device-control costs, noise sensing, access to environmental outputs, or protocol classes that prohibit such encoding.

The exact zero/positive boundary is a corollary. With $\Lambda_r(S)$ the largest eigenvalue of

```math
L_r(u,v)=e^{-(u+v)/2}[e^{-r|u-v|}+e^{-r(u+v)}],
```

positive Q is possible exactly when $\beta_M\Lambda_r(S)>1$. The qubit damping–dephasing criterion and logarithmic-entropy argument behind that sign classification are established [SAJS24](references.md#sajs24). The current result neither reevaluates $\Lambda_r(S)$ nor solves the positive capacity value. In particular, the [finite pulse](finite-time.md) with an achievable lower bound near 0.0296 qubits/use does not become an exact capacity by this theorem.

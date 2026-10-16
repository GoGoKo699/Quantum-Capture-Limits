# From quantum input-output theory to capture limits

This tutorial connects one external teaching source to the results in this repository: Combes, Kerckhoff and Sarovar, **The SLH framework for modeling quantum input-output networks** [CKS17](references.md#cks17). Familiarity with state vectors, density operators, partial traces and ordinary differential equations is sufficient to begin. The review supplies the field–system dynamics; the steps below turn those dynamics into a channel comparison and a communication limit.

The physical question is simple: if the sender can choose the temporal encoding, can changing a noisy receiver's coupling outperform holding it at its maximum? Answering it requires following both excitation probability and coherence with vacuum. An absorption calculation alone does not describe transfer of an unknown qubit.

Consider a qubit whose logical zero is vacuum and whose logical one is a photon with a chosen temporal envelope. The receiver starts in its ground state. Capturing the logical one can excite it, while transferring a superposition also requires coherence between its ground and excited states. These two responses become the population parameter $`p`$ and coherence amplitude $`m`$ below. The final worked example evaluates both for a finite pulse.

## 1. A focused route through the review

Use the original PDF's section and equation numbers; the experimental HTML numbers some equations differently. These selections form one reading route within the same source.

| Read in CKS17 | What to bring into this repository |
|---|---|
| §III.A; §IV | Input-output fields, the Markov approximation and quantum stochastic increments |
| §V.A; §V.D | The joint unitary specified by an SLH triple and the reduced master equation |
| §VII.A.3 | Temporal photon wavepackets and superpositions or mixtures with vacuum |
| §VII.B.1, Eqs. (153)–(155) | The auxiliary Fock-state operators that retain number-sector coherence |
| §VII.C; Appendix A, item 12, Eq. (A13) | Inaccessible loss ports and a waveguide-coupled two-level atom |

The review's broader network applications are useful context but unnecessary for this route. Quantum capacity is introduced below and developed in the repository's proofs; another external tutorial is not a prerequisite.

## 2. Translate the receiver into an SLH model

An SLH triple specifies a field scattering matrix $`\mathbf S`$, a vector of system–field coupling operators $`\mathbf L`$, and a system Hamiltonian $`H`$. For this receiver, $`g(t)`$ is the incident-port population-decay rate, $`\kappa`$ is the inaccessible radiative-loss rate, $`\gamma_\phi`$ is extra coherence decay, and $`\Delta(t)`$ is known detuning. The rates $`\kappa`$ and $`\gamma_\phi`$ stay fixed under control.

Write $`P_e=|e\rangle\langle e|`$ and $`\sigma_-=|g\rangle\langle e|`$. The receiver's triple is

```math
\mathbf S=I_3,\qquad
\mathbf L(t)=
\begin{pmatrix}
\sqrt{g(t)}\,\sigma_-\\
\sqrt\kappa\,\sigma_-\\
\sqrt{2\gamma_\phi}\,P_e
\end{pmatrix},\qquad
H(t)=\Delta(t)P_e.
```

The first port carries the signal, which may contain vacuum and one photon. The second and third input ports are in vacuum and are inaccessible. Squaring a coupling amplitude gives its dissipator coefficient, so the population-decay couplings contain $`\sqrt{g(t)}`$ and $`\sqrt\kappa`$. The third coupling produces pure dephasing rather than an additional population-decay channel. The scalar $`S=\Gamma_M T`$ used later is a dimensionless window length, distinct from the scattering matrix $`\mathbf S`$.

| Review language | Receiver notation and meaning |
|---|---|
| Coupling operator $`L`$ | Incident coupling $`\ell(t)=\sqrt{g(t)}\sigma_-`$ |
| Additional vacuum inputs | Fixed radiative loss $`\kappa`$ and extra coherence decay $`\gamma_\phi`$ |
| Temporal envelope $`\xi(t)`$ | A normalized photon on $`[-T,0]`$; $`f(x)=\xi(-x)`$ in lookback time |
| Local system state | Receiver qubit, initially $`\lvert g\rangle`$ and evaluated at time zero |
| Joint unitary and traced fields | A channel from the input field to the final receiver |

Define $`\mathcal D[L]X=LXL^\dagger-\{L^\dagger L,X\}/2`$. The vacuum generator is

```math
\mathcal L_t(X)=-i[\Delta(t)P_e,X]
 +[g(t)+\kappa]\mathcal D[\sigma_-]X
 +2\gamma_\phi\mathcal D[P_e]X.
```

Thus population decays at $`\Gamma(t)=g(t)+\kappa`$, whereas optical coherence decays at $`\Gamma(t)/2+\gamma_\phi`$. The factor of two in the dephasing dissipator matters.

The Markov assumption specifies the system's coupling to white-noise fields for every allowed schedule. A photon wavepacket nevertheless correlates different arrival times: it is a coherent superposition over them. Those input correlations explain the auxiliary hierarchy below; they do not introduce a colored bath into the model. Likewise, an uncontrolled receiver's measured decay constants alone do not specify its dynamics under arbitrary control. Read the [model's resource assumptions](model.md) with the SLH dilation in mind.

## 3. Follow four operators, not one driven density matrix

Specialize the review's Fock hierarchy to photon numbers zero and one. In the convention implemented by [`matrix_fock`](../capture/receiver.py), the auxiliary operators satisfy

```math
\begin{aligned}
\dot q_{00}&=\mathcal L_t(q_{00}),\\
\dot q_{10}&=\mathcal L_t(q_{10})+
 \xi(t)[q_{00},\ell^\dagger(t)],\\
\dot q_{01}&=\mathcal L_t(q_{01})+
 \xi(t)^*[\ell(t),q_{00}],\\
\dot q_{11}&=\mathcal L_t(q_{11})+
 \xi(t)[q_{01},\ell^\dagger(t)]
 +\xi(t)^*[\ell(t),q_{10}].
\end{aligned}
```

At $`t=-T`$,

```math
q_{00}=q_{11}=|g\rangle\langle g|,
\qquad q_{10}=q_{01}=0.
```

For vacuum input the receiver state is $`q_{00}`$; for a photon it is $`q_{11}`$. The other two operators describe the responses to $`|1_\xi\rangle\langle\mathrm{vac}|`$ and its adjoint. They are not normalized states and must not be interpreted as extra populations. For input $`\alpha|\mathrm{vac}\rangle+\zeta|1_\xi\rangle`$, the physical output is

```math
|\alpha|^2q_{00}+\zeta\alpha^*q_{10}
 +\alpha\zeta^*q_{01}+|\zeta|^2q_{11}.
```

In the physical basis $`(|g\rangle,|e\rangle)`$, put $`q_{10}=-a(t)|e\rangle\langle g|`$ and $`P(t)=\langle e|q_{11}|e\rangle`$. The hierarchy reduces to

```math
\dot a=-[\Gamma/2+\gamma_\phi+i\Delta]a+\sqrt g\,\xi,
\qquad
\dot P=-\Gamma P+2\sqrt g\,\mathrm{Re}(\xi^*a),
\qquad a(-T)=P(-T)=0.
```

The code reports $`p=P(0)`$ and $`m=a(0)=-q_{10}(0)_{eg}`$. Its minus sign is a fixed output phase convention. Applying the known phase correction $`|e\rangle\mapsto-|e\rangle`$, or equivalently using $`|1\rangle_{\rm out}=-|e\rangle`$, gives the channel convention used below, with $`m\rho_{10}`$ in its lower off-diagonal entry. The amplitude $`m`$ can be complex; this calibration does not change excitation probabilities or capacity.

## 4. Solve once, then recognize the whole channel

Change to $`x=-t`$ and write $`g(x)=g(-x)`$ and $`\Delta(x)=\Delta(-x)`$. Define

```math
R(x)=\int_0^x[g(s)+\kappa]ds,
\qquad D(x)=\int_0^x\Delta(s)ds.
```

Solving the first scalar equation gives

```math
m=\int_0^T\sqrt{g(x)}\,e^{-R(x)/2-\gamma_\phi x-iD(x)}f(x)\,dx
 =\langle h_g,f\rangle,
```

where the inner product is conjugate-linear in its first argument and

```math
h_g(x)=\sqrt{g(x)}\,e^{-R(x)/2-\gamma_\phi x+iD(x)}.
```

Substituting $`a`$ into the population equation gives a time-ordered double integral. Combining its two conjugate triangles produces

```math
p=\langle f,K_gf\rangle,\qquad
K_g(x,y)=\sqrt{g(x)g(y)}e^{-[R(x)+R(y)]/2}
 e^{-\gamma_\phi|x-y|+i[D(x)-D(y)]}.
```

Dephasing acts differently in these two expressions: $`m`$ compares a photon amplitude with vacuum, while $`p`$ compares amplitudes with different arrival times. In particular, $`K_g=h_gh_g^\dagger`$ when $`\gamma_\phi=0`$; that rank-one simplification generally disappears with dephasing.

For one fixed mode, in the calibrated output basis,

```math
\mathcal E_f(\rho)=
\begin{pmatrix}
\rho_{00}+(1-p)\rho_{11}&m^*\rho_{01}\\
m\rho_{10}&p\rho_{11}
\end{pmatrix},\qquad 0\le|m|^2\le p\le1.
```

To cover all encodings, enlarge the input to $`\mathbb C|\mathrm{vac}\rangle\oplus L^2([0,T])`$. Its block $`\rho_{11}`$ is now an operator on temporal modes, not a scalar occupation. Polarization of the response gives

```math
\mathcal E_g(\rho)=
\begin{pmatrix}
\rho_{00}+\mathrm{Tr}[(I-K_g)\rho_{11}]&\rho_{01}|h_g\rangle\\
\langle h_g|\rho_{10}&\mathrm{Tr}(K_g\rho_{11})
\end{pmatrix}.
```

The inequalities $`0\preceq|h_g\rangle\langle h_g|\preceq K_g\preceq I`$ ensure complete positivity and trace preservation. Missing excitation returns to the ground state; it is never removed by conditional normalization. Exactly-one-photon inputs produce diagonal receiver states, a binary measure–prepare channel. Positive quantum communication in this model uses coherence with vacuum, with a shared phase reference.

## 5. Why two clocks yield a channel comparison

Impose $`0\le g(t)\le g_{\max}`$, and put $`\Gamma_M=g_{\max}+\kappa`$ and $`\beta=g_{\max}/\Gamma_M`$. Accumulated radiative time is $`u=R(x)`$; physical time is $`x`$. Their difference at maximum rate is

```math
d(u)=x(u)-u/\Gamma_M\ge0.
```

The rate ceiling makes $`d`$ nondecreasing. Consequently

```math
|x(u)-x(v)|=|u-v|/\Gamma_M+|d(u)-d(v)|.
```

This identity separates the reference receiver's dephasing from additional dephasing. It becomes an actual input channel through three operations:

1. Reshape the photon with the norm-preserving Jacobian $`f(x(u))/\sqrt{\Gamma(x(u))}`$, and remove the known phase $`D`$.
2. Attenuate its amplitude by $`c(u)=\sqrt{[g(x(u))/\Gamma(x(u))]/\beta}\le1`$, sending lost weight to vacuum.
3. Average the phases $`e^{i\sqrt{2\gamma_\phi}B_{d(u)}}`$, where $`B_s`$ is a standard Wiener process; leave vacuum unchanged.

A standard Wiener process starts at $`B_0=0`$ and has centered Gaussian increments with variance $`\mathrm{Var}(B_s-B_t)=|s-t|`$. The Gaussian phase average therefore multiplies photon-sector coherences by $`e^{-\gamma_\phi|d(u)-d(v)|}`$ and vacuum–photon coherences by $`e^{-\gamma_\phi d(u)}`$. Every realization is unitary, so their average is completely positive. This simultaneous match of both kernels is stronger than a pointwise comparison of excitation probabilities.

A generalized inverse handles sets with zero total radiative rate $`\Gamma=0`$: their photon weight becomes vacuum, while elapsed time remains in $`d`$. If $`g=0`$ but $`\kappa>0`$, the coordinate map retains the sector and the attenuation $`c=0`$ sends its photon weight to vacuum. The [full proof](channel-simulation.md) establishes

```math
\mathcal E_{g,\Delta}=\mathcal E_{\max}\circ\mathcal P_{g,\Delta}.
```

The preprocessor preserves the input window and never increases photon number. It may turn a pure pulse into a mixed state. Thus the theorem assumes unrestricted temporal source channels, including discarded source ancillas, rather than only deterministic waveform shaping. It does not assume free waveform conversion in a bandwidth-limited experiment.

## 6. From exact simulation to communication

Quantum capacity is the asymptotic number of qubits transmitted per independently reset use with vanishing error. It permits encoding across uses and entanglement with a reference. Since the identity above holds for the full channel, an encoder for any controlled receiver can be followed by $`\mathcal P_{g,\Delta}`$ and used with the maximum-rate receiver. The same decoder then receives the same joint state. Hence maximum coupling attains the supremum over allowed schedules, without any additivity assumption.

For a fixed mode, the established damping–dephasing criterion is $`Q>0`$ exactly when $`p+|m|^2>1`$. For all temporal encodings, write $`A=K_M`$ and $`b=h_M`$, the kernels above specialized to $`g=g_{\max}`$. The [capacity proof](capacity-boundary.md) constructs a symmetric two-output extension when $`A+|b\rangle\langle b|\preceq I`$. Both outputs reproduce the noisy receiver, which makes the channel antidegradable: its environment can be processed to simulate its output. Such a channel has zero unassisted capacity, because reliable quantum transmission would also give the environment a copy of the unknown state. If the inequality fails, a single mode has $`p+|m|^2>1`$ and witnesses positive capacity. Therefore

```math
Q(\mathcal E_{\max})>0
\quad\Longleftrightarrow\quad
\mu=\|A+|b\rangle\langle b|\|>1.
```

This determines zero versus positive capacity, not its positive value. The [upper-bound proof](capacity-upper-bound.md) instead attenuates the output to an antidegradable comparison channel and uses continuity to bound every code's rate. A small achievable-rate example alone supplies no such ceiling.

## 7. An analytical checkpoint

Take maximum coupling, zero detuning and zero extra dephasing, with $`u=\Gamma_M x`$ and $`S=\Gamma_M T`$. A physical-time envelope $`f_x`$ becomes $`f_u(u)=f_x(u/\Gamma_M)/\sqrt{\Gamma_M}`$, so its normalization is unchanged. Write $`f=f_u`$ for this calculation. The kernels in this coordinate are

```math
h_M(u)=\sqrt\beta\,e^{-u/2},\qquad
K_M(u,v)=\beta e^{-(u+v)/2},\qquad 0\le u,v\le S.
```

Normalize $`f(u)=Ne^{-u/2}`$. Calculate $`m`$ and $`p`$, then determine the window in which this mode has positive capacity. Only the integral of $`e^{-u}`$ is needed.

**Solution.** $`N=(1-e^{-S})^{-1/2}`$, so

```math
m=\sqrt{\beta(1-e^{-S})},\qquad
p=|m|^2=\beta(1-e^{-S}).
```

The fixed-mode channel is amplitude damping. Its positive-capacity condition is $`2p>1`$, giving

```math
\beta>\tfrac12,\qquad
S>\log\frac{2\beta}{2\beta-1}.
```

Equality still gives zero capacity. At $`\beta=1`$, the infimum window is $`S=\log2`$. Cauchy–Schwarz shows this pulse maximizes $`p`$, and the rank-one operator also gives $`\mu=2\beta(1-e^{-S})`$: the fixed-mode and all-mode boundaries coincide here. This is the elementary limit of the [finite-time spectral result](finite-time.md), whose dephasing case replaces the rank-one calculation with an integral eigenvalue problem.

Continue with the [model](model.md), [channel proof](channel-simulation.md), and [capacity boundary](capacity-boundary.md). The [counterexamples](counterexamples.md) show what changes with fixed source pulses or colored noise; the [verification guide](../verification/README.md) connects the formulas to independent computational checks.

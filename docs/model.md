# Receiver model and conventions

The model specifies the controlled dynamics, input resources, and capture instant. These conditions are part of the theorem, not implementation details.

## Physical and operational scope

The receiver is one two-level system, initially in $|g\rangle$. Its incident-port population-decay rate $g(t)$ is controllable within

```math
0\le g(t)\le g_{\max},\qquad 0<g_{\max}<\infty.
```

An inaccessible vacuum port has fixed population-decay rate $\kappa\ge0$. A fixed Markov dephasing process adds optical-coherence decay rate $\gamma_\phi\ge0$. The dissipators are $g(t)\mathcal D[\sigma_-]$, $\kappa\mathcal D[\sigma_-]$, and $2\gamma_\phi\mathcal D[|e\rangle\langle e|]$. A known detuning $\Delta(t)|e\rangle\langle e|$ and known incident-coupling phase are allowed. Known phase factors can be absorbed into the input temporal mode.

These are white-noise, Markov, rotating-wave dynamics for **every allowed schedule**, not merely a fit of the uncontrolled receiver's $T_1$ and $T_2$. The photon-sector equations are standard Fock-wavepacket dynamics [B12](references.md#b12). A reduced free-evolution semigroup alone does not generally specify controlled dynamics; the dilation caveat is important [GN17](references.md#gn17) and [BFH23](references.md#bfh23).

A use has input support $[-T,0]$, with $0<T<\infty$ and the output qubit evaluated at time 0. Inputs belong to

```math
\mathcal H_T=\mathbb C|\mathrm{vac}\rangle\oplus L^2([0,T],dx),
```

where $x=-t$ is lookback time. Codes may use arbitrary temporal modes and joint states across independently reset uses, with **at most one photon in each use**, including coherent superpositions with vacuum. Standard encoders may be channels, not only isometries into a single deterministic wavepacket. A shared phase reference is assumed for vacuum–one-photon coherence. Restricting every use to exactly one photon gives a different, zero-capacity input restriction, as shown below.

The comparison does not charge a source bandwidth, a temporal-mode conversion cost, a peak photon-flux limit, a coupling slew-rate limit, or a total integrated coupling budget. Maximum coupling throughout $T$ is an allowed schedule. The preprocessing in [the proof](channel-simulation.md) can compress a waveform and increase its bandwidth; it must not be called a free operation when a bandwidth-limited or fixed source is the real resource model. Temporal conversion itself is a physical operation with implementation conditions, as studied in [DMR15](references.md#dmr15); no laboratory implementation of the simulator is asserted here.

No transverse qubit drives, extra storage levels, additional emitters, multiphoton inputs per use, feedback, heralding, bath monitoring, or access to outgoing light is added. Fixed rates do not change when the receiver is controlled. Post-capture storage, output transfer, decoding and reset costs are not included; capacity is per independent capture-map use, not demonstrated throughput per second.

## Controlled kernels and the full channel

Write $g(x)=g(-x)$ and define

```math
\Gamma(x)=g(x)+\kappa,\qquad
R(x)=\int_0^x\Gamma(s)\,ds,\qquad
D(x)=\int_0^x\Delta(-s)\,ds.
```

The kernels are

```math
h_g(x)=\sqrt{g(x)}\,e^{-R(x)/2-\gamma_\phi x+iD(x)},
```

```math
K_g(x,y)=\sqrt{g(x)g(y)}\,e^{-[R(x)+R(y)]/2}
 e^{-\gamma_\phi|x-y|+i[D(x)-D(y)]}.
```

A normalized pulse f gives $p=\langle f,K_gf\rangle$ and $m=\langle h_g,f\rangle$. More importantly, for an arbitrary block density operator with photon block $\rho_{11}$,

```math
\mathcal E_g(\rho)=
\begin{pmatrix}
\rho_{00}+\operatorname{Tr}[(I-K_g)\rho_{11}]&\rho_{01}|h_g\rangle\\
\langle h_g|\rho_{10}&\operatorname{Tr}(K_g\rho_{11})
\end{pmatrix}.
```

The physical kernels obey $0\preceq|h_g\rangle\langle h_g|\preceq K_g\preceq I$. Thus this is a CPTP map. The entire map, rather than just the selected-mode damping–dephasing parameters, is used in the simulation proof.

### Why the vacuum sector matters

For a normalized input $\rho_{11}$ supported entirely in the one-photon temporal space, the restricted map is

```math
\mathcal E_g^{(1)}(\rho_{11})=
\operatorname{Tr}[(I-K_g)\rho_{11}]|g\rangle\langle g|
+\operatorname{Tr}(K_g\rho_{11})|e\rangle\langle e|.
```

This measures the binary POVM $\{I-K_g,K_g\}$ and prepares the corresponding energy eigenstate. It is entanglement breaking and has zero unassisted quantum capacity, even with arbitrary temporal-mode encodings entangled across independently reset uses that each contain exactly one photon. Temporal interference can change the excitation probability, but cannot create off-diagonal coherence in this two-level output.

The positive-capacity statements in this repository permit vacuum–photon coherence; their explicit witnesses encode in $\operatorname{span}\{|\mathrm{vac}\rangle,|1_f\rangle\}$. The fixed-one-photon-per-use restriction above does not describe a photon coherently shared across several uses, or global phase averaging of such a code. It is not a no-capacity claim for all protocols without a shared phase reference.

### Maximum-rate reference channel

For maximum coupling define

```math
\Gamma_M=g_{\max}+\kappa,\qquad
\beta_M=\frac{g_{\max}}{\Gamma_M},\qquad
r=\frac{\gamma_\phi}{\Gamma_M},\qquad S=\Gamma_M T.
```

In its dimensionless lookback coordinate u, the reference kernels are

```math
h_M(u)=\sqrt{\beta_M}e^{-u/2-r u},\qquad
K_M(u,v)=\beta_M e^{-(u+v)/2-r|u-v|},\quad 0\le u,v\le S.
```

The notation $\mathcal E_{\max}$ includes the number-preserving conversion between this coordinate and physical time $u/\Gamma_M$.

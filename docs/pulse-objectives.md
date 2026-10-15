# Absorption, coherence, and qubit fidelity

These closed-form pulse optimizations distinguish receiver excitation from coherent single-rail transfer. They support the channel analysis; the lead result is the [controlled-receiver simulation](channel-simulation.md). For this page, $\Gamma$ is a constant total radiative rate, $r=\gamma_\phi/\Gamma$ and $\beta=\Gamma_{\rm in}/\Gamma$.

## Exact channel from two pulse functionals

Define the positive integral operator K_r and the real function h_r by

```math
K_r(u,v)=e^{-(u+v)/2-r|u-v|},\qquad
h_r(u)=e^{-(r+1/2)u}.
```

After correcting the known overall optical phase, the receiver channel is

```math
\mathcal E_f(\rho)=
\begin{pmatrix}
\rho_{00}+(1-p)\rho_{11}&m^*\rho_{01}\\
m\rho_{10}&p\rho_{11}
\end{pmatrix},
\quad
p=\beta\langle f,K_rf\rangle,\quad
m=\sqrt\beta\langle h_r,f\rangle.
```

Here p is excitation conditional on sending one photon. The off-diagonal multiplier m determines coherence with the vacuum part of the logical input. No averaging over a photon-number mixture is substituted for the encoded superposition.

For real f the exact scalar Fock-state equations are, in forward dimensionless time s=Gamma t<=0,

```math
\dot c=-(1/2+r)c+\sqrt\beta f(-s),\qquad
\dot p=-p+2\sqrt\beta f(-s)c.
```

At s=0, c=m. These equations follow from the one-photon sector of the established wavepacket master-equation framework [B12](references.md#b12). The checker also integrates all four 2x2 generalized Fock density-matrix blocks independently; it does not only compare two formulas built from the same scalar equations.

The kernel has trace one. Its square is integrable, and direct integration gives Tr(K_r^2)=1/(1+2r). Under this particular Markov model, K_r is also the time-reversed temporal density kernel of a spontaneously emitted photon. This is a model-specific relation, not an assertion that time-reversing a reduced photon reverses an unobserved bath.

Complete positivity is visible from the exact Gram decomposition

```math
K_r(u,v)-h_r(u)h_r(v)
=2r\,e^{-(u+v)/2}
\int_0^{\min(u,v)}e^{-r(u+v-2s)}ds.
```

Thus p>=|m|^2, and p<=beta. The r=0 kernel is rank one and the two objectives coincide.

For a uniform unknown input qubit, the entanglement and average fidelities are

```math
F_e=\frac{1+p+2|m|}{4},\qquad
\overline F=\frac{3+p+2|m|}{6}.
```

These use the standard relation between average and entanglement fidelity [N02](references.md#n02); both are independently checked from the Choi matrix and the six Pauli eigenstates, a qubit state 2-design. A fixed known phase makes m nonnegative. Since the kernel and h are pointwise positive, replacing f by |f| cannot decrease either relevant term; the positive-weight optima can be sought among real nonnegative pulses without restricting the attainable optimum.

The two-thirds classical benchmark is for deterministic measure–prepare transfer of one uniformly unknown qubit [MP95](references.md#mp95), not for heralded or differently weighted ensembles. It can also be checked directly: a qubit POVM element E_y=a_y(I+v_y.sigma), followed by preparation with Bloch vector s_y, contributes a_y[1/2+(v_y.s_y)/6] to the uniform average. Summing and using sum a_y=1 and |v_y.s_y|<=1 gives 2/3. A projective measurement followed by preparation attains it.

A fidelity below 2/3 does **not** imply that this channel is entanglement-breaking or that its quantum capacity vanishes. Its Choi state is non-positive under partial transpose whenever m is nonzero. The separate [capacity analysis](capacity-boundary.md) supplies the zero/positive criterion; a fidelity threshold alone does not.

## Exact population-optimal pulse

For r>0 set nu=2r, and let j_{nu-1,n} be the nth positive zero of J_{nu-1}. The complete nonzero spectrum of K_r is

```math
\mu_n=\frac{8r}{j_{\nu-1,n}^2}.
```

The first eigenvalue is simple and its eigenfunction is positive. The maximal population and normalized optimizing pulse are

```math
p_{\max}=\beta\mu_1,\qquad
f_{\rm pop}(u)=
\frac{e^{-u/2}J_\nu(j_{\nu-1,1}e^{-u/2})}
{J_\nu(j_{\nu-1,1})}.
```

For -1<nu-1<0 the first positive root is the one tending to zero as r decreases; using integer-order root routines or discarding this root gives an incorrect small-noise limit. At r=0 the answer is obtained continuously as f=e^{-u/2}, p=beta, m=sqrt(beta).

### Derivation and boundary conditions

For K_r f=mu f, let

```math
g(u)=\int_0^\infty e^{-r|u-v|}e^{-v/2}f(v)dv.
```

The Green-function jump of the exponential yields

```math
g''-r^2g=-\frac{2r}{\mu}e^{-u}g,
\qquad g'(0)=r g(0).
```

With x=sqrt(8r/mu)e^{-u/2}, this is Bessel's equation of order nu=2r. The decaying solution at u=infinity selects J_nu and excludes Y_nu. At u=0 the Robin condition becomes

```math
xJ_\nu'(x)+\nu J_\nu(x)=xJ_{\nu-1}(x)=0.
```

The normalization identity is

```math
2\int_0^1 yJ_\nu(xy)^2dy
=J_\nu(x)^2-J_{\nu-1}(x)J_{\nu+1}(x).
```

At x=j_{nu-1,1}, the denominator therefore reduces to J_nu(x). The normalized pulse has f_pop(0)=1. Evaluating the eigenvalue equation at u=0 then gives

```math
\langle h_r,f_{\rm pop}\rangle=\mu_1,
\qquad
m_{\rm pop}=\sqrt\beta\,\mu_1.
```

For beta=1 this is the particularly simple relation **m_pop=p_max**. Consequently its mean qubit fidelity is 1/2+mu_1/2. This equality is a property of the excitation-optimal pulse, not a universal channel identity.

## Exact optimum for any positive population/coherence weighting

Maximize

```math
Q_w(f)=p+2wm,\qquad w>0,\quad\|f\|=1.
```

The case w=1 is the qubit-fidelity objective. Let A=beta K_r and b=w sqrt(beta)h_r. There is a unique z>beta mu_1 for which

```math
f_*=(zI-A)^{-1}b,\qquad \|f_*\|=1.
```

Existence and uniqueness follow because the squared norm is continuous and strictly decreasing in z, diverges at the largest eigenvalue (b has nonzero overlap with its positive eigenfunction), and tends to zero at infinity.

This solution is globally optimal. For every other unit pulse f,

```math
\boxed{
Q_w(f_*)-Q_w(f)
=\langle f-f_*,(zI-\beta K_r)(f-f_*)\rangle\ge0.
}
```

The statement extends to complex pulses after the known overall phase is aligned. It is not merely a stationarity equation or a finite pulse-parameter search.

### Closed Bessel form

Define, for 0<x<j_{nu-1,1},

```math
\Delta_\nu(x)=J_\nu(x)^2-J_{\nu-1}(x)J_{\nu+1}(x),
\qquad
f_x(u)=\frac{e^{-u/2}J_\nu(xe^{-u/2})}{\sqrt{\Delta_\nu(x)}}.
```

The unique x associated with weight w solves

```math
w^2x^2\Delta_\nu(x)=4\beta J_{\nu-1}(x)^2.
```

For compact expressions set

```math
d(x)=\frac{2J_{\nu+1}(x)}{x\sqrt{\Delta_\nu(x)}},
\qquad c(x)=\frac{2J_{\nu-1}(x)}{x\sqrt{\Delta_\nu(x)}}.
```

Then

```math
m(x)=\sqrt\beta\,d(x),\qquad
p(x)=\beta\left[\frac{8r}{x^2}-c(x)d(x)\right].
```

At the optimum c=w/sqrt(beta), z=8r beta/x^2. This gives the exact supported upper tradeoff for every positive weight, without asserting that every suboptimal (p,m) pair has been classified. Taking x to its first-root endpoint gives the population optimum. Taking x to zero gives the coherence optimum

```math
f_{\rm coh}(u)=\sqrt{1+2r}\,e^{-(r+1/2)u},
\quad m_{\max}=\sqrt{\frac\beta{1+2r}},
\quad p_{\rm coh}=\frac\beta{1+r}.
```

## The operational difference

For beta=1 and r=3:

| Optimized objective | Population p | Vacuum–one coherence amplitude m | Average qubit fidelity |
|---|---:|---:|---:|
| Excitation | 0.311935720 | 0.311935720 | 0.655967860 |
| Unknown-qubit fidelity | 0.281616606 | 0.371214185 | 0.670674163 |

The fidelity-optimal pulse sacrifices population but improves the off-diagonal transfer enough to cross the deterministic 2/3 benchmark. These are unconditional channel values. No photon-detection event or successful-capture subensemble is selected.

The population-optimized fidelity crosses 2/3 at r=2.6836838712. The fidelity-optimized result crosses at r=3.1435035139. These are numerical evaluations of the exact model formulas, not experimentally measured thresholds.

A useful explanation is that population depends on interference between two possible absorption times, whose relative dephasing is set by their separation. Coherence with the vacuum branch instead compares an excitation history with the branch in which the receiver was never excited. Long histories are penalized differently. The unknown-qubit objective therefore favors a different temporal concentration. No assumption that the bath literally performs a projective measurement of absorption time is needed.

### Large dephasing exposes a parametric distinction

For fixed beta>0 and r tending to infinity, the Bessel-root asymptotic gives mu_1~2/r. Hence

```math
\overline F_{\rm pop}-\frac12
\sim\frac{\beta+2\sqrt\beta}{3r}.
```

In contrast,

```math
\overline F_{\rm opt}-\frac12
\sim\frac{\sqrt\beta}{3\sqrt{2r}}.
```

The second relation does not require a delicate expansion of the optimal root. The coherence-only exponential supplies a lower bound with that leading term, while p<=beta mu_1 and m<=sqrt(beta/(1+2r)) supply a matching upper bound. These are improvements relative to the **reset-channel baseline 1/2**. Both fidelities eventually fall below 2/3; this is not a diverging practical quantum advantage in the strong-noise limit. The white-noise Markov approximation must also remain valid over the increasing pulse bandwidth for an actual device interpretation.

## Strong comparator: optimize an ordinary exponential for the correct task

A normalized rising exponential has f_a(u)=sqrt(a)e^{-a u/2}, a>0. Exact integration yields

```math
p(a)=\frac{4\beta a}{(1+a)(1+a+2r)},\qquad
m(a)=\frac{2\sqrt{\beta a}}{1+a+2r}.
```

Optimizing population within this family gives $a=\sqrt{1+2r}$. For $\beta=1$, this population formula and its exponential-family optimum are established in [L22](references.md#l22), final supplementary Sec. S8, Eq. (S22). Their pulse $\xi(t)=\sqrt{2/\tau}e^{t/\tau}\Theta(-t)$ gives the present formula under $a=2/(\Gamma\tau)$ and $\gamma=\Gamma(1/2+r)$; their $\tau_{\rm opt}=\sqrt{2/(\gamma\Gamma)}$ becomes the same $a$. This identification uses their pure-dephasing limit: their combined linewidth parameter can also include nonradiative relaxation, which must be distinguished from pure dephasing as explained in [Physical interpretation](physical-scope.md#rates-must-describe-the-same-dynamics).

Optimizing fidelity is a separate one-variable problem. At r=3, beta=1 its answer is a=4.86350965, with

```math
p=0.27966553,\quad m=0.37178486,\quad
\overline F=0.670539208.
```

The exact unrestricted optimum is only 0.000134955 higher. At r=1 the gap is 0.000335130; at r=10 it is 0.0000082673. These are pointwise checks, not a proof of a uniform bound over every r and beta.

This comparison prevents a misleading “Bessel pulses beat conventional pulses” claim. The clear effect is the change from optimizing population to optimizing coherent transfer. Complex pulse synthesis supplies only a small further improvement in the tested examples. Even the strong-noise leading fidelity scaling is already achieved by the coherence-optimized exponential.

This page optimizes a fixed receiver. The [simulation theorem](channel-simulation.md) handles bounded variable coupling and known detuning with free temporal encoding. Extra storage levels, bath monitoring, heralded recovery, multiple excited-state encodings and non-Markovian reservoir engineering are different resource models; no universal quantum-memory bound is claimed. The hypothetical ideal beta=1 removes radiative port mismatch in the example; it does not claim such a port is available for every emitter.

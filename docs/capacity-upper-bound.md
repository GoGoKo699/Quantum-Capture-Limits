# A capacity upper bound from the boundary margin

The [exact sign criterion](capacity-boundary.md) also yields a quantitative upper bound on the full, regularized unassisted capacity. It applies to every temporal encoding and every block length under the [model's resources](model.md). It complements the coherent-information lower bounds without asserting their optimality.

For a vacuum-preserving capture map with photon-sector effect $A$ and vacuum-coherence vector $b$, set

```math
\mu=\|A+|b\rangle\langle b|\|,\qquad
\epsilon=\begin{cases}0,&\mu\le1,\\1-\mu^{-1},&\mu>1.\end{cases}
```

Write $g_2(x)=(1+x)\log_2(1+x)-x\log_2x$, with $g_2(0)=0$. Then

```math
\boxed{Q(\mathcal E)\le U(\mu)
\equiv\min\{1,\;2\epsilon+g_2(\epsilon)\}
\quad\text{qubits per independent use}.}
```

At maximum coupling, $\mu=\beta\Lambda_r(S)$; the infinite-window value is $\beta\Lambda_r(\infty)$. By [channel simulation](channel-simulation.md), the same maximum-coupling upper bound applies to every allowed receiver schedule. The entropy-continuity ingredient is established: [SH17, Proposition 3A and Proposition 11, Eq. (87)](references.md#sh17) refines the channel-capacity method of [LS09](references.md#ls09). The receiver-specific step below constructs the nearby antidegradable channel directly from $\mu$.

## An explicit antidegradable comparison channel

If $\mu\le1$, the original map already has the broadcast extension in [Capacity boundary](capacity-boundary.md). Suppose $\mu>1$, and let $\mathcal D_\tau$ be qubit amplitude damping with excitation survival $\tau=\mu^{-1}$. Its Kraus operators are

```math
D_0=|0\rangle\langle0|+\sqrt\tau|1\rangle\langle1|,
\qquad D_1=\sqrt{1-\tau}|0\rangle\langle1|.
```

The comparison channel $\mathcal F=\mathcal D_\tau\circ\mathcal E$ has

```math
A_{\mathcal F}=\tau A,\qquad b_{\mathcal F}=\sqrt\tau\,b,
\qquad A_{\mathcal F}+|b_{\mathcal F}\rangle\langle b_{\mathcal F}|
=\tau(A+|b\rangle\langle b|)\preceq I.
```

Thus $\mathcal F$ is antidegradable. This auxiliary output attenuation is used to bound the original channel's capacity. It is separate from the input preprocessing that compares different receiver schedules.

## Distance bound, including arbitrary references

Qubit amplitude damping satisfies

```math
\frac12\|\mathrm{id}-\mathcal D_\tau\|_\diamond=1-\tau.
```

Here is a direct proof. By convexity it suffices to use a pure qubit–reference state, written

```math
|\psi\rangle=\sqrt{1-q}|0\rangle|a\rangle
+\sqrt q|1\rangle|b\rangle,
```

where $a,b$ are unit vectors, need not be orthogonal, and $0\le q\le1$. Let $v=(D_0\otimes I)\psi$, $w=(D_1\otimes I)\psi$, $\epsilon=1-\tau$, and $c=1-\sqrt\tau$. The trace norm of the difference between the two rank-one operators $|\psi\rangle\langle\psi|$ and $|v\rangle\langle v|$, followed by the triangle inequality for $|w\rangle\langle w|$, gives

```math
\big\||\psi\rangle\langle\psi|-(\mathcal D_\tau\otimes\mathrm{id})(|\psi\rangle\langle\psi|)\big\|_1
\le\epsilon q+\sqrt{\epsilon^2q^2+4c^2q(1-q)}\le2\epsilon.
```

The last inequality follows from $c\le\epsilon$ and

```math
\epsilon^2(2-q)^2-[\epsilon^2q^2+4c^2q(1-q)]
=4(1-q)(\epsilon^2-c^2q)\ge0.
```

An excited-state input attains $2\epsilon$. Contractivity under composition therefore gives

```math
\frac12\|\mathcal E-\mathcal F\|_\diamond
\le\frac12\|\mathrm{id}-\mathcal D_\tau\|_\diamond
=\epsilon.
```

This bound has no temporal-input dimension factor.

## From channel distance to arbitrary quantum codes

Use the finite-output mutual-information estimate [SH17, Proposition 3A](references.md#sh17). For any state on a finite reference $R$ and $n$ temporal inputs, the two output states obey

```math
|I(R:B^n)_{\mathcal E^{\otimes n}}
-I(R:B^n)_{\mathcal F^{\otimes n}}|
\le n[2\epsilon\log_2d_B+g_2(\epsilon)],\qquad d_B=2.
```

The proof replaces one channel use at a time; the other-output marginal stays fixed. Its constant depends only on the finite output dimension, so it covers the separable infinite-dimensional temporal input space.

Both channels leave $\rho_R$ unchanged, so the same difference bound holds for $I_c(R\rangle B^n)=I(R:B^n)-H(R)$. The broadcast extension of $\mathcal F^{\otimes n}$ and weak monotonicity of conditional entropy give

```math
I_c(R\rangle B^n)_{\mathcal F^{\otimes n}}\le0,
\qquad
I_c(R\rangle B^n)_{\mathcal E^{\otimes n}}
\le n[2\epsilon+g_2(\epsilon)].
```

This also gives an operational converse for general encoders, including encoders that discard an ancillary system. If a code transmits half of a maximally entangled state of dimension $K_n$, with decoded trace distance $\eta_n=\tfrac12\|\rho_{R\widehat R}-\Phi_{K_n}\|_1$ from the target, data processing and conditional-entropy continuity imply

```math
(1-2\eta_n)\log_2K_n-g_2(\eta_n)
\le I_c(R\rangle B^n)
\le n[2\epsilon+g_2(\epsilon)].
```

For vanishing decoding error, divide by $n$ and take the limit. The output-dimension bound supplies $Q\le1$. This proof uses finite reference and output systems for each code, without assuming an infinite-input regularized capacity formula.

## What the upper bound adds near the threshold

With positive spectral margin $\delta=\mu-1$, one has $\epsilon=\delta/(1+\delta)$ and

```math
U(1+\delta)=\delta\log_2(1/\delta)+O(\delta),
\qquad \delta\downarrow0.
```

Thus the full capacity tends to zero as the sign boundary is approached from the positive side, even after optimizing temporal encoding and allowing entangled multi-use codes. A tiny coherent-information witness alone would not establish this upper limit.

At $\beta=1$ with an unlimited input window, evaluations of this upper bound are:

| Extra-dephasing ratio $r$ | Capacity upper bound $U(\Lambda_r)$, qubits/use |
|---:|---:|
| 0.72 | 0.1015959055 |
| 0.73 | 0.0426072045 |
| 0.735 | 0.0073368283 |
| 0.7357 | 0.0011154523 |
| 0.7358 | 0.0000123617 |

These are floating-point evaluations of the analytic bound, rather than interval-certified numerical enclosures. They upper-bound the entire temporal-input channel, whereas the [positive-rate witnesses](capacity-boundary.md#capacity-positivity-is-not-a-useful-rate-guarantee) use particular pulses and inputs.

The bound is generally loose. It neither matches the small-input lower bounds nor determines the true critical scaling of the capacity. Combined with the [logarithmic capture-window divergence](finite-time.md#the-window-diverges-logarithmically-near-the-boundary), it distinguishes two statements: longer pulses can recover a strictly positive channel margin near the boundary, while every achievable rate still vanishes as that margin vanishes.

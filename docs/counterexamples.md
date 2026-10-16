# Why source freedom and Markov noise matter

These examples delimit the [simulation theorem](channel-simulation.md). A prescribed source and a different noise dilation are different operational tasks.

## A fixed-source counterexample makes the resource boundary concrete

Take a noiseless single-port receiver, $\kappa=\gamma_\phi=0$, with $g_{\max}=1$, capture duration $T=4$, and a **fixed square** incoming pulse $\xi(t)=1/\sqrt T$ for $0\le t\le T$. The excitation amplitude satisfies

```math
\dot a(t)=-\frac{g(t)}2a(t)+\sqrt{g(t)}\,\xi(t),\qquad a(0)=0.
```

At constant maximum coupling,

```math
p_{\rm square,max}=\frac4T(1-e^{-T/2})^2=0.7476450724\ldots .
```

Choose the bounded schedule

```math
g(t)=\begin{cases}
1,&0\le t\le t_s,\\
[1+t-t_s]^{-1},&t_s<t\le T,
\end{cases}
\qquad t_s=2\log2.
```

It yields

```math
p_{\rm square,tuned}=1-\frac{2\log2-1}{T}=0.9034264097\ldots .
```

Thus tuning raises capture from 74.76% to 90.34% for this prescribed photon. The formula is checked against a separate amplitude ODE. This is an illustrative instance of familiar waveform matching, not a new receiver protocol [NJY16](references.md#njy16). No optimality claim for every fixed-waveform problem is needed.

With source shape free, a maximum-rate receiver instead accepts the normalized rising exponential and attains

```math
p_{\rm free-source,max}=1-e^{-T}=0.9816843611\ldots .
```

The three numbers answer different optimization problems. They are model calculations, not measured efficiencies. They demonstrate why the theorem must not be advertised as “tunable coupling never helps.” It says that tunability does not add communication power **over a source-optimized maximum-rate receiver**, given the stated resources.

## Markov behavior is a substantive assumption, not a fitted T2 label

The simulation relied on exponential phase correlations splitting into an extra independent phase process. It does not automatically extend to colored detuning noise.

For a simple check, take reference lookback labels $t=(0,1,2)$, with label zero denoting vacuum, and actual labels $x=(0,1,3)$. They are compatible with a maximum rate $1$, followed farther back by a rate-$1/2$ interval. The extra-time labels are $d=(0,0,1)$.

Under Markov phase diffusion the residual correlation matrix is $e^{-\gamma_\phi|d_i-d_j|}$, which is positive semidefinite. Under a quasistatic Gaussian detuning with variance 1 and coherence $C(t)=e^{-t^2/2}$, the analogous proposed ratio would be

```math
R=\begin{pmatrix}
1&1&e^{-5/2}\\
1&1&e^{-3/2}\\
e^{-5/2}&e^{-3/2}&1
\end{pmatrix}.
```

Its determinant is

```math
\det R=-(e^{-5/2}-e^{-3/2})^2<0,
```

and its minimum eigenvalue is approximately -0.0100788. Therefore that ratio cannot define a phase-noise channel. This is a failure of this particular simulator construction in a changed noise model; it is not a proof that no other simulator exists, and not a computed capacity advantage for modulation under colored noise.

Related dynamical-decoupling literature already warns that a specified quantum stochastic Markov dilation and a merely fitted semigroup are not interchangeable [GN17](references.md#gn17) and [BFH23](references.md#bfh23). Those results concern noise suppression by system rotations and do not directly establish this source/receiver channel factorization. Conversely, this factorization is not a general theorem that all physical dynamics fitted by a Lindblad equation remain uncontrollable under every available intervention.

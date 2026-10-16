# Finite capture time and achievable-rate witnesses

Use $`\Gamma_M=g_{\max}+\kappa`$, $`\beta=g_{\max}/\Gamma_M`$, $`r=\gamma_\phi/\Gamma_M`$ and $`S=\Gamma_M T`$. The [simulation theorem](channel-simulation.md) reduces the all-control question to the maximum-rate receiver. Capacity is positive precisely when $`\beta\Lambda_r(S)>1`$; equality belongs to the zero-capacity side.

The notation $`\beta_M`$ in the threshold formulas denotes this same maximum-rate branching ratio $`\beta`$.

## Exact finite-window spectrum

Let $`L_rf=\Lambda f`$ on $`[0,S]`$, and define $`g`$ by factoring out the outer $`e^{-u/2}`$. Direct differentiation of the Green kernel gives

```math
g''-r^2g=-\frac{2r}{\Lambda}e^{-u}g,
\qquad g'(0)=0,\qquad g'(S)=-r g(S).
```

For $`r>0`$, set $`\nu=2r`$, $`k=\sqrt{8r/\Lambda}`$, and $`x=k e^{-u/2}`$. The equation becomes Bessel's equation of order $`\nu`$. Its solution is

```math
g(u)=A J_\nu(x)+B Y_\nu(x).
```

Write $`y=k e^{-S/2}`$. The endpoint conditions imply

```math
\boxed{J_\nu'(k)Y_{\nu+1}(y)-Y_\nu'(k)J_{\nu+1}(y)=0.}
```

The smallest appropriate positive root yields the largest eigenvalue $`\Lambda_r(S)=8r/k^2`$. Equivalently, take

```math
f(u)=\mathcal N e^{-u/2}
\left[J_\nu(k e^{-u/2})+cY_\nu(k e^{-u/2})\right],
\quad c=-\frac{J_{\nu+1}(y)}{Y_{\nu+1}(y)},
```

and choose the positive, normalized leading solution. This is a continuum optimum. The numerical code also checks the integral equation directly and compares nested Galerkin lower bounds obtained from exact cell integrals.

At $`r=0`$ the kernel is rank one and no Bessel root is needed:

```math
\Lambda_0(S)=2(1-e^{-S}),\qquad
f(u)=\frac{e^{-u/2}}{\sqrt{1-e^{-S}}}.
```

The zero-noise finite-coupling restriction is elementary and consistent with prior perfect-capture constructions. Nurdin, James and Yamamoto derive ideal arbitrary-wavepacket capture with a tunable coupling that is singular at the start, then examine its truncation [NJY16](references.md#njy16). Those results concern receiving a prescribed waveform in the ideal model; they are not contradicted by a maximum-rate noisy, jointly optimized comparison.

## Minimum window for positive unassisted capacity

The leading eigenvalue increases strictly with $`S`$, from zero to $`\Lambda_r(\infty)`$. Positive kernels and their strictly positive maximizing eigenfunctions establish the strictness. If $`\beta_M\Lambda_r(\infty)\le1`$, no window succeeds. Otherwise define $`S_Q`$ by

```math
\beta_M\Lambda_r(S_Q)=1.
```

The equality point still has $`Q=0`$. Every strictly longer window admits $`Q>0`$; every shorter window, for every allowed control, has $`Q=0`$. This is an infimum, not an attained minimum with positive rate.

At zero dephasing,

```math
\boxed{S_Q(0,\beta_M)=\log\frac{2\beta_M}{2\beta_M-1},\qquad\beta_M>1/2.}
```

For $`r>0`$, the threshold determinant can be evaluated at $`k=\sqrt{8\beta_M r}`$, leaving a scalar root for $`S`$. At perfect incident coupling:

| $`r`$ = extra coherence decay / maximum population decay | Infimum $`S_Q=\Gamma_M T`$ |
|---|---:|
| 0 | 0.69314718 |
| 0.1 | 0.73906771 |
| 0.3 | 0.86403415 |
| 0.5 | 1.08345528 |
| 0.6 | 1.29351974 |
| 0.7 | 1.82985541 |
| 0.72 | 2.16784210 |
| 0.73 | 2.58372062 |
| 0.735 | 3.40157853 |
| 0.7357 | 4.24968578 |
| $`0.735800796744\ldots`$ or larger | No successful window |

For $`\beta=0.9`$ and $`r=0.5`$ the infimum is $`1.61272277`$, versus $`1.08345528`$ for $`\beta=1`$. For $`\beta\le1/2`$ there is no positive-capacity regime, even without dephasing.

At maximum coupling define $`T_1=1/\Gamma_M`$ and $`T_2=1/(\Gamma_M/2+\gamma_\phi)`$ in this Markov model. For $`\beta=1`$, the unlimited-window boundary is equivalently

```math
\frac{T_2}{T_1}>\frac1{1/2+r_c}=0.8091919043\ldots.
```

This is a translated model criterion, not a universal hardware threshold inferred from arbitrary experimental $`T_1/T_2`$ fits or noise with memory.

### The window diverges logarithmically near the boundary

For fixed $`r>0`$, let $`j=j'_{2r,1}`$, $`\nu=2r`$,

```math
D_\nu=J_\nu(j)^2-J_{\nu-1}(j)J_{\nu+1}(j).
```

The small-$`y`$ expansion of the finite-endpoint coefficient and the simple derivative root give

```math
\Lambda_r(\infty)-\Lambda_r(S)
\sim A(r)e^{-(1+2r)S},
```

```math
A(r)=\frac{\Lambda_r(\infty)(j/2)^{4r}}
{(1+2r)\Gamma(2r+1)^2D_{2r}}>0.
```

For example, $`-J_{\nu+1}(y)/Y_{\nu+1}(y)\sim\pi y^{2\nu+2}/[2^{2\nu+2}\Gamma(\nu+1)\Gamma(\nu+2)]`$. Expanding the Neumann condition around $`j`$ and using the Bessel Wronskian produces the stated coefficient. The $`r=0`$ limit is exactly $`2e^{-S}`$.

With $`\beta=1`$ and $`r`$ approaching $`r_c`$ from below,

```math
\boxed{S_Q=\frac{1}{1+2r_c}\log\frac1{r_c-r}+O(1).}
```

The logarithmic coefficient is 0.4045959521. This is not a closing dynamical Liouvillian gap or a dissipative phase transition: it is the duration needed for a truncated input to retain enough of the optimal temporal mode to cross the channel criterion. The radiative/dephasing rates themselves remain finite. Near this boundary the achievable rate may still be extremely small, as shown in [Capacity boundary](capacity-boundary.md).

The [capacity upper bound](capacity-upper-bound.md) also applies at each finite $`S`$. For fixed $`r`$ and $`\beta`$ admitting a threshold $`S_Q`$, smoothness of the leading eigenvalue gives $`Q=O((S-S_Q)\log(1/(S-S_Q)))`$ as $`S\downarrow S_Q`$ from above. This bounds every unassisted code; it does not assert an exact rate law.

## A finite-pulse achievable rate

To separate the exact sign boundary from infinitesimal-rate examples, choose $`\beta=1`$, $`r=1/2`$ and the ordinary rising exponential, written in lookback coordinates,

```math
f(u)=\sqrt{\frac{a}{1-e^{-aS}}}\,e^{-au/2},
\qquad 0\le u\le S,
\qquad a=\sqrt3.
```

This is a compactly supported input, not an infinite precursor tail. With $`c=(a+1)/2`$,

```math
p=\frac{2a}{(1-e^{-aS})(c-r)}
\left[\frac{1-e^{-(c+r)S}}{c+r}
-\frac{1-e^{-2cS}}{2c}\right],
```

```math
m=\sqrt{\frac{a}{1-e^{-aS}}}\frac{1-e^{-(c+r)S}}{c+r}.
```

At $`S=4`$,

```math
p=0.6789544878092364\ldots,\qquad
m^2=0.4973397041579534\ldots.
```

For $`\rho_q=\mathrm{diag}(1-q,q)`$, $`q=1/8`$, evaluate the exact complementary-channel entropy. The environment eigenvalues are

```math
e_0=(1-p)q,\qquad e_\pm=\frac{t\pm\sqrt{t^2-4(p-m^2)q(1-q)}}2,
\quad t=1-(1-p)q.
```

Therefore

```math
\boxed{I_c=H_2(pq)-H(e_0,e_+,e_-)
=0.02957015518955046\ldots\ \text{qubits per use}.}
```

Positive coherent information gives an achievable unassisted rate; this coding fact and the damping–dephasing application are established [SAJS24](references.md#sajs24). No selection of successful absorptions is made. The average one-photon probability for this input is $`1/8`$, not exponentially small. Explicit values for the same pulse family and input are:

| Dimensionless window $`S`$ | Coherent-information lower bound per use |
|---:|---:|
| 2 | 0.01577406 |
| 3 | 0.02784609 |
| 4 | 0.02957016 |
| 6 | 0.02981631 |

These are achievable rates for the stated pulse and input. The [capacity upper bound](capacity-upper-bound.md) constrains the regularized rate over all codes; damping–dephasing coherent information can be nonadditive [SAJS24](references.md#sajs24) and [BL25](references.md#bl25). Rates are per independently reset capture use, with the operational assumptions specified in the [model](model.md).

Any strict positive margin and positive coherent-information witness persist under a sufficiently small $`L^2`$ smoothing of the waveform's support endpoint.

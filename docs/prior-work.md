# Relation to established capture and channel results

The contribution is a comparison of imperfect capture channels under specified resources. For a two-level receiver with bounded incident coupling, fixed vacuum loss, and fixed Markov dephasing, the [construction](channel-simulation.md) gives

```math
\mathcal E_{g,\Delta}=\mathcal E_{\max}\circ\mathcal P_{g,\Delta}.
```

The normal CPTP preprocessor acts on vacuum plus the entire one-photon temporal space, preserves the input window, and never increases photon number. It depends on the known schedule, not the transmitted state. With unrestricted temporal encoding, every controlled-receiver code can therefore be reproduced by the maximum-rate receiver. The [model assumptions](model.md), including fixed controlled Markov dynamics and uncharged source bandwidth, are part of this claim.

## Capture, time reversal, and source shaping

| Antecedent and source location | Established result and relation to this work |
|---|---|
| [G07L](references.md#g07l), arXiv:quant-ph/0604037v3, Eqs. (6)–(10) and the time-reversal argument following Eq. (8) | Optimized storage in a Lambda ensemble follows from retrieval; the optimum can be attained by adjusting the receiver control or the photon waveform. A spectral efficiency problem and accumulated-control variable already appear. Source/control interchange and time reversal are therefore established ideas. The result is not the present capped, noisy two-level channel factorization. |
| [G07C](references.md#g07c), arXiv:quant-ph/0612082v2, Sec. V.A Eq. (15), Sec. V.B Eqs. (23)–(26), Appendix B | In the adiabatic cavity model, rescaling amplitudes and changing time to integrated control intensity removes control dependence. Optimal storage is a mode-overlap problem. Equation (22) also treats spin-wave decay by modifying the desired mode and efficiency. This is a direct antecedent for changing the clock; that substitution alone is not a contribution here. |
| [WMHS12](references.md#wmhs12), arXiv:1111.4261, Sec. II.C Eqs. (19)–(24) | A moving-mirror two-level receiver has a bounded radiative rate. Its general finite-window capture amplitude and optimized efficiency $`1-e^{-\int\gamma_z(t)dt}`$ are derived after specializing inaccessible loss to zero. The noiseless maximum-rate absorption conclusion follows immediately. Fixed extra dephasing and an all-input noisy preprocessor are not established by these formulas. |
| [YJ14](references.md#yj14), arXiv:1403.1698v2, Sec. 4.2 Eqs. (18)–(20), Sec. 5.1 Eqs. (23)–(25) | Zero output determines rising-exponential vector modes in passive linear networks. Equation (25) transfers unknown coherent temporal-mode superpositions into network excitations. This is quantum-state transfer, not merely a population calculation. The comparison here instead concerns all inputs of imperfect noisy channels and all allowed schedules. |
| [NJY16](references.md#njy16), arXiv:1609.05643v1, Secs. 3.1–3.3 | Prescribed-wavepacket capture uses a coupling amplitude proportional to $`\xi(t)/\sqrt{\int_0^t\lvert\xi(s)\rvert^2ds}`$. The initial singularity and its truncation are explicitly considered. This tailors the receiver to a given source; it does not order the full noisy channel family under free source encoding. Their coupling amplitude is the square root of a population-decay rate. |

For constant coupling without extra dephasing, [RG17](references.md#rg17), Sec. II.A.1 Eqs. (1)–(8), already gives the loss branching ratio and matched rising exponential. Its 2019 erratum concerns large-photon-number Eqs. (38)–(39), not these single-photon formulas. The recent [LSYMM25](references.md#lsymm25), v1 Sec. II.2 Eqs. (9)–(12), explicitly reviews integrated coupling and ideal single-excitation matching before extending to multiphoton processes.

## Noisy pulse objectives and quantum information

[L22](references.md#l22), final supplementary information Sec. S8 Eqs. (S18), (S22), already gives the noisy rising-exponential population optimum. With zero inaccessible loss, its parameters satisfy $`\gamma=\Gamma/2+\gamma_\phi`$. Substituting $`a=2/(\Gamma\tau)`$ yields exactly

```math
p(a)=\frac{4a}{(1+a)(1+a+2r)},\qquad a_{\rm opt}=\sqrt{1+2r}.
```

This parameter mapping is our comparison; the optimum itself is prior work. The final article's measured coherent loading and calculated Fock-state loading must not be relabeled as arbitrary-qubit process fidelity.

[AO07](references.md#ao07), arXiv:quant-ph/0703203v1, Eqs. (2), (4) and pp. 1–3, already distinguishes particle survival from coherence with vacuum, including a dephasing two-level example. Thus absorption versus quantum coherence is established motivation.

For a fixed temporal mode, the capacity-sign criterion is also established. [MSG24](references.md#msg24), arXiv:2401.15634v3, Theorem 2 and Appendix B.2.2, Theorem 38, Proof 3, gives qubit-restriction antidegradability. Setting $`\lambda=p`$ and $`e^{-\gamma}=|m|^2/p`$ gives $`p+|m|^2\le1`$; $`p=0`$ is trivial. [SAJS24](references.md#sajs24), v2 Sec. 4.1 Lemma 1, Eqs. (24)–(25), gives positive coherent information on the other side through entropy singularities. Its Sec. 4.2 also supplies nonadditivity examples. Neither a fixed-mode witness nor the present channel ordering determines the positive capacity value.

[J21](references.md#j21), v2 Sec. 3.4, Theorem 6, supplies the channel-preprocessing comparison framework. Encoder composition is standard; our continuum result uses an explicit normal-channel construction rather than importing a finite-dimensional theorem unchanged. [WLJ22](references.md#wlj22), published Eqs. (1)–(6), already uses capacity to optimize bounded-coupling bosonic transducers. Its stationary, unrestricted-energy field-to-field rates are in qubits/second; our reset capture capacity is in qubits/use.

## The receiver-specific construction

The coupling ceiling makes the excess physical time $`d(u)=x(u)-u/\Gamma_M`$ nonnegative and nondecreasing. Its Wiener phase average simultaneously produces photon-block factors $`e^{-\gamma_\phi|d(u)-d(v)|}`$ and vacuum-coherence factors $`e^{-\gamma_\phi d(u)}`$. Together with attenuation and retained vacuum probability, these define one CPTP preprocessor for every state and reference system. Scalar efficiency bounds, pointwise kernel inequalities, and ideal waveform matching do not by themselves establish that statement.

The receiver-specific result combines these ingredients into an explicit fixed-noise channel identity on the entire vacuum-plus-one-photon temporal space. The cited section and equation numbers refer to the source versions identified above.

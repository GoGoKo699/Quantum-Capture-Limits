# Physical interpretation and experimental scope

The [simulation theorem](channel-simulation.md) answers a receiver-design question: with freely chosen temporal source encoding, does varying a bounded receiving rate add communication capability? For the [specified two-level Markov model](model.md), every allowed schedule is reproduced by maximum coupling preceded by a source channel. This statement covers reference entanglement and arbitrary codes; it identifies an optimal receiving schedule without determining the positive capacity.

## A directly coupled receiver

A two-level emitter coupled to a single incident waveguide port is the relevant idealization. Lin et al. study a transmon near the termination of a semi-infinite waveguide [L22](references.md#l22), an experimentally recognizable direct-coupling geometry. Their final paper, Fig. 1 and Eqs. (3)–(4), motivates this setting; it does not validate every controlled-noise assumption in our theorem. Its weak-coherent loading measurements and supplementary single-photon calculation are distinct from reconstructing the full vacuum-plus-one-photon capture channel.

The comparison fixes the same coupling ceiling, loss and dephasing rates for every receiver schedule. Within that class, modulation cannot recover positive unassisted capacity on the zero-capacity side of the [exact boundary](capacity-boundary.md), where $`\beta\Lambda_r(S)\le1`$. Changing source restrictions, the coupling ceiling, the bath or the receiving architecture changes the problem. In particular, cavity-assisted Raman interfaces such as [K18](references.md#k18) require their own effective-model justification before applying this bound.

The output is the qubit at the capture instant. Coupling control can still be useful afterward to suppress re-emission or transfer the excitation into storage. Neither long-term storage nor output transfer and reset costs are included. Capacity here is per independently reset capture use; dividing by the pulse duration does not establish a device's communication throughput.

## Rates must describe the same dynamics

The conventions at maximum coupling are

```math
\Gamma_M=g_{\max}+\kappa,\qquad
T_1^{-1}=\Gamma_M,\qquad
T_2^{-1}=\Gamma_M/2+\gamma_\phi.
```

Here $`g_{\max}`$ and $`\kappa`$ are population-decay rates into the incident and inaccessible vacuum ports; $`\gamma_\phi`$ is additional coherence decay from pure dephasing. They enter different dissipators and cannot be interchanged.

For example, [L22](references.md#l22), final supplementary Sec. S2, groups pure dephasing and half the nonradiative relaxation rate into $`\Gamma_{\phi,n}=\Gamma_\phi+\Gamma_n/2`$. Identifying this fitted combination with our $`\gamma_\phi`$ would omit the associated population loss when $`\Gamma_n`$ is appreciable. A parameter translation must account for that loss in the population dynamics as well. The [exponential-pulse comparison](pulse-objectives.md#strong-comparator-optimize-an-ordinary-exponential-for-the-correct-task) uses their pure-dephasing limit explicitly.

The theorem assumes the stated white-noise dilation remains valid under every allowed schedule, with fixed noise rates. An uncontrolled exponential decay fit alone does not establish this. Different dilations of the same free semigroup can respond differently to control [GN17](references.md#gn17), [BFH23](references.md#bfh23), particularly BFH23 Corollary 4.6. The numerical $`T_2/T_1`$ boundary is therefore a translated model criterion, not a universal hardware threshold. The [colored-noise example](counterexamples.md#markov-behavior-is-a-substantive-assumption-not-a-fitted-t2-label) exhibits the failure of the proposed phase simulator under a different correlation law.

## What free source encoding supplies

The preprocessor is a CPTP map composed with the sender's encoder. It can reshape temporal modes, attenuate photons and mix phases, without increasing photon number, postselecting, or sending an environmental record to the receiver. It need not produce one deterministic pure wavepacket. This is an exact communication comparison, not a claimed inline converter for an arbitrary unknown arriving photon.

No source bandwidth, peak flux or conversion cost is charged. Actual waveform-conversion schemes have phasematching, power and efficiency conditions [DMR15](references.md#dmr15), Secs. II and IV.C; they are not implementations of every preprocessor constructed here. For a prescribed source, receiving modulation can improve capture, as the [square-pulse example](counterexamples.md#a-fixed-source-counterexample-makes-the-resource-boundary-concrete) shows.

A shared phase reference is supplied for vacuum–one-photon coherence; its establishment and maintenance are uncharged. Independent phase randomization at each use and a common unknown phase across a block are different channels: the latter can preserve relational encodings. The collective-reference constructions of [BRS03](references.md#brs03) give a broader precedent for this distinction. A phase-uncertainty model must therefore specify correlations across uses.

## Observables that distinguish the claims

For a chosen temporal mode, excited-state population determines $`p`$, while transfer of coherence between vacuum and that mode determines $`m`$. The [qubit-fidelity formula](pulse-objectives.md) depends on both. Energy loading, coherent emitted energy and quantum process fidelity are different observables; experimental values for them cannot be substituted into one another. Even a complete single-use characterization supplies a channel, not automatically its regularized capacity. Here the all-code conclusions come from the channel simulation, antidegradability proof and explicit achievable witnesses, under the stated source and reset assumptions.

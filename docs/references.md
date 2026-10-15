# References and attribution

This source map connects the receiver-specific results to established capture and channel theory. The stochastic dilation and permitted source operations determine which control results apply.

| Subject | Established predecessor | Role here |
|---|---|---|
| Single tutorial anchor | [CKS17](#cks17) | Input–output dynamics, Fock-state coherences and loss; the local tutorial supplies the channel and capacity bridge |
| Photon-sector dynamics | [B12](#b12) | Equations used for the independent receiver calculation |
| Source/control interchange and effective time | [G07L](#g07l), [G07C](#g07c) | Earlier optimized-memory equivalence and integrated-control clock |
| Bounded two-level absorption | [WMHS12](#wmhs12), [RG17](#rg17) | Finite-window noiseless optimum and constant-rate loss limit |
| Prescribed-waveform capture | [NJY16](#njy16) | Tunable waveform matching, distinct from free-source channel ordering |
| Coherent temporal-state transfer | [YJ14](#yj14) | Earlier transfer of unknown superpositions into a linear memory network |
| Channel comparison | [J21](#j21) | Operational consequences of a constructed CPTP preprocessor |
| Controlled Markov noise | [GN17](#gn17), [BFH23](#bfh23) | Why uncontrolled T1/T2 do not specify all controlled dynamics |
| Loss versus vacuum coherence | [AO07](#ao07), [LC23](#lc23) | Earlier channel structure and distinction between survival and coherence |
| Shared reference resources | [BRS03](#brs03) | Collective reference uncertainty can preserve encoded quantum information |
| Damping–dephasing capacity sign | [MSG24](#msg24), [SAJS24](#sajs24), [CJ14](#cj14), [S21](#s21) | Established qubit criterion, extension and entropy argument |
| Capacity upper bound near the boundary | [LS09](#ls09), [SH17](#sh17) | Finite-output continuity applied to an explicitly attenuated antidegradable comparison channel |
| Capacity as a device objective | [WLJ22](#wlj22) | Stationary bosonic transducers; different input resources and rate units |
| Dephasing and exponential pulses | [L22](#l22) | Existing noisy exponential-family population optimum and direct artificial-atom setting |
| Temporal shaping | [RG17](#rg17), [W11](#w11), [DMR15](#dmr15) | Excitation optima and physical costs of source conversion |

The [prior-work comparison](prior-work.md) explains the equation-level relationships to the full input-preprocessing identity. The [physical interpretation](physical-scope.md) connects the resource assumptions to device measurements. Equation numbers below refer to the linked source versions.

<a id="ao07"></a>
## AO07

J. Åberg, D. K. L. Oi, *Generalized spectroscopy; coherence, superposition, and loss*, arXiv:quant-ph/0703203v1 (2007). [Preprint](https://arxiv.org/abs/quant-ph/0703203v1). Equations (2), (4) distinguish vacuum coherence from particle survival, including equal-loss channels with different retained coherence. The two-level-atom example includes dephasing. The broad distinction between absorption and coherent transfer is established.

<a id="b12"></a>
## B12

B. Q. Baragiola, R. L. Cook, A. M. Branczyk, J. Combes, *N-Photon wave packets interacting with an arbitrary quantum system*, Physical Review A **86**, 013811 (2012). DOI 10.1103/PhysRevA.86.013811. https://arxiv.org/abs/1202.3430 .

<a id="bfh23"></a>
## BFH23

D. Burgarth, P. Facchi, R. Hillier, *Control of Quantum Noise: On the Role of Dilations*, Annales Henri Poincaré **24**, 325–347 (2023), published online July 24, 2022. DOI 10.1007/s00023-022-01211-y. https://link.springer.com/article/10.1007/s00023-022-01211-y . Distinct dilations of the same free Markov semigroup can differ in their controllability.

<a id="bl25"></a>
## BL25

S. Bhalerao, F. Leditzky, *Improving quantum communication rates with permutation-invariant codes*, arXiv:2508.09978v1 (2025). https://arxiv.org/html/2508.09978v1 . Relevant to why a diagonal-input coherent-information lower bound is not the capacity.

<a id="brs03"></a>
## BRS03

S. D. Bartlett, T. Rudolph, R. W. Spekkens, *Classical and quantum communication without a shared reference frame*, Physical Review Letters **91**, 027901 (2003). DOI [10.1103/PhysRevLett.91.027901](https://doi.org/10.1103/PhysRevLett.91.027901). [arXiv:quant-ph/0302111v3](https://arxiv.org/abs/quant-ph/0302111v3). Equations (1)–(3) distinguish individual from collective reference uncertainty and motivate invariant encodings. This supports distinguishing one common unknown reference across a block from independent randomization of every use; it does not evaluate the capture capacity without a phase reference.

<a id="cj14"></a>
## CJ14

J. Chen, Z. Ji, D. Kribs, N. Lutkenhaus, B. Zeng, *Symmetric extension of two-qubit states*, Physical Review A 90, 032318 (2014). https://arxiv.org/abs/1310.3530 . DOI 10.1103/PhysRevA.90.032318.

<a id="cks17"></a>
## CKS17

J. Combes, J. Kerckhoff, M. Sarovar, *The SLH framework for modeling quantum input-output networks*, Advances in Physics: X **2**, 784–888 (2017). DOI [10.1080/23746149.2017.1343097](https://doi.org/10.1080/23746149.2017.1343097). [arXiv:1611.00375v2](https://arxiv.org/abs/1611.00375v2); [PDF](https://arxiv.org/pdf/1611.00375v2). This is the single external teaching anchor for the [tutorial bridge](tutorial.md). The selected route uses PDF Secs. III–IV, V.A and V.D, VII.A.3, VII.B.1, VII.C, and Appendix A item 12. Equations (153)–(155) describe Fock-state coherences and their master equations; Eq. (A13) gives a two-level emitter coupled to guided and nonguided modes. Section and equation numbers refer to the PDF; the experimental HTML conversion numbers them differently. The review supplies the input–output framework, while the repository develops its specific noisy preprocessor and communication consequences.

<a id="dmr15"></a>
## DMR15

J. M. Donohue, M. D. Mazurek, K. J. Resch, *Theory of high-efficiency sum-frequency generation for single-photon waveform conversion*, Physical Review A **91**, 033809 (2015). DOI 10.1103/PhysRevA.91.033809. https://journals.aps.org/pra/abstract/10.1103/PhysRevA.91.033809 . Included to identify temporal reshaping as a resource with physical implementation conditions, not to claim that it implements every simulator above at finite cost.

<a id="g07l"></a>
## G07L

A. V. Gorshkov, A. André, M. Fleischhauer, A. S. Sørensen, M. D. Lukin, *Universal Approach to Optimal Photon Storage in Atomic Media*, Physical Review Letters **98**, 123601 (2007). DOI [10.1103/PhysRevLett.98.123601](https://doi.org/10.1103/PhysRevLett.98.123601). [arXiv:quant-ph/0604037v3](https://arxiv.org/abs/quant-ph/0604037v3). The introduction and discussion following Eq. (8) state that optimal efficiency can be reached by changing the control or the input waveform. Equations (6)–(10) supply the spectral/time-reversal optimization for a Lambda ensemble. Source/control interchange is an antecedent, not a new general principle here.

<a id="g07c"></a>
## G07C

A. V. Gorshkov, A. André, M. D. Lukin, A. S. Sørensen, *Photon storage in Lambda-type optically dense atomic media. I. Cavity model*, Physical Review A **76**, 033804 (2007). DOI [10.1103/PhysRevA.76.033804](https://doi.org/10.1103/PhysRevA.76.033804). [arXiv:quant-ph/0612082v2](https://arxiv.org/abs/quant-ph/0612082v2). Section V.A, Eq. (15), introduces integrated control intensity as an effective time that removes the control from the adiabatic equations. Section V.B, Eqs. (23)–(26), derives optimal mode matching; Eq. (22) and nearby text address spin-wave decay. The extra physical-time phase channel in the present construction is a separate step.

<a id="gn17"></a>
## GN17

J. E. Gough, H. I. Nurdin, *Can Quantum Markov Evolutions Ever Be Dynamically Decoupled?*, IEEE CDC 2017, 6155–6160. DOI 10.1109/CDC.2017.8264587. https://arxiv.org/abs/1710.11133 . The no-decoupling statement concerns a specified quantum stochastic Markov dilation, not every microscopic realization with the same free semigroup.

<a id="j21"></a>
## J21

A. Jenčová, *A General Theory of Comparison of Quantum Channels (and Beyond)*, IEEE Transactions on Information Theory **67**, 3945–3964 (2021). DOI 10.1109/TIT.2021.3070120. https://arxiv.org/abs/2002.04240 . General simulation/preprocessing criteria and their operational interpretation are prior art.

<a id="k18"></a>
## K18

P. Kurpiers et al., *Deterministic quantum state transfer and remote entanglement using microwave photons*, Nature **558**, 264–267 (2018), arXiv:1712.08593. DOI: 10.1038/s41586-018-0195-y. The absorption and process-fidelity values cited here come from its primary abstract and are not reinterpreted as this model's average fidelity. https://arxiv.org/abs/1712.08593 .

<a id="l22"></a>
## L22

W.-J. Lin et al., *Deterministic Loading of Microwaves onto an Artificial Atom Using a Time-Reversed Waveform*, Nano Letters **22**, 8137–8142 (2022). DOI [10.1021/acs.nanolett.2c02578](https://doi.org/10.1021/acs.nanolett.2c02578). [Published paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9615994/); [final supplementary information](https://nori-physics.org/wp-content/uploads/NanoLett_22_pp8137_6_Oct_2022_SM.pdf). Supplement Sec. S8, Eq. (S22) and the following optimum give exactly the population efficiency and optimal duration of the lossless exponential family used here, with $a=2/(\Gamma\tau)$ and $r=\Gamma_{\phi,n}/\Gamma$ when nonradiative loss is negligible. Sec. S2 distinguishes $\Gamma_{\phi,n}=\Gamma_\phi+\Gamma_n/2$ from pure dephasing alone. Experimental coherent-pulse loading and calculated Fock-state loading are different observables; neither is an arbitrary-input quantum-capacity measurement. The preprint title differs from the published title.

<a id="lc23"></a>
## LC23

D. Lonigro, D. Chruściński, *Excitation-damping quantum channels*, Journal of Physics A: Mathematical and Theoretical **56**, 255301 (2023). DOI [10.1088/1751-8121/acd734](https://doi.org/10.1088/1751-8121/acd734). [arXiv:2206.04623v1](https://arxiv.org/html/2206.04623v1). Section 2.1, Eqs. (2)–(5), and Theorem 2.1 give a general excitation/ground block structure and complete-positivity criterion. Applied to finite temporal truncations with a one-dimensional excited output, this reduces to $0\preceq|b\rangle\langle b|\preceq A\preceq I$.

<a id="ls09"></a>
## LS09

D. Leung, G. Smith, *Continuity of quantum channel capacities*, Communications in Mathematical Physics **292**, 201–215 (2009). DOI 10.1007/s00220-009-0833-1. https://arxiv.org/abs/0810.4931v2 . Theorem 11 supplies the channel-by-channel telescoping method; the sharper constant used here is from [SH17](#sh17).

<a id="lsymm25"></a>
## LSYMM25

H. Luo, P. S. Shah, F. Yang, M. Mirhosseini, S. Mahmoodian, *Dynamic stimulated emission for deterministic addition and subtraction of propagating photons*, arXiv:2512.09711v1, submitted December 10, 2025. https://arxiv.org/html/2512.09711v1 . Section II.2 reviews the single-excitation tunable-coupling construction.

<a id="mp95"></a>
## MP95

S. Massar and S. Popescu, *Optimal Extraction of Information from Finite Quantum Ensembles*, Physical Review Letters **74**, 1259 (1995). DOI: 10.1103/PhysRevLett.74.1259. The single uniformly unknown qubit gives the deterministic 2/3 measure–prepare benchmark, also derived directly in [Pulse objectives](pulse-objectives.md). https://doi.org/10.1103/PhysRevLett.74.1259 .

<a id="msg24"></a>
## MSG24

F. A. Mele, F. Salek, V. Giovannetti, L. Lami, *Quantum communication on the bosonic loss-dephasing channel*, Physical Review A **110**, 012460 (2024). DOI [10.1103/PhysRevA.110.012460](https://doi.org/10.1103/PhysRevA.110.012460). [arXiv:2401.15634v3](https://arxiv.org/html/2401.15634v3). Theorem 2 and Appendix B.2.2, Theorem 38, Proof 3, establish the antidegradability boundary for its qubit restriction. With $\lambda=p$ and $e^{-\gamma}=|m|^2/p$, the criterion is $p+|m|^2\le1$. Equation (62) displays its Choi state. This is the direct antecedent credited by [SAJS24](#sajs24); the full bosonic channel has different input resources and must be distinguished from the qubit restriction.

<a id="n02"></a>
## N02

M. A. Nielsen, *A simple formula for the average gate fidelity of a quantum dynamical operation*, Physics Letters A **303**, 249–252 (2002), arXiv:quant-ph/0205035. DOI: 10.1016/S0375-9601(02)01272-0. Gives the average-fidelity identity used in the pulse comparison. https://arxiv.org/abs/quant-ph/0205035 .

<a id="njy16"></a>
## NJY16

H. I. Nurdin, M. R. James, N. Yamamoto, *Perfectly capturing traveling single photons of arbitrary temporal wavepackets with a single tunable device*, arXiv:1609.05643 (2016); associated CDC paper DOI 10.1109/CDC.2016.7798639. https://arxiv.org/abs/1609.05643 . The inspected preprint explicitly derives the prescribed-pulse coupling and discusses initial singularity truncation.

<a id="rg17"></a>
## RG17

H. S. Rag, J. Gea-Banacloche, *Two-level-atom excitation probability for single- and N-photon wave packets*, Physical Review A **96**, 033817 (2017). [arXiv:1707.02351](https://arxiv.org/abs/1707.02351). DOI [10.1103/PhysRevA.96.033817](https://doi.org/10.1103/PhysRevA.96.033817). Section II.A.1, Eqs. (1)–(8), gives the constant-rate single-photon optimum with inaccessible radiative loss but no extra pure dephasing. [Erratum, Physical Review A **99**, 029904 (2019)](https://journals.aps.org/pra/pdf/10.1103/PhysRevA.99.029904), corrects the large-photon-number Sec. III.A.2, Eqs. (38)–(39), and leaves Eq. (40) unchanged; it does not change the cited single-photon equations.

<a id="s21"></a>
## S21

V. Siddhu, *Entropic singularities give rise to quantum transmission*, Nature Communications 12, 5750 (2021). DOI 10.1038/s41467-021-25954-0. https://doi.org/10.1038/s41467-021-25954-0 .

<a id="sajs24"></a>
## SAJS24

V. Siddhu, D. Abdelhadi, T. Jochym-O'Connor, J. Smolin, *Entanglement sharing across a damping-dephasing channel*, [arXiv:2405.06231v2](https://arxiv.org/html/2405.06231v2), July 28, 2025 (initial preprint 2024). Section 4.1, Lemma 1, Eqs. (24)–(25), proves positive coherent information outside $p+|m|^2\le1$ by an entropy-singularity argument; the same section credits [MSG24](#msg24) for the exact antidegradability criterion and constructs an explicit antidegrading map. Section 4.2 supplies nonadditivity context. No new fixed-mode capacity-sign classification is asserted here.

<a id="sh17"></a>
## SH17

M. E. Shirokov, *Tight uniform continuity bounds for the quantum conditional mutual information, for the Holevo quantity, and for capacities of quantum channels*, Journal of Mathematical Physics **58**, 102202 (2017). DOI 10.1063/1.4987135. https://arxiv.org/abs/1512.09047v7 . Proposition 3A, Eq. (45), bounds output mutual-information variation by $n[2\epsilon\log d_B+g_2(\epsilon)]$, with $\epsilon$ the half diamond distance; Proposition 11, Eq. (87), gives the capacity statement. The [receiver upper bound](capacity-upper-bound.md) applies Proposition 3A directly to finite-reference codes, so it does not require the infinite-input capacity-formula assumption stated in that preprint's footnote 10. Numbering here refers to arXiv v7.

<a id="w11"></a>
## W11

Y. Wang, J. Minar, L. Sheridan, and V. Scarani, *Efficient excitation of a two-level atom by a single photon in a propagating mode*, Physical Review A **83**, 063842 (2011). DOI: 10.1103/PhysRevA.83.063842. https://doi.org/10.1103/PhysRevA.83.063842 .

<a id="wmhs12"></a>
## WMHS12

Y. Wang, J. Minář, G. Hétet, V. Scarani, *Quantum memory with a single two-level atom in a half cavity*, Physical Review A **85**, 013823 (2012). DOI [10.1103/PhysRevA.85.013823](https://doi.org/10.1103/PhysRevA.85.013823). [arXiv:1111.4261](https://arxiv.org/abs/1111.4261). Section II.C, Eqs. (19)–(24), gives finite-window capture under time-dependent radiative coupling, an optimized efficiency $1-\exp(-\int\gamma_z)$, and an explicit rate ceiling. Environmental loss is set to zero before these equations and no extra dephasing is included. Maximum-rate source-optimized absorption in that limit follows from these established formulas.

<a id="wlj22"></a>
## WLJ22

C.-H. Wang, F. Li, L. Jiang, *Quantum capacities of transducers*, Nature Communications **13**, 6698 (2022). DOI [10.1038/s41467-022-34373-8](https://doi.org/10.1038/s41467-022-34373-8). [Published article](https://pmc.ncbi.nlm.nih.gov/articles/PMC9637183/); [arXiv:2203.00012v2](https://arxiv.org/abs/2203.00012v2). Equations (1)–(5) describe stationary bosonic transducers and integrate frequency-resolved capacities into qubits per second under unrestricted input energy. The bounded coherent intermode couplings differ from the incident population-decay rate here. Capacity as a device-design objective is established; its throughput cannot be identified with the present qubits per reset use.

<a id="yj14"></a>
## YJ14

N. Yamamoto, M. R. James, *Zero-dynamics principle for perfect quantum memory in linear networks*, New Journal of Physics **16**, 073032 (2014). DOI [10.1088/1367-2630/16/7/073032](https://doi.org/10.1088/1367-2630/16/7/073032). [arXiv:1403.1698v2](https://arxiv.org/abs/1403.1698v2). Sections 4.1–4.2 derive perfect capture by zero output; Sec. 5.1, Eqs. (23)–(25), transfers unknown coefficients of temporal-mode superpositions to network excitations. This is coherent quantum-state transfer, not only an excitation-probability result. Its linear network and protected storage subsystem differ from the present noisy two-level receiver.

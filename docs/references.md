# References and attribution

The established ingredients are credited separately from the receiver-specific results. This is a source map, not a priority certificate. The exact stochastic dilation and permitted source operations matter when comparing control theorems.

| Subject | Established predecessor | Role here |
|---|---|---|
| Photon-sector dynamics | [B12](#b12) | Equations used for the independent receiver calculation |
| Prescribed-waveform capture | [NJY16](#njy16) | Tunable waveform matching, distinct from free-source channel ordering |
| Channel comparison | [J21](#j21) | Operational consequences of a constructed CPTP preprocessor |
| Controlled Markov noise | [GN17](#gn17), [BFH23](#bfh23) | Why uncontrolled T1/T2 do not specify all controlled dynamics |
| Damping–dephasing capacity sign | [SAJS24](#sajs24), [CJ14](#cj14), [S21](#s21) | Established channel criterion, extension and entropy argument |
| Temporal shaping | [RG17](#rg17), [W11](#w11), [DMR15](#dmr15) | Excitation optima and physical costs of source conversion |

The precise object to compare with these works is the full map $\mathcal E_{g,\Delta}=\mathcal E_{\max}\circ\mathcal P_{g,\Delta}$ on vacuum plus the entire one-photon temporal Hilbert space, with the coupling ceiling, source freedom, and fixed Markov noise retained. General data processing alone does not construct that map. Conversely, no new channel-comparison principle or qubit-capacity classification is asserted here.

A focused primary-source recheck on 15 October 2026 confirmed the scope of the waveform-matching, channel-comparison, and damping–dephasing antecedents below. It was not an exhaustive reference or citation-network search.

<a id="b12"></a>
## B12

B. Q. Baragiola, R. L. Cook, A. M. Branczyk, J. Combes, *N-Photon wave packets interacting with an arbitrary quantum system*, Physical Review A **86**, 013811 (2012). DOI 10.1103/PhysRevA.86.013811. https://arxiv.org/abs/1202.3430 .

<a id="bfh23"></a>
## BFH23

D. Burgarth, P. Facchi, R. Hillier, *Control of Quantum Noise: On the Role of Dilations*, Annales Henri Poincaré **24**, 325–347 (2023), published online July 24, 2022. DOI 10.1007/s00023-022-01211-y. https://link.springer.com/article/10.1007/s00023-022-01211-y . Distinct dilations of the same free Markov semigroup can differ in their controllability.

<a id="bl25"></a>
## BL25

S. Bhalerao, F. Leditzky, *Improving quantum communication rates with permutation-invariant codes*, arXiv:2508.09978v1 (2025). https://arxiv.org/html/2508.09978v1 . Relevant to why a diagonal-input coherent-information lower bound is not the capacity.

<a id="cj14"></a>
## CJ14

J. Chen, Z. Ji, D. Kribs, N. Lutkenhaus, B. Zeng, *Symmetric extension of two-qubit states*, Physical Review A 90, 032318 (2014). https://arxiv.org/abs/1310.3530 . DOI 10.1103/PhysRevA.90.032318.

<a id="dmr15"></a>
## DMR15

J. M. Donohue, M. D. Mazurek, K. J. Resch, *Theory of high-efficiency sum-frequency generation for single-photon waveform conversion*, Physical Review A **91**, 033809 (2015). DOI 10.1103/PhysRevA.91.033809. https://journals.aps.org/pra/abstract/10.1103/PhysRevA.91.033809 . Included to identify temporal reshaping as a resource with physical implementation conditions, not to claim that it implements every simulator above at finite cost.

<a id="gn17"></a>
## GN17

J. E. Gough, H. I. Nurdin, *Can Quantum Markov Evolutions Ever Be Dynamically Decoupled?*, IEEE CDC 2017, 6155–6160. DOI 10.1109/CDC.2017.8264587. https://arxiv.org/abs/1710.11133 . The no-decoupling statement concerns a specified quantum stochastic Markov dilation, not every microscopic realization with the same free semigroup.

<a id="j21"></a>
## J21

A. Jenčová, *A General Theory of Comparison of Quantum Channels (and Beyond)*, IEEE Transactions on Information Theory **67**, 3945–3964 (2021). DOI 10.1109/TIT.2021.3070120. https://arxiv.org/abs/2002.04240 . General simulation/preprocessing criteria and their operational interpretation are prior art.

<a id="k18"></a>
## K18

P. Kurpiers et al., *Deterministic quantum state transfer and remote entanglement using microwave photons*, Nature **558**, 264–267 (2018), arXiv:1712.08593. DOI: 10.1038/s41586-018-0195-y. The absorption and process-fidelity values cited here come from its primary abstract and are not reinterpreted as this model's average fidelity. https://arxiv.org/abs/1712.08593 .

<a id="lsymm25"></a>
## LSYMM25

H. Luo, P. S. Shah, F. Yang, M. Mirhosseini, S. Mahmoodian, *Dynamic stimulated emission for deterministic addition and subtraction of propagating photons*, arXiv:2512.09711v1, submitted December 10, 2025. https://arxiv.org/html/2512.09711v1 . Section II.2 reviews the single-excitation tunable-coupling construction. No latest-version claim is made.

<a id="mp95"></a>
## MP95

S. Massar and S. Popescu, *Optimal Extraction of Information from Finite Quantum Ensembles*, Physical Review Letters **74**, 1259 (1995). DOI: 10.1103/PhysRevLett.74.1259. The single uniformly unknown qubit gives the deterministic 2/3 measure–prepare benchmark, also derived directly in [Pulse objectives](pulse-objectives.md). https://doi.org/10.1103/PhysRevLett.74.1259 .

<a id="n02"></a>
## N02

M. A. Nielsen, *A simple formula for the average gate fidelity of a quantum dynamical operation*, Physics Letters A **303**, 249–252 (2002), arXiv:quant-ph/0205035. DOI: 10.1016/S0375-9601(02)01272-0. This paper also attributes earlier fidelity identities; no new average-fidelity formula is claimed here. https://arxiv.org/abs/quant-ph/0205035 .

<a id="njy16"></a>
## NJY16

H. I. Nurdin, M. R. James, N. Yamamoto, *Perfectly capturing traveling single photons of arbitrary temporal wavepackets with a single tunable device*, arXiv:1609.05643 (2016); associated CDC paper DOI 10.1109/CDC.2016.7798639. https://arxiv.org/abs/1609.05643 . The inspected preprint explicitly derives the prescribed-pulse coupling and discusses initial singularity truncation.

<a id="rg17"></a>
## RG17

H. S. Rag, J. Gea-Banacloche, *Two-level atom excitation probability for single- and N-photon wavepackets*, Physical Review A 96, 033817 (2017). https://arxiv.org/abs/1707.02351 . DOI 10.1103/PhysRevA.96.033817. Erratum: Physical Review A 99, 029904 (2019); include it in a formula-level comparison.

<a id="s21"></a>
## S21

V. Siddhu, *Entropic singularities give rise to quantum transmission*, Nature Communications 12, 5750 (2021). DOI 10.1038/s41467-021-25954-0. https://doi.org/10.1038/s41467-021-25954-0 .

<a id="sajs24"></a>
## SAJS24

V. Siddhu, D. Abdelhadi, T. Jochym-O'Connor, J. Smolin, *Entanglement sharing across a damping-dephasing channel*, arXiv:2405.06231v2. https://arxiv.org/html/2405.06231v2 . Established capacity-sign criterion, entropic argument, and nonadditive-information context; no new general qubit-capacity classification is claimed.

<a id="w11"></a>
## W11

Y. Wang, J. Minar, L. Sheridan, and V. Scarani, *Efficient excitation of a two-level atom by a single photon in a propagating mode*, Physical Review A **83**, 063842 (2011). DOI: 10.1103/PhysRevA.83.063842. https://doi.org/10.1103/PhysRevA.83.063842 .

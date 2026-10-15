# Repository guide

This repository concerns bounded-coupling, Markov-noisy, single-photon capture. Read README.md, docs/model.md, then docs/channel-simulation.md. The lead result is a full input-preprocessing identity, not just an absorption inequality.

## Scientific invariants

- Keep radiative population decay distinct from EXTRA coherence decay. Define g, kappa, gamma_phi, Gamma_M, beta, r, S consistently.
- State free temporal source encoding, independent reset, shared phase reference and one photon PER USE. Do not turn the result into a fixed-source, finite-bandwidth, assisted, multiphoton, or long-term memory claim.
- Preserve vacuum coherence, trace preservation, reference entanglement and lost-photon vacuum weight. Never renormalize a retained photon subensemble.
- The simulator is input preprocessing, not a claimed output degradation. Pointwise kernel inequalities are not automatically Loewner inequalities.
- Channel ordering identifies an optimal schedule, not the positive capacity value. A coherent-information witness is a lower bound, not an optimized or regularized capacity. The p+|m|^2 criterion is attributed to established damping–dephasing results.
- The antidegradable-comparison capacity upper bound covers regularized capacity but is not claimed tight. Its epsilon is HALF the diamond norm. Use exact spectra or spectral upper estimates; never use a Galerkin lower bound as an upper certificate. Attribute the continuity inequality to its established source.
- A failed proposed simulator under colored noise does not prove that every simulator fails or that modulation has positive-capacity advantage there.

## Maintenance

Keep the reading path centered on the present theorems, proofs and examples. Do not import unrelated projects, exploratory conversation transcripts, journal targets, progress bars or an omnibus unfinished-work list. Necessary scientific assumptions belong next to their claims.

Run `python tools/verify.py` from the repository root after any scientific or computational edit. Assertions must remain enabled. Generated raw reports belong in artifacts/; they are not proofs. Update verification/validated.json only from an actual successful run using `python tools/verify.py --record`. Never hand-edit a PASS or claim that GitHub Actions succeeded without checking the workflow result.

Use `$...$` for inline mathematics and fenced `math` blocks for display equations. Keep relative links valid. Source code is deduplicated in capture/; independent checks are in tests/. Retain independent numerical formulations instead of making a test call only the formula it purportedly verifies. Preserve LICENSE.

## Discovery and citation

Relevant terms: quantum capture, single-photon absorption, two-level receiver, temporal modes, channel simulation, bounded coupling, Markov dephasing, antidegradability, quantum capacity. Cite the exact theorem together with its resource assumptions and cite its antecedents from docs/references.md. This is a research repository; passing local checks is not independent review or a priority certificate.

# Reviews

[`independent-review.md`](independent-review.md) is an independent, AI-assisted
adversarial review of the core uniqueness proof (Sections 2–3 of the
[paper](../papers/injectivity-of-the-tent-transform.md)) and of the auxiliary periodic
claims. The reviewer was instructed to find a counterexample, an invalid implication, or a
missing hypothesis, and to re-derive the central identities independently.

**Verdict of the review:** no counterexample, fatal implication, or missing hypothesis was
found in the core proof. It recommended three changes, all incorporated in the paper:
an explicit finite-side energy estimate, an explicit use of the integrability of $f^{-2}$
when proving $f\to\infty$, and a precise statement of the competitor class.

Scope: this is not peer review, formal verification, or a claim of publication. Numerical
checks in [`code/`](../code), in particular `kernel_lemma_checks.py`, test the kernel
lemmas used in the proof; they are evidence, not proof.

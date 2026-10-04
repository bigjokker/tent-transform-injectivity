# Appendices

These appendices contain the proofs of the results quoted in the
[paper](../papers/injectivity-of-the-tent-transform.md) outside its main theorem. The main
theorem (Theorem 3.1) is proved in full in the paper itself.

| Appendix | Contents | Used in |
|---|---|---|
| [A. Wave identity, local information, and mass detection](wave-identity-and-mass-detection.md) | The identity $f''=2gf-4I$ for arbitrary locally integrable inputs; support, local-mass and pointwise bounds; detection of finite mass on each half-line; identification of integrable densities | Sections 2, 4, 5 |
| [B. Periodic local inversion, stability, and regularity](periodic-inversion-stability-regularity.md) | Local $C^1$ inversion $L^\infty\to W^{2,\infty}$ at every bounded nonnegative periodic density; Lipschitz stability on Hölder-bounded classes; $f\in C^{k+2,\alpha}$ iff $g\in C^{k,\alpha}$ | Section 5 |
| [C. Controlled tails](controlled-tails.md) | The approximation $g=\pi/(4f^2)+O(x^{-2})$ on controlled tails; a forward $L^1$ estimate; injectivity on the controlled-tail class | Sections 4, 5 |
| [D. Periodic reconstruction with an unknown mean](periodic-reconstruction.md) | Certified terminating search for the mean; geometric convergence of the primitive iteration | Section 5 |

Numerical checks corresponding to each appendix are in [`code/`](../code); run them all
with `python scripts/run_checks.py`.

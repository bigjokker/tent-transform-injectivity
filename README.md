# Injectivity of an Exponential Tent Transform

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23138548.svg)](https://zenodo.org/doi/10.5281/zenodo.23138548)

Proofs and reproducible numerical checks for the uniqueness of the inverse problem for

$$
T(g)(x)=\int_0^\infty\exp\Big(-\int_0^t\!\!\int_{x-\tau}^{x+\tau}g(y)\,dy\,d\tau\Big)\,dt ,
$$

defined for nonnegative, locally integrable densities $g$ on $\mathbb R$ that are not a.e.
zero. In the one-dimensional Kolmogorov–Johnson–Mehl–Avrami model, $T(g)(x)$ is the
expected transformation time at $x$ for nucleation rate $g$ and unit growth speed.

**Main theorem.** If $T(g_1)=T(g_2)$ and $\sup_a\int_a^{a+1}|g_1-g_2|<\infty$, then
$g_1=g_2$ almost everywhere.

Consequences:

- $T$ is injective on nonnegative densities in $L^p(\mathbb R)$ for every $1\le p\le\infty$,
  including $L^2$ densities of infinite mass, with no smoothness or decay assumption.
- Equal transforms force equality whenever the difference is in $L^p$, bounded, or
  periodic, including periodic densities with different periods.
- An integrable density is determined by its transform among all nonnegative locally
  integrable densities.

The work also includes the exact identity $T(g)''=2g\,T(g)-4I$, information about $g$
read directly from $T(g)$, an injectivity theorem for controlled tails, and, for periodic
densities, local inversion, Lipschitz stability on Hölder-bounded classes, a regularity
equivalence, and a certified reconstruction procedure.

No closed-form inverse, range characterization, or real-line reconstruction scheme is
claimed. Pairs whose difference has unbounded mass on unit intervals (outside the
controlled-tail class) and signed densities are not covered. No literature-priority claim
is made.

## Read the results

- [Paper: Injectivity of an Exponential Tent Transform](papers/injectivity-of-the-tent-transform.md)
- [Appendices](docs/README.md)
- [Independent review and its scope](reviews/README.md)

## Reproduce

Use **Python 3.11 or later** with NumPy and SciPy. From the repository root:

```sh
pip install -r requirements.txt
python scripts/run_checks.py
```

The full suite takes under a minute. `python scripts/run_checks.py --quick` skips the
kernel-lemma stress test. Each check exits with a nonzero status on failure; output files
are written to `scratch/`, which is not tracked.

The checks are numerical illustrations of the analytic arguments. They are not part of
the proofs and are not interval-arithmetic certificates.

## Repository layout

| Directory | Contents |
|---|---|
| `papers/` | The main manuscript, with the complete proof of the main theorem |
| `docs/` | Appendices A–D: wave identity and mass detection, periodic results, controlled tails, periodic reconstruction |
| `code/` | Forward map and numerical checks |
| `scripts/` | Runner for all checks |
| `reviews/` | Independent AI-assisted review of the core proof |

## Status

The results were developed with AI assistance and checked by an independent AI-assisted
adversarial review, which found no counterexample or fatal defect in the core proof. They
have not been peer reviewed or formally verified.

## Citation

Archived on Zenodo. Cite the concept DOI to refer to the work in general, or
the version DOI to pin a specific release.

* All versions: [10.5281/zenodo.23138548](https://zenodo.org/doi/10.5281/zenodo.23138548)
* v1.0.0: [10.5281/zenodo.23138549](https://zenodo.org/doi/10.5281/zenodo.23138549)

```bibtex
@misc{tent_transform_injectivity,
  title  = {Injectivity of an Exponential Tent Transform},
  author = {bigjokker},
  year   = {2026},
  doi    = {10.5281/zenodo.23138548},
  url    = {https://github.com/bigjokker/tent-transform-injectivity}
}
```

See also [CITATION.cff](CITATION.cff).

This repository is released under the MIT License. See [LICENSE](LICENSE).

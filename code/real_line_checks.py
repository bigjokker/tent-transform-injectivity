"""Calculation checks for the real-line identities (papers/, Sections 2-3).

Checks the autoconvolution identity and an example whose extremum escapes to infinity.
These are numerical illustrations, not a uniqueness proof.
"""

import numpy as np
from scipy.integrate import quad

from transform import H, forward, gaussian_G2


def check_autoconvolution():
    components = [(0.5, 0.4, -2.0), (0.9, 1.2, 1.0), (0.3, 0.7, 3.0)]
    primitives = [(gaussian_G2(a, s), center) for a, s, center in components]

    def G(z):
        return sum(primitive(np.asarray(z) - center)
                   for primitive, center in primitives)

    locations = np.array([-6.0, -2.0, 0.0, 1.3, 4.0, 9.0])
    numeric = forward(G, locations)
    errors = []
    for x, value in zip(locations, numeric):
        # This is (W*W)(2x)/(2W(x)^2), evaluated on the whole real line.
        # Scale inside the exponential to avoid underflow in tiny convolutions.
        convolution = 0.5 * quad(
            lambda y: np.exp(2.0*G(x) - G(y) - G(2.0*x-y)),
            -np.inf, np.inf, epsabs=1e-11, epsrel=1e-11, limit=300,
        )[0]
        errors.append(abs(convolution - value))
    maximum = max(errors)
    assert maximum < 1e-8, maximum
    print(f"Autoconvolution identity: max absolute forward error = {maximum:.3e}")


def check_escaping_extremum():
    a, s, shift = 1.0, 1.0, 1.0
    centered = gaussian_G2(a, s)
    half_mass = a*s*np.sqrt(np.pi/2.0)

    def normalized_G(z):
        return centered(z) + half_mass*np.asarray(z)

    def Q(z):
        return normalized_G(np.asarray(z)-shift) - normalized_G(z)

    # Equal masses and a translated mean: Q lies strictly between its
    # endpoint limits. Its infimum is approached only at positive infinity.
    first_moment_difference = 2.0*half_mass*shift
    sample = np.linspace(-5.0, 5.0, 101)
    values = Q(sample)
    assert np.all(values > -first_moment_difference)
    assert np.all(values < 0.0)
    locations = np.array([2.0, 4.0, 6.0])
    f1 = forward(lambda z: centered(np.asarray(z)-shift), locations)
    f2 = forward(centered, locations)
    assert np.all(f1 < f2), f1-f2
    print("Escaping-extremum example, shifted Gaussian minus original:")
    for x, difference in zip(locations, f1-f2):
        print(f"  x={x:g}: Q(x)={Q(x):.9f}, f1-f2={difference:.9f}")

    # Check the penalty's algebraic sign over the interval used in the proof.
    x = 20.0
    tau = np.sqrt(2.0*x*x-2.0)
    t = np.linspace(0.0, tau, 1001)
    second_difference = np.log1p((x+t)**2) + np.log1p((x-t)**2) - 2*np.log1p(x*x)
    assert np.max(second_difference) < 1e-13
    print(f"Logarithmic penalty: max second difference on [0,tau] = "
          f"{np.max(second_difference):.3e}")


if __name__ == "__main__":
    check_autoconvolution()
    check_escaping_extremum()
    print("All real-line calculation checks passed.")

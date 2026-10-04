"""Forward map g -> f = T(g) for the exponential tent transform.

    f(x) = ∫_0^∞ exp(-H(t,x)) dt,   H(t,x) = ∫_{-t}^{t} (t-|u|) g(x+u) du.

Key identity: if G2 is any second primitive of g (G2'' = g), then

    H(t,x) = G2(x+t) + G2(x-t) - 2 G2(x),

because both sides vanish with their t-derivative at t = 0 and both have
∂_tt = g(x+t) + g(x-t). So H costs three lookups instead of a double integral.
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson
from scipy.special import erf


class SecondPrimitive:
    """Numerical G2 (with G2'' = g), tabulated on [-L, L], linear outside.

    The linear extension is exact when g vanishes outside [-L, L], so pick L
    large enough that g is negligible there.
    """

    def __init__(self, g, L, h=1e-3):
        self.L = L
        self.y = np.linspace(-L, L, int(round(2 * L / h)) + 1)
        gy = g(self.y)
        self.G1 = cumulative_trapezoid(gy, self.y, initial=0.0)
        self.G2 = cumulative_trapezoid(self.G1, self.y, initial=0.0)

    def __call__(self, z):
        z = np.asarray(z, dtype=float)
        out = np.interp(z, self.y, self.G2)
        out = np.where(z > self.L, self.G2[-1] + self.G1[-1] * (z - self.L), out)
        out = np.where(z < -self.L, self.G2[0] + self.G1[0] * (z + self.L), out)
        return out


def constant_G2(c):
    """Exact G2 for g = c, so H(t,x) = c t^2."""
    return lambda z: 0.5 * c * np.asarray(z, dtype=float) ** 2


def gaussian_G2(A, s):
    """Exact G2 for g(y) = A exp(-y^2 / (2 s^2))."""
    a = s * np.sqrt(2.0)
    k = A * s * np.sqrt(np.pi / 2.0)

    def G2(z):
        z = np.asarray(z, dtype=float)
        return k * (z * erf(z / a) + a / np.sqrt(np.pi) * np.exp(-(z / a) ** 2))

    return G2


def H(G2, t, x):
    return G2(x + t) + G2(x - t) - 2.0 * G2(x)


def forward(G2, x, dt=2e-3, H_stop=40.0, T0=4.0, T_max=1e5):
    """Compute f(x) for each x by Simpson's rule in t.

    H is nondecreasing and convex in t, so the t-range is doubled until
    H(T,x) > H_stop. Convexity bounds the dropped tail by
    T*exp(-H(T,x))/H(T,x), not just exp(-H_stop). This routine does not
    return a certified quadrature error bound.
    dt must resolve the decay of e^{-H}; for g ~ c that scale is 1/sqrt(c).
    """
    x = np.atleast_1d(np.asarray(x, dtype=float))
    f = np.empty_like(x)
    for i, xi in enumerate(x):
        T = T0
        while H(G2, T, xi) < H_stop:
            T *= 2.0
            if T > T_max:
                raise RuntimeError(f"H(t, {xi}) stays below {H_stop} up to t = {T_max}; is g zero?")
        t = np.linspace(0.0, T, int(np.ceil(T / dt)) + 1)
        f[i] = simpson(np.exp(-H(G2, t, xi)), x=t)
    return f


def inverse_leading(f):
    """Leading-order inverse g ≈ π / (4 f^2), exact for constant g."""
    return np.pi / (4.0 * np.asarray(f) ** 2)


def inverse_corrected(f, x):
    """First correction for slowly varying g: g ≈ g0 (1 - g0'' / (8 g0^2)), g0 = π/(4f^2).

    From H ≈ g t^2 + g'' t^4 / 12, which gives f ≈ (√π/2) g^{-1/2} (1 - g''/(16 g^2)).
    g0'' is taken by finite differences on the x grid, so x must be uniform.
    """
    g0 = inverse_leading(f)
    g0_xx = np.gradient(np.gradient(g0, x), x)
    return g0 * (1.0 - g0_xx / (8.0 * g0 ** 2))

"""Numerical checks of the periodic linearization (papers/, Section 5).

Uses finite Fourier sums for the periodic tent operator, with scalar
quadrature as an independent check of the analytic derivative multipliers.
These checks do not prove nonlinear injectivity.
"""

import numpy as np
from scipy.integrate import quad, simpson


def multiplier(c, n):
    if n == 0:
        return -np.sqrt(np.pi) / (4.0 * c**1.5)
    return -np.sqrt(np.pi / c) * (-np.expm1(-n*n / (4.0*c))) / (n*n)


def tent_mode(t, n):
    if n == 0:
        return t*t
    # 4 sin^2(nt/2) avoids cancellation for small t.
    return 4.0 * np.sin(n*t / 2.0)**2 / (n*n)


def forward_modes(c, modes, x, dt=0.002, cutoff=45.0):
    """Compute T for c + sum amplitude*cos(nx) + amplitude*sin(nx).

    Each mode is (frequency, cosine amplitude, sine amplitude). The density
    lower bound is conservative. If it is m>0, K_t g >= m t^2, and the
    discarded integral is <= exp(-m T^2)/(2 m T).
    """
    lower = c - sum(np.hypot(a, b) for _, a, b in modes)
    if lower <= 0:
        raise ValueError("A positive certified lower bound is required")
    end = np.sqrt(cutoff / lower)
    count = int(np.ceil(end / dt))
    count += count % 2  # even number of subintervals for Simpson
    t = np.linspace(0.0, end, count + 1)
    exponent = np.broadcast_to(c*t*t, (len(x), len(t))).copy()
    for n, a, b in modes:
        shape = a*np.cos(n*x) + b*np.sin(n*x)
        exponent += shape[:, None] * tent_mode(t, n)[None, :]
    return simpson(np.exp(-exponent), x=t, axis=1)


def sobolev_norm(values, s):
    n = np.fft.fftfreq(len(values), d=1.0 / len(values))
    coeff = np.fft.fft(values) / len(values)
    return np.sqrt(np.sum((1.0 + n*n)**s * np.abs(coeff)**2))


def check_scalar_integrals():
    errors = []
    for c in (0.25, 1.0, 4.0):
        for n in (0, 1, 4, 16, 64):
            if n == 0:
                numeric = -quad(lambda t: t*t*np.exp(-c*t*t), 0, np.inf,
                                epsabs=1e-12, epsrel=1e-12)[0]
            else:
                gaussian = np.sqrt(np.pi / c) / 2.0
                cosine = quad(lambda t: np.exp(-c*t*t), 0, np.inf,
                              weight="cos", wvar=n,
                              epsabs=1e-12, limit=200)[0]
                numeric = -2.0 * (gaussian - cosine) / (n*n)
            exact = multiplier(c, n)
            errors.append(abs(numeric - exact) / abs(exact))
    maximum = max(errors)
    assert maximum < 1e-9, maximum
    print(f"Scalar quadrature: max relative multiplier error = {maximum:.3e}")


def check_nonlinear_derivative(x):
    eps = 1e-4
    errors = []
    for c in (0.25, 1.0, 4.0):
        for n in (0, 1, 4, 16, 64):
            direction = np.cos(n*x)
            plus = forward_modes(c, [(n, eps, 0.0)], x)
            minus = forward_modes(c, [(n, -eps, 0.0)], x)
            numeric = (plus - minus) / (2.0*eps)
            exact = multiplier(c, n) * direction
            errors.append(np.max(np.abs(numeric - exact)) / abs(multiplier(c, n)))
    maximum = max(errors)
    assert maximum < 2e-6, maximum
    print(f"Nonlinear central difference: max scaled error = {maximum:.3e}")


def check_remainder(x):
    c = 1.0
    modes = [(0, 0.3, 0.0), (1, 1.0, 0.0), (3, 0.0, 0.2)]
    direction = sum(a*np.cos(n*x) + b*np.sin(n*x) for n, a, b in modes)
    linear = sum(multiplier(c, n)*(a*np.cos(n*x) + b*np.sin(n*x))
                 for n, a, b in modes)
    baseline = np.sqrt(np.pi / c) / 2.0
    input_norm = sobolev_norm(direction, 1)
    norms = []
    print("Quadratic remainder for a perturbation including the mean:")
    print("  epsilon       ||remainder||_H3       remainder/(epsilon^2 ||h||_H1^2)")
    for eps in (0.08, 0.04, 0.02, 0.01):
        data = forward_modes(c, [(n, eps*a, eps*b) for n, a, b in modes], x)
        remainder = data - baseline - eps*linear
        norm = sobolev_norm(remainder, 3)
        norms.append(norm)
        print(f"  {eps:7.3f}       {norm:17.8e}       {norm/(eps*input_norm)**2:12.8f}")
    orders = np.log2(np.array(norms[:-1]) / np.array(norms[1:]))
    assert np.all((orders > 1.8) & (orders < 2.2)), orders
    print("  Observed orders:", ", ".join(f"{order:.4f}" for order in orders))


def check_frequency_scaling():
    n = np.arange(1, 10001)
    m = -np.sqrt(np.pi) * (-np.expm1(-n*n / 4.0)) / (n*n)
    weights = (1.0 + n*n)*np.abs(m)
    assert np.min(weights) > 0.39
    assert abs(n[-1]**2 * abs(m[-1]) - np.sqrt(np.pi)) < 1e-12
    print(f"At c=1, modes 1--10000: (1+n^2)|m_n| in "
          f"[{weights.min():.6f}, {weights.max():.6f}]")
    print(f"At n=10000: n^2 |m_n| = {n[-1]**2 * abs(m[-1]):.9f}")


if __name__ == "__main__":
    x = np.arange(512) * (2.0*np.pi / 512)
    check_scalar_integrals()
    check_nonlinear_derivative(x)
    check_remainder(x)
    check_frequency_scaling()
    print("All periodic calculation checks passed.")

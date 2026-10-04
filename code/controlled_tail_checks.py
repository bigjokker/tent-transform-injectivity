"""Checks for docs/controlled-tails.md; run python code/controlled_tail_checks.py.

These floating-point calculations support, but do not prove, the results.
The density is evaluated on the whole line. Time-tail bounds are derived
from convexity; their printed numerical values are not interval arithmetic.
"""

import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre


def transform_at(g, x, order=96, cutoff=48.0):
    """Integrate in units of 1/sqrt(g(x)), avoiding primitive cancellation."""
    nodes, weights = roots_legendre(order)
    v, weights = (nodes + 1.0) / 2.0, weights / 2.0
    a = float(g(x))
    scale = 1.0 / np.sqrt(a)

    def exponent(s):
        t = s * scale
        samples = g(x + t * v) + g(x - t * v)
        return t * t * np.dot(weights * (1.0 - v), samples)

    end = 2.0
    while exponent(end) < cutoff:
        end *= 2.0
    value, quad_error = quad(lambda s: np.exp(-exponent(s)), 0.0, end,
                            epsabs=2e-12, epsrel=2e-12, limit=200)
    t = end * scale
    ht = t * np.dot(weights, g(x + t * v) + g(x - t * v))
    tail_bound = np.exp(-exponent(end)) / ht
    return value * scale, quad_error * scale, tail_bound


def check_power_tails():
    records = []
    for beta in (0.7, 1.0, 1.6):
        g = lambda y: (1.0 + np.asarray(y)**2)**(-beta / 2.0)
        expected = beta * (beta + 1.0) / 8.0
        for x in (30.0, 100.0, 300.0, 1000.0):
            f, quad_error, tail_bound = transform_at(g, x)
            f_fine, _, _ = transform_at(g, x, order=192)
            estimate = np.pi / (4.0 * f**2)
            scaled_error = x*x*(estimate - g(x))
            refinement = abs(f-f_fine) / f_fine
            assert refinement < 2e-10, (beta, x, refinement)
            assert tail_bound < 1e-14 * f, (beta, x, tail_bound)
            assert 0.0 < scaled_error < 1.0, (beta, x, scaled_error)
            record = dict(beta=beta, x=x, f=f, scaled_error=scaled_error,
                          predicted_limit=expected, quadrature_error=quad_error,
                          time_tail_bound=tail_bound, relative_refinement=refinement)
            records.append(record)
            print(f"beta={beta:.1f}, x={x:6.0f}: x^2(A(g)-g)={scaled_error:.8f}; "
                  f"limit={expected:.8f}; refinement={refinement:.1e}")
        # beta=1.6 approaches its limit more slowly (relative correction x^-0.4).
        tolerance = 0.025 if beta < 1.5 else 0.25
        assert abs(records[-1]['scaled_error'] - expected) < tolerance
    return records


def check_log_oscillation():
    def g(y):
        y = np.asarray(y)
        return (1.0+y*y)**(-0.35) * (1.0+0.3*np.sin(0.5*np.log1p(y*y)))

    records = []
    for x in (30.0, 100.0, 300.0, 1000.0):
        f, _, tail = transform_at(g, x)
        f_fine, _, _ = transform_at(g, x, order=192)
        scaled = x*x*(np.pi/(4.0*f*f)-g(x))
        assert abs(f-f_fine) / f_fine < 2e-10
        assert abs(scaled) < 2.0
        records.append(dict(x=x, scaled_error=scaled, time_tail_bound=tail))
        print(f"log-oscillating tail, x={x:6.0f}: x^2(A(g)-g)={scaled:.8f}")
    return records


def check_tent_increment():
    # Directly test the pointwise inequality behind the global L1 estimate.
    rng = np.random.default_rng(814)
    x, y, z = rng.normal(size=(3, 20000)) * 20.0
    s = rng.exponential(10.0, size=len(x))
    d = abs(x-y)
    lhs = np.maximum(d+s-abs(z-x), 0.0)-np.maximum(d-abs(z-x), 0.0)
    rhs = np.maximum(s-abs(z-y), 0.0)
    residual = float(np.min(lhs-rhs))
    assert residual > -5e-14, residual
    print(f"Tent-increment inequality: minimum numerical residual {residual:.2e}")
    return residual


def check_small_local_mass():
    # Here x^2*g(x) -> 0, so a local Gaussian approximation to f is
    # inappropriate. An exact primitive avoids underresolving the central
    # mass while testing the second branch of the proof.
    primitive = lambda y: 0.5*y*np.arctan(y)
    records = []
    for x in (30.0, 100.0, 300.0, 1000.0):
        def exponent(t):
            return primitive(x+t)+primitive(x-t)-2*primitive(x)

        first = quad(lambda t: np.exp(-exponent(t)), 0.0, x,
                     epsabs=1e-9, epsrel=1e-11, limit=200)[0]
        rest = quad(lambda s: np.exp(-exponent(x+s)), 0.0, 40.0,
                    epsabs=1e-11, epsrel=1e-11)[0]
        f = first+rest
        scaled = x*x*(np.pi/(4*f*f)-(1+x*x)**-2)
        assert 0.0 < scaled < 1.0
        records.append(dict(x=x, f=f, scaled_error=scaled))
        print(f'x^2*g -> 0, x={x:6.0f}: x^2(A(g)-g)={scaled:.8f}; '
              f'limit={np.pi/4:.8f}')
    assert abs(records[-1]['scaled_error']-np.pi/4) < 0.005
    return records


if __name__ == '__main__':
    result = dict(power_tails=check_power_tails(),
                  log_oscillation=check_log_oscillation(),
                  small_local_mass=check_small_local_mass(),
                  tent_increment_residual=check_tent_increment())
    path = Path(__file__).resolve().parents[1] / 'scratch' / 'controlled_tail_checks.json'
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(f"All checks passed. Saved {path}")

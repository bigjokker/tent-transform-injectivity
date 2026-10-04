"""Periodic inversion with two trial means and explicit unresolved status.

The positive quadrature map is a finite model of the continuum transform.
Floating-point residual signs are numerical evidence, not validated
continuum certificates. In particular, exhausting an iteration budget
never collapses the mean bracket or claims that a root was found.

Run python code/periodic_reconstruction.py for regression checks and a
reconstruction from data produced by independent continuum quadrature.
"""

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp


@dataclass
class MeanResult:
    status: str
    bracket: tuple
    primitive: np.ndarray
    refinements: int
    trial_steps: int

    @property
    def mean(self):
        return sum(self.bracket) / 2.0


class PeriodicMap:
    def __init__(self, points=128, period=2*np.pi, time_end=16.0):
        self.points = points
        self.period = period
        self.step = period / points
        count = int(np.ceil(time_end / self.step))
        self.time = np.arange(count+1) * self.step
        weights = np.full(count+1, self.step)
        weights[[0, -1]] *= 0.5
        self.log_weights = np.log(weights)
        i = np.arange(points)[:, None]
        k = np.arange(count+1)[None, :]
        self.plus = (i+k) % points
        self.minus = (i-k) % points
        self.frequency = 2*np.pi*np.fft.fftfreq(points, d=self.step)

    def phi(self, c, p, logf):
        exponent = (self.log_weights-c*self.time**2
                    -p[self.plus]-p[self.minus])
        return 0.5*(logf-logsumexp(exponent, axis=1))

    def gaussian_integral(self, c):
        return float(np.exp(logsumexp(self.log_weights-c*self.time**2)))

    def gaussian_inverse(self, target):
        # This brackets the finite quadrature model, not a silently
        # substituted continuum Gaussian integral.
        minimum = self.step / 2.0
        maximum = self.time[-1]
        if not minimum < target < maximum:
            raise ValueError('Data outside this quadrature model; refine the grid/time range')
        high = max(1.0, np.pi/(4*target*target))
        while self.gaussian_integral(high) > target:
            high *= 2.0
        return brentq(lambda c: self.gaussian_integral(c)-target, 0.0, high,
                      xtol=1e-14, rtol=1e-14)

    def initial_bracket(self, f):
        return (self.gaussian_inverse(float(np.max(f))),
                self.gaussian_inverse(float(np.min(f))))

    def density(self, c, p):
        return c+np.fft.ifft(-self.frequency**2*np.fft.fft(p)).real


class Trial:
    def __init__(self, c, initial):
        self.c = c
        self.initial = initial.copy()
        self.p = initial.copy()
        self.offset = 0.0
        self.steps = 0

    def advance(self, model, logf, sign_guard):
        q = model.phi(self.c, self.p, logf)
        one_step = q-self.p
        shift = float(np.mean(q))
        self.offset += shift
        self.p = q-shift
        self.steps += 1
        n_step = self.p+self.offset-self.initial
        # Both one-step and accumulated n-step tests are sound for the
        # exact finite map, assuming in-range data. The guard is a
        # numerical tolerance, not a rigorous roundoff enclosure.
        for residual, guard in ((one_step, sign_guard),
                                (n_step, sign_guard*self.steps)):
            if float(np.min(residual)) > guard:
                return 1
            if float(np.max(residual)) < -guard:
                return -1
        return 0


def find_mean(model, f, tolerance=1e-8, max_trial_steps=1000,
              sign_guard=1e-12, bracket=None, initial=None):
    f = np.asarray(f, dtype=float)
    if f.shape != (model.points,) or not np.all(np.isfinite(f)) or np.any(f <= 0):
        raise ValueError('Expected one finite positive datum per grid point')
    if tolerance <= 0 or sign_guard < 0 or max_trial_steps < 0:
        raise ValueError('Invalid tolerance or iteration budget')
    a, b = model.initial_bracket(f) if bracket is None else bracket
    if not 0 < a <= b:
        raise ValueError('Expected a positive mean bracket')
    p = np.zeros(model.points) if initial is None else np.asarray(initial, dtype=float).copy()
    logf = np.log(f)
    refinements = steps = 0
    while b-a > tolerance:
        trials = [Trial(a+(b-a)/3.0, p), Trial(a+2.0*(b-a)/3.0, p)]
        found = False
        for _ in range(max_trial_steps):
            for trial in trials:
                sign = trial.advance(model, logf, sign_guard)
                steps += 1
                if sign:
                    if sign > 0:
                        b = trial.c
                    else:
                        a = trial.c
                    p = trial.p
                    refinements += 1
                    found = True
                    break
            if found:
                break
        if not found:
            # Crucial: retain the entire last valid bracket.
            return MeanResult('unresolved', (a, b), p, refinements, steps)
    return MeanResult('mean_tolerance_reached', (a, b), p, refinements, steps)


def reconstruct(model, f, **kwargs):
    result = find_mean(model, f, **kwargs)
    if result.status != 'mean_tolerance_reached':
        return result, None
    logf, p = np.log(f), result.primitive
    for _ in range(2000):
        q = model.phi(result.mean, p, logf)
        q -= np.mean(q)
        change = np.ptp(q-p)
        p = q
        if change < 1e-13:
            result.primitive = p
            return result, model.density(result.mean, p)
    result.status = 'primitive_unresolved'
    result.primitive = p
    return result, None


def regression_checks():
    model = PeriodicMap(points=64)
    f = np.full(model.points, model.gaussian_integral(1.0))
    # The first trisection point is exactly the root, yet the other point
    # decides a refinement. Constant data keep this test independent of
    # errors in approximating an unknown primitive.
    trial = Trial(1.0, np.zeros(model.points))
    assert all(trial.advance(model, np.log(f), 1e-12) == 0 for _ in range(5))
    result = find_mean(model, f, bracket=(0.5, 2.0), tolerance=1e-7)
    assert result.status == 'mean_tolerance_reached', result
    assert result.bracket[0] <= 1.0 <= result.bracket[1], result.bracket
    assert result.bracket[1]-result.bracket[0] <= 1e-7
    # An exhausted budget is not an equality certificate.
    stopped = find_mean(model, f, bracket=(0.5, 2.0), max_trial_steps=0)
    assert stopped.status == 'unresolved' and stopped.bracket == (0.5, 2.0)
    # A tolerance can make every sign indecisive. This also must not
    # collapse a bracket, even after actual iterations were performed.
    guarded = find_mean(model, f, bracket=(0.5, 2.0), max_trial_steps=3, sign_guard=10.0)
    assert guarded.status == 'unresolved' and guarded.bracket == (0.5, 2.0)
    constant, density = reconstruct(model, f)
    assert constant.status == 'mean_tolerance_reached'
    assert np.max(abs(density-1.0)) < 1e-11
    print('Exact-root trisection, constant data, and unresolved-budget checks passed.')


def independent_data_checks():
    from periodic_linearization_checks import forward_modes

    results = []
    modes = [(1, 0.35, 0.0), (3, 0.0, 0.15), (5, 0.08, 0.0)]
    for points in (64, 128):
        model = PeriodicMap(points=points)
        x = np.arange(points)*model.step
        true = 1.0+sum(a*np.cos(n*x)+b*np.sin(n*x) for n, a, b in modes)
        # Independent Simpson integration of the exact Fourier tent
        # formula; reconstruction uses a different time grid and rule.
        f = forward_modes(1.0, modes, x, dt=0.001, cutoff=48.0)
        result, recovered = reconstruct(model, f, tolerance=1e-9)
        assert result.status == 'mean_tolerance_reached', result
        assert result.bracket[0]-1e-12 <= 1.0 <= result.bracket[1]+1e-12
        error = float(np.max(abs(recovered-true)))
        assert error < 2e-7, (points, error)
        record = dict(points=points, status=result.status, bracket=result.bracket,
                      max_density_error=error, refinements=result.refinements,
                      trial_steps=result.trial_steps)
        results.append(record)
        print(f'N={points}: bracket={result.bracket}, max density error={error:.3e}, '
              f'trial steps={result.trial_steps}')
    return results


if __name__ == '__main__':
    regression_checks()
    results = independent_data_checks()
    path = Path(__file__).resolve().parents[1] / 'scratch' / 'periodic_reconstruction_results.json'
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(results, indent=2)+'\n', encoding='utf-8')
    print(f'All checks passed. Saved {path}')

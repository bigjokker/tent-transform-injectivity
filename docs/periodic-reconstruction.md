# Appendix D. Periodic reconstruction with an unknown mean

Assume exact continuum data f=T(gstar) from a nonzero nonnegative bounded
periodic density of period P. Let cstar be its positive mean and pstar
its periodic second primitive. Define

$$
\Phi_c(p)(x)=\tfrac12\log f(x)-\tfrac12\log
\int_0^\infty e^{-ct^2-p(x+t)-p(x-t)}dt.
$$

This appendix gives a terminating search for the unknown mean cstar and a
convergent iteration for the periodic second primitive, for exact data.

## Bracket

At a maximum of pstar, Delta_t pstar<=0, giving f>=sqrt(pi/cstar)/2.
At a minimum the inequalities reverse. Hence

$$
\frac\pi{4(\max f)^2}\le cstar\le
\frac\pi{4(\min f)^2}.
$$

The bracket collapses for constant data, determining the mean immediately.

## Uniform gap and sound certificates

For each x, phi(t)=exp(-p(x+t)-p(x-t)) is positive and P-periodic. Its
Gaussian integral is

$$
\int_0^P\phi(s)w_c(s)ds,\qquad
w_c(s)=\sum_{k\ge0}e^{-c(s+kP)^2}.
$$

For c<cstar, w_c/w_cstar is continuous and strictly greater than 1 on
[0,P]. At s=0, the positive k>=1 terms make the inequality strict too.
Thus

$$
\Phi_{cstar}(p)-\Phi_c(p)\ge
\delta(c)=\tfrac12\log\min_{[0,P]}(w_c/w_{cstar})>0
$$

uniformly in p and x. The symmetric statement holds for c>cstar.
Monotonicity and constant equivariance give by induction

$$
\Phi_c^n p_0\le\Phi_{cstar}^n p_0-n\delta(c)
\quad(c<cstar).
$$

The orbit at cstar is bounded between pstar+min(p0-pstar) and
pstar+max(p0-pstar). Therefore a strict negative certificate
max(Phi_c^n p0-p0)<0 appears for n>osc(p0-pstar)/delta(c).
A strict positive certificate appears similarly for c>cstar.

Conversely, a strict negative n-step certificate implies c<cstar. If
c>=cstar, repeated application of that n-step map would drive its orbit
linearly down, contradicting comparison with the bounded true-mean orbit.
The positive certificate is analogous. Thus the signs are sound, and a
certificate appears after finitely many steps whenever c differs from cstar.

## Trials at the true mean

At c=cstar neither strict certificate can ever appear, so a plain
bisection could stall if a trial mean equals the solution.

A two-point interval search avoids this without needing an equality test.
Given a valid bracket [a,b] with a<b, use u=a+(b-a)/3 and
v=a+2(b-a)/3. Advance the n-step certificates at both trial means until
either has a strict sign. At least one trial differs from cstar, so this
takes finitely many steps. Use that sign to replace the bracket by the
appropriate subinterval:

- a negative certificate at a trial z means z<cstar: replace a by z;
- a positive certificate means z>cstar: replace b by z.

Either update reduces the interval length to at most two thirds of its
previous length. Thus any requested positive accuracy for the mean is
reached after finitely many refinements, even if one trial exactly equals
the true mean. This proves a terminating mean-search procedure for exact
data in the range. It does not require knowing pstar or the gap constant.

## Geometric convergence at the true mean

Phi_c is order-preserving and commutes with adding constants. Hence it is
nonexpansive in the oscillation seminorm:
osc(Phi_c(p)-Phi_c(q))<=osc(p-q).

At c=cstar convergence modulo constants is globally geometric. Let
R=osc(pstar)+osc(p0-pstar). Nonexpansiveness bounds the oscillation of every
iterate, and of every point of the segment between an iterate and pstar, by R.
The difference Phi(p)-Phi(q) equals a Markov kernel on the circle applied to
p-q (the mean value of the derivative along the segment). Normalizing a segment
point to have minimum zero, its time-integral denominator is at most
sqrt(pi/c)/2. For a landing point on the circle choose one positive jump of
length in [0,P]; this gives kernel density at least

$$
\alpha=\sqrt{c/\pi}\,e^{-cP^2-4R}>0 .
$$

Splitting the kernel into alpha times Lebesgue measure plus a residual
probability kernel proves the contraction factor 1-alpha*P in oscillation.
Subtracting the spatial mean at each iteration then gives uniform convergence
of the iterates to the mean-zero pstar.

## Recovering the primitive with approximate means

The geometric iteration at cstar above is in the oscillation seminorm. Approximate means can also be coupled to it: on any compact
positive interval of means, log w_c has a uniformly bounded c-derivative.
The folded-weight formula then gives

$$
\|\Phi_c(p)-\Phi_d(p)\|_\infty\le C|c-d|
$$

independently of p. Each Phi is nonexpansive in the uniform norm, so n
raw iterations from the same p0 differ by at most n*C*|c-cstar|. Find a
mean within n^{-2}, run n iterations, and subtract the final spatial mean.
The primitive converges uniformly to the mean-zero pstar, with error
bounded by the true-mean geometric error plus 2C/n.

This yields recovery of the primitive and distributional recovery of
c+p''. Uniform recovery of the density by differentiating approximate
primitives requires additional error estimates or regularization.

## Computational qualification

These are exact continuum certificates. A floating-point sign of a
quadrature approximation is not a mathematical certificate without
validated error bounds.

`code/periodic_reconstruction.py` implements the two-point search for a
positive finite quadrature model. It obtains its initial bracket by
inverting that model's constant-input integral, advances both trisection
trials, and preserves the current bracket if its budget or sign tolerance
prevents a decision. It reports `unresolved` in that case. It also reports
failure of the subsequent primitive iteration separately.

Its regression checks include a trial exactly at the true mean, constant
data, exhausted budgets, and intentionally indecisive signs. The data for
the nonconstant recovery checks use independent continuum Fourier/Simpson
quadrature, rather than the inverse's trapezoidal time grid. Errors on the
64- and 128-point examples are below 4e-11. These are floating-point
checks, not validated continuum certificates or a general accuracy bound.

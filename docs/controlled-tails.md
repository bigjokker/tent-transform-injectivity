# Appendix C. Controlled tails

This appendix proves the controlled-tail approximation quoted in Section 5 of
the [paper](../papers/injectivity-of-the-tent-transform.md), a forward L1 estimate, and an
injectivity theorem for a class defined by tail derivative bounds. Those bounds allow, for example,
tails of polynomial growth, so pairs in this class need not satisfy the
hypothesis of Theorem 3.1. The proof below is analytic; the accompanying
numerical checks (`code/controlled_tail_checks.py`) are illustrations only.

Write T(g)=f and A(g)=pi/(4f^2). All densities are nonnegative, nonzero,
and locally integrable. Equality of densities means equality a.e.

## 1. Statement

Call a half-line tail **controlled** if it is positive and C2 there and,
for some finite constants A,B>=0 and R>0,

$$
 |g'(x)|\le A\frac{g(x)}{|x|},\qquad
 |g''(x)|\le B\frac{g(x)}{x^2},\qquad |x|\ge R.
\tag{1}
$$

The constants may differ on the two half-lines and between densities.

**Tail estimate.** On each controlled tail,

$$
 \boxed{\left|g(x)-\frac{\pi}{4T(g)(x)^2}\right|
       \le\frac{C(A,B)}{x^2}\quad (|x|\ge2R).}
\tag{2}
$$

The other tail and the density on compact sets are unrestricted. Neither
finite mass nor a lower bound on x^2 g(x) is required.

**Injectivity theorem.** Let S consist of the densities for which, on each
half-line, either the mass is finite or (1) holds. Let

$$
 \widetilde S=\{g\ge0:g\not\equiv0,\quad
     g-b\in L^1(\mathbb R)\text{ for some }b\in S\}.
$$

Then T is injective on $\widetilde S$. In particular, it is injective
on S, including its infinite-mass L2 members.

This is uniqueness **within this class**. Unlike the identification of
integrable densities in Appendix A,
it does not identify an infinite-mass member among every nonnegative
locally integrable competitor.

## 2. Proof of the tail estimate

It suffices to treat the right tail. Fix x>=2R, put a=g(x)>0 and r=x/2,
and let b=2^{-A}. Integration of |(log g)'|<=A/y gives

$$
 ba\le g(y)\le b^{-1}a\quad(x/2\le y\le3x/2),
 \qquad |g''(y)|\le D a/x^2,\quad D=4B/b.
\tag{3}
$$

For 0<=t<=r, Taylor's theorem and symmetry of the tent give

$$
 ba t^2\le H(t,x)\le b^{-1}a t^2,
 \qquad |H(t,x)-at^2|\le\frac{D a t^4}{12x^2}.
\tag{4}
$$

In detail, integrate
$|g(x+u)+g(x-u)-2a|\le Da u^2/x^2$
against t-u on [0,t]. Since b<=1,

$$
 \int_0^r|e^{-H(t,x)}-e^{-at^2}|dt
 \le\frac{Da}{12x^2}\int_0^\infty t^4e^{-bat^2}dt
 =\frac{D\sqrt\pi}{32b^{5/2}a^{3/2}x^2}.
\tag{5}
$$

There is also a bound on the entire time tail, without a spatial
truncation of g. Convexity in t and (3) give

$$
 H(r,x)\ge ba x^2/4,\quad H_t(r,x)\ge ba x,
 \quad\int_r^\infty e^{-H(t,x)}dt
       \le\frac{e^{-ba x^2/4}}{ba x}.
\tag{6}
$$

The Gaussian tail is at most $e^{-a x^2/4}/(a x)$. As
$\sup_{z>0}\sqrt z e^{-bz/4}=\sqrt{2/(be)}$, (5)--(6) imply

$$
 \left|f(x)-\frac{\sqrt\pi}{2\sqrt a}\right|
 \le\frac{C_0}{a^{3/2}x^2},\qquad
 C_0=\frac{D\sqrt\pi}{32b^{5/2}}
          +(1+b^{-1})\sqrt{\frac2{be}}.
\tag{7}
$$

Converting this into (2) requires care when a is very small.

If ax^2>=1, the interval [0,1/(2sqrt(a))] lies in [0,r], so

$$
 f(x)\ge\frac{c_b}{\sqrt a},\qquad
 c_b=\int_0^{1/2}e^{-s^2/b}ds>0.
$$

Both f and sqrt(pi)/(2sqrt(a)) obey this lower bound. The derivative
of z -> pi/(4z^2) now turns (7) into

$$
 |A(g)(x)-a|\le\frac{\pi C_0}{2c_b^3x^2}.
\tag{8}
$$

If ax^2<1, (4) instead gives

$$
 f(x)\ge (x/2)e^{-1/(4b)},\qquad
 |A(g)(x)-a|\le\frac{1+\pi e^{1/(2b)}}{x^2}.
\tag{9}
$$

Taking the larger constant in (8)--(9) proves (2), including regions
where the local Gaussian scale is larger than x. Reflection proves the
left-tail assertion.

## 3. Stability of the data estimate under L1 changes

The following forward estimate allows arbitrary integrable perturbations
in the injectivity class:

$$
 \|T(g_1)^{-2}-T(g_2)^{-2}\|_{L^1}
       \le C_*\|g_1-g_2\|_{L^1}.
\tag{10}
$$

Here C_* is universal, and the individual masses need not be finite.
Only the difference on the right is assumed integrable.

First fix one density and put w(t,x)=exp(-H(t,x)). With
lambda=(log 2)/2, convexity of H and
f(x)>=2f(x) exp(-H(2f(x),x)) imply

$$
 w(t,x)\le2e^{-\lambda t/f(x)}\quad(t\ge0).
\tag{11}
$$

For d=|x-y| and s>=0, the tent kernels satisfy

$$
 k_{d+s}(z-x)-k_d(z-x)\ge k_s(z-y).
$$

For |z-x|<=d the left side is s. Otherwise it is
(d+s-|z-x|)_+, and the triangle inequality proves the assertion.
Consequently,

$$
 w(d+s,x)\le w(d,x)w(s,y).
\tag{12}
$$

The nonnegative sensitivity kernel of T(g)^{-2} is

$$
 J_g(x,y)=\frac2{f(x)^3}
       \int_d^\infty(t-d)w(t,x)dt.
$$

Write F=f(y). When d<=F/2, Lipschitz continuity gives f(x)>=F/2;
(11) therefore yields

$$
 J_g(x,y)\le 8/(\lambda^2 F).
\tag{13}
$$

When d>F/2, (11)--(12) yield

$$
 J_g(x,y)\le\frac{8F^2}{\lambda^2f(x)^3}
                  e^{-\lambda d/f(x)}
       \le\frac{216 e^{-3}F^2}{\lambda^5d^3}.
\tag{14}
$$

The last inequality uses max_{z>0} z^3 exp(-lambda z)
=27 exp(-3)/lambda^3. Integrating (13)--(14) over x proves

$$
 \sup_y\int_{\mathbb R}J_g(x,y)dx\le
 C_*:=8/\lambda^2+864e^{-3}/\lambda^5.
\tag{15}
$$

For h=g1-g2 in L1 use the nonnegative segment g_theta=g2+theta h.
Differentiation under the time integral is justified by
|K_t h|<=t||h||_1 and a common exponential tail on this segment.
It gives

$$
 \left|\partial_\theta T(g_\theta)(x)^{-2}\right|
 \le\int J_{g_\theta}(x,y)|h(y)|dy.
$$

Integrating over theta and x, using Tonelli and (15), proves (10).
This is an absolute L1 estimate for inverse squares of the data, not
an inverse stability estimate for recovering g.

## 4. Complete proof of the injectivity assertion

For b in S, (2) handles every controlled tail. On a finite-mass tail,
the mass-detection theorem of Appendix A makes T(b)^{-2} integrable
there. On compact sets, b is integrable and T(b) is positive continuous.
Thus

$$
 b-A(b)\in L^1(\mathbb R).
$$

If g-b is integrable, (10) gives A(g)-A(b) in L1 and hence

$$
 g-A(g)\in L^1(\mathbb R)\qquad(g\in\widetilde S).
\tag{16}
$$

If two such densities have the same transform, their A values agree,
and (16) gives g1-g2 in L1. An L1 difference has bounded mass on unit
intervals, so [Theorem 3.1 of the paper](../papers/injectivity-of-the-tent-transform.md)
gives g1=g2. Every step also allows infinite total mass.

## 5. Examples

The class includes positive smooth tails proportional to |x|^{-beta},
with derivative bounds (1). For 1/2<beta<=1 these are in L2 and have
infinite mass. The constants and exponents can differ at the two ends.
Logarithmic factors and positive, smooth periodic factors in log|x|
also satisfy (1) when their first two derivatives are bounded and the
factors are bounded away from zero. Integrable changes need satisfy no
tail derivative conditions.

The class also includes tails that oscillate at multiplicative spatial scales, rather
than requiring an asymptotic limit of g(x)|x|^beta.

Smoothness and L2 membership do not imply (1); for example,
(1+x^2)^(-0.35)(2+sin(x^2)) violates these derivative bounds. Such
densities are covered instead by Theorem 3.1 of the paper.

The local approximation pi/(4f^2) is not asserted to be an exact inverse.
For power tails with beta<2, the leading correction is
g''/(8g); for g=(1+x^2)^(-beta/2), this predicts
x^2(A(g)-g) -> beta(beta+1)/8. The checks in
`code/controlled_tail_checks.py` test this prediction without replacing an
infinite-mass tail by a compactly supported density.

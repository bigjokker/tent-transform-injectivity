# Appendix B. Periodic local inversion, stability, and regularity

Fix a period P. The transform extends to the open subset
U={g in real L-infinity(T_P): mean(g)>0}, allowing signed inputs.
For a nonnegative base density that vanishes somewhere, an open L-infinity
neighborhood necessarily contains signed inputs; this extension is part of
the local-diffeomorphism statement. It is not a claim that every nearby
datum comes from a nonnegative density.

## 1. Mapping properties

Write g=c+p'', with p periodic and of mean zero. The map g->(c,p) is bounded
from L-infinity to R times W2,infinity. Then

$$
H=ct^2+\Delta_tp,\quad
M_L=ct+p'(x)-p'(x-t),\quad
M_R=ct+p'(x+t)-p'(x).
$$

Translations preserve all the periodic norms in use. The W2,infinity norm
of Delta_t p is bounded independently of t. The W1,infinity norms of M_L
and M_R grow at most like 1+t. On a neighborhood with c>=c0/2>0 and bounded
norm of p, the exponential and its spatial derivatives have Gaussian
majorants times polynomials in t. W2,infinity and W1,infinity are algebras.
Consequently differentiation of the exponential series under the integral
shows T is C1 (in fact smooth) from U to W2,infinity, and I is C1 from U
to W1,infinity.

These integrals are defined pointwise, with spatial derivatives taken in
distributions and estimated by the integrable essential bounds. No
Bochner measurability of translations in the nonseparable W2,infinity
norm is assumed. The remainder estimates below apply directly to those
pointwise integrals and weak derivatives.

For completeness, for a perturbation h define N_L,N_R as the corresponding
one-sided h integrals and J=K_t h. Then

$$
DI[g]h=\int_0^\infty e^{-H}
\big(N_LM_R+M_LN_R-JM_LM_R\big)dt.
\tag{1}
$$

The W1,infinity norms of N_L,N_R are <=C(1+t)||h||_infinity, and that of
J is <=C(1+t^2)||h||_infinity. Formula (1) is therefore bounded in
W1,infinity by C||h||_infinity. The same estimates with an extra polynomial
factor dominate the quadratic remainders and prove continuity of the
Frechet derivatives in operator norm.

The identity T(g)''=2gT(g)-4I[g] is valid on all of U, including signed
inputs: their mean gives Gaussian integrability, and all the bounded
periodic derivatives justify the finite-time calculation and boundary limit.

## 2. Local inversion near every bounded nonnegative periodic density

Let g0>=0 be in L-infinity, with positive mean, and f0=T(g0). Set

$$
A=\partial_x^2-2g_0:W^{2,\infty}\to L^\infty,\qquad
B=2f_0\operatorname{Id}-4DI[g_0]:L^\infty\to L^\infty.
$$

Differentiating the identity gives A DT[g0]=B.

The map partial_x^2-1 is an isomorphism W2,infinity->L-infinity on the
circle (its periodic Green function, or the elementary periodic ODE,
provides the bounded inverse). A is its compact perturbation, because
W2,infinity embeds compactly into L-infinity. It is therefore Fredholm of
index zero. Its kernel is zero: if u''=2g0 u, periodic integration by parts
gives -integral |u'|^2=2 integral g0 u^2. Hence u is constant and positive
mean of g0 makes that constant zero. Thus A is an isomorphism.

Likewise, DI is compact on L-infinity by its W1,infinity bound and
Arzela--Ascoli. Multiplication by 2f0 is an isomorphism, since f0 is
continuous, strictly positive and periodic. Thus B is Fredholm of index zero.
The relevant abstract facts are in [Bremer's PDE notes, Section 2.3](https://www.math.toronto.edu/bremer/pdenotes.pdf#page=19).

To check injectivity of L=DT[g0], write h=a+r'', with r periodic in
W2,infinity. Then

$$
Lh(x)=-\int_0^\infty e^{-H_0(t,x)}\big(at^2+\Delta_tr(x)\big)dt.
$$

If a is nonzero, replace h by -h to make a>0 and evaluate at a minimum of
r. The integrand in parentheses is strictly positive for t>0, so Lh
cannot vanish. If a=0 and r is nonconstant, at that minimum Delta_t r is
nonnegative and strictly positive on an interval, again contradicting
Lh=0. Thus h=0. If Bh=0, the isomorphism A implies Lh=0, so B is injective
and hence an isomorphism by its index. Therefore L=A^{-1}B is an isomorphism.

The Banach inverse function theorem now proves that T is a local C1
diffeomorphism from L-infinity to W2,infinity at g0, for the extension on U.
In a sufficiently small neighborhood,

$$
\|g_1-g_2\|_\infty\le C\|T(g_1)-T(g_2)\|_{W^{2,\infty}}.
$$

## 3. Uniform stability on compact prior classes

Fix 0<alpha<1, M,m0>0, and restrict to nonnegative periodic densities
with C0,alpha norm at most M and mean at least m0. This set K is compact
in L-infinity: uniform limits preserve the mean, nonnegativity and the
Hölder bound, and equicontinuity supplies compactness.

There is C(K) such that the preceding Lipschitz estimate holds for all
pairs in K. Otherwise choose pairs with unbounded ratio. Since their
L-infinity distance is bounded, their W2,infinity data distance tends to
zero. Subsequence compactness and continuity of T give two limits with
the same transform. Theorem 3.1 of the [paper](../papers/injectivity-of-the-tent-transform.md) makes those limits
equal, since periodic densities have differences of bounded mass on unit intervals. Both
pairs eventually lie in a neighborhood where the local inverse is
Lipschitz, a contradiction. This proves the global claim on K.

It does not establish global stability without the prior bounds or in
a common norm. In particular the data norm still controls two derivatives.

## 4. Regularity equivalence

Suppose g is periodic and belongs to Cj,alpha. Then p belongs to
C(j+2),alpha and p' to C(j+1),alpha. The displayed expressions in Section
1 give

$$
\|M_L(\cdot,t)\|_{C^{j+1,\alpha}}+
\|M_R(\cdot,t)\|_{C^{j+1,\alpha}}\le C(1+t),
$$

while exp(-Delta_t p) has a uniformly bounded C(j+2),alpha norm. Therefore
the C(j+1),alpha norm of the I integrand is bounded by
C(1+t)^2 exp(-ct^2), an integrable majorant. It follows that

$$
g\in C^{j,\alpha}\Longrightarrow I\in C^{j+1,\alpha}.
\tag{2}
$$

No derivatives of g beyond j are used: differentiating the mass integrals
j+1 times produces only differences of g^{(j)} at their endpoints. For
Hölder spaces one may estimate the derivatives and difference quotients
directly, rather than assume strong continuity of translations in the
full Hölder norm.

If g is initially only L-infinity, Section 1 still gives I in
W1,infinity, hence I in C0,alpha. If f belongs to C(k+2),alpha, the
identity g=(f''+4I)/(2f) first gives a C0,alpha representative of g.
Apply (2) repeatedly to bootstrap that representative to Ck,alpha.
Conversely g in Ck,alpha makes p in C(k+2),alpha, and the defining
integral gives f in C(k+2),alpha by Gaussian bounds. Thus

$$
\boxed{f\in C^{k+2,\alpha}\iff
g\text{ has a }C^{k,\alpha}\text{ representative}.}
$$

For continuous inputs the representative qualification is unnecessary.
Smoothness of f is equivalent to smoothness of g in this periodic class.
The claim is global and periodic; it is not a local regularity theorem
for unbounded densities on the real line.

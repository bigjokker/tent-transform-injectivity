# Appendix A. Wave identity, local information, and mass detection

This appendix proves the identities and estimates used in Sections 2, 4 and 5 of the
[paper](../papers/injectivity-of-the-tent-transform.md). Densities below
are nonnegative, locally integrable, and not zero almost everywhere: their
mass on some bounded interval is positive. Equality of
measurable densities means equality almost everywhere.

## 1. Identity with the correct interpretation of boundary terms

Write w=e^{-H}, M_L=integral from x-t to x of g, and M_R=integral from x to
x+t of g. Second primitives give

$$
H_t=M_L+M_R,\quad H_x=M_R-M_L,\quad H_{tt}-H_{xx}=2g(x).
$$

The last identity holds in distributions, and is the wave equation with
zero initial value and zero initial time derivative.

The function f is C1, with

$$
\int_0^\infty wH_tdt=1,\qquad
f'=\int_0^\infty w(M_L-M_R)dt.
\tag{1}
$$

Here is a justification without any growth assumption on g at infinity.
Finite-time integrals f_T are C1, since H_x is continuous. On each compact
set of x the tail of f_T is uniformly exponentially small, using positive
mass in a fixed bounded interval. The derivative tails satisfy

$$
\int_T^\infty w|H_x|dt\le\int_T^\infty wH_tdt=e^{-H(T,x)},
$$

which also tends to zero uniformly on compact sets. Thus f_T and f_T'
converge locally uniformly and (1) follows.

Set I_T=integral from 0 to T of w M_L M_R, and
B_T(x)=e^{-H(T,x)}H_t(T,x)>=0. The exact finite-time identity is

$$
f_T''=2gf_T-4I_T-B_T.
\tag{2}
$$

For merely locally integrable g, this is a distributional identity:
chain rules on bounded x,t sets and Fubini justify it. For smooth bounded
inputs it also follows directly by integrating the second time derivative
of w. In particular the finite-time boundary term has a minus sign.

Fix a nonnegative smooth compactly supported test function phi. Equation
(2), nonnegativity of I_T and B_T, and f_T<=f imply a bound independent of T:

$$
4\int\phi I_T\le2\int\phi gf+\int|\phi''|f<\infty.
$$

Monotone convergence proves I=lim I_T belongs to L1_loc and is finite
almost everywhere. Equation (2) then shows integral of phi B_T has a
finite nonnegative limit. On the other hand, Tonelli and (1) give

$$
\int_0^\infty\int\phi(x)B_T(x)dx\,dT=\int\phi(x)dx<\infty.
$$

A nonnegative function with a finite integral and an existing limit has
limit zero. Hence B_T tends to zero in distributions, and

$$
\boxed{f''=2gf-4I,\qquad I=\int_0^\infty e^{-H}M_LM_Rdt\ge0}
\tag{3}
$$

holds in distributions and almost everywhere. Since the right-hand side
is L1_loc, f belongs to W2,1_loc and f' is locally absolutely continuous.

**Remark.** For arbitrary g in L1_loc the boundary term need
not vanish at every individual x, and I need not be finite at every x.
One can construct a smooth even density with successively very massive,
very narrow bumps at radii n: choose the nth bump's mass at least
exp(H(n,0)) and its width so small that H increases by less than 1 across
it. At its outer edge, e^{-H}H_t stays bounded below along a sequence.
The ensuing intervals of length 1/H_t each contribute a fixed positive
amount to integral of e^{-H}H_t^2, so I(0)=infinity. The weak argument above
is what makes the almost-everywhere theorem valid despite such examples.
For bounded periodic densities, Gaussian tails justify the ordinary
boundary calculation and I is Lipschitz and finite everywhere.

## 2. Direct information from f

Equations (1) give integral of w M_L=(1+f')/2 and integral of w M_R=(1-f')/2.
Since the factors are nonnegative and the weight is strictly positive:

- f'(x)=1 exactly when g vanishes almost everywhere to the right of x;
- f'(x)=-1 exactly when g vanishes almost everywhere to the left of x.

For 0<s<f(x), monotonicity of the one-sided masses and
integral from s to infinity of w >= f(x)-s give

$$
\int_x^{x+s}g\le\frac{1-f'(x)}{2(f(x)-s)},\qquad
\int_{x-s}^xg\le\frac{1+f'(x)}{2(f(x)-s)}.
\tag{4}
$$

Also f>=2f exp(-H(2f,x)), and H(2f,x)<=2f integral of g over [x-2f,x+2f].
Therefore

$$
\int_{x-2f(x)}^{x+2f(x)}g\ge\frac{\log2}{2f(x)}.
\tag{5}
$$

Under the probability measure w dt/f, both one-sided masses are
nondecreasing in t. Chebyshev's integral inequality, or its double-integral
covariance proof, gives

$$
If\ge(1-f'^2)/4.
$$

Truncating the nonnegative masses first justifies the inequality even
where I is infinite. Combining it with (3) yields, almost everywhere,

$$
\boxed{g\ge\frac{f''}{2f}+\frac{1-f'^2}{2f^2}.}
\tag{6}
$$

Where g vanishes on an interval, (3) implies f''<=0 almost everywhere,
so f is concave there.

## 3. Mass identity and detection on each half-line

All integrations by parts below are valid because f is positive, C1 and
W2,1_loc. Equation (3) gives

$$
\int_A^B g=\left[\frac{f'}{2f}\right]_A^B
+\frac12\int_A^B(f'/f)^2+2\int_A^B I/f.
\tag{7}
$$

Using (6) and |f'|<=1,

$$
\boxed{\int_A^B f^{-2}\le2\int_A^B g+1/f(A)+1/f(B).}
\tag{8}
$$

For the reverse direction, take W_x=[x-f(x)/2,x+f(x)/2]. Equation (4)
gives integral over W_x of g <=2/f(x). The Lipschitz property gives
f(y)<=3f(x)/2 there, so

$$
\int_{W_x}g\le\tfrac92\int_{W_x}f^{-2}.
$$

Select a covering of [A,B] from these centered intervals with multiplicity
at most two, using the one-dimensional Besicovitch covering lemma.
A primary reference is [Lemma 3.3.1 in Quint's course notes](https://www.math.u-bordeaux.fr/~jquint/publications/RiggaCourse.pdf#page=39).
It follows that

$$
\int_A^B g\le9\int_{A'}^{B'}f^{-2},\quad
A'=\min_{[A,B]}(x-f(x)/2),\quad B'=\max_{[A,B]}(x+f(x)/2).
\tag{9}
$$

For centers in [0,B], A'>=-f(0)/2 by Lipschitz continuity. Thus if f^{-2}
is integrable on the right, (9) proves that g is too.

Conversely, if g has finite right mass, (5) proves f(x)->infinity on the
right: either f(x)>x/4, or its window in (5) lies in [x/2,infinity), whose
mass tends to zero, forcing f(x) to infinity. Equation (8) now bounds the
right integral of f^{-2}. Reflection gives the left-hand statement.

We have proved

$$
\boxed{\int_0^\infty g<\infty\iff\int_0^\infty f^{-2}<\infty,}
$$

and the corresponding equivalence on the left. For finite total mass,
letting both endpoints escape in (8) gives integral of f^{-2} <=2 integral
of g. The coefficient 2 is sharp: a mass M concentrated at the origin has
f(x)=|x|+1/M and integral of f^{-2}=2M. Smooth approximations achieve the
same limit. Indeed the tent kernel is 1-Lipschitz in its spatial argument,
so concentrating mass in [-epsilon,epsilon] changes H from the point-mass
exponent by at most M epsilon, which supplies a dominated limit for f^{-2}.

## 4. Identification of integrable densities

If a nonnegative locally integrable density has the same transform as an
integrable density, the equivalences of Section 3 force it to be integrable
too. Their difference is then integrable, so Theorem 3.1 of the
[paper](../papers/injectivity-of-the-tent-transform.md) gives equality almost everywhere. Thus an integrable density
is identifiable among all nonnegative locally integrable densities.

More generally, equal transforms imply that on each half-line both masses
are finite or both are infinite.

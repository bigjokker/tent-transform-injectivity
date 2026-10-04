# Injectivity of an Exponential Tent Transform

**bigjokker** · Version 1.0.0 · October 2026

**Abstract.** For a nonnegative, locally integrable density $g$ on $\mathbb R$ that is not
almost everywhere zero, let

$$
T(g)(x)=\int_0^\infty\exp\Big(-\int_0^t\!\!\int_{x-\tau}^{x+\tau}g(y)\,dy\,d\tau\Big)\,dt .
$$

We prove that $T(g_1)=T(g_2)$ implies $g_1=g_2$ almost everywhere whenever
$\sup_a\int_a^{a+1}|g_1-g_2|<\infty$. In particular, $T$ is injective on nonnegative
densities in $L^p(\mathbb R)$ for every $1\le p\le\infty$, with no smoothness, decay or
finite-mass assumption. The proof combines an averaging identity along the segment between
the two densities, a nonlocal zero-flux equation, an energy kernel that dominates the
associated flux kernel, and a growth argument in the intrinsic distance $\int dx/T(g)$.
We also record an exact second-order identity $T(g)''=2g\,T(g)-4I$, information about $g$
that can be read off $T(g)$ directly, and, for periodic densities, local inversion,
Lipschitz stability on Hölder-bounded classes, a regularity equivalence, and a certified
reconstruction procedure. No closed-form inverse is given.

## 1. Introduction

For $g$ as above, write

$$
H(t,x)=\int_0^t\!\!\int_{x-\tau}^{x+\tau}g(y)\,dy\,d\tau=\int_{-t}^{t}(t-|z|)\,g(x+z)\,dz,
\qquad f(x)=T(g)(x)=\int_0^\infty e^{-H(t,x)}\,dt .
$$

So $H(\cdot,x)$ is a tent-weighted local mass of $g$ around $x$. The transform is finite
for every such $g$ (Section 2).

**Motivation.** In the one-dimensional Kolmogorov–Johnson–Mehl–Avrami (KJMA) model of
phase transformation [1–3], nuclei appear as a Poisson process in space-time with
intensity $g(y)\,dy\,ds$, and each grows at unit speed in both directions. A point $x$ is
still untransformed at time $t$ exactly when no nucleus lies in the backward cone
$\{(y,s):0\le s\le t,\ |y-x|\le t-s\}$, whose expected number of nuclei is $H(t,x)$. Hence
$e^{-H(t,x)}$ is the probability that $x$ is untransformed at time $t$, and $T(g)(x)$ is the
expected transformation time at $x$. The inverse problem asks whether the map of mean
transformation times determines the nucleation rate.

**Main results.** Theorem 3.1 states that if $T(g_1)=T(g_2)$ and
$\sup_a\int_a^{a+1}|g_1-g_2|<\infty$, then $g_1=g_2$ almost everywhere. Section 4 derives:

- $T$ is injective on nonnegative, not a.e. zero densities in $L^p(\mathbb R)$ for every
  $1\le p\le\infty$; in particular on $L^2(\mathbb R)$, including densities of infinite mass;
- equal transforms force equality whenever the difference lies in $L^p$, is bounded, or is
  locally integrable and periodic; this includes periodic densities with different periods;
- an integrable density is determined by its transform among all nonnegative locally
  integrable densities.

Appendix C adds injectivity on a class of controlled tails that need not satisfy the
hypothesis of Theorem 3.1. Sections 2 and 5 record the identity $f''=2gf-4I$, information
about $g$ that can be read off $f$, and the periodic results.

**Idea of the proof.** Interpolating $g_\theta=(1-\theta)g_2+\theta g_1$ and differentiating
$\log T(g_\theta)$ shows that the difference $Q$ of second primitives is invariant under an
average $P$ of symmetric jump kernels (Step 1). Equivalently, $u=Q'$ satisfies a nonlocal
zero-flux equation $\Lambda u=0$ (Step 2). An energy kernel $k$, obtained by differentiating
the jump distribution functions, dominates the flux kernel; this is estimate (K), proved
with Chebyshev's inequality (Step 3). It yields an energy flux $B$ with $B'=2E\ge0$ and
$|B|\le fB'$ (Step 4). A nonzero $B$ would grow exponentially in the intrinsic distance
$\int dx/f$. That forces finite mass on one side, and then contradicts the decay of the
energy there (Step 5).

**What is not claimed.** No closed-form inverse, characterization of the range, or
real-line reconstruction scheme is given. Pairs whose difference has unbounded mass on unit
intervals (outside Appendix C) and signed densities are not covered. No claim of priority
in the literature is made.

**Notation.** In addition to $H$ and $f$:

- $M_L(t,x)=\int_{x-t}^{x}g$ and $M_R(t,x)=\int_x^{x+t}g$.

## 2. Preliminaries

### 2.1 Basic properties

* **Growth for $L^2$ inputs.** If $g\in L^2$, then $f(x)\to\infty$ as $|x|\to\infty$; in
  particular $f\notin L^2$, so $g$ and $T(g)$ cannot both lie in $L^2$. For fixed $A$,
  Cauchy–Schwarz gives $H(A,x)\le C_A\|g\|_{L^2(x-A,x+A)}\to0$ as $|x|\to\infty$. So
  $f(x)\ge Ae^{-H(A,x)}$ has $\liminf_{|x|\to\infty}f\ge A$ for every $A$.

* **$f$ is 1-Lipschitz.** If $d=|x-z|$, the tent of radius $t+d$ at $z$
  dominates the tent of radius $t$ at $x$. So $H(t+d,z)\ge H(t,x)$, hence
  $f(z)=\int_0^de^{-H(s,z)}ds+\int_0^\infty e^{-H(t+d,z)}dt\le d+f(x)$.

* $f$ is finite for every nonnegative, locally integrable $g\not\equiv0$: if
  $\int_{a-1}^{a+1}g=m>0$, then $\partial_tH\ge m$ for $t\ge|x-a|+1$.

### 2.2 Second primitives, d'Alembert's formula, and an implicit inverse

Let $G''=g$. Then
$$
H(t,x)=G(x+t)+G(x-t)-2G(x).
$$
That is, $H$ is d'Alembert's solution of $H_{tt}-H_{xx}=2g(x)$ with zero initial data.
Two consequences:

1. With $W=e^{-G}$,
   $$
   2f(x)\,W(x)^2=(W*W)(2x),\qquad\text{equivalently}\qquad G(x)=\tfrac12\log f(x)-\tfrac12\log\int_0^\infty e^{-G(x+t)-G(x-t)}dt .
   $$
   Adding an affine function to $G$ changes nothing, and $g=G''$.
2. Since $H_t=M_L+M_R$ and $H_x=M_R-M_L$,
   $$
   \int_0^\infty e^{-H}(M_L+M_R)\,dt=1,\qquad f'=\int_0^\infty e^{-H}(M_L-M_R)\,dt .
   $$
   Using finite-time integrals and passing to a weak limit gives, in distributions and
   almost everywhere,
   $$
   \boxed{\,f''=2g\,f-4I,\qquad I(x)=\int_0^\infty e^{-H}M_LM_R\,dt\ \ge 0\,}
   $$
   Here $f\in C^1\cap W^{2,1}_{\mathrm{loc}}$ and $I\in L^1_{\mathrm{loc}}$. For arbitrary
   unbounded inputs, the time-boundary term need not vanish pointwise; this is why the
   weak limit is needed. In the periodic setting, $g\in C^{j,\alpha}$ implies
   $I\in C^{j+1,\alpha}$ for $0<\alpha<1$; bounded inputs also give $I\in W^{1,\infty}$.
   These gains give the regularity statement in §5
   ([Appendix B](../docs/periodic-inversion-stability-regularity.md)). Thus $g=(f''+4I)/(2f)$
   almost everywhere.

   For completeness, the finite-time identity is $f_T''=2gf_T-4I_T-B_T$, where
   $B_T=e^{-H(T,x)}H_t(T,x)\ge0$. Testing with a nonnegative compactly supported smooth
   function bounds $I_T$ locally in $L^1$. Monotone convergence gives $I$, and the
   identity gives an existing limit for each integral $\int\varphi B_T$.
   But $\int_0^\infty\int\varphi B_T\,dx\,dT=\int\varphi\,dx<\infty$, so that limit
   is zero. This justifies the weak passage to the displayed identity. Full details are in
   [Appendix A](../docs/wave-identity-and-mass-detection.md).

## 3. The uniqueness theorem

**Theorem 3.1.** Let $g_1,g_2\ge0$ be locally integrable and not a.e. zero, and put

$$
h=g_1-g_2,\qquad K=\sup_{a\in\mathbb R}\int_a^{a+1}|h|<\infty.
$$

If $T(g_1)=T(g_2)$, then $g_1=g_2$ almost everywhere.

Two elementary facts hold for any admissible $g$ with transform $f$, writing
$M_s(x)=\int_{x-s}^{x+s}g$.

**(F1)** For $t\ge s$, $H(t,x)\ge(t-s)M_s(x)$. Hence $f(x)\le s+1/M_s(x)$, interpreting
$1/M_s=\infty$ when $M_s=0$. Taking
$s=f(x)/2$ gives $M_{f(x)/2}(x)\le2/f(x)$.

**(F2)** If $\int_0^\infty f^{-2}<\infty$, then $\int_0^\infty g<\infty$ and $f(x)\to\infty$.

*Proof of (F2).* Integrate $M_{f(x)/2}(x)/f(x)\le2f(x)^{-2}$ over $x\ge0$ and use
Tonelli. For $y\ge f(0)/2$, the interval $|x-y|\le f(y)/3$ lies in $[0,\infty)$. On it,
$\frac23f(y)\le f(x)\le\frac43f(y)$, so $|x-y|\le f(x)/2$, and the interval contributes at
least $\frac12$ to $\int dx/f(x)$. Hence $\int_{f(0)/2}^\infty g\le4\int_0^\infty f^{-2}$, and $g$ is integrable on
$[0,f(0)/2]$ by local integrability. Finally $f\to\infty$: otherwise $f(x_n)\le M$ along a sequence $x_n\to\infty$ with
$x_{n+1}\ge x_n+1$, and by the Lipschitz bound each interval $[x_n,x_n+1]$ contributes at
least $(M+1)^{-2}$ to $\int_0^\infty f^{-2}<\infty$. The left half-line is the same. ∎

Now fix $g_1,g_2$ as in the theorem, with common transform $f$. If $K=0$, the conclusion
is immediate. Otherwise put $Q=G_1-G_2$ and $u=Q'$. Covering an interval by unit intervals gives

$$
|u(y)-u(x)|\le K(1+|y-x|),\qquad |\Delta_tQ(x)|\le Kt(1+2t).
$$

Thus $u$ and $Q$ grow at most polynomially.

**Step 1 (averaging along the segment).** For $\theta\in[0,1]$ put
$g_\theta=(1-\theta)g_2+\theta g_1$, $w_\theta=e^{-H_\theta}$ and
$f_\theta=\int_0^\infty w_\theta\,dt$. Let $T_\theta$ have density
$w_\theta(t,x)/f_\theta(x)$ on $(0,\infty)$. Let $P_\theta$ be the symmetric jump
$x\mapsto x\pm T_\theta$, i.e. $p_\theta(x,y)=w_\theta(|y-x|,x)/(2f_\theta(x))$. Then
$$
\partial_\theta\log f_\theta(x)=-\mathbb E\big[\Delta_{T_\theta}Q(x)\big]=2\big(Q-P_\theta Q\big)(x).
$$
Differentiation is justified by the polynomial tent bound and a common exponential time
tail: one bounded interval carries positive mass for each endpoint and hence for every
convex combination. Integrate over $\theta\in[0,1]$ and use $f_1=f_0$:
$$
Q=PQ,\qquad P=\int_0^1P_\theta\,d\theta .
$$
Also $f/(2e)\le f_\theta\le f$. The upper bound is Hölder. For the lower bound, $H(t,x)\le t\,M_t(x)$ and (F1) give
$H_j(f/2,x)\le\frac f2\cdot\frac2f=1$ for $j=1,2$ (both have transform $f$). Since
$H_\theta\le\max(H_1,H_2)$ is nondecreasing in $t$, $f_\theta\ge\int_0^{f/2}e^{-1}dt=f/(2e)$.

**Step 2 (zero flux).** Let $A_\theta(r,x)=\mathbb P(T_\theta>r)$ and
$$
\Lambda v(x)=\int_0^1\!\!\int_0^\infty A_\theta(s,x)\,\big(v(x+s)-v(x-s)\big)\,ds\,d\theta .
$$
Since $\Delta_tQ(x)=\int_0^t\big(u(x+s)-u(x-s)\big)ds$, Fubini gives $\Lambda u=2(PQ-Q)$.
So for **every** $x$,
$$
\Lambda u(x)=0,\qquad \Lambda 1=0. \tag{Z}
$$
Differentiating $Q=PQ$ would only show that $\Lambda u$ is constant. The vanishing of that
constant is what the proof needs.

**Step 3 (an energy kernel).** Let $F_\theta(x,y)=P_\theta\big(x,(-\infty,y]\big)$ and
$k_\theta=-\partial_xF_\theta$. Let $k=\int_0^1k_\theta\,d\theta$ and
$\ell(x,y)=\operatorname{sgn}(y-x)\int_0^1A_\theta(|y-x|,x)\,d\theta$, so that
$\Lambda v=\int\ell(\cdot,y)v(y)\,dy$. Off the diagonal
$F_\theta=\mathbf 1_{x<y}-\frac12\ell_\theta$. $F_\theta(\cdot,y)$ is $C^1$ for $x\ne y$, since
$\partial_xH_\theta=M_R-M_L$ is continuous, and continuous across $x=y$ (both one-sided
limits are $\frac12$). So $\partial_x\ell=2k$ off the diagonal, while $\ell(\cdot,y)$ jumps
from $+1$ to $-1$ at $x=y$. Differentiating $\Lambda v(x)=\int\ell(x,y)v(y)\,dy$, the jump
contributes $-2v(x)$. This gives the row identity below, classically, for every
continuous $v$ of polynomial growth. Hence:

* rows: $\int k(x,y)\,dy=1$, and $\int k(x,y)v(y)\,dy-v(x)=\tfrac12(\Lambda v)'(x)$;
* columns: $\int k(x,y)\varphi(x)\,dx=\varphi(y)-\tfrac12\int\ell(x,y)\varphi'(x)\,dx$ for
  $\varphi\in C^1_c$.

The kernel can also be written, off the diagonal, as

$$
k_\theta(x,y)=p_\theta(x,y)
+\tfrac12\operatorname{sgn}(y-x)
\left.\partial_x A_\theta(r,x)\right|_{r=|y-x|},
$$

where the derivative holds $r$ fixed. The derivative-tail bound
$\int_T^\infty w_\theta|H_{\theta,x}|\,dt\le w_\theta(T,x)$ also
justifies spatial differentiation locally uniformly. The key estimate is
$$
\frac{A_\theta(|y-x|,x)}{2f_\theta(x)}\ \le\ k_\theta(x,y)\ \le\ 3\,p_\theta(x,y).\tag{K}
$$

*Lower bound.* Take $y=x+r$. Since $F_\theta=1-\frac12A_\theta(y-x,x)$ there,
$k_\theta=p_\theta+\frac12\partial_xA_\theta(r,x)$. With $r$ fixed,
$f_\theta\partial_xA_\theta=-\int_r^\infty w_\theta\,\partial_xH_\theta\,dt-A_\theta f_\theta'$.
Also $w_\theta(r)=\int_r^\infty w_\theta\,\partial_tH_\theta\,dt$ and
$\partial_tH_\theta-\partial_xH_\theta=2M_L$. Hence
$2f_\theta k_\theta=2\int_r^\infty w_\theta M_L\,dt-A_\theta f_\theta'$. Under the law of
$T_\theta$, both $M_L(t)$ and $\mathbf 1_{t>r}$ are nondecreasing in $t$, so by
Chebyshev's inequality $\mathbb E[M_L\mathbf 1_{T_\theta>r}]\ge\mathbb E[M_L]\,\mathbb P(T_\theta>r)$. Combined with
$f_\theta'=f_\theta\,\mathbb E[M_L-M_R]$ and $\mathbb E[M_L+M_R]=1/f_\theta$, this gives
$2f_\theta k_\theta\ge A_\theta$. For $y<x$, use $M_R$.

*Upper bound.* $H_\theta$ is convex in $t$ with $H_\theta(0)=0$, hence superadditive, hence
$A_\theta(r)\le w_\theta(r)$. Since $|\partial_xH_\theta|\le\partial_tH_\theta$ and
$|f_\theta'|\le1$ ($f_\theta$ is itself a transform),
$|f_\theta\partial_xA_\theta|\le w_\theta(r)+A_\theta\le2w_\theta(r)$. So
$|\partial_xA_\theta|\le4p_\theta$ and $k_\theta\le3p_\theta$.

Consequences of (K):

- $k>0$ everywhere;
- $|\ell|\le2f\,k$;
- the two values $k_\theta(x,x\pm r)$ add up to $w_\theta(r)/f_\theta$, so
  $\int|y-x|\,k_\theta\,dy=\mathbb E\,T_\theta\le f_\theta\le f$.
  Also, since $A_\theta\le w_\theta$,

  $$
  \mathbb E T_\theta^2=2\int_0^\infty rA_\theta(r,x)\,dr
  \le2f_\theta\mathbb E T_\theta\le2f_\theta^2.
  $$

**Step 4 (the energy flux).** Set
$$
E(x)=\int k(x,y)\big(u(y)-u(x)\big)^2dy,\qquad B(x)=\Lambda(u^2)(x),
$$
and $Sv=\int k(\cdot,y)v(y)\,dy$. By the row identity and (Z), $Su-u=\frac12(\Lambda u)'=0$.
Apply the row identity to $v=u^2$: $B$ is $C^1$ and $B'=2(S(u^2)-u^2)$. Since $Su=u$,

$$
S(u^2)(x)-u(x)^2=\int k(x,y)\big(u(y)-u(x)\big)^2dy+2u(x)\,(Su-u)(x)=E(x),
$$

so $B'=2E$. Polynomial growth and locally uniform exponential tails justify
all integrals. Also $B=\Lambda((u-u(x))^2)$ pointwise by (Z). Hence

$$
B'=2E\ge0,\qquad |B|\le2fE=fB',\qquad
E\le K^2(1+2f+2f^2).\tag{B}
$$

The energy bound uses $|u(y)-u(x)|\le K(1+|y-x|)$ and the two moment bounds above.
In particular $B$ is nondecreasing. Also, by $|\ell|\le2fk$,
$|B(x)|\le\int|\ell(x,y)|(u(y)-u(x))^2dy\le2f(x)E(x)$.

**Step 5 (a growth contradiction).** Let $D(x)=\int_0^xdz/f(z)$. Since
$f(x)\le f(0)+|x|$, we have $D(x)\ge\log\big(1+x/f(0)\big)\to\infty$ on the right.

Suppose $B(x_0)>0$. Since $B$ is nondecreasing we may take $x_0\ge0$. By (B),
$B'=2E\ge|B|/f=B/f$, i.e. $(\log B)'\ge1/f$ on $[x_0,\infty)$, so
$$
B(x)\ \ge\ B(x_0)\,e^{D(x)-D(x_0)}\ \ge\ c\,x .
$$
Since also $B\le2K^2f(1+2f+2f^2)$, this forces $f\to\infty$. Once $f\ge1$, $B\le10K^2f^3$, so
$f\ge c'e^{D/3}$ for sufficiently large $x$. Hence
$$
\int_{x_0}^\infty f^{-2}dx=\int^{\infty}\frac{dD}{f}<\infty .
$$
Applying (F2) to $g_1$ and to $g_2$ separately, both have finite mass on $(0,\infty)$, and
$f\to\infty$ there.

But then $E(x)\to0$. Let $\tau(r)=\int_r^\infty(g_1+g_2)$; since $|h|\le g_1+g_2$,
$|u(a)-u(b)|\le\tau(r)$ for $a,b\ge r$. Since $w_\theta\le1$ and $f_\theta\ge f/(2e)$, we have $p_\theta\le\frac1{2f_\theta}\le\frac ef$,
so (K) gives $k(x,y)\le3e/f(x)$.

* **$y\ge0$:** fix $R>0$ and let $x\ge R$. Split at $y=R$ and use that rows of $k$ have
  mass $1$:
  $$
  \int_0^\infty k(x,y)\big(u(y)-u(x)\big)^2dy\le\frac{3eR\,\tau(0)^2}{f(x)}+\tau(R)^2 .
  $$
  Let $x\to\infty$, then $R\to\infty$.
* **$y=-s<0$:** choose $A,m>0$ with $\int_{-A}^Ag_j\ge m$ for $j=1,2$, hence for every
  $g_\theta$. For $x\ge0$ and $s\ge A$, the window $[x-t,x+t]$ at $t=x+s$ contains $[-A,A]$,
  so $w_\theta(x+s,x)\le e^{-m(s-A)}$. For $0\le s\le A$, use $w_\theta\le1$. With
  $\rho(s)=\min\big(1,e^{-m(s-A)}\big)$, this gives $k(x,-s)\le\frac{3e}{f(x)}\rho(s)$,
  uniformly in $\theta$ and $x\ge0$. Since $|u(-s)-u(x)|\le K(1+s)+\tau(0)$,
  $$
  \int_{-\infty}^0k(x,y)\big(u(y)-u(x)\big)^2dy\le\frac{6e}{f(x)}\int_0^\infty\rho(s)\big(K^2(1+s)^2+\tau(0)^2\big)\,ds\longrightarrow0 .
  $$

So $B\le2fE=o(x)$ (as $f(x)\le f(0)+x$), contradicting $B\ge cx$. Hence $B\le0$ everywhere.

Reflecting both densities ($x\mapsto-x$) gives a pair with equal transforms and the same
$K$. Its $u$ is $-u(-x)$ and its flux is $-B(-x)$, so likewise $B\ge0$. Thus $B\equiv0$, $E=\frac12B'\equiv0$, and since $k>0$,
$u$ is constant. Hence $h=u'=0$. ∎

*Remarks.*

- Equality of the transforms enters in $f_1=f_0$ in Step 1, the
  common $f$ in $f/(2e)\le f_\theta\le f$, and applying (F2) to each density.

- Nonnegativity ensures convexity and superadditivity of $H_\theta$ in $t$, monotonicity
  of $M_L,M_R$, and $|H_{\theta,x}|\le H_{\theta,t}$. Signed densities are outside this proof.

- The proof is a Liouville theorem for the zero-flux equation (Z). $D=\int dx/f$ measures
  distance in units of the local jump length $f(x)$.

## 4. Consequences

For any $1\le p\le\infty$, Hölder's inequality on a unit interval gives
$K\le\|g_1-g_2\|_p$. Thus the theorem covers every $L^p$ difference, and in particular
all nonnegative nonzero $L^p$ inputs. A locally integrable periodic function also has
uniformly bounded absolute mass on unit intervals. Consequently periodic differences
are covered, and so are pairs of periodic densities even with different periods.

An integrable density is identifiable among **all nonnegative locally integrable**
competitors. Divide $f''=2gf-4I$ by $2f$ and integrate by parts, using
$f\in W^{2,1}_{\rm loc}$:
$\int_A^Bg=\big[\frac{f'}{2f}\big]_A^B+\frac12\int_A^B\frac{f'^2}{f^2}+2\int_A^B\frac If$.
With $2I/f\ge(1-f'^2)/(2f^2)$ (Chebyshev, §5) and $|f'|\le1$, this gives

$$
\int_A^Bf^{-2}\le2\int_A^Bg+\frac1{f(A)}+\frac1{f(B)}.
$$

For an integrable input, $f\to\infty$ at both ends. Otherwise, by the lower mass bound
$\int_{x-2f}^{x+2f}g\ge\log2/(2f)$ of §5, infinitely many disjoint bounded windows would
each carry mass bounded below. Hence $f^{-2}$ is integrable. Applying (F2) to a competitor with the same data makes
it integrable too. Its difference from the original input is in $L^1$, so the theorem
applies.

[Appendix C](../docs/controlled-tails.md) proves a further injectivity theorem for densities whose tails satisfy
the derivative bounds $|g'|\le Ag/|x|$ and $|g''|\le Bg/x^2$ (or have finite mass), together
with integrable perturbations of them. Such pairs, for example with polynomially growing
tails, need not satisfy the unit-interval hypothesis of Theorem 3.1.

## 5. Information about $g$ contained in $f$, and stability

Proofs of the first three items are in [Appendix A](../docs/wave-identity-and-mass-detection.md).

* **Support:** $f'(x)=1$ iff $g=0$ a.e. on $(x,\infty)$, and $f'(x)=-1$ iff $g=0$ a.e. on
  $(-\infty,x)$.

* **Local masses:** for $0<s<f(x)$,
  $$
  \int_x^{x+s}g\le\frac{1-f'(x)}{2(f(x)-s)},\qquad\int_{x-s}^{x}g\le\frac{1+f'(x)}{2(f(x)-s)},\qquad\int_{x-2f(x)}^{x+2f(x)}g\ge\frac{\log2}{2f(x)}.
  $$

* **Pointwise:** a.e., $\displaystyle g\ \ge\ \frac{f''}{2f}+\frac{1-f'^2}{2f^2}$. This is
  Chebyshev again: $I\,f\ge\frac{(1+f')(1-f')}{4}$. In particular, $f$ is concave on any
  interval where $g=0$ a.e.

* **Controlled-tail approximation:** on positive $C^2$ tails satisfying
  $|g'|\le Ag/|x|$ and $|g''|\le Bg/x^2$,
  $g=\pi/(4f^2)+O(x^{-2})$. The formula is exact for constants. It is not an asserted
  approximation for every unrestricted input ([Appendix C](../docs/controlled-tails.md)).

* **Regularity and stability** (for bounded nonnegative nonzero densities with a fixed period;
  [Appendix B](../docs/periodic-inversion-stability-regularity.md)):

  - For $k\ge0$ and $0<\alpha<1$, $f\in C^{k+2,\alpha}$ iff $g$ has a $C^{k,\alpha}$
    representative.
  - $T$ is a local $C^1$-diffeomorphism $L^\infty\to W^{2,\infty}$ near every such density,
    using its extension to signed periodic inputs of positive mean. This does not mean
    that every nearby datum comes from a nonnegative input.
  - This gives $\|g_1-g_2\|_\infty\le C\,\|Tg_1-Tg_2\|_{W^{2,\infty}}$ on classes bounded in
    $C^\alpha$ with mean bounded below.
  - For period $2\pi$ at a constant $c>0$, the linearization has Fourier multipliers
    $m_n=-\frac{\sqrt\pi}{\sqrt c}\frac{1-e^{-n^2/(4c)}}{n^2}$ for $n\ne0$, and
    $m_0=-\frac{\sqrt\pi}{4c^{3/2}}$. They are nonzero, and $|m_n|$ decays like $n^{-2}$.

* **Reconstruction** (periodic, exact data): the unknown mean is located by a terminating
  certified search, and the periodic second primitive is recovered by an iteration that
  converges geometrically at the true mean
  ([Appendix D](../docs/periodic-reconstruction.md)). Uniform recovery of the density itself
  from approximate primitives needs additional regularization.

The regularity, stability and reconstruction statements are periodic. The multiplier decay
rules out an inverse Lipschitz estimate in the same Sobolev norm near a constant; the
displayed stability estimate uses two additional derivatives of the data. General
real-line reconstruction convergence has not been established here.

## 6. Scope and open questions

* Uniqueness is relative to the competitor class. A nonnegative nonzero locally integrable
  input is determined among nonnegative nonzero locally integrable inputs whose difference
  from it has uniformly bounded absolute mass on unit intervals; for example, any two $L^2$
  inputs. An integrable input is determined among all nonnegative locally integrable
  inputs (§4).

* For an $L^2$ input, an even transform implies an even density: reflection commutes
  with $T$, and injectivity applies to the reflected pair.

* The proof uses nonnegativity of the densities (convexity of $H$ in $t$, monotonicity of
  the one-sided masses, and $|H_x|\le H_t$). Signed densities are not covered.

* Open questions:
  - uniqueness for pairs whose difference has unbounded mass on unit intervals, beyond the
    controlled-tail class of Appendix C (for example, densities with growing local mass);
  - a closed-form or otherwise explicit inverse, and a characterization of the range of $T$;
  - a convergent reconstruction scheme and stability estimates on the whole real line.

## Appendices

- A. [Wave identity, local information, and mass detection](../docs/wave-identity-and-mass-detection.md)
- B. [Periodic local inversion, stability, and regularity](../docs/periodic-inversion-stability-regularity.md)
- C. [Controlled tails](../docs/controlled-tails.md)
- D. [Periodic reconstruction with an unknown mean](../docs/periodic-reconstruction.md)

## References

1. A. N. Kolmogorov, On the statistical theory of the crystallization of metals,
   *Izv. Akad. Nauk SSSR Ser. Mat.* 1 (1937), no. 3, 355–359.
2. W. A. Johnson and R. F. Mehl, Reaction kinetics in processes of nucleation and growth,
   *Trans. AIME* 135 (1939), 416–442.
3. M. Avrami, Kinetics of phase change. I. General theory, *J. Chem. Phys.* 7 (1939),
   1103–1112.

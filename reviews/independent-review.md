# Independent adversarial review of the uniqueness proof

Reviewed on 2026-10-04.

> **Editorial note.** This review was carried out on a self-contained review copy of the
> proof (called "the packet" below), containing the material now in Sections 2–6 of the
> [paper](../papers/injectivity-of-the-tent-transform.md). Line numbers refer to that copy,
> which is not included. The review's three recommended changes are incorporated in the
> paper: the explicit finite-side energy estimate (Sections 7.3–7.4 below), the explicit
> use of the integrability of $f^{-2}$ when proving $f\to\infty$ in (F2), and the precise
> statement of the competitor class. Apart from this note and the title, the review is
> reproduced unchanged.

## Verdict and strongest objection

**Verdict: no defect found in the core uniqueness proof after independently
checking its central identities and supplying the estimates detailed below.**
I did not find a counterexample, a fatal implication, or a missing hypothesis
needed for the theorem as stated.

The strongest candidate gap was the negative-landing estimate in Step 5,
packet lines 295–299:

> “the probability of landing in $y<0$, and its first two overshoot moments,
> are $O(1/f(x))\to0$.”

An exponential bound starting only at $t=x+A$ does not, by itself, cover
the strip $-A<y<0$. Also, “overshoot moments” must mean the **unconditioned,
truncated moments**
$$
 \int_{-\infty}^0 |y|^j p_\theta(x,y)\,dy,\qquad j=0,1,2,
$$
not moments conditional on landing on the negative half-line.

Both issues are repairable without strengthening the hypotheses. The strip
has fixed length and density at most $e/f(x)$; beyond it there is a
uniform exponential envelope. Section 7 below gives explicit constants
and proves the required energy limit. Conditional overshoot moments need
not tend to zero: an exact one-sided example in Section 9 demonstrates this.
The proof only needs the unconditioned version, which is correct.

Other points worth making explicit:

- The row identity in Step 3 does hold classically for continuous functions
  of polynomial growth. Its justification needs both the diagonal jump
  and the kernel domination, not merely the assertion that $H_x$ is
  continuous. I give that justification below.
- The last sentence of (F2), “$f\to\infty$ follows from the Lipschitz
  bound,” also uses the premise $\int_0^\infty f^{-2}<\infty$.
  Lipschitz continuity alone would not suffice.
- The theorem is a statement about a **pair whose difference has uniformly
  bounded absolute mass on unit intervals**. It does not prove uniqueness
  among all locally integrable competitors merely because one input is
  in $L^2$. It does prove injectivity on the whole nonnegative $L^2$
  class, including its infinite-mass members.

These are explanatory and scope issues, not demonstrated counterexamples
to the stated core theorem. After filling in the displayed estimates, I
have no unresolved analytic step in that core argument.

I used the packet's assumptions and re-derived the mathematics. Previous
audit verdicts and stress-test outputs were not used as evidence. This
review is not formal verification, and it makes no claim about originality.

## 1. Basic regularity and exactly which bounds are uniform

Write
$$
 H_g(t,x)=\int_{-t}^t(t-|z|)g(x+z)\,dz,\qquad w_g=e^{-H_g}.
$$
For nonnegative $g\in L^1_{\rm loc}$, choose a second primitive
$G\in C^1\cap W^{2,1}_{\rm loc}$. Then
$$
 H_g(t,x)=G(x+t)+G(x-t)-2G(x),
$$
and
$$
 H_t=M_L+M_R,\qquad H_x=M_R-M_L,\qquad |H_x|\le H_t.
$$
Both first derivatives are continuous; no pointwise boundedness of $g$
is required.

For the segment $g_\theta=(1-\theta)g_2+\theta g_1$, choose $A>0$ and
$m>0$ such that
$$
 \int_{-A}^A g_j\ge m\qquad(j=1,2).
$$
This is possible even when the supports of the two densities are disjoint.
For every $\theta\in[0,1]$,
$$
 H_{\theta,t}(t,x)\ge m\quad\text{if }t\ge |x|+A,
$$
and consequently
$$
 w_\theta(t,x)\le e^{-m(t-|x|-A)}
 \quad(t\ge |x|+A).
\tag{1}
$$
The constants $A,m$ are independent of $\theta$. On a fixed compact
set of $x$'s, (1) is a common exponential envelope.

The envelope is **not** uniform over all $x$ with a fixed time origin.
The proof must not use it that way. In the finite-side estimate it is
instead applied at $t=x+s$, which produces an envelope in the overshoot
$s$ independent of $x\ge0$. That use is legitimate.

For spatial differentiation, the essential estimate is
$$
 \int_T^\infty w_\theta(t,x)|H_{\theta,x}(t,x)|\,dt
 \le \int_T^\infty w_\theta H_{\theta,t}\,dt
 =w_\theta(T,x).
\tag{2}
$$
Finite-time integrals are $C^1$. Equations (1)–(2) imply locally uniform
convergence of the integrals and their first derivatives, uniformly in
$\theta$. Therefore
$$
 f_\theta\in C^1,\qquad
 f_\theta'=\int_0^\infty w_\theta(M_L-M_R)\,dt,
 \qquad |f_\theta'|\le1.
\tag{3}
$$
Also $\int_0^\infty w_\theta H_{\theta,t}\,dt=1$.

This argument survives arbitrarily high, narrow spikes. Trying instead
to bound $H_x$ pointwise by a fixed polynomial would be unjustified
under the stated assumptions.

The unit-interval bound on $h=g_1-g_2$ gives
$$
 |u(y)-u(x)|\le K(1+|y-x|),\qquad u=Q',\quad Q''=h.
\tag{4}
$$
In particular $u=O(1+|x|)$ and $Q=O(1+x^2)$, after fixing their affine
normalizations. A slightly sharper bound than the packet needs is
$$
 |\Delta_tQ(x)|
 \le \int_0^t K(1+2s)\,ds=K(t+t^2).
\tag{5}
$$
Thus its weaker bound $Kt(1+2t)$ is valid.

For later reference:

| Estimate | Uniformity actually available |
|---|---|
| Common exponential time tail | Uniform in $\theta$, locally uniform in $x$ |
| Polynomial bounds for $u$ and $Q$ | Global in space; constants depend on $K$ and normalization |
| $f/(2e)\le f_\theta\le f$ | All $x,\theta$, under equal endpoint data |
| Radial moments of $k_\theta$ | Controlled by $f(x)$ and $f(x)^2$, not by absolute constants |
| Negative-landing moments | Uniform in $\theta,x\ge0$ after multiplication by $f(x)$ |

No later step needs a stronger uniformity statement.

## 2. Independent check of (F1), (F2), and absence of circularity

For $t\ge s$,
$$
 H(t,x)\ge(t-s)\int_{x-s}^{x+s}g.
$$
Writing the mass on the right as $M_s(x)$, integration gives
$$
 f(x)\le s+\frac1{M_s(x)}.
$$
If $M_s=0$ this carries no information, as acknowledged in the packet.
Taking $s=f(x)/2$ otherwise gives
$$
 M_{f(x)/2}(x)\le\frac2{f(x)}.
\tag{6}
$$
The direction of this inequality is correct.

Now suppose $J=\int_0^\infty f(x)^{-2}\,dx<\infty$. Tonelli and (6) give
$$
 \int_{\mathbb R}g(y)
 \left(\int_{\substack{x\ge0\\|x-y|\le f(x)/2}}\frac{dx}{f(x)}\right)dy
 \le2J.
\tag{7}
$$
For $y\ge f(0)/2$, define
$$
 I_y=[y-f(y)/3,\ y+f(y)/3].
$$
Since $f(y)\le f(0)+y$, the left endpoint is nonnegative. On $I_y$,
the Lipschitz bound implies
$$
 \frac23 f(y)\le f(x)\le\frac43 f(y),\qquad
 |x-y|\le f(x)/2.
$$
The inner integral in (7) is therefore at least
$$
 |I_y|\frac3{4f(y)}=\frac12.
$$
Consequently
$$
 \int_{f(0)/2}^\infty g\le4J,
$$
and the omitted compact interval has finite mass by local integrability.
This proves the finite-mass conclusion of (F2).

The assertion $f(x)\to\infty$ follows from $J<\infty$ and Lipschitz
continuity: if $f(x_n)\le M$ along $x_n\to\infty$, choose disjoint
intervals $[x_n,x_n+1]$. On each, $f\le M+1$, giving a contribution
at least $(M+1)^{-2}$ to $J$, a contradiction.

This proof of (F2) uses neither the wave identity nor the uniqueness
theorem. Its use in Step 5 is not circular.

## 3. The segment equation gives zero flux, not an unspecified constant

By (1) and (5), differentiation under the time integral is justified
uniformly in $\theta$, for each fixed $x$, and locally uniformly in
$x$:
$$
 \partial_\theta\log f_\theta(x)
 =-\frac{\int_0^\infty w_\theta(t,x)\Delta_tQ(x)\,dt}{f_\theta(x)}
 =2\bigl(Q-P_\theta Q\bigr)(x).
$$
The endpoint equality $f_1=f_0$ therefore implies
$$
 Q=PQ,\qquad P=\int_0^1P_\theta\,d\theta.
\tag{8}
$$
This is the mixture of the individually normalized kernels. Replacing it
by an independently normalized average of their unnormalized weights
would be a different operator. The packet makes the correct choice.

Hölder gives $f_\theta\le f$. For the opposite bound, (6) applied to
each endpoint density gives
$$
 H_j(f/2,x)\le1.
$$
The same is true of $H_\theta$, and monotonicity in time gives
$$
 f_\theta(x)\ge\int_0^{f(x)/2}e^{-1}\,dt=\frac{f(x)}{2e}.
\tag{9}
$$
Thus the claimed comparison is valid without a positive lower bound on
either density.

Let
$$
 A_\theta(r,x)=\frac1{f_\theta(x)}\int_r^\infty w_\theta(t,x)\,dt.
$$
Since
$$
 \Delta_tQ(x)=\int_0^t[u(x+s)-u(x-s)]\,ds,
$$
absolute Fubini, using (4) and the exponential time tail, yields
$$
 \Lambda u(x)=\int_0^1\mathbb E_\theta[\Delta_TQ(x)]\,d\theta
            =2(PQ-Q)(x)=0.
\tag{10}
$$
Also $\Lambda1=0$ by symmetry. This proves the actual zero-flux equation
at every $x$.

This step is indispensable. In a constant background $g=c$, take the
test function $u(x)=x$, without claiming it comes from equal data.
Then $k=p$, $Su=u$, but
$$
 \Lambda u=\mathbb E T^2=\frac1{2c}\ne0.
$$
Moreover
$$
 B(x)=\Lambda(u^2)(x)=\frac{x}{c},\qquad
 E(x)=\frac1{2c}.
$$
The desired bound $|B|\le fB'$ fails for large $|x|$. Thus a proof
using only $Su=u$ really would be wrong. The packet avoids this error
by deriving (10) before using the energy argument.

## 4. Independent derivation of the kernel and its comparisons

In this section suppress $\theta$ temporarily.

### 4.1 The diagonal jump and its sign

For $y>x$,
$$
 F(x,y)=1-\frac12A(y-x,x),
$$
and for $y<x$,
$$
 F(x,y)=\frac12A(x-y,x).
$$
Define
$$
 \ell(x,y)=\operatorname{sgn}(y-x)A(|y-x|,x).
$$
Then, off the diagonal,
$$
 F(x,y)=\mathbf1_{\{x<y\}}-\frac12\ell(x,y).
$$
As $x$ passes through fixed $y$, $\ell$ jumps from $+1$ to $-1$.
Thus its distributional derivative is
$$
 \partial_x\ell=2k-2\delta_{x=y},\qquad k=-\partial_xF.
\tag{11}
$$
The delta terms cancel in $\partial_xF$, because $F$ itself is
continuous at the diagonal. Equivalently, for $r=|y-x|>0$,
$$
 k(x,y)=p(x,y)+\frac12\operatorname{sgn}(y-x)A_x(r,x),
\tag{12}
$$
where $A_x$ holds $r$ fixed. The sign in (12), and the negative
diagonal term in (11), agree with the packet.

### 4.2 The upper comparison

Convexity of $H(\cdot,x)$, with $H(0,x)=0$, gives
$$
 H(r+s,x)\ge H(r,x)+H(s,x).
$$
Hence
$$
 A(r,x)\le w(r,x).
\tag{13}
$$
At fixed $r$, (2)–(3) justify
$$
 fA_x=-\int_r^\infty wH_x\,dt-Af'.
\tag{14}
$$
Using $|H_x|\le H_t$, $|f'|\le1$, and (13),
$$
 |fA_x|\le w(r)+A(r)\le2w(r).
$$
Substitution into (12) gives
$$
 k(x,y)\le3p(x,y).
\tag{15}
$$
This upper bound does not assume that $k$ is positive; positivity is
established separately next.

### 4.3 The lower comparison

For $y=x+r$, equations (12) and (14) give
$$
 2fk(x,x+r)
 =w(r)-\int_r^\infty wH_x\,dt-A(r)f'
 =2\int_r^\infty wM_L\,dt-A(r)f'.
$$
Both $M_L(t)$ and $\mathbf1_{\{t>r\}}$ are nondecreasing functions of
$t$. Their covariance under the probability law $w(t)\,dt/f$ is
nonnegative. The necessary first moment exists because
$\int wM_L\le\int wH_t=1$. Therefore
$$
 \int_r^\infty wM_L\,dt
 \ge A(r)\int_0^\infty wM_L\,dt
 =A(r)\frac{1+f'}2,
$$
so
$$
 2fk(x,x+r)\ge A(r).
$$
For $y=x-r$, the corresponding formula is
$$
 2fk(x,x-r)=2\int_r^\infty wM_R\,dt+A(r)f'\ge A(r).
$$
Thus
$$
 \boxed{\frac{A(|y-x|,x)}{2f(x)}
        \le k(x,y)\le3p(x,y).}
\tag{16}
$$
In particular $k>0$ at every off-diagonal point.

There is no singularity that affects the argument at $y=x$.
Indeed $A(0,x)=1$, $A_x(0,x)=0$, and the continuous extension is
$k(x,x)=1/(2f(x))>0$.

### 4.4 Row identities, regularity, and the mixture

Equation (11) gives
$$
 (\Lambda v)'=2Sv-2v,\qquad Sv(x)=\int k(x,y)v(y)\,dy.
\tag{17}
$$
Here is a sufficient justification for the claim that this is classical.
The fixed-$r$ derivative $A_x$ is continuous by (2), (14), and locally
uniform convergence of the tail integrals. Its limit at $r=0$ is zero,
so $k$ has the continuous extension just described. For $x$ in a
compact set, (15) and (1) dominate $k(x,y)v(y)$ by an integrable function
of $y$ whenever $v$ is continuous with polynomial growth.
The same holds for $\ell v$, using (13).
Consequently $Sv$ and $\Lambda v$ are continuous; the distributional
identity (17) has a continuous right side and implies $\Lambda v\in C^1$.
All these dominations are uniform in $\theta$.

This argument does not differentiate $v$, which matters because $u$
is only locally absolutely continuous.

Directly from (12),
$$
 k_\theta(x,x+r)+k_\theta(x,x-r)
 =\frac{w_\theta(r,x)}{f_\theta(x)}.
\tag{18}
$$
Integrating shows $\int k_\theta(x,y)\,dy=1$. For a compactly supported
$C^1$ test function $\varphi$, (11) also gives exactly
$$
 \int k(x,y)\varphi(x)\,dx
 =\varphi(y)-\frac12\int\ell(x,y)\varphi'(x)\,dx.
$$
The column formula has the correct sign. The proof does not need to
upgrade it to any unproved statement about an infinite invariant measure.

Averaging preserves positivity, row normalization, and (17). Since all
$\ell_\theta(x,y)$ have the same sign for fixed $x,y$, (9), (16), and
$f_\theta\le f$ give
$$
 |\ell(x,y)|=\int_0^1A_\theta(|y-x|,x)\,d\theta
 \le2f(x)\int_0^1k_\theta(x,y)\,d\theta=2f(x)k(x,y).
\tag{19}
$$

The kernel $k$ need not be symmetric around $x$, and $S$ is not
asserted to be a martingale operator. Only the paired radial identity
(18) is used. From it and $A_\theta\le w_\theta$,
$$
 \int |y-x|k_\theta(x,y)\,dy
 =\mathbb E T_\theta
 =\int_0^\infty A_\theta(r,x)\,dr\le f_\theta,
$$
$$
 \int |y-x|^2k_\theta(x,y)\,dy
 =\mathbb E T_\theta^2
 =2\int_0^\infty rA_\theta(r,x)\,dr
 \le2f_\theta\mathbb E T_\theta\le2f_\theta^2.
\tag{20}
$$
These identities justify the claimed moment estimates for the mixture.

## 5. Independent derivation of the energy-flux identity

Apply (17) to $u$. Equation (10) implies $Su=u$. Apply it also to
$u^2$, which is continuous and has at most quadratic growth:
$$
 B=\Lambda(u^2)\in C^1,\qquad
 B'=2(S(u^2)-u^2).
$$
The row normalization and $Su=u$ give
$$
 S(u^2)(x)-u(x)^2
 =\int k(x,y)(u(y)-u(x))^2\,dy=E(x).
$$
Therefore
$$
 \boxed{B'=2E\ge0.}
\tag{21}
$$

For each fixed $x$, expand the square in the following integral.
Both the cross term and the constant term vanish by
$\Lambda u=\Lambda1=0$:
$$
 B(x)=\int\ell(x,y)(u(y)-u(x))^2\,dy.
\tag{22}
$$
This is a pointwise algebraic identity; it does not illicitly differentiate
the parameter $u(x)$ while applying (17). Equations (19), (21), and
(22) give
$$
 |B(x)|\le2f(x)E(x)=f(x)B'(x).
\tag{23}
$$
Finally, (4), row normalization, and (20) give
$$
 E(x)\le K^2\bigl(1+2f(x)+2f(x)^2\bigr).
\tag{24}
$$
Every integral in this derivation is absolutely convergent by the
polynomial growth and local exponential dominations already established.
No boundedness of $u$, finite total mass, or derivative of $h$ was used.

## 6. The positive-flux implication really forces finite right mass

Suppose $B(x_0)>0$. Since $B$ is nondecreasing, increase $x_0$ if
necessary to make $x_0\ge0$. On that half-line $B>0$, so (23) implies
$$
 (\log B)'\ge1/f.
$$
With $D(x)=\int_0^x f(z)^{-1}\,dz$,
$$
 B(x)\ge B(x_0)e^{D(x)-D(x_0)}.
\tag{25}
$$
Lipschitz continuity implies $f(x)\le f(0)+x$ for $x\ge0$, hence
$$
 D(x)\ge\log(1+x/f(0)).
$$
In particular $D\to\infty$ and (25) gives $B(x)\ge c x$ for large $x$.

Equations (23)–(24) also imply
$$
 B(x)\le2K^2f(x)(1+2f(x)+2f(x)^2).
\tag{26}
$$
The right side stays bounded if $f(x)$ stays bounded. Since $B(x)\to
\infty$, (26) implies the full limit $f(x)\to\infty$, not merely an
unbounded subsequence. Eventually $f\ge1$, and
$$
 B(x)\le10K^2 f(x)^3.
$$
Combining with (25),
$$
 f(x)\ge c_1 e^{D(x)/3}.
\tag{27}
$$
The substitution $dD=dx/f(x)$ is valid because $f>0$ is continuous
and $D$ is strictly increasing onto an unbounded interval. Thus
$$
 \int^\infty f(x)^{-2}\,dx
 =\int^\infty\frac{dD}{f(x(D))}
 \le c_1^{-1}\int^\infty e^{-D/3}\,dD<\infty.
\tag{28}
$$
Compact initial intervals contribute a finite amount. Applying the
independently proved (F2) to each endpoint density now gives
$$
 \int_0^\infty(g_1+g_2)<\infty.
\tag{29}
$$
There is no circular use of the uniqueness theorem, and no assumption
of (29) before it is established.

## 7. The finite-side energy limit, including explicit overshoot moments

Assume the conclusions just obtained: (29) and $f(x)\to\infty$.
Set
$$
 \tau(r)=\int_r^\infty(g_1+g_2),\qquad r\ge0.
$$
Then $\tau(0)<\infty$, $\tau(r)\downarrow0$, and for $a,b\ge r$,
$$
 |u(a)-u(b)|\le\tau(r).
\tag{30}
$$

### 7.1 Landings with $y\ge x/2$

For $x>0$, (30) and the fact that a row of $k$ has mass one imply
$$
 \int_{x/2}^\infty k(x,y)(u(y)-u(x))^2\,dy
 \le\tau(x/2)^2\longrightarrow0.
\tag{31}
$$

### 7.2 Landings with $0\le y<x/2$

Let $\lambda=(\log2)/2$. For any one of the transforms $f_\theta$,
$$
 f_\theta\ge2f_\theta e^{-H_\theta(2f_\theta,x)},
$$
so $H_\theta(2f_\theta,x)\ge\log2$. Convexity in time, and the trivial
bound $w_\theta\le1$ at shorter times, give the all-time estimate
$$
 w_\theta(t,x)\le2e^{-\lambda t/f_\theta(x)}.
$$
Using (9) and $f_\theta\le f$,
$$
 p_\theta(x,y)\le\frac{2e}{f(x)}
                  e^{-\lambda |y-x|/f(x)},\qquad
 k(x,y)\le\frac{6e}{f(x)}
                  e^{-\lambda |y-x|/f(x)}.
\tag{32}
$$
On the middle region, $|u(y)-u(x)|\le\tau(y)$. Its energy is bounded by
$$
 \frac{6e}{f(x)}e^{-\lambda x/(2f(x))}
       \int_0^{x/2}\tau(y)^2\,dy
 =
 3e\,z e^{-\lambda z/2}
       \left(\frac2x\int_0^{x/2}\tau(y)^2\,dy\right),
 \quad z=\frac{x}{f(x)}.
\tag{33}
$$
The first factor is bounded for all $z>0$. The factor in parentheses
tends to zero by Cesàro averaging of the bounded function $\tau^2$
that tends to zero. This verifies the packet's middle-region limit.

In particular, there is no hidden assumption $f=o(x)$. The case
$f(x)\asymp x$, which occurs for one-sided support, is included.

### 7.3 Negative landings: a uniform second-moment bound

Use the fixed $A,m$ from Section 1, and take $x\ge0$. At time $x+A$,
the interval $[x-t,x+t]$ contains $[-A,A]$, so for $s\ge A$,
$$
 w_\theta(x+s,x)
 \le w_\theta(x+A,x)e^{-m(s-A)}
 \le e^{-m(s-A)}.
$$
For $0\le s\le A$ use $w_\theta\le1$. Define the fixed envelope
$$
 \rho(s)=
 \begin{cases}
  1,&0\le s\le A,\\
  e^{-m(s-A)},&s>A.
 \end{cases}
$$
Then uniformly in $\theta\in[0,1]$ and $x\ge0$,
$$
 p_\theta(x,-s)\le\frac{e}{f(x)}\rho(s),\qquad
 k(x,-s)\le\frac{3e}{f(x)}\rho(s).
\tag{34}
$$
For the $p_\theta$ kernel this gives explicit truncated moments
$$
 \int_{-\infty}^0 |y|^j p_\theta(x,y)\,dy
 \le\frac{e}{f(x)}C_j,\quad j=0,1,2,
\tag{35}
$$
where
$$
 C_0=A+\frac1m,
$$
$$
 C_1=\frac{A^2}{2}+\frac{A}{m}+\frac1{m^2},
$$
$$
 C_2=\frac{A^3}{3}+\frac{A^2}{m}
                    +\frac{2A}{m^2}+\frac2{m^3}.
$$
For $k$, replace $e$ by $3e$. These constants do not depend on
$x$ or $\theta$. This supplies the explicit uniform second overshoot
moment requested in the review instructions.

The increment bound needed here follows by splitting at zero:
$$
 |u(-s)-u(x)|
 \le |u(-s)-u(0)|+|u(0)-u(x)|
 \le K(1+s)+\tau(0).
$$
Therefore the negative-landing energy is bounded by
$$
 E_-(x)\le\frac{6e}{f(x)}
   \left[K^2(C_0+2C_1+C_2)+\tau(0)^2C_0\right]
 \longrightarrow0.
\tag{36}
$$
No control on the mass of either density on the negative half-line was
assumed. Only their difference is controlled there by $K$; their
common background can be arbitrarily large.

Equations (31), (33), and (36) prove $E(x)\to0$ with the stated assumptions.

### 7.4 A shorter replacement for the two nonnegative regions

There is a simpler argument that avoids (32)–(33) entirely.
Since $w_\theta\le1$ and $f_\theta\ge f/(2e)$,
$$
 k(x,y)\le3e/f(x).
\tag{37}
$$
Fix $R>0$, and take $x\ge R$. On $0\le y<R$, the increment is at
most $\tau(0)$, while on $y\ge R$ it is at most $\tau(R)$. Thus
$$
 \int_0^\infty k(x,y)(u(y)-u(x))^2\,dy
 \le\frac{3eR\,\tau(0)^2}{f(x)}+\tau(R)^2.
\tag{38}
$$
First let $x\to\infty$, then $R\to\infty$. Combining (38) with (36)
proves the same energy limit.

I recommend (34)–(38) as a replacement for the packet's three-region
paragraph. It is shorter and makes the uniformity transparent. This is
a simplification, not a new hypothesis or a repair of a false conclusion.

## 8. The final contradiction, reflection, and zero energy

In the positive-flux branch, (23), $E(x)\to0$, and
$f(x)\le f(0)+x$ imply
$$
 0<B(x)\le2f(x)E(x)=o(x).
$$
This contradicts (25), which gave $B(x)\ge cx$. Therefore $B\le0$
at every point.

For reflection choose $\widetilde G_j(x)=G_j(-x)$. Then
$$
 \widetilde Q(x)=Q(-x),\quad
 \widetilde u(x)=-u(-x),\quad
 \widetilde f(x)=f(-x),\quad
 \widetilde A_\theta(r,x)=A_\theta(r,-x).
$$
The unit-interval bound is unchanged, and the reflected pair still has
equal transforms. Substituting in the definition of the flux gives
$$
 \widetilde B(x)
 =\int_0^1\int_0^\infty A_\theta(s,-x)
       [u(-x-s)^2-u(-x+s)^2]\,ds\,d\theta
 =-B(-x).
$$
Applying the already proved nonpositivity to the reflected pair yields
$B\ge0$. Hence $B\equiv0$ and $E\equiv0$.

For any fixed $x$, the function
$k(x,y)(u(y)-u(x))^2$ is nonnegative, and $k(x,y)>0$. Its integral
vanishes. Since $u$ is continuous, any nonzero increment would persist
on an open interval and give a strictly positive integral. Thus $u$
is constant and $h=u'=0$ almost everywhere.

This proves the claimed conclusion without assuming that an extremum of
$Q$ is attained or invoking optional stopping for an unbounded function.

## 9. Adversarial examples and what they actually test

### 9.1 One-sided support and conditional overshoot

Take the admissible density $g=\mathbf1_{(-\infty,0)}$. For $x\ge0$,
$$
 H(t,x)=0\quad(t\le x),\qquad
 H(x+s,x)=s^2/2\quad(s\ge0),
$$
and
$$
 f(x)=x+c_0,\qquad c_0=\sqrt{\pi/2}.
$$
For $y=-s<0$,
$$
 p(x,-s)=\frac{e^{-s^2/2}}{2(x+c_0)}.
$$
Therefore
$$
 \mathbb P_x(Y<0)=\frac{c_0}{2(x+c_0)},\qquad
 \mathbb E_x[(-Y)^2\mathbf1_{\{Y<0\}}]
       =\frac{c_0}{2(x+c_0)}.
$$
Both are $O(1/f(x))$, as needed. But
$$
 \mathbb E_x[(-Y)^2\mid Y<0]=1,
$$
which does not tend to zero. This is the exact reason the packet should
specify “unconditioned truncated moments.”

This example also tests $f(x)\asymp x$, infinite mass on the opposite
half-line, and a large interval with no density near the observation
point. None breaks the corrected estimate.

### 9.2 Narrow, very high spikes

The unit-interval assumption bounds the mass of the *difference*, not its
pointwise height. It allows spikes of height $n^2$ and width $n^{-2}$,
placed more than two units apart. It also permits arbitrarily high common
spikes with no uniform bound on their mass.

The argument uses $\int |h|$ to bound increments of $u$, so tall
difference spikes do not evade (4). Common spikes do not enter that
bound at all. They make some weights small, but the exact inequalities
$\int wH_t=1$, (2), and (16) still hold. No replacement of an a.e.
bound by a pointwise bound on $g$ was found.

The possibility that $w(T,x)H_t(T,x)$ has nonvanishing spikes at large
$T$ is real for unrestricted backgrounds. The packet avoids relying on
its pointwise convergence. Its weak-limit treatment is checked in Section
10.1.

### 9.3 Arbitrarily long empty intervals

Separated bumps can leave gaps whose lengths grow arbitrarily fast, and
the local transform scale can grow with a gap length. No proof step
requires a uniform positive density, bounded jump length, or $f=o(x)$.
The coefficient $z e^{-\lambda z/2}$ in (33) remains bounded for every
possible $z=x/f(x)$; alternatively (38) bypasses that ratio entirely.

### 9.4 Unbounded common background with bounded difference

For example,
$$
 g_2(x)=e^{x^4},\qquad g_1(x)=e^{x^4}+\tfrac1{10}\sin x
$$
are nonnegative, locally integrable, and have $K\le1/10$.
They are an admissible stress family for all the analytic constructions,
although they are not claimed to have equal data.

The proof does not need a global upper bound on either $g_j$, a
polynomial bound on their derivatives, or a lower bound on $f$ at
infinity. Local differentiation uses (2). In the hypothetical positive
flux branch, growth of $B$ itself forces $f\to\infty$, rather than
assuming it.

### 9.5 A nonzero constant flux

The explicit $g=c,\ u=x$ calculation in Section 3 breaks a plausible
weakened proof based only on $Su=u$. It does **not** break this packet,
because equal transform data give the stronger equation $\Lambda u=0$.

### 9.6 What was not found

None of these tests produces two distinct admissible densities with
exactly equal transforms. Examples that violate a discarded intermediate
claim, or merely produce close transforms, are not counterexamples to
the theorem. Point masses themselves are outside the stated density
class and were not used as alleged counterexamples.

## 10. Separate audit of the auxiliary assertions

### 10.1 Wave identity for arbitrary locally integrable inputs

Independently,
$$
 H_{tt}-H_{xx}=2g(x),\qquad
 w_{xx}-w_{tt}=2g(x)w-4wM_LM_R.
$$
Integrating over $0\le t\le T$, with $w_t(0,x)=0$, gives
$$
 f_T''=2gf_T-4I_T-w(T,x)H_t(T,x).
\tag{39}
$$
The sign of the boundary term is negative.

Finite-time quantities have the necessary local weak derivatives by
Fubini on compact rectangles. For nonnegative
$\varphi\in C_c^\infty$, (39) gives
$$
 4\int\varphi I_T+\int\varphi B_T
 =2\int\varphi gf_T-\int\varphi''f_T.
$$
The right side is uniformly bounded, using local uniform convergence of
$f_T$ and local integrability of $g$. Since $I_T$ increases,
$I\in L^1_{\rm loc}$. The same identity makes
$\int\varphi B_T$ converge. Its integral over $T$ is
$$
 \int_0^\infty\int\varphi B_T\,dx\,dT=\int\varphi\,dx<\infty.
$$
A nonnegative function with a finite time integral and an existing limit
must have limit zero. This justifies the weak boundary passage. General
test functions follow by domination by a nonnegative compactly supported
test function.

Together with $f\in C^1$, this proves
$$
 f\in W^{2,1}_{\rm loc},\qquad f''=2gf-4I
$$
in distributions and a.e. The packet does not improperly assert that
$I(x)$ is finite at every exceptional point.

### 10.2 Support, local masses, and the pointwise lower bound

The identities
$$
 \int_0^\infty wM_L\,dt=(1+f')/2,\qquad
 \int_0^\infty wM_R\,dt=(1-f')/2
$$
prove the support criterion. If either integral vanishes, positivity of
$w$, nonnegativity, and continuity of the corresponding window mass
force that window mass to vanish for every time.

For $0<s<f(x)$, monotonicity of $M_R$ gives
$$
 \frac{1-f'}2
 \ge M_R(s,x)\int_s^\infty w\,dt
 \ge M_R(s,x)(f(x)-s).
$$
The left bound is identical with $L,R$ exchanged. Also
$$
 f\ge2f\,e^{-H(2f,x)},\qquad H(2f,x)\le2f\,M_{2f}(x)
$$
prove the displayed lower local-mass bound.

Chebyshev under $w\,dt/f$ gives
$$
 If\ge\frac{1-f'^2}{4},
$$
first with bounded truncations if needed. Combining this with the
wave identity proves the pointwise lower bound a.e. The concavity claim
on an interval where $g=0$ follows in the distributional sense, hence
for the continuous representative of $f$.

### 10.3 Identification of an integrable density among all competitors

Dividing the wave identity by $2f$ and integrating by parts is valid
on compact intervals because $f>0$, $f\in C^1\cap W^{2,1}_{\rm loc}$.
It gives the claimed inequality
$$
 \int_A^B f^{-2}
 \le2\int_A^B g+\frac1{f(A)}+\frac1{f(B)}.
$$
For $g\in L^1$, one can show $f(x)\to\infty$ directly: for any fixed
$L>0$,
$$
 H(L,x)\le L\int_{x-L}^{x+L}g\longrightarrow0,
 \qquad f(x)\ge L e^{-H(L,x)}.
$$
Letting $L$ be arbitrary proves the limit. Taking the endpoints to
infinity in the inequality makes $f^{-2}$ integrable on the whole line.
Applying (F2) on both sides to any competitor with the same data makes
that competitor integrable too. The core theorem then applies.

This stronger “among all competitors” conclusion is justified for an
integrable input, in contrast to the general $L^2$ scope issue.

### 10.4 Periodic regularity: independently checked

Fix a period $P$. Write
$$
 g=c+p'',\qquad \int_0^Pp=0,\quad c>0.
$$
For bounded $g$, $p\in W^{2,\infty}$, and
$$
 H(t,x)=ct^2+p(x+t)+p(x-t)-2p(x).
\tag{40}
$$
The periodic part is uniformly bounded in $t$, giving a Gaussian
envelope. If $g\in C^{j,\alpha}$, then $p\in C^{j+2,\alpha}$.
The window masses satisfy
$$
 M_L=ct+p'(x)-p'(x-t),\qquad
 M_R=ct+p'(x+t)-p'(x),
$$
so they belong to $C^{j+1,\alpha}$, with norms growing at most
linearly in $t$. Gaussian integration gives
$$
 I\in C^{j+1,\alpha}.
$$
For merely bounded $g$, the same formulas in weak derivatives give
$I\in W^{1,\infty}$. In particular no derivative $g^{(j+1)}$
is needed for the stated gain.

The forward implication $g\in C^{k,\alpha}\Rightarrow f\in
C^{k+2,\alpha}$ follows directly from (40). For the converse, start
with $I\in W^{1,\infty}\subset C^{0,\alpha}$ and
$$
 g=(f''+4I)/(2f).
$$
This first yields $g\in C^{0,\alpha}$, then the gain for $I$ allows
iteration up to $C^{k,\alpha}$. These are periodic assertions for
bounded inputs, exactly as qualified in the packet.

### 10.5 Periodic local inverse theorem: extra argument supplied

The packet states the result without its Fredholm argument. I checked
that the argument can be completed from the stated hypotheses.

Let $U$ be the open subset of periodic $L^\infty$ functions with
positive mean, allowing signs. Formula (40) and Gaussian bounds for the
exponential and its variations show that
$$
 T:U\longrightarrow W^{2,\infty}
$$
is $C^1$. These bounds can be applied to scalar integrals and weak
spatial derivatives; one need not assume strong continuity of translations
in $L^\infty$.

At a nonnegative, nonzero $g_0$, differentiate the wave identity:
$$
 (\partial_x^2-2g_0)\,DT[g_0]h
       =2T(g_0)h-4DI[g_0]h.
\tag{41}
$$
The periodic operator
$\partial_x^2-2g_0:W^{2,\infty}\to L^\infty$ is Fredholm of index
zero, as a compact perturbation of $\partial_x^2-1$.
Its kernel is zero by
$$
 -\int|v'|^2=2\int g_0v^2.
$$
It is therefore an isomorphism.

Differentiating $I=\int wM_LM_R$ shows that
$DI[g_0]:L^\infty\to W^{1,\infty}$ is bounded: first spatial
derivatives of the varying window masses involve endpoint values of
$h$, not $h'$, and the remaining time factors are polynomials under
a Gaussian envelope. Thus $DI[g_0]$, regarded as an operator into
$L^\infty$, is compact. Since $T(g_0)>0$ on the compact period,
$2T(g_0)\operatorname{Id}-4DI[g_0]$ has index zero.

One must separately show derivative injectivity; nonlinear injectivity
alone would not establish it. Write $h=a+q''$, with $q$ periodic.
Then
$$
 DT[g_0]h(x)=-\int_0^\infty e^{-H_0(t,x)}
                   [at^2+\Delta_tq(x)]\,dt.
$$
If $a>0$, evaluate at a minimum of $q$: the result is strictly
negative. If $a<0$, evaluate at a maximum: it is strictly positive.
If $a=0$ and $q$ is nonconstant, a maximum of $q$ gives a
nonnegative, nonzero integral after the minus sign. Thus a vanishing
derivative forces $a=0$ and $q$ constant, so $h=0$.

Equation (41), index zero, and the inverse function theorem prove the
local diffeomorphism. The signed extension is essential at densities
that vanish somewhere; the cone of nonnegative inputs is not an open
neighborhood there. The packet states this qualification correctly.

### 10.6 Uniform periodic stability on Hölder-bounded classes

For fixed period, Hölder exponent, norm bound, and positive mean lower
bound, the nonnegative class is compact in $L^\infty$. Suppose no
uniform Lipschitz inverse estimate held. Pairs with divergent ratios
would have subsequences converging in $L^\infty$. Their data differences
would tend to zero, so continuity and injectivity would make both limits
the same density. The local Lipschitz inverse at that common limit then
contradicts the divergent ratios.

Thus the stated estimate follows, with the constant depending on the
class parameters. The argument supplies no explicit useful numerical
value for that constant, and the packet does not claim one.

### 10.7 Fourier multipliers and the loss of two derivatives

At $g=c$, the tent multiplier for frequency $n\ne0$ is
$2(1-\cos nt)/n^2$. Hence
$$
 DT[c](e^{inx})
 =-\frac2{n^2}\int_0^\infty e^{-ct^2}(1-\cos nt)\,dt\,e^{inx}
 =-\frac{\sqrt\pi}{\sqrt c}
       \frac{1-e^{-n^2/(4c)}}{n^2}\,e^{inx}.
$$
For the zero mode,
$$
 -\int_0^\infty t^2e^{-ct^2}\,dt
 =-\frac{\sqrt\pi}{4c^{3/2}}.
$$
The constants, signs, and zero-mode value are correct.

A uniform inverse Lipschitz bound in the same Sobolev norm would, by
taking arbitrarily small perturbations in each smooth Fourier direction,
imply $1\le C|m_n|$ for every $n$. This contradicts $m_n\to0$.
It does not contradict inversion with two additional derivatives of
data, and it says nothing by itself about nonlinear non-injectivity.

### 10.8 Controlled-tail approximation

I also checked the qualification on this optional assertion. For a right
tail with the displayed derivative bounds, fix $x$ large and put
$a=g(x)$. On $[x/2,3x/2]$, the logarithmic derivative bound gives
comparability with $a$, and the second derivative bound gives
$$
 |H(t,x)-at^2|\le C at^4/x^2,\qquad
 H(t,x)\ge c at^2,\qquad 0\le t\le x/2.
$$
Integration of the exponential difference gives a local error at most
$C/(a^{3/2}x^2)$. At $r=x/2$,
$$
 H(r,x)\ge cax^2,\quad H_t(r,x)\ge cax,
$$
so the entire remaining time integral is bounded by
$C e^{-cax^2}/(ax)$. Comparing also the Gaussian tail yields
$$
 \left|f(x)-\frac{\sqrt\pi}{2\sqrt a}\right|
 \le \frac{C}{a^{3/2}x^2}.
$$
If $ax^2\ge1$, a local integral also gives $f(x)\ge c/\sqrt a$,
which turns this into
$|\pi/(4f(x)^2)-a|\le C/x^2$.
If $ax^2<1$, instead $f(x)\ge cx$, and both $a$ and
$\pi/(4f(x)^2)$ are $O(x^{-2})$. Reflection handles the left side.

Thus the statement is valid with its derivative hypotheses. Merely
writing a local Gaussian expansion without the second case would leave
a gap. It is not an estimate for arbitrary $L^2$ inputs.

## 11. Scope, conclusions, and what remains unchecked

The core theorem, as checked here, says:
$$
 \begin{gathered}
 g_1,g_2\ge0,\quad g_j\in L^1_{\rm loc},\quad g_j\not\equiv0,\\
 \sup_a\int_a^{a+1}|g_1-g_2|<\infty,\quad T(g_1)=T(g_2)
 \quad\Longrightarrow\quad g_1=g_2\ \text{a.e.}
 \end{gathered}
$$

For $1\le p\le\infty$, Hölder on intervals of length one gives the
required bound for an $L^p$ difference. Thus **the proof does cover the
full nonnegative $L^2$ input class, with no finite-mass or controlled-tail
restriction**. Smoothness is not needed for that conclusion.

The original simultaneous assumptions $f,g\in L^2$ remain inconsistent:
the fixed-window estimate in the packet makes $f(x)\to\infty$ when
nonnegative nonzero $g\in L^2$. This obstruction is distinct from
injectivity after dropping the condition on $f$.

The following limits of scope should remain explicit:

1. Both candidates must satisfy the pairwise unit-interval difference
   condition. One $L^2$ input alone does not establish that condition
   against every arbitrary nonnegative locally integrable competitor.
2. An integrable input does have uniqueness among all such competitors,
   by the separate mass-detection argument.
3. Periodic inputs with different periods are included: each has bounded
   absolute mass on intervals of length one, so their difference does too.
4. Reflection preserves the $L^2$ class and commutes with $T$.
   Consequently an even transform of an $L^2$ input forces an even
   density.
5. The theorem does not cover signed densities, all arbitrary
   $L^1_{\rm loc}$ pairs, or arbitrary measures with atoms as stated.
6. The periodic regularity and stability assertions must retain their
   fixed-period and bounded-input hypotheses.
7. Injectivity does not supply a closed-form inverse, a range
   characterization, or convergence of a real-line reconstruction
   algorithm. None was established by this review.

**Checked:** the spatial and segment differentiations, all Fubini uses
in the core proof, the diagonal distribution, both kernel comparisons,
the normalized mixture, radial moments, zero flux, energy algebra,
growth implication, elementary mass detection, all three landing-region
limits, uniform second overshoot control, reflection, and the final
zero-energy conclusion. The auxiliary wave, mass, periodic regularity,
periodic inverse, stability, and multiplier assertions were checked
separately above.

**Not checked or claimed:** priority in the literature, previous audit
verdicts, correctness of accompanying numerical implementations, validated
floating-point certificates, existence for arbitrary proposed data,
real-line reconstruction convergence, or machine-checked formal validity.
No numerical agreement was used to establish a limit or exact equality.

The most useful change before circulating the proof is to replace the
compressed finite-side paragraph with (34)–(38), and to make the final
uniqueness wording explicitly relative to the allowed competitor class.
Those changes improve auditability without changing the theorem.

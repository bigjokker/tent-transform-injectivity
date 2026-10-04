"""Adversarial numerical stress test for the lemmas used in the uniqueness proof
(papers/injectivity-of-the-tent-transform.md, Section 3, Steps 1-4).  Two DIFFERENT densities g1, g2 and the
theta-mixture are used, so nothing relies on equal transforms except where noted.

Checked, with assertions:
 (K-)  k_theta(x,y) >= A_theta(|y-x|,x) / (2 f_theta(x))       [k from finite differences of F]
 (K+)  k_theta(x,y) <= 3 p_theta(x,y)
 (Kf)  k_theta(x,x+r) = p_theta + (1/2) d/dx A_theta(r,x)       [formula vs finite differences]
 (R1)  rows of k_theta integrate to 1
 (M1)  E T_theta <= f_theta,  E T_theta^2 <= 2 f_theta^2,  A_theta <= w_theta
 (Hol) f_theta <= f_2^(1-theta) f_1^theta                       [Hoelder]
 (Row) (Lambda v)'(x) = 2 (S v - v)(x)  for an arbitrary v      [mixture operator]
 (Sq)  S(v^2) - v^2 = E_v + 2 v (S v - v)                        [algebra used for B' = 2E]
 (Cmp) |Lambda((v - v(x))^2)(x)| <= 2 max_theta f_theta(x) * E_v(x)
Density families include zeros on intervals, narrow spikes, and one-sided support.
"""
import sys
import numpy as np
from scipy.integrate import cumulative_trapezoid

L, hy = 60.0, 1e-3
yy = np.arange(-L, L + hy/2, hy)
t = np.linspace(0.0, 40.0, 8001)
dt = t[1] - t[0]
X = np.linspace(-3.0, 3.0, 121)           # centres
Ymat = np.arange(-20.0, 20.0 + 1e-9, 0.05)   # wide landing grid (rows must not be truncated)
dx = X[1] - X[0]
dy = Ymat[1] - Ymat[0]
THETAS = np.linspace(0.0, 1.0, 9)
EPS = 1e-4


def prims(g):
    G1 = cumulative_trapezoid(g(yy), yy, initial=0.0)
    G2 = cumulative_trapezoid(G1, yy, initial=0.0)
    return G1, G2


def family(seed, kind):
    r = np.random.default_rng(seed)
    if kind == "smooth":
        c, w, a = r.uniform(-3, 3, 4), r.uniform(.3, 1.2, 4), r.uniform(.2, 2, 4)
        return lambda y: np.sum(a*np.exp(-((np.asarray(y)[..., None]-c)/w)**2/2), -1)
    if kind == "zeros":       # vanishes on (-0.5, 1.5)
        base = family(seed, "smooth")
        return lambda y: base(y)*((np.asarray(y) < -0.5) | (np.asarray(y) > 1.5))
    if kind == "spiky":       # narrow tall spikes on a weak background
        c = r.uniform(-2, 2, 3)
        return lambda y: 0.02 + np.sum(5.0*np.exp(-((np.asarray(y)[..., None]-c)/0.03)**2/2), -1)
    if kind == "onesided":    # support in y < -1 only
        base = family(seed, "smooth")
        return lambda y: (base(y) + 0.3)*(np.asarray(y) < -1.0)
    raise ValueError(kind)


def objects(G1a, G2a, G1b, G2b, theta, x):
    """theta-mixture of densities a (weight 1-theta) and b (weight theta) at centre x."""
    G1 = (1-theta)*G1a + theta*G1b
    G2 = (1-theta)*G2a + theta*G2b
    P1 = lambda z: np.interp(z, yy, G1)
    P2 = lambda z: np.interp(z, yy, G2)
    H = P2(x+t) + P2(x-t) - 2*P2(x)
    w = np.exp(-H)
    ML = P1(x) - P1(x-t)
    MR = P1(x+t) - P1(x)
    tail = np.concatenate([np.cumsum((w[::-1][1:]+w[::-1][:-1])*dt/2)[::-1], [0.0]])
    f = tail[0]
    return dict(f=f, w=w, A=tail/f, ML=ML, MR=MR)


def F_row(obj, x):
    r = np.abs(Ymat - x)
    A = np.interp(r, t, obj["A"])
    return np.where(Ymat >= x, 1 - 0.5*A, 0.5*A)


def run(kind1, kind2, seed):
    ga, gb = family(seed, kind1), family(seed + 100, kind2)
    G1a, G2a = prims(ga)
    G1b, G2b = prims(gb)
    worst = dict(Kminus=0, Kplus=0, Kf=0, R1=0, M1=0, Hol=0, Row=0, Sq=0, Cmp=0)
    vfun = lambda z: np.sin(1.7*z) + 0.3*z + 0.2*np.cos(4.1*z + 0.3)
    v = vfun(Ymat)
    Amix = np.zeros((len(X), len(t)))
    f0 = np.array([objects(G1a, G2a, G1b, G2b, 0.0, x)["f"] for x in X])
    f1 = np.array([objects(G1a, G2a, G1b, G2b, 1.0, x)["f"] for x in X])
    Kmix = np.zeros((len(X), len(Ymat)))
    Lmix = np.zeros_like(Kmix)
    fmax = np.zeros(len(X))
    wq = np.full(len(THETAS), THETAS[1]-THETAS[0]); wq[0] = wq[-1] = wq[0]/2
    for th, q in zip(THETAS, wq):
        for i, x in enumerate(X):
            o = objects(G1a, G2a, G1b, G2b, th, x)
            op, om = objects(G1a, G2a, G1b, G2b, th, x+EPS), objects(G1a, G2a, G1b, G2b, th, x-EPS)
            k = -(F_row(op, x+EPS) - F_row(om, x-EPS))/(2*EPS)   # note: F evaluated at moving x
            r = np.abs(Ymat - x)
            A = np.interp(r, t, o["A"])
            p = np.interp(r, t, o["w"])/(2*o["f"])
            off = (r > 2*dy) & (p > 1e-7) & (np.abs(Ymat) < 2.5)
            worst["Kminus"] = max(worst["Kminus"], np.max(((A/(2*o["f"]) - k)/p)[off]))
            worst["Kplus"] = max(worst["Kplus"], np.max((k/p)[off]) - 3)
            # formula k = p + 1/2 sgn d/dx A(r,x) at fixed r
            dA = (np.interp(r, t, op["A"]) - np.interp(r, t, om["A"]))/(2*EPS)
            kform = p + 0.5*np.sign(Ymat - x)*dA
            worst["Kf"] = max(worst["Kf"], np.max(np.abs(k - kform)[off]/p[off]))
            worst["M1"] = max(worst["M1"],
                              np.trapezoid(o["A"], t)/o["f"] - 1,
                              2*np.trapezoid(t*o["A"], t)/o["f"]**2 - 2,
                              np.max(o["A"] - o["w"]))
            worst["Hol"] = max(worst["Hol"], o["f"] - f0[i]**(1-th)*f1[i]**th)
            Kmix[i] += q*k
            Lmix[i] += q*np.sign(Ymat - x)*A
            Amix[i] += q*o["A"]
            fmax[i] = max(fmax[i], o["f"])
    inner = np.abs(X) < 1.5
    rows = Kmix.sum(1)*dy
    worst["R1"] = np.max(np.abs(rows[inner] - 1))
    Sv = Kmix @ v*dy
    Lv = np.array([np.trapezoid(Amix[i]*(vfun(X[i]+t) - vfun(X[i]-t)), t) for i in range(len(X))])
    dLv = np.gradient(Lv, dx)
    vx = np.sin(1.7*X) + 0.3*X + 0.2*np.cos(4.1*X + 0.3)
    worst["Row"] = np.max(np.abs(dLv - 2*(Sv - vx))[inner])
    Ev = np.array([np.sum(Kmix[i]*(v - vx[i])**2)*dy for i in range(len(X))])
    Sv2 = Kmix @ v**2*dy
    worst["Sq"] = np.max(np.abs(Sv2 - vx**2 - Ev - 2*vx*(Sv - vx))[inner])
    Bloc = np.array([np.trapezoid(Amix[i]*((vfun(X[i]+t)-vx[i])**2 - (vfun(X[i]-t)-vx[i])**2), t) for i in range(len(X))])
    worst["Cmp"] = np.max((np.abs(Bloc) - 2*fmax*Ev)[inner])
    return worst


tol = dict(Kminus=2e-3, Kplus=1e-3, Kf=5e-3, R1=5e-3, M1=2e-3, Hol=1e-9, Row=2e-2, Sq=2e-3, Cmp=1e-3)
fails = 0
for kinds in [("smooth", "smooth"), ("zeros", "smooth"), ("spiky", "smooth"), ("onesided", "spiky")]:
    for seed in (1, 2):
        wv = run(*kinds, seed)
        bad = [k for k in wv if wv[k] > tol[k]]
        fails += len(bad)
        print(f"{kinds[0]:>8}/{kinds[1]:<8} seed {seed}: " +
              "  ".join(f"{k}={wv[k]:+.1e}" for k in wv) + ("   FAIL " + str(bad) if bad else "   ok"))
print("ALL PASSED" if fails == 0 else f"{fails} FAILURES")
sys.exit(1 if fails else 0)

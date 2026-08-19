#!/usr/bin/env python3
"""
Maximum-likelihood fit of  delta_q = c_q * q^(-alpha),  c_q ~ N(0, sigma^2),
with the (known, q-independent) binomial measurement noise in the model.
Run:  python3 src/fit_alpha.py results/delta_q.json
"""
import sys, json, math
import numpy as np

d = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "results/delta_q.json"))
qs = np.array(sorted(int(k) for k in d), dtype=float)
dl = np.array([d[str(int(q))][0] for q in qs])
Z  = np.array([d[str(int(q))][1] for q in qs])
s = float(np.median(np.abs(dl/Z)))
print(f"n = {len(qs)} primes, q in [{int(qs[0])}, {int(qs[-1])}], measurement SE = {s:.5f}\n")

def negll(alpha, sigma):
    var = sigma**2 * qs**(-2*alpha) + s**2
    return 0.5*np.sum(np.log(var) + dl**2/var)

def profile(alpha):
    lo, hi = 1e-4, 5.0
    for _ in range(300):
        m1, m2 = lo + (hi-lo)/3, hi - (hi-lo)/3
        if negll(alpha, m1) < negll(alpha, m2): hi = m2
        else: lo = m1
    return negll(alpha, (lo+hi)/2), (lo+hi)/2

al = np.linspace(0.15, 0.95, 321)
ll = np.array([profile(x)[0] for x in al])
sg = np.array([profile(x)[1] for x in al])
i = int(ll.argmin())
print(f"alpha_hat = {al[i]:.3f}   sigma_hat = {sg[i]:.4f}")
for lvl, lab in [(0.5, "68%"), (2.0, "95%")]:
    ok = al[ll <= ll[i]+lvl]
    print(f"  {lab} interval: [{ok.min():.3f}, {ok.max():.3f}]")
print("\ncandidate exponents (chi^2 with 1 dof; reject above 3.84):")
for cand, lab in [(0.5, "1/2  Euler-factor / Sato-Tate"), (1/3, "1/3"),
                  (0.25, "1/4"), (3/8, "3/8"), (0.4, "2/5"), (5/12, "5/12"),
                  (2/3, "2/3"), (0.75, "3/4")]:
    # profile at the EXACT candidate. Snapping to the alpha grid put 1/3, 5/12
    # and 2/3 up to 0.00083 off their true values, which visibly shifted their
    # d(-2logL); the grid is for locating the optimum, not for evaluating tests.
    llc, _ = profile(cand)
    d = 2*(llc-ll[i])
    print(f"  {cand:.4f}  {lab:<30}  d(-2logL) = {d:7.2f}"
          f"  {'consistent' if d < 3.84 else 'REJECTED'}")
print("\nnoise-deconvolved signal sd by q-bin (constant <=> alpha = 1/2):")
x = dl*np.sqrt(qs)
edges = [5, 60, 200, 600, 1200, 2400, 5500, 20000, 100000]
for lo, hi in zip(edges[:-1], edges[1:]):
    msk = (qs >= lo) & (qs < hi); n = int(msk.sum())
    if n < 8: continue
    var = x[msk].var() - (s**2*qs[msk]).mean()
    sd = math.sqrt(max(var, 0))
    print(f"  q in [{lo:>6},{hi:>6}): n={n:>4}  sd = {sd:.4f} +- {sd/math.sqrt(2*n):.4f}")

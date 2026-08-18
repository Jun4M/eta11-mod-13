#!/usr/bin/env python3
"""
Reproduces the two measured laws:
   E(t)  = 13*P(a(t)=0 | t squarefree) - 1        -- decays like t^-0.28
   f_q   = delta_q / E                            -- decays like q^-0.41
Run:  python3 src/analyze.py a13.bin [--qmax 5500]
Writes results/E_of_t.json and results/delta_q.json
"""
import sys, math, json, os, time
import numpy as np

path = sys.argv[1] if len(sys.argv) > 1 else "a13.bin"
QMAX = 5500
if "--qmax" in sys.argv: QMAX = int(sys.argv[sys.argv.index("--qmax")+1])
os.makedirs("results", exist_ok=True)
t0 = time.time()

a = np.fromfile(path, dtype=np.int8); NC = len(a)
i = np.arange(NC, dtype=np.int64); m = 24*i + 11; MMAX = int(m[-1])
print(f"{NC:,} coefficients, m up to {MMAX:,}", flush=True)

# squarefree kernel of m  (m odd and coprime to 3, so only p>=5 matter)
LIM = int(MMAX**0.5) + 1
sv = np.ones(LIM+1, bool); sv[:2] = False
for j in range(2, int(LIM**0.5)+1):
    if sv[j]: sv[j*j::j] = False
sq = np.ones(NC, dtype=np.int64)
for p in np.nonzero(sv)[0]:
    p = int(p)
    if p < 5: continue
    q = p*p
    if q > MMAX: break
    while q <= MMAX:
        k0 = ((-11)*pow(24, -1, q)) % q
        idx = np.arange(k0, NC, q)
        if len(idx) == 0: break
        sq[idx] *= p
        q *= p*p
t = m // (sq*sq)
sf = (sq == 1)
print("kernel sieve done", round(time.time()-t0, 1), flush=True)

z = (a == 0)
QR13 = np.array([1, 3, 4, 9, 10, 12])
leg13 = np.where(m % 13 == 0, 0, np.where(np.isin(m % 13, QR13), 1, -1))

# ---- E(t) by half-decade window -------------------------------------------
print("\n=== E(t): normalised excess of vanishing, squarefree kernels ===")
print("  t window                    B        E        n")
Ecurve = []
e = 4.0
while 10**(e+0.5) <= MMAX + 24:   # +24: m runs in steps of 24, so a window
                                  # whose top is within one step is fully covered
    lo, hi = 10**e, 10**(e+0.5)
    w = (t >= lo) & (t < hi) & (leg13 != 0)
    if w.sum() >= 2000:
        B = 13*z[w].mean()
        Ecurve.append((math.sqrt(lo*hi), B-1, int(w.sum())))
        print(f" [{lo:>13,.0f},{hi:>14,.0f})  {B:.4f}  {B-1:+.4f}  {w.sum():>11,}")
    e += 0.5
json.dump(Ecurve, open("results/E_of_t.json", "w"))
# E(t) is NOT well described by a single power law: the smallest windows are
# essentially flat (E ~ 0.75) and only the large-t ones decay, so the exponent
# is a statement about which windows were fitted.  Both are printed; neither is
# "the" exponent.  The stable quantity to quote is the N_exc(X) exponent below,
# which is a fit to a cumulative count rather than to windowed differences.
if len(Ecurve) >= 4:
    g = np.array([r[0] for r in Ecurve]); v = np.array([r[1] for r in Ecurve])
    print()
    for lab, msk in [("all windows", v > 0), ("t >= 1e6 only", (v > 0) & (g > 1e6))]:
        if msk.sum() >= 4:
            s, c = np.polyfit(np.log(g[msk]), np.log(v[msk]), 1)
            print(f"  fit [{lab:<13}, {int(msk.sum())} windows]: "
                  f"E(t) = {math.exp(c):.3f} * t^({s:+.4f})")

# ---- cumulative excess count ----------------------------------------------
print("\n=== cumulative excess count  N(X) = #{t<=X : a(t)=0} - #{t<=X}/13 ===")
rows = []
for X in [10**x for x in np.arange(6.0, math.log10(MMAX)+0.001, 0.5)]:
    cut = sf & (m <= X) & (leg13 != 0)
    n = int(cut.sum())
    if n < 5000: continue
    N = z[cut].sum() - n/13
    rows.append((X, N))
    print(f"  X = {X:>15,.0f}   N_exc = {N:>12,.0f}")
if len(rows) >= 4:
    X = np.array([r[0] for r in rows]); Y = np.array([r[1] for r in rows])
    ok = Y > 0
    b, c = np.polyfit(np.log(X[ok]), np.log(Y[ok]), 1)
    print(f"\n  fit: N_exc(X) = {math.exp(c):.4f} * X^({b:.4f})")

# ---- delta_q --------------------------------------------------------------
# numerator and denominator must share one population: the factorisation
# delta_q = E * f_q presupposes uniform E, and the 13|m stratum has its own,
# much larger E (reported below).  (t|13) = 0 is structurally excluded for
# q = 13, so excluding 13|m everywhere keeps all 723 primes comparable.
w = sf & (m >= MMAX//10) & (leg13 != 0)   # top decade, 13 does not divide m
mw = m[w]; zw = z[w]
Ew = 13*zw.mean() - 1
print(f"\n=== delta_q on the top decade (E = {Ew:+.4f}, n = {len(mw):,}) ===")
s13 = sf & (m >= MMAX//10) & (leg13 == 0)
print(f"  excluded 13|m stratum: E = {13*z[s13].mean()-1:+.4f}, n = {int(s13.sum()):,}")
def isprime(x):
    if x < 2: return False
    d = 2
    while d*d <= x:
        if x % d == 0: return False
        d += 1
    return True
out = {}
for q in range(5, QMAX):
    if not isprime(q): continue
    flag = np.zeros(q, dtype=np.int8)
    for x in range(1, q): flag[(x*x) % q] = 1
    flag[0] = 0
    r = mw % q
    lab = flag[r]
    isq = lab == 1
    isn = (r != 0) & (~isq)
    n1, n2 = int(isq.sum()), int(isn.sum())
    z1, z2 = int(zw[isq].sum()), int(zw[isn].sum())
    p1, p2 = z1/n1, z2/n2
    pp = (z1+z2)/(n1+n2)
    Z = (p1-p2)/math.sqrt(pp*(1-pp)*(1/n1+1/n2))
    out[q] = (13*(p1-p2), Z, 13*(p1-p2)/Ew)
json.dump(out, open("results/delta_q.json", "w"))
sig = sum(1 for v in out.values() if abs(v[1]) >= 3)
print(f"  {len(out)} primes; {sig} with |Z| >= 3  (chance: {0.0027*len(out):.1f})")
print(f"  written to results/delta_q.json  [delta_q, Z, f_q]")
print("\ntotal", round(time.time()-t0, 1), "s")

#!/usr/bin/env python3
"""
Reproduces the three tables of section 4 from res_<l>.bin.

  Table 1  purity of the ratio A(t p^2)/A(t) on each Legendre class of (t|p)
           -- 100% for 13 <= l <= 31, about 1/l for l >= 37
  Table 2  the fitted Hecke exponent e, predicted to be k-1 with k = (l-3)/2
  Table 3  the character eps_l(p), predicted to be -(p|3)*(-1|p)^(k+1)

Run:  python3 src/controls.py            (expects res_<l>.bin in the cwd)
      python3 src/controls.py --tmax 150000 --fitmax 40
Time: a few minutes.

The purity percentages of Table 1 for l >= 37 are finite-sample estimates and
move with the kernel bound. Every run prints the bounds it used; do not quote a
purity figure without them.
"""
import sys, math, os
import numpy as np
from collections import Counter

PRIMES = [13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
AUX    = [5, 7, 11, 17, 19]          # auxiliary p used for the purity test

# Kernel bound. The dichotomy of Table 1 (100% vs about 1/l) and the whole of
# Tables 2 and 3 are insensitive to it. The purity PERCENTAGES for l >= 37 are
# not: they are finite-sample estimates of an equidistribution and they move.
TMAX    = 400000
FIT_MAX = 60                         # auxiliary p < FIT_MAX for exponent/character
for _i, _a in enumerate(sys.argv):
    if _a == "--tmax":   TMAX = int(sys.argv[_i+1])
    if _a == "--fitmax": FIT_MAX = int(sys.argv[_i+1])
FIT_P  = [p for p in range(5, FIT_MAX)]

def isprime(x):
    if x < 2: return False
    d = 2
    while d*d <= x:
        if x % d == 0: return False
        d += 1
    return True

def squarefree_sieve(N):
    """SF[t] is True iff t is squarefree.

    Computed once at import. The trial-division test this replaces was re-run
    for every (l, p) pair over the same t values and dominated the runtime."""
    sf = np.ones(N + 1, dtype=bool)
    i = 2
    while i * i <= N:
        sf[i*i :: i*i] = False
        i += 1
    return sf

SF = squarefree_sieve(TMAX + 1)

def load(l):
    """A_l(m) = p((l*m+1)/24) mod l, for m in the progression m0 + 24*Z."""
    path = f"res_{l}.bin"
    if not os.path.exists(path): return None
    res = np.fromfile(path, dtype=np.uint8).astype(np.int64)
    d = next(x for x in range(l) if (24*x) % l == 1)
    m0 = (24*d - 1)//l
    mmax = m0 + 24*(len(res)-1)
    return res, d, m0, mmax

def A_of(res, m0, mmax, m):
    if m % 24 != m0 % 24 or m > mmax or m < m0: return None
    return int(res[(m - m0)//24])

# ---------------------------------------------------------------- table 1
print("=== Table 1: is A(t p^2)/A(t) constant on each class of (t|p)? ===")
print("Worst case over p in", AUX, "and over the classes (t|p) = +1, -1.")
print(f"Kernel bound TMAX = {TMAX:,}.  Purity for l >= 37 is a finite-sample")
print("estimate at this bound and moves with it.\n")
print("    l | delta | m mod 24 | worst purity |   1/l   | verdict")
print("  " + "-"*62)
for l in PRIMES:
    got = load(l)
    if got is None:
        print(f"  {l:>4} | res_{l}.bin missing"); continue
    res, d, m0, mmax = got
    inv = [0] + [pow(i, -1, l) for i in range(1, l)]
    cls = m0 % 24
    worst = 1.0
    for p in AUX:
        if p == l: continue
        g = {1: Counter(), -1: Counter()}
        t = cls if cls > 0 else 24
        while t*p*p <= mmax and t < TMAX:
            if SF[t]:
                a0 = A_of(res, m0, mmax, t)
                if a0:
                    v = A_of(res, m0, mmax, t*p*p)
                    if v is not None:
                        r = t % p
                        if r:
                            lg = 1 if pow(r, (p-1)//2, p) == 1 else -1
                            g[lg][(v*inv[a0]) % l] += 1
            t += 24
        for k in (1, -1):
            tot = sum(g[k].values())
            if tot >= 200:
                worst = min(worst, g[k].most_common(1)[0][1]/tot)
    print(f"  {l:>4} |  {d:>3}  |    {cls:>2}    |   {100*worst:6.2f}%    | {100/l:5.2f}% |"
          f" {'EIGENFORM' if worst > 0.999 else 'no'}")

# ------------------------------------------------------- tables 2 and 3
print("\n=== Table 2/3: Hecke exponent and character, for the primes where it holds ===")
print("Prediction: exponent e = k-1 and eps_l(p) = -(p|3)*(-1|p)^(k+1), k = (l-3)/2.")
print(f"Kernel bound TMAX = {TMAX:,}, auxiliary primes 5 <= p < {FIT_MAX}.\n")
print("    l |  k | predicted e | fitted e | character rule holds? | primes used")
print("  " + "-"*74)
for l in [13, 17, 19, 23, 29, 31]:
    got = load(l)
    if got is None: continue
    res, d, m0, mmax = got
    inv = [0] + [pow(i, -1, l) for i in range(1, l)]
    cls = m0 % 24
    k = (l-3)//2
    cc = {}
    for p in FIT_P:
        if not isprime(p) or p == l: continue
        g = {1: Counter(), -1: Counter()}
        t = cls if cls > 0 else 24
        while t*p*p <= mmax and t < TMAX:
            if SF[t]:
                a0 = A_of(res, m0, mmax, t)
                if a0:
                    v = A_of(res, m0, mmax, t*p*p)
                    if v is not None:
                        r = t % p
                        if r:
                            lg = 1 if pow(r, (p-1)//2, p) == 1 else -1
                            g[lg][(v*inv[a0]) % l] += 1
            t += 24
        if min(sum(g[1].values()), sum(g[-1].values())) < 50: continue
        rp = g[1].most_common(1)[0][0]
        rm = g[-1].most_common(1)[0][0]
        cc[p] = ((rp - rm) * inv[2]) % l

    # exponent: the e for which c_p / p^e is always +-1
    fitted = None
    for e in range(l-1):
        ok = True
        for p, cp in cc.items():
            pe = pow(p, e, l)
            if pe == 0 or (cp*inv[pe]) % l not in (1, l-1): ok = False; break
        if ok: fitted = e; break

    # character, using the fitted exponent
    holds, used = True, 0
    if fitted is not None:
        for p, cp in cc.items():
            s = (cp * inv[pow(p, fitted, l)]) % l
            eps = 1 if s == 1 else -1
            leg3 = 1 if p % 3 == 1 else -1          # (p|3)
            m1   = 1 if p % 4 == 1 else -1          # (-1|p)
            used += 1
            if eps != -leg3 * (m1**(k+1)): holds = False
    print(f"  {l:>4} | {k:>2} |     {k-1:>2}      |    {str(fitted):>4}  |"
          f"        {'YES' if (holds and fitted == k-1) else 'no':<4}           |   {used}")
print(f"\nDone.  TMAX = {TMAX:,}, FIT_MAX = {FIT_MAX}.")

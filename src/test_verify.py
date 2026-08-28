#!/usr/bin/env python3
"""
Independent verification of everything the project relies on.
Run:  python3 src/test_verify.py [path/to/a13.bin]

Every check here uses a DIFFERENT method from the one that produced a13.bin,
so agreement is real cross-validation rather than a repeat of the same code.
"""
import sys, math, numpy as np

FAIL = []
def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok: FAIL.append(name)

# ---------------------------------------------------------------- exact p(n)
def exact_partitions(N):
    """p(n) as exact Python integers, pentagonal recurrence."""
    p = [0]*(N+1); p[0] = 1
    pen = []; k = 1
    while True:
        g1 = k*(3*k-1)//2
        if g1 > N: break
        s = 1 if k % 2 else -1
        pen.append((g1, s))
        g2 = k*(3*k+1)//2
        if g2 <= N: pen.append((g2, s))
        k += 1
    pen.sort()
    for n in range(1, N+1):
        t = 0
        for o, s in pen:
            if o > n: break
            t += p[n-o] if s > 0 else -p[n-o]
        p[n] = t
    return p

N_EXACT = 40000
print(f"computing exact p(n) for n <= {N_EXACT} ...")
P = exact_partitions(N_EXACT)
check("reference values",
      P[10] == 42 and P[100] == 190569292 and P[54] == 386155,
      f"p(100)={P[100]}, {len(str(P[N_EXACT]))} digits at n={N_EXACT}")

# ------------------------------------------------- load a13.bin if available
path = sys.argv[1] if len(sys.argv) > 1 else "a13.bin"
try:
    a = np.fromfile(path, dtype=np.int8)
except FileNotFoundError:
    a = None
    print(f"(!) {path} not found -- skipping checks that need it")

def A(m):
    """coefficient of q^m in prod (1-q^k)^11, mod 13"""
    if a is None or m % 24 != 11: return None
    i = (m-11)//24
    return int(a[i]) if 0 <= i < len(a) else None

# ------------------------------------------ 1. the reduction p(n) = 11 a(m)
if a is not None:
    bad = []
    for n in range(6, N_EXACT+1, 13):
        m = (24*n-1)//13
        v = A(m)
        if v is None: break
        if P[n] % 13 != (11*v) % 13: bad.append(n)
    check("reduction  p(n) = 11*a((24n-1)/13) mod 13",
          not bad, f"tested n <= {N_EXACT}, n = 6 mod 13; {len(bad)} mismatches")

# --------------------------- 2. eta^11 recomputed independently (small range)
def eta_power_mod(e, N, mod):
    """prod (1-q^k)^e mod `mod`, coefficients 0..N-1, via numpy (different code path)"""
    pen = []; k = 1
    while True:
        g1 = k*(3*k-1)//2
        if g1 >= N: break
        s = -1 if k % 2 else 1
        pen.append((g1, s))
        g2 = k*(3*k+1)//2
        if g2 < N: pen.append((g2, s))
        k += 1
    f = np.zeros(N, dtype=np.int64); f[0] = 1
    for _ in range(e):
        g = f.copy()
        for off, s in pen:
            if s > 0: g[off:] += f[:N-off]
            else:     g[off:] -= f[:N-off]
        f = g % mod
    return f

NS = 200000
ref = eta_power_mod(11, NS, 13)
if a is not None:
    n = min(NS, len(a))
    check("a13.bin agrees with independent numpy recomputation",
          np.array_equal(ref[:n], a[:n].astype(np.int64)), f"{n} coefficients")

# ------------------- 3. no EXACT zeros (vanishing is genuinely 13-adic)
big = eta_power_mod(11, NS, (1 << 31) - 1)
check("no exact zeros of eta^11 in the tested range",
      int((big == 0).sum()) == 0, f"checked {NS} coefficients mod 2^31-1")

# ------------------- 4. Hecke multiplicative law (Folsom-Kent-Ono, l=13)
# lambda_13 = 6 included: at p = 13 both relations degenerate, since 13^4 and
# 13^9 are 0 mod 13, so the Legendre term and the second term of the recursion
# drop out and the formulas below reduce to a(t*13^(2j)) = 6^j * a(t). The same
# code therefore tests it without a special case.
LAM = {5:10, 7:8, 11:5, 13:6, 17:1, 19:8, 23:8, 29:4, 31:4, 37:5, 41:9, 43:12, 47:6}
def sqfree(x):
    y, d = x, 2
    while d*d <= y:
        if y % (d*d) == 0: return False
        if y % d == 0: y //= d
        d += 1
    return True
def leg(t, p):
    r = t % p
    return 0 if r == 0 else (1 if pow(r, (p-1)//2, p) == 1 else -1)

if a is not None:
    MMAX = 24*(len(a)-1)+11
    bad = tot = 0
    for p, lam in LAM.items():
        eps = 1 if p % 3 == 2 else -1
        t = 11
        while t*p*p <= MMAX and t < 200000:
            if sqfree(t):
                a0 = A(t); a1 = A(t*p*p)
                if a0 is not None and a1 is not None:
                    pred = ((lam + eps*leg(t, p)*pow(p, 4, 13)) * a0) % 13
                    tot += 1
                    if pred != a1 % 13: bad += 1
            t += 24
    check("Hecke law  a(t p^2) = (lam_p + eps(p)(t|p)p^4) a(t)",
          bad == 0, f"{tot} congruences tested, {bad} failures")

    bad2 = tot2 = 0
    for p, lam in LAM.items():
        t = 11
        while t*p**4 <= MMAX and t < 200000:
            if sqfree(t):
                a0, a1, a2 = A(t), A(t*p*p), A(t*p**4)
                if None not in (a0, a1, a2):
                    tot2 += 1
                    if (lam*a1 - pow(p, 9, 13)*a0) % 13 != a2 % 13: bad2 += 1
            t += 24
    check("three-term recursion  a(t p^4) = lam_p a(t p^2) - p^9 a(t)",
          bad2 == 0, f"{tot2} congruences tested, {bad2} failures")

# ------------------- 5. square-class propagation (Corollary of the above)
if a is not None:
    t0, bad3, n3 = 155, 0, 0
    s = 1
    while t0*s*s <= MMAX:
        if s % 2 == 1 and s % 3 != 0:
            v = A(t0*s*s)
            if v is not None:
                n3 += 1
                if v != 0: bad3 += 1
        s += 1
    check("square class of t=155 vanishes entirely",
          bad3 == 0 and n3 > 5, f"{n3} coefficients, {bad3} nonzero")

# ------------------- 6. regenerated delta_q reproduces the reference
# Runs only if analyze.py has already been run; test_verify.py is normally run
# first, so a missing results/delta_q.json is not a failure.
import json, os
REF, CUR = "data/delta_q_consistent.json", "results/delta_q.json"
if os.path.exists(REF) and os.path.exists(CUR):
    ref, cur = json.load(open(REF)), json.load(open(CUR))
    if set(ref) != set(cur):
        check("delta_q reproduces data/delta_q_consistent.json", False,
              f"prime sets differ: {len(ref)} reference, {len(cur)} current")
    else:
        ks = sorted(ref, key=int)
        dd = max(abs(ref[k][0]-cur[k][0]) for k in ks)
        dz = max(abs(ref[k][1]-cur[k][1]) for k in ks)
        exact = sum(1 for k in ks if ref[k][0] == cur[k][0] and ref[k][1] == cur[k][1])
        # tolerance is far below any real effect (delta_q ~ 1e-2, SE ~ 1.2e-3);
        # bit-exact agreement is what a same-platform rerun actually gives
        check("delta_q reproduces data/delta_q_consistent.json",
              dd < 1e-12 and dz < 1e-9,
              f"{len(ks)} primes, {exact} bit-exact, "
              f"max|d delta|={dd:.2e}, max|d Z|={dz:.2e}")
elif os.path.exists(REF):
    print(f"(!) {CUR} not found -- run analyze.py, then rerun this check")

print()
if FAIL:
    print("FAILURES:", ", ".join(FAIL)); sys.exit(1)
print("all checks passed")

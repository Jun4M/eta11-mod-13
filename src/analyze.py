#!/usr/bin/env python3
"""
Reproduces the two measured laws:
   E(t)  = 13*P(a(t)=0 | t squarefree) - 1
   f_q   = delta_q / E
Run:  python3 src/analyze.py a13.bin [--qmax 5500]
Writes results/E_of_t.json and results/delta_q.json

Memory: the earlier version materialised i, m, t, leg13 and sq as int64 over the
whole range, which is 21 GB at MMAX = 1e10 and does not fit. This version keeps
only sq (uint32) and two bool arrays full-length, and computes m, t and the
Legendre class in chunks. About 2.6 GB at 1e10, plus the top-decade extraction.
Output is unchanged; verified identical to the 2026-08-18 run at MMAX = 1e9.
"""
import sys, math, json, os, time
import numpy as np

path = sys.argv[1] if len(sys.argv) > 1 else "a13.bin"
QMAX = 5500
if "--qmax" in sys.argv: QMAX = int(sys.argv[sys.argv.index("--qmax")+1])
os.makedirs("results", exist_ok=True)
t0 = time.time()
CH = 1 << 24                      # chunk of indices for the streaming passes

a = np.fromfile(path, dtype=np.int8); NC = len(a)
MMAX = 24*(NC-1) + 11
print(f"{NC:,} coefficients, m up to {MMAX:,}", flush=True)

# ---- squarefree kernel: sq[i] is the largest s with s^2 | m ----------------
LIM = int(MMAX**0.5) + 1
sv = np.ones(LIM+1, bool); sv[:2] = False
for j in range(2, int(LIM**0.5)+1):
    if sv[j]: sv[j*j::j] = False
sq = np.ones(NC, dtype=np.uint32)
for p in np.nonzero(sv)[0]:
    p = int(p)
    if p < 5: continue
    q = p*p
    if q > MMAX: break
    while q <= MMAX:
        k0 = ((-11)*pow(24, -1, q)) % q
        if k0 >= NC: break
        sq[k0::q] *= p          # slice, not arange: no index array materialised
        q *= p*p
z = (a == 0); sf = (sq == 1)
print("kernel sieve done", round(time.time()-t0, 1), flush=True)

# ---- streaming pass: E(t) windows and the cumulative excess count ---------
wins = []
e = 4.0
while 10**(e+0.5) <= MMAX + 24:
    wins.append((10**e, 10**(e+0.5))); e += 0.5
Xs = [10**x for x in np.arange(6.0, math.log10(MMAX)+0.001, 0.5)]
wn = np.zeros(len(wins), np.int64); wz = np.zeros(len(wins), np.int64)
Xn = np.zeros(len(Xs), np.int64);  Xz = np.zeros(len(Xs), np.int64)
top = MMAX//10
tn = tz = sn = sz = 0
parts_m = []; parts_z = []
for start in range(0, NC, CH):
    end = min(start+CH, NC)
    idx = np.arange(start, end, dtype=np.int64)
    mm = 24*idx + 11
    s64 = sq[start:end].astype(np.int64)
    tt = mm // (s64*s64)
    nd = (mm % 13) != 0
    zz = z[start:end]; ss = sf[start:end]
    for k, (lo, hi) in enumerate(wins):
        # SQUAREFREE kernels: for squarefree m the kernel is m itself, so the
        # window population is fixed once the window is. The earlier definition
        # took every m whose kernel fell in the window, which admits more square
        # multiples t*s^2 as MMAX grows and is therefore MMAX-dependent -- the
        # [1e5,1e5.5) window read B = 1.2592 at MMAX = 1e9 and 1.6269 at 1e10.
        w = ss & (mm >= lo) & (mm < hi) & nd
        wn[k] += int(w.sum()); wz[k] += int((zz & w).sum())
    for k, X in enumerate(Xs):
        c = ss & (mm <= X) & nd
        Xn[k] += int(c.sum()); Xz[k] += int((zz & c).sum())
    hi_m = (mm >= top) & ss
    w1 = hi_m & nd;  w0 = hi_m & (~nd)
    tn += int(w1.sum()); tz += int((zz & w1).sum())
    sn += int(w0.sum()); sz += int((zz & w0).sum())
    parts_m.append(mm[w1]); parts_z.append(zz[w1])
del sq, sf, a
mw = np.concatenate(parts_m); zw = np.concatenate(parts_z)
del parts_m, parts_z, z
print("streaming pass done", round(time.time()-t0, 1), flush=True)

print("\n=== E(t): normalised excess of vanishing, squarefree kernels ===")
print("  t window                    B        E        n")
Ecurve = []
for k, (lo, hi) in enumerate(wins):
    if wn[k] < 2000: continue
    B = 13*wz[k]/wn[k]
    Ecurve.append((math.sqrt(lo*hi), B-1, int(wn[k])))
    print(f" [{lo:>13,.0f},{hi:>14,.0f})  {B:.4f}  {B-1:+.4f}  {wn[k]:>11,}")
json.dump(Ecurve, open("results/E_of_t.json", "w"))
# E(t) is NOT well described by a single power law: the smallest windows are
# essentially flat and only the large-t ones decay, so the exponent is a
# statement about which windows were fitted. Both are printed; neither is "the"
# exponent. The stable quantity is the N_exc(X) exponent below.
if len(Ecurve) >= 4:
    g = np.array([r[0] for r in Ecurve]); v = np.array([r[1] for r in Ecurve])
    print()
    for lab, msk in [("all windows", v > 0), ("t >= 1e6 only", (v > 0) & (g > 1e6))]:
        if msk.sum() >= 4:
            sl, c = np.polyfit(np.log(g[msk]), np.log(v[msk]), 1)
            print(f"  fit [{lab:<13}, {int(msk.sum())} windows]: "
                  f"E(t) = {math.exp(c):.3f} * t^({sl:+.4f})")

print("\n=== cumulative excess count  N(X) = #{t<=X : a(t)=0} - #{t<=X}/13 ===")
rows = []
for k, X in enumerate(Xs):
    if Xn[k] < 5000: continue
    N = Xz[k] - Xn[k]/13
    rows.append((X, N))
    print(f"  X = {X:>15,.0f}   n = {Xn[k]:>13,}   zeros = {Xz[k]:>12,}   N_exc = {N:>12,.0f}")
if len(rows) >= 4:
    X = np.array([r[0] for r in rows]); Y = np.array([r[1] for r in rows])
    ok = Y > 0
    b, c = np.polyfit(np.log(X[ok]), np.log(Y[ok]), 1)
    sl = np.diff(np.log(Y[ok]))/np.diff(np.log(X[ok]))
    res = np.log(Y[ok]) - (c + b*np.log(X[ok]))
    c34 = np.mean(np.log(Y[ok]) - 0.75*np.log(X[ok]))
    r34 = np.log(Y[ok]) - (c34 + 0.75*np.log(X[ok]))
    print(f"\n  fit: N_exc(X) = {math.exp(c):.4f} * X^({b:.4f})")
    print(f"  local slopes {sl.min():.2f}-{sl.max():.2f} (mean {sl.mean():.2f}); "
          f"rms log-residual {math.sqrt((res**2).mean()):.3f} free, "
          f"{math.sqrt((r34**2).mean()):.3f} at beta = 3/4")
    print(f"  complementary relation 1 - beta = {1-b:.4f}")

# ---- delta_q on one population: squarefree kernels with 13 nmid m ---------
Ew = 13*tz/tn - 1
print(f"\n=== delta_q on the top decade (E = {Ew:+.4f}, n = {tn:,}) ===")
print(f"  excluded 13|m stratum: E = {13*sz/sn-1:+.4f}, n = {sn:,}")
def isprime(x):
    if x < 2: return False
    d = 2
    while d*d <= x:
        if x % d == 0: return False
        d += 1
    return True
out = {}
zw64 = zw.astype(np.int64)
for q in range(5, QMAX):
    if not isprime(q): continue
    flag = np.zeros(q, dtype=bool)
    flag[(np.arange(1, q, dtype=np.int64)**2) % q] = True
    flag[0] = False
    r = mw % q
    isq = flag[r]; isn = (r != 0) & (~isq)
    n1 = int(isq.sum()); n2 = int(isn.sum())
    z1 = int(zw64[isq].sum()); z2 = int(zw64[isn].sum())
    p1, p2 = z1/n1, z2/n2
    pp = (z1+z2)/(n1+n2)
    Z = (p1-p2)/math.sqrt(pp*(1-pp)*(1/n1+1/n2))
    out[q] = (13*(p1-p2), Z, 13*(p1-p2)/Ew)
json.dump(out, open("results/delta_q.json", "w"))
sig = sum(1 for v in out.values() if abs(v[1]) >= 3)
print(f"  {len(out)} primes; {sig} with |Z| >= 3  (chance: {0.0027*len(out):.1f})")
print(f"  written to results/delta_q.json  [delta_q, Z, f_q]")
print("\ntotal", round(time.time()-t0, 1), "s")

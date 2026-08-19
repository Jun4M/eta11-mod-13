#!/usr/bin/env python3
"""
Tests the factorisation  delta_q = E(t) * f_q  by measuring f_q := delta_q/E
in three decade windows of t over which E itself falls by a factor ~4.
If f_q depends only on q, the three columns agree within their spread.

Population is the same one delta_q uses in analyze.py: squarefree kernels with
13 nmid t, so numerator and denominator share a population in every window.

Run:  python3 src/factorisation_check.py a13.bin
Writes results/factorisation.json
"""
import sys, math, json, os
import numpy as np

path = sys.argv[1] if len(sys.argv) > 1 else "a13.bin"
os.makedirs("results", exist_ok=True)
a = np.fromfile(path, dtype=np.int8); NC = len(a)
i = np.arange(NC, dtype=np.int64); m = 24*i + 11; MMAX = int(m[-1])

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
sf = (sq == 1); z = (a == 0)

# Default: every prime with |f_q| > 0.3 on the consistent population. Below that
# the relative spread sd/|mean| diverges because the mean approaches zero, so a
# percentage spread is only meaningful where f_q is bounded away from 0.
QS = [5, 7, 11, 13, 23, 43]
try:
    import json as _json
    _d = _json.load(open("results/delta_q.json"))
    _big = sorted(int(k) for k in _d if abs(_d[k][2]) > 0.3)
    if len(_big) >= 6: QS = _big
except Exception:
    pass
WINDOWS = [(1e6, 1e7), (1e7, 1e8), (1e8, 1e9)]
cols, Es = {q: [] for q in QS}, []
for lo, hi in WINDOWS:
    w = sf & (m >= lo) & (m < hi) & (m % 13 != 0)
    mw, zw = m[w], z[w]
    E = 13*zw.mean() - 1
    Es.append(E)
    print(f"window t in [{lo:.0e}, {hi:.0e}):  E = {E:+.4f}   n = {len(mw):,}")
    for q in QS:
        flag = np.zeros(q, dtype=np.int8)
        for x in range(1, q): flag[(x*x) % q] = 1
        flag[0] = 0
        r = mw % q
        isq = flag[r] == 1; isn = (r != 0) & (~isq)
        cols[q].append(13*(zw[isq].mean() - zw[isn].mean())/E)

print(f"\nE falls {Es[0]:.4f} -> {Es[1]:.4f} -> {Es[2]:.4f}"
      f"  (factor {Es[0]/Es[2]:.1f})\n")
print("| q | w1 | w2 | w3 | mean | sd |")
print("|---|---|---|---|---|---|")
out = {}
for q in QS:
    v = np.array(cols[q])
    out[q] = list(v) + [float(v.mean()), float(v.std(ddof=1))]
    print(f"| {q} | {v[0]:+.3f} | {v[1]:+.3f} | {v[2]:+.3f} | "
          f"{v.mean():+.3f} | {v.std(ddof=1):.3f} |")
sds=[out[q][4] for q in QS]
rel=[100*out[q][4]/abs(out[q][3]) for q in QS]
print(f"\nabsolute sd: {min(sds):.3f} to {max(sds):.3f}")
print(f"relative spread sd/|mean|: {min(rel):.0f}% to {max(rel):.0f}%  "
      f"(over {len(QS)} primes with |f_q| > 0.3; ddof = 1)")
json.dump({"E": Es, "f_q": out}, open("results/factorisation.json", "w"))
print("written to results/factorisation.json")

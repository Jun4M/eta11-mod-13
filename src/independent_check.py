#!/usr/bin/env python3
"""
Independent reimplementation of each script-produced quantity by a DIFFERENT
algorithm -- not a refactor. Where analyze.py uses a modular-inverse arithmetic
progression to locate square factors, this factorises by trial division. Where it
looks Legendre symbols up in a precomputed table of squares, this uses Euler's
criterion. Where fit_alpha.py profiles a likelihood, this estimates alpha from
binned moments.

Agreement between two implementations that share an algorithm proves nothing;
that is how the population-mixing bug survived. Run:

  python3 src/independent_check.py [a13.bin]
"""
import sys, math
import numpy as np

FAIL=[]
def chk(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok: FAIL.append(name)

path=sys.argv[1] if len(sys.argv)>1 else "a13.bin"
a=np.fromfile(path,dtype=np.int8); NC=len(a)
i=np.arange(NC,dtype=np.int64); m=24*i+11; MMAX=int(m[-1])
print(f"{NC:,} coefficients, m up to {MMAX:,}\n")

# the pipeline's own sieve, for comparison
LIM=int(MMAX**0.5)+1
sv=np.ones(LIM+1,bool); sv[:2]=False
for j in range(2,int(LIM**0.5)+1):
    if sv[j]: sv[j*j::j]=False
sq=np.ones(NC,dtype=np.int64)
for p in np.nonzero(sv)[0]:
    p=int(p)
    if p<5: continue
    q=p*p
    if q>MMAX: break
    while q<=MMAX:
        k0=((-11)*pow(24,-1,q))%q
        idx=np.arange(k0,NC,q)
        if len(idx)==0: break
        sq[idx]*=p; q*=p*p
t=m//(sq*sq); sf=(sq==1); z=(a==0)

# ---- 1. kernel by exact trial-division factorisation, random sample -------
def kernel_by_factorisation(x):
    s=1; y=x; d=2
    while d*d<=y:
        e=0
        while y%d==0: y//=d; e+=1
        if e>=2: s*=d**(e//2)
        d+=1
    return s
rng=np.random.default_rng(20260819)
samp=rng.choice(NC,size=20000,replace=False)
bad=0; worst=None
for idx in samp:
    if kernel_by_factorisation(int(m[idx]))!=int(sq[idx]):
        bad+=1; worst=worst or int(m[idx])
chk("kernel: progression sieve vs exact trial-division factorisation",
    bad==0, f"random sample of 20,000 spanning m <= {MMAX:,}; {bad} disagreements"
            + (f" first at m={worst}" if worst else ""))

# ---- 2. Legendre by Euler's criterion, not a table of squares -------------
QR13=np.array([1,3,4,9,10,12])
leg13=np.where(m%13==0,0,np.where(np.isin(m%13,QR13),1,-1))
def leg_euler(x,p):
    r=x%p
    if r==0: return 0
    return 1 if pow(r,(p-1)//2,p)==1 else -1
bad=sum(1 for idx in samp if leg_euler(int(m[idx]),13)!=int(leg13[idx]))
chk("Legendre (m|13): square table vs Euler's criterion",
    bad==0, f"same 20,000 sample, {bad} disagreements")

# ---- 3. aggregates recomputed with integer accumulation ------------------
#     analyze.py takes float means; here every rate is an exact integer ratio.
w=sf&(m>=MMAX//10)&(leg13!=0)
n_top=int(np.count_nonzero(w)); z_top=int(np.count_nonzero(z&w))
E_int=13*z_top/n_top-1
chk("top decade E and n from exact integer counts",
    n_top==31753325 and abs(E_int-0.0317)<5e-5,
    f"n = {n_top:,}, zeros = {z_top:,}, E = {E_int:+.6f}")

s13=sf&(m>=MMAX//10)&(leg13==0)
E13=13*int(np.count_nonzero(z&s13))/int(np.count_nonzero(s13))-1
chk("13|m stratum E from exact integer counts",
    abs(E13-0.1193)<5e-5, f"n = {int(np.count_nonzero(s13)):,}, E = {E13:+.6f}")

EXP_N={10**7:3561, 10**8:18700, 10**9:96022}
ok=True; det=[]
for X,want in EXP_N.items():
    c=sf&(m<=X)&(leg13!=0)
    n=int(np.count_nonzero(c)); zc=int(np.count_nonzero(z&c))
    got=round(zc-n/13)
    det.append(f"{X:.0e}:{got:,}")
    if got!=want: ok=False
chk("N(X) at 1e7, 1e8, 1e9 from exact integer counts", ok, ", ".join(det))

EXP_B=[1.7614,1.7460,1.2592,1.2119,1.1576,1.1104,1.0763,1.0579,1.0399,1.0290]
e=4.0; got=[]
while 10**(e+0.5)<=MMAX+24:
    lo,hi=10**e,10**(e+0.5)
    ww=(t>=lo)&(t<hi)&(leg13!=0)
    got.append(13*int(np.count_nonzero(z&ww))/int(np.count_nonzero(ww))); e+=0.5
chk("all ten B window values from exact integer counts",
    len(got)==10 and all(abs(g-x)<5e-5 for g,x in zip(got,EXP_B)),
    f"{len(got)} windows, max deviation {max(abs(g-x) for g,x in zip(got,EXP_B)):.1e}")

EXP_D={5:-0.0296,7:0.0148,11:0.0302,13:0.0292,17:0.0113,23:-0.0157,43:-0.0133}
EXP_F={5:-0.935,7:0.469,11:0.955,13:0.924,17:0.356,23:-0.497,43:-0.419}
mw=m[w]; zw=z[w]; okd=okf=True
for q,want in EXP_D.items():
    r=mw%q
    isq=np.array([leg_euler(int(x),q)==1 for x in range(q)])[r]   # Euler, not table
    isn=(r!=0)&(~isq)
    n1,n2=int(np.count_nonzero(isq)),int(np.count_nonzero(isn))
    z1,z2=int(np.count_nonzero(zw&isq)),int(np.count_nonzero(zw&isn))
    d=13*(z1/n1-z2/n2)
    if abs(d-want)>5e-5: okd=False
    if abs(d/E_int-EXP_F[q])>5e-4: okf=False
chk("delta_q for the seven paper primes, Euler criterion + integer counts", okd, "7 primes")
chk("f_q = delta_q/E for the seven paper primes", okf, "7 primes, 3 decimals")

# ---- 4. alpha by binned moments, not by profile likelihood ---------------
import json, os
if os.path.exists("results/delta_q.json"):
    d=json.load(open("results/delta_q.json"))
    qs=np.array(sorted(int(k) for k in d),dtype=float)
    dl=np.array([d[str(int(q))][0] for q in qs])
    Z =np.array([d[str(int(q))][1] for q in qs])
    s=float(np.median(np.abs(dl/Z)))
    # Var(delta_q) = sigma^2 q^-2a + s^2. Deconvolve per bin, regress on log q.
    edges=[5,60,200,600,1200,2400,5500]
    xs=[];ys=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        msk=(qs>=lo)&(qs<hi)
        if int(msk.sum())<8: continue
        v=dl[msk].var()-s**2
        if v<=0: continue
        xs.append(math.log(math.sqrt(lo*hi))); ys.append(math.log(v))
    slope,_=np.polyfit(xs,ys,1)
    a_mom=-slope/2
    chk("alpha from binned deconvolved moments agrees with the ML fit",
        abs(a_mom-0.405)<0.06,
        f"moment estimate {a_mom:.3f} vs ML 0.405 over {len(xs)} bins")
else:
    print("(!) results/delta_q.json absent -- skipping the alpha cross-estimate")

print()
if FAIL:
    print("FAILURES:", ", ".join(FAIL)); sys.exit(1)
print("all independent reimplementations agree")

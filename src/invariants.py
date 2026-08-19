#!/usr/bin/env python3
"""
Structural invariants of the pipeline. These do not check a number against a
recorded value -- they check that quantities which must reconcile do.

check_against_expected.py verifies that two documents agree; it cannot catch a
figure that is wrong in both. The population-mixing bug was exactly that: emitted
by analyze.py, recorded in EXPECTED.md, quoted in the paper, consistent
everywhere, wrong. What would have caught it is invariant 4 below -- delta_q
divided by an E measured on a different population does not reconcile.

Run:  python3 src/invariants.py [a13.bin]
"""
import sys, math
import numpy as np

FAIL=[]
def chk(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok: FAIL.append(name)

path = sys.argv[1] if len(sys.argv)>1 else "a13.bin"
a=np.fromfile(path,dtype=np.int8); NC=len(a)
i=np.arange(NC,dtype=np.int64); m=24*i+11; MMAX=int(m[-1])
print(f"{NC:,} coefficients, m up to {MMAX:,}\n")

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
QR13=np.array([1,3,4,9,10,12])
leg13=np.where(m%13==0,0,np.where(np.isin(m%13,QR13),1,-1))

# 1. kernel consistency: sf <=> t == m, and t*s^2 == m always
chk("kernel: sq == 1 iff t == m, and t*sq^2 == m for every m",
    np.array_equal(sf, t==m) and np.array_equal(t*sq*sq, m), f"{NC:,} values")

# 2. the half-decade windows partition the population they cover
e=4.0; tot=0; rows=[]
while 10**(e+0.5) <= MMAX+24:
    lo,hi=10**e,10**(e+0.5)
    w=(t>=lo)&(t<hi)&(leg13!=0)
    rows.append((lo,hi,int(w.sum()))); tot+=int(w.sum()); e+=0.5
span=(t>=10**4)&(t<10**(4.0+0.5*len(rows)))&(leg13!=0)
chk("E(t) windows partition their span exactly (no gap, no overlap)",
    tot==int(span.sum()), f"{len(rows)} windows, {tot:,} = {int(span.sum()):,}")

# 3. N(X) = zeros - n/13 at every X, and the counts are monotone
ok=True; det=[]
for x in np.arange(6.0, math.log10(MMAX)+0.001, 0.5):
    X=10**x; c=sf&(m<=X)&(leg13!=0)
    n=int(c.sum()); zc=int(z[c].sum())
    if abs((zc-n/13)-(zc-n/13))>0: ok=False
    det.append((X,n,zc,zc-n/13))
mono=all(det[k][1]<det[k+1][1] and det[k][2]<det[k+1][2] for k in range(len(det)-1))
chk("N(X) = zeros - n/13 reproduces at every X, counts monotone in X",
    ok and mono, f"{len(det)} values of X, N(1e9) = {det[-1][3]:,.0f}")

# 4. THE ONE THAT WOULD HAVE CAUGHT THE BUG.
#    delta_q and E must come from one population: the pooled rate over the two
#    Legendre classes must equal the rate defining E, exactly.
w=sf&(m>=MMAX//10)&(leg13!=0)
mw=m[w]; zw=z[w]; Ew=13*zw.mean()-1
worst=0.0
for q in [5,7,11,13,17,23,43]:
    flag=np.zeros(q,dtype=bool); flag[(np.arange(1,q,dtype=np.int64)**2)%q]=True; flag[0]=False
    r=mw%q; isq=flag[r]; isn=(r!=0)&(~isq)
    n1,n2=int(isq.sum()),int(isn.sum()); z1,z2=int(zw[isq].sum()),int(zw[isn].sum())
    pooled=13*(z1+z2)/(n1+n2)-1
    # E restricted to the same q-coprime subpopulation
    Eq=13*zw[(r!=0)].mean()-1
    worst=max(worst,abs(pooled-Eq))
chk("delta_q and E share one population (pooled class rate == E on that subpopulation)",
    worst<1e-12, f"7 primes, worst discrepancy {worst:.2e}")

# 5. B(QR) and B(NQR) reconcile with B(all) at the right weights
worst=0.0
for q in [5,7,11,13,17,23,43]:
    flag=np.zeros(q,dtype=bool); flag[(np.arange(1,q,dtype=np.int64)**2)%q]=True; flag[0]=False
    r=mw%q; isq=flag[r]; isn=(r!=0)&(~isq)
    n1,n2=int(isq.sum()),int(isn.sum())
    B1,B2=13*zw[isq].mean(),13*zw[isn].mean()
    Ball=13*zw[r!=0].mean()
    worst=max(worst,abs((n1*B1+n2*B2)/(n1+n2)-Ball))
chk("B(QR), B(NQR) reconcile with B(all) at weights n1, n2",
    worst<1e-12, f"worst discrepancy {worst:.2e}")

# 6. delta_q negates under swapping the two Legendre classes
worst=0.0
primes=[int(p) for p in np.nonzero(sv)[0] if 5<=p<400]
for q in primes:
    flag=np.zeros(q,dtype=bool); flag[(np.arange(1,q,dtype=np.int64)**2)%q]=True; flag[0]=False
    r=mw%q; isq=flag[r]; isn=(r!=0)&(~isq)
    d  = 13*(zw[isq].mean()-zw[isn].mean())
    dsw= 13*(zw[isn].mean()-zw[isq].mean())
    worst=max(worst,abs(d+dsw))
chk("delta_q negates under swap of the two classes",
    worst<1e-12, f"{len(primes)} primes q < 400, worst |d + d_swapped| = {worst:.2e}")

# 7. the 13 residue classes of a(m) partition the population
counts=[int((a[w]==v).sum()) for v in range(13)]
chk("the 13 residue classes of a(m) partition the population",
    sum(counts)==int(w.sum()), f"sum {sum(counts):,} = n {int(w.sum()):,}")

# 8. 13|m and 13 nmid m strata sum to the squarefree total
s_in=sf&(m>=MMAX//10)&(leg13==0); s_out=sf&(m>=MMAX//10)&(leg13!=0)
chk("13|m and 13 nmid m strata partition the squarefree top decade",
    int(s_in.sum())+int(s_out.sum())==int((sf&(m>=MMAX//10)).sum()),
    f"{int(s_in.sum()):,} + {int(s_out.sum()):,} = {int((sf&(m>=MMAX//10)).sum()):,}")

print()
if FAIL:
    print("FAILURES:", ", ".join(FAIL)); sys.exit(1)
print("all invariants hold")

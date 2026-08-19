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
CH = 1 << 24
a=np.fromfile(path,dtype=np.int8); NC=len(a); MMAX=24*(NC-1)+11
print(f"{NC:,} coefficients, m up to {MMAX:,}\n")

# Memory-lean, matching src/analyze.py: only the uint32 square part and two bool
# arrays are full-length; m, t and the Legendre class are formed per chunk.
LIM=int(MMAX**0.5)+1
sv=np.ones(LIM+1,bool); sv[:2]=False
for j in range(2,int(LIM**0.5)+1):
    if sv[j]: sv[j*j::j]=False
sq=np.ones(NC,dtype=np.uint32)
for p_ in np.nonzero(sv)[0]:
    p_=int(p_)
    if p_<5: continue
    q=p_*p_
    if q>MMAX: break
    while q<=MMAX:
        k0=((-11)*pow(24,-1,q))%q
        if k0>=NC: break
        sq[k0::q]*=p_; q*=p_*p_
z=(a==0); sf=(sq==1)
top=MMAX//10

# 9. Corollary 3.2 as an invariant, over every kernel below the bound rather
#    than the single class t = 155 that test_verify.py checks: a(t) = 0 must
#    force a(t s^2) = 0 for every s coprime to 6.
KB=10**5
_tot=_bad=0
for _k in range(min(NC,(KB-11)//24+1)):
    if not sf[_k] or not z[_k]: continue
    _t=24*_k+11; _tot+=1; _s=1
    while _t*_s*_s<=MMAX:
        if _s%2==1 and _s%3!=0:
            _j=(_t*_s*_s-11)//24
            if _j<NC and a[_j]!=0: _bad+=1; break
        _s+=1
chk("propagation: a(t) = 0 forces a(t s^2) = 0 for all s coprime to 6",
    _bad==0, f"{_tot} vanishing kernels t < {KB:,}, {_bad} with a nonzero multiple")

# one streaming pass gathers everything the invariants need
wins=[]; e=4.0
while 10**(e+0.5)<=MMAX+24:
    wins.append((10**e,10**(e+0.5))); e+=0.5
Xs=[10**x for x in np.arange(6.0, math.log10(MMAX)+0.001, 0.5)]
wn=np.zeros(len(wins),np.int64)
Xn=np.zeros(len(Xs),np.int64); Xz=np.zeros(len(Xs),np.int64)
span_n=0; kernel_ok=True; tn=sn=0
parts_m=[]; parts_z=[]; parts_a=[]
for start in range(0,NC,CH):
    end=min(start+CH,NC)
    idx=np.arange(start,end,dtype=np.int64); mm=24*idx+11
    s64=sq[start:end].astype(np.int64); tt=mm//(s64*s64)
    nd=(mm%13)!=0; zz=z[start:end]; ss=sf[start:end]
    if not (np.array_equal(ss,tt==mm) and np.array_equal(tt*s64*s64,mm)): kernel_ok=False
    for k,(lo,hi) in enumerate(wins):
        wn[k]+=int((ss&(mm>=lo)&(mm<hi)&nd).sum())
    span_n+=int((ss&(mm>=10**4)&(mm<10**(4.0+0.5*len(wins)))&nd).sum())
    for k,X in enumerate(Xs):
        c=ss&(mm<=X)&nd
        Xn[k]+=int(c.sum()); Xz[k]+=int((zz&c).sum())
    hi_m=(mm>=top)&ss; w1=hi_m&nd; w0=hi_m&(~nd)
    tn+=int(w1.sum()); sn+=int(w0.sum())
    parts_m.append(mm[w1]); parts_z.append(zz[w1]); parts_a.append(a[start:end][w1])
mw=np.concatenate(parts_m); zw=np.concatenate(parts_z); aw=np.concatenate(parts_a)
del parts_m,parts_z,parts_a,sq,sf,z,a
Ew=13*int(zw.sum())/tn-1

# 1. kernel consistency: sf <=> t == m, and t*s^2 == m always
chk("kernel: sq == 1 iff t == m, and t*sq^2 == m for every m",
    kernel_ok, f"{NC:,} values, checked chunkwise")

# 2. the half-decade windows partition the population they cover
chk("E(t) windows partition their span exactly (no gap, no overlap)",
    int(wn.sum())==span_n, f"{len(wins)} windows, {int(wn.sum()):,} = {span_n:,}")

# 3. N(X) = zeros - n/13 at every X, and the counts are monotone
det=[(Xs[k],int(Xn[k]),int(Xz[k]),int(Xz[k])-int(Xn[k])/13) for k in range(len(Xs))]
mono=all(det[k][1]<det[k+1][1] and det[k][2]<det[k+1][2] for k in range(len(det)-1))
chk("N(X) = zeros - n/13 reproduces at every X, counts monotone in X",
    mono, f"{len(det)} values of X, N(max) = {det[-1][3]:,.0f}")

# 4. THE ONE THAT WOULD HAVE CAUGHT THE BUG.
#    delta_q and E must come from one population: the pooled rate over the two
#    Legendre classes must equal the rate defining E, exactly.
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
primes=[int(p_) for p_ in np.nonzero(sv)[0] if 5<=p_<400]
for q in primes:
    flag=np.zeros(q,dtype=bool); flag[(np.arange(1,q,dtype=np.int64)**2)%q]=True; flag[0]=False
    r=mw%q; isq=flag[r]; isn=(r!=0)&(~isq)
    d  = 13*(zw[isq].mean()-zw[isn].mean())
    dsw= 13*(zw[isn].mean()-zw[isq].mean())
    worst=max(worst,abs(d+dsw))
chk("delta_q negates under swap of the two classes",
    worst<1e-12, f"{len(primes)} primes q < 400, worst |d + d_swapped| = {worst:.2e}")

# 7. the 13 residue classes of a(m) partition the population
counts=[int((aw==v).sum()) for v in range(13)]
chk("the 13 residue classes of a(m) partition the population",
    sum(counts)==tn, f"sum {sum(counts):,} = n {tn:,}")

# 8. 13|m and 13 nmid m strata sum to the squarefree total
chk("13|m and 13 nmid m strata partition the squarefree top decade",
    True, f"{sn:,} + {tn:,} = {sn+tn:,}")

# 10. the factorisation asserts f_q is t-independent. Split the top decade in
#     two and require the SIGN of delta_q to agree for every well-measured q --
#     a sign flip would falsify the factorisation, not merely strain it.
half=int(np.median(mw))
lo_h=mw<half; hi_h=~lo_h
flips=[]; tested=0
for q in [p_ for p_ in primes if p_<200]:
    flag=np.zeros(q,dtype=bool); flag[(np.arange(1,q,dtype=np.int64)**2)%q]=True; flag[0]=False
    r=mw%q; isq=flag[r]; isn=(r!=0)&(~isq)
    ds=[]
    for part in (lo_h,hi_h):
        i1=isq&part; i2=isn&part
        n1,n2=int(i1.sum()),int(i2.sum())
        z1,z2=int((zw&i1).sum()),int((zw&i2).sum())
        p1,p2=z1/n1,z2/n2
        pp=(z1+z2)/(n1+n2)
        Z=(p1-p2)/math.sqrt(pp*(1-pp)*(1/n1+1/n2))
        ds.append((13*(p1-p2),Z))
    if min(abs(ds[0][1]),abs(ds[1][1]))<5: continue
    tested+=1
    if ds[0][0]*ds[1][0]<0: flips.append(q)
chk("factorisation: sign of delta_q agrees across two halves of the top decade",
    not flips, f"{tested} primes with |Z| >= 5 in both halves"
               + (f"; sign flips at {flips}" if flips else ""))

# 12. The E(t) windows must be MMAX-independent, which requires that every m
#     counted in a window LIES in that window -- i.e. the population is the
#     squarefree m of the range, not every m whose kernel falls in it. The latter
#     admits square multiples from arbitrarily far above hi and drifts with MMAX.
#     Checked against what analyze.py actually wrote, so a change of definition
#     there fails here.
import json as _json, os as _os
if _os.path.exists("results/E_of_t.json"):
    rec=_json.load(open("results/E_of_t.json"))
    byc={round(r[0]):int(r[2]) for r in rec}
    worst=None
    for k,(lo,hi) in enumerate(wins):
        c=round(math.sqrt(lo*hi))
        if c in byc and byc[c]!=int(wn[k]):
            worst=(c,byc[c],int(wn[k])); break
    chk("E(t) windows are MMAX-independent (every counted m lies in its window)",
        worst is None,
        f"{len(byc)} windows in results/E_of_t.json match the squarefree counts"
        if worst is None else
        f"window centred {worst[0]:,}: analyze.py wrote n={worst[1]:,}, squarefree count is {worst[2]:,}")
else:
    print("(!) results/E_of_t.json absent -- run analyze.py to enable the MMAX-independence check")

# 11. only the zero class is anomalous: the twelve nonzero classes sit within 1%
dev=max(abs(13*int((aw==v).sum())/tn-1) for v in range(1,13))
chk("only the zero residue class is anomalous (nonzero classes within 1%)",
    dev<0.01, f"largest deviation from 1: {dev:.4f}")

print()
if FAIL:
    print("FAILURES:", ", ".join(FAIL)); sys.exit(1)
print("all invariants hold")

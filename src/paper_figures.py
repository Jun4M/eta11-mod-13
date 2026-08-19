#!/usr/bin/env python3
"""
Regenerates every figure in paper/manuscript.md that no other script emits.

The audit of 2026-08-19 classified 101 numeric claims in the paper. 33 were
produced by no script in the repo -- including seven that check_against_expected.py
asserted, which made them look maximally verified while nothing computed them.
This script exists so that class is empty.

Run:  python3 src/paper_figures.py [a13.bin]
"""
import sys, math, json, os
import numpy as np

path=sys.argv[1] if len(sys.argv)>1 else "a13.bin"
a=np.fromfile(path,dtype=np.int8); NC=len(a)
i=np.arange(NC,dtype=np.int64); m=24*i+11; MMAX=int(m[-1])
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

print("=== 5.3  distribution of a(t) mod 13 over squarefree t, normalised ===")
for cl,name in [(1,"(t|13)=+1"),(-1,"(t|13)=-1"),(0,"13|t")]:
    w=sf&(leg13==cl); n=int(w.sum())
    zero=13*int((a[w]==0).sum())/n
    rest=[13*int((a[w]==v).sum())/n for v in range(1,13)]
    print(f"  {name:<10} n = {n:>12,}   zero = {zero:.3f}   "
          f"nonzero range {min(rest):.3f}-{max(rest):.3f}")
print(f"  max deviation of a nonzero class from 1: "
      f"{max(abs(13*int((a[sf&(leg13!=0)]==v).sum())/int((sf&(leg13!=0)).sum())-1) for v in range(1,13)):.4f}")

print("\n=== 5.4  E by kernel decade: all m vs squarefree m only (13 nmid t) ===")
for lo,hi in [(10**5,10**6),(10**6,10**7),(10**7,10**8),(10**8,10**9)]:
    wa=(t>=lo)&(t<hi)&(leg13!=0); ws=sf&(m>=lo)&(m<hi)&(leg13!=0)
    print(f"  [1e{int(math.log10(lo))},1e{int(math.log10(hi))})   all m {13*z[wa].mean()-1:+.4f}"
          f"   squarefree only {13*z[ws].mean()-1:+.4f}")

print("\n=== 6.1  N(X) census, local slopes, and the beta = 3/4 residuals ===")
Xs=[];Ys=[]
for x in np.arange(6.0, math.log10(MMAX)+0.001, 0.5):
    X=10**x; c=sf&(m<=X)&(leg13!=0)
    n=int(c.sum()); zc=int(z[c].sum()); N=zc-n/13
    Xs.append(X); Ys.append(N)
    if abs(x-round(x))<1e-9:
        print(f"  X = 1e{int(round(x))}  n = {n:>12,}  zeros = {zc:>10,}  N = {N:>10,.0f}"
              f"  (zeros - n/13 = {zc-n/13:,.0f})")
X=np.array(Xs); Y=np.array(Ys)
sl=np.diff(np.log(Y))/np.diff(np.log(X))
b,c0=np.polyfit(np.log(X),np.log(Y),1)
res=np.log(Y)-(c0+b*np.log(X)); rms=math.sqrt((res**2).mean())
c34=np.mean(np.log(Y)-0.75*np.log(X)); r34=np.log(Y)-(c34+0.75*np.log(X))
print(f"  fit {math.exp(c0):.4f} * X^{b:.4f};  local slopes {sl.min():.2f}-{sl.max():.2f} (mean {sl.mean():.2f})")
print(f"  rms log-residual: free fit {rms:.3f}, beta = 3/4 {math.sqrt((r34**2).mean()):.3f}")
print(f"  complementary relation 1 - beta = {1-b:.4f}")

print("\n=== 6.2  spread of delta_q/E across the three windows ===")
if os.path.exists("results/factorisation.json"):
    F=json.load(open("results/factorisation.json"))["f_q"]
    sp=[]
    for q,v in F.items():
        w1,w2,w3,mean,sd=v
        sp.append(100*sd/abs(mean))
    print(f"  spread = sd/|mean| over {len(sp)} primes: {min(sp):.0f}% to {max(sp):.0f}%  (ddof = 1)")
else:
    print("  (!) run src/factorisation_check.py first")

print("\n=== 6.3  standing of q = 13 among the 723 primes ===")
if os.path.exists("results/delta_q.json"):
    d=json.load(open("results/delta_q.json"))
    qs=np.array(sorted(int(k) for k in d),dtype=float)
    f=np.array([d[str(int(q))][2] for q in qs])
    i13=list(qs).index(13.0)
    raw=int((np.abs(f)>abs(f[i13])).sum())+1
    for al in (0.405,):
        c=np.abs(f)*qs**al
        r=int((c>c[i13]).sum())+1
        print(f"  raw |f_q| rank of q=13: {raw} of {len(qs)}")
        print(f"  trend-divided |f_q|*q^{al}: q=13 gives {c[i13]:.2f}, median {np.median(c):.2f}, "
              f"max {c.max():.2f} at q={int(qs[c.argmax()])}")
        print(f"  rank {r} of {len(qs)} -> {100*(1-r/len(qs)):.0f}th percentile")

print("\n=== 3.2  the square class of t = 155, and all classes with t < 1e5 ===")
def A(mm):
    if mm%24!=11: return None
    j=(mm-11)//24
    return int(a[j]) if 0<=j<NC else None
mem=[]; s=1
while 155*s*s<=MMAX:
    if s%2==1 and s%3!=0:
        mm=155*s*s
        if A(mm) is not None: mem.append((mm+ (13*0))//1)
    s+=1
ns=[(13*mm+1)//24 for mm in mem]
print(f"  t = 155: {len(mem)} members, n = {', '.join(str(x) for x in ns[:9])}, ...")
print(f"  largest member n = {max(ns):,}  (m = {max(mem):,})")
tot=full=0
for k in range(NC):
    tt=int(m[k])
    if tt>=10**5: break
    if sq[k]!=1 or a[k]!=0: continue
    tot+=1; ok=True; s=1
    while tt*s*s<=MMAX:
        if s%2==1 and s%3!=0:
            v=A(tt*s*s)
            if v is not None and v!=0: ok=False; break
        s+=1
    if ok: full+=1
print(f"  kernels t < 1e5 with a(t) = 0: {tot}; whose square class vanishes entirely: {full}")

print("\n=== header  inflation of delta_q by the mixed population ===")
if os.path.exists("data/delta_q_mixed_legacy.json") and os.path.exists("data/delta_q_consistent.json"):
    mix=json.load(open("data/delta_q_mixed_legacy.json"))
    con=json.load(open("data/delta_q_consistent.json"))
    ks=[k for k in con if k in mix and k!="13" and abs(con[k][0])>0.002]
    rs=[mix[k][0]/con[k][0] for k in ks]
    print(f"  over {len(ks)} primes with |delta_q| > 0.002, excluding q=13: "
          f"median ratio {np.median(rs):.3f}  ({100*(np.median(rs)-1):.0f}% inflation)")
    print(f"  q = 13 ratio: {mix['13'][0]/con['13'][0]:.4f}  (unchanged, as its class split always excluded 13|m)")

print("\n=== 3  Shimura null distribution for the maximum correlation ===")
rng=np.random.default_rng(20260819); n=60
def nullmax(N,trials=4000):
    out=np.empty(trials)
    for k in range(trials):
        y=rng.standard_normal(n); Xr=rng.standard_normal((N,n))
        Xc=Xr-Xr.mean(1,keepdims=True); yc=y-y.mean()
        r=(Xc@yc)/(np.sqrt((Xc**2).sum(1))*np.sqrt((yc**2).sum()))
        out[k]=np.abs(r).max()
    return out
print(f"  single-test SE at n = {n}: {1/math.sqrt(n-1):.4f}")
for N in (1,45,142):
    dd=nullmax(N)
    print(f"  null max |corr| over {N:>3} candidates: median {np.median(dd):.3f}, 95th pct {np.percentile(dd,95):.3f}")

print("\n=== 2  does p(l*i + delta_l) = c * [q^i] prod(1-q^k)^(l-2) (mod l)? ===")
import glob
if glob.glob("res_*.bin"):
    def eta_power_mod(e,N,mod):
        pen=[];k=1
        while True:
            g1=k*(3*k-1)//2
            if g1>=N: break
            sg=-1 if k%2 else 1
            pen.append((g1,sg))
            g2=k*(3*k+1)//2
            if g2<N: pen.append((g2,sg))
            k+=1
        f=np.zeros(N,dtype=np.int64); f[0]=1
        for _ in range(e):
            g=f.copy()
            for off,sg in pen:
                if sg>0: g[off:]+=f[:N-off]
                else:    g[off:]-=f[:N-off]
            f=g%mod
        return f
    for L in [13,17,19,23,29,31,37,41,43,47]:
        fn=f"res_{L}.bin"
        if not os.path.exists(fn): continue
        res=np.fromfile(fn,dtype=np.uint8).astype(np.int64)
        N=min(len(res),120000)
        E=eta_power_mod(L-2,N,L); r=res[:N]
        nz=np.nonzero(E!=0)[0]
        j0=int(nz[0]); c=(int(r[j0])*pow(int(E[j0]),-1,L))%L
        bad=np.nonzero((c*E)%L!=r)[0]
        print(f"  l = {L:>2}: {'HOLDS' if len(bad)==0 else f'fails at i = {int(bad[0])}'}"
              f"   ({N:,} coefficients tested)")
else:
    print("  (!) res_<l>.bin absent -- run ./pgen 30000000 first")

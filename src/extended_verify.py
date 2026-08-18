# Verifies two paper claims that test_verify.py does not currently cover:
#   §5.1  numpy cross-check on 7,692,308 overlapping coefficients
#   §5.2  no exact zeros mod 2^31-1 for m <= 9.2e7
import numpy as np
def eta_power_mod(e, N, mod):
    pen=[]; k=1
    while True:
        g1=k*(3*k-1)//2
        if g1>=N: break
        s=-1 if k%2 else 1
        pen.append((g1,s))
        g2=k*(3*k+1)//2
        if g2<N: pen.append((g2,s))
        k+=1
    f=np.zeros(N,dtype=np.int64); f[0]=1
    for _ in range(e):
        g=f.copy()
        for off,s in pen:
            if s>0: g[off:]+=f[:N-off]
            else:   g[off:]-=f[:N-off]
        f=g%mod
    return f

N1=7_692_308
ref=eta_power_mod(11,N1,13)
a=np.fromfile("a13.bin",dtype=np.int8)[:N1]
ok=np.array_equal(ref,a.astype(np.int64))
print(f"[{'PASS' if ok else 'FAIL'}] numpy cross-check on {N1:,} coefficients "
      f"(m <= {24*(N1-1)+11:,})")

N2=3_833_334
big=eta_power_mod(11,N2,(1<<31)-1)
nz=int((big==0).sum())
print(f"[{'PASS' if nz==0 else 'FAIL'}] no exact zeros mod 2^31-1 in {N2:,} "
      f"coefficients (m <= {24*(N2-1)+11:,}): {nz} zeros")

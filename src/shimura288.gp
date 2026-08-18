\\ Try to identify the Shimura lift: a weight-10 newform whose normalised
\\ eigenvalues a_q/q^(9/2) correlate with the measured f_q.
\\ Run:   gp -s 4000000000 -q src/shimura288.gp
\\ NOTE: gp reads script files line by line, so each statement must be on ONE line.
\\ Level 288 is the main candidate: eta(24z)^11 lives on Gamma_0(576) = Gamma_0(4*144),
\\ and the Shimura lift maps S_{k+1/2}(4N) -> S_{2k}(2N), so the lift has level | 288.
read("results/fq_for_pari.txt");
n = length(QS); mv = sum(i=1,n,FQ[i])/n; sv = sqrt(sum(i=1,n,(FQ[i]-mv)^2)); mq = QS[n];
best = 0; bestdesc = "none";
for(di=1, length(LEVELS), lev = LEVELS[di]; print("--- level ", lev, " ---"); MM = mfinit([lev,10],0); if(mfdim(MM)==0, next); BB = mfeigenbasis(MM); for(j=1, length(BB), F = BB[j]; co = mfcoefs(F, mq); EM = mfembed(F, co); if(type(EM[1])!="t_VEC", EM=[EM]); for(e=1, length(EM), v = EM[e]; u = vector(n, i, real(v[QS[i]+1])/QS[i]^4.5); mu = sum(i=1,n,u[i])/n; su = sqrt(sum(i=1,n,(u[i]-mu)^2)); if(su>1e-12, cr = sum(i=1,n,(u[i]-mu)*(FQ[i]-mv))/(su*sv); printf("  level %4d form %d emb %d  corr = %+.4f\n", lev, j, e, cr); if(abs(cr)>abs(best), best = cr; bestdesc = Str("level ",lev," form ",j," emb ",e)))))));
print(); print("best correlation: ", best, "  at ", bestdesc);
print("a genuine match should give |corr| > 0.9; anything below ~0.5 is noise.");
quit;

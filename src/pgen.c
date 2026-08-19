/* p(n) mod l for l in {13,17,19,23,29,31,37,41,43,47}, for n up to N.
 *
 * Writes one file per prime, res_<l>.bin, containing p(n) mod l as a single
 * unsigned char for each n in the arithmetic progression n = delta_l (mod l),
 * where 24*delta_l = 1 (mod l).  Those are exactly the n on the Ramanujan
 * branch for that prime, i.e. the n for which 24n-1 is divisible by l.
 *
 * Method: pentagonal-number recurrence, evaluated modulo the product of all
 * ten primes at once, then reduced.  The inner loop is blocked so that the
 * accumulator stays in cache; the large-offset pass reads p[] sequentially,
 * which is what makes this tractable at n = 3e7.
 *
 * Build:  gcc -O3 -march=native -funroll-loops -o pgen src/pgen.c
 * Usage:  ./pgen 30000000
 * Memory: 8 bytes per n, so about 240 MB at n = 3e7.
 * Time:   roughly 2 minutes at n = 3e7 on one core.
 *
 * Overflow: M = 13*17*...*47 = 266,186,053,068,611 (2.66e14). At N = 3e7 there
 * are 8943 pentagonal terms, so the accumulator stays under 8943*M = 2.38e18,
 * inside int64 (9.22e18) with a factor 3.9 to spare. That margin is NOT large:
 * the term count grows like 2*sqrt(2N/3), so this overflows silently for
 * N > 4.5e8. Reduce M (drop primes) or accumulate in __int128 before going
 * past that.
 */
#include <stdio.h>
#include <stdlib.h>
typedef long long ll;
#define BSH 18
#define BS (1<<BSH)

int main(int argc, char **argv){
    if(argc < 2){ fprintf(stderr,"usage: %s N\n", argv[0]); return 1; }
    ll N = atoll(argv[1]);
    int primes[10] = {13,17,19,23,29,31,37,41,43,47};
    ll M = 1;
    for(int i=0;i<10;i++) M *= primes[i];

    ll *p   = malloc((N+1)*sizeof(ll));
    ll *acc = malloc(BS*sizeof(ll));
    if(!p || !acc){ fprintf(stderr,"alloc failed (need %.2f GB)\n", 8.0*N/1e9); return 1; }

    ll cap = 300000;
    ll *off = malloc(cap*sizeof(ll));
    int *sg = malloc(cap*sizeof(int));
    ll nt = 0;
    for(ll k=1;;k++){
        ll g1 = k*(3*k-1)/2; if(g1 > N) break;
        int s = (k & 1) ? 1 : -1;
        off[nt]=g1; sg[nt]=s; nt++;
        ll g2 = k*(3*k+1)/2; if(g2 <= N){ off[nt]=g2; sg[nt]=s; nt++; }
    }
    ll nsmall = 0; while(nsmall < nt && off[nsmall] < BS) nsmall++;
    fprintf(stderr,"N=%lld  pentagonal terms=%lld  small=%lld  M=%lld\n", N, nt, nsmall, M);

    p[0] = 1;
    for(ll A=0; A<=N; A+=BS){
        ll hi = A+BS-1; if(hi > N) hi = N;
        ll len = hi-A+1;
        for(ll i=0;i<len;i++) acc[i]=0;
        /* large offsets: every value read is already final, and read in order */
        for(ll i=nsmall;i<nt;i++){
            ll o = off[i]; if(o > hi) break;
            ll lo = (o > A) ? o : A;
            int s = sg[i];
            const ll *src = p - o;
            if(s > 0) for(ll n=lo;n<=hi;n++) acc[n-A] += src[n];
            else      for(ll n=lo;n<=hi;n++) acc[n-A] -= src[n];
        }
        /* small offsets: sequential within the block, stays in cache */
        for(ll n=(A==0?1:A); n<=hi; n++){
            ll t = acc[n-A];
            for(ll i=0;i<nsmall;i++){
                ll o = off[i]; if(o > n) break;
                t += sg[i]*p[n-o];
            }
            t %= M; if(t < 0) t += M;
            p[n] = t;
        }
        if(((A>>BSH) & 31) == 0) fprintf(stderr,"\r%lld", A);
    }
    fprintf(stderr,"\nrecurrence done\n");

    for(int a=0;a<10;a++){
        int L = primes[a], d = 0;
        for(int x=0;x<L;x++) if((24*x) % L == 1){ d = x; break; }
        char fn[128]; sprintf(fn,"res_%d.bin",L);
        FILE *f = fopen(fn,"wb");
        ll cnt = 0;
        unsigned char *buf = malloc(N/L + 2);
        for(ll n=d; n<=N; n+=L) buf[cnt++] = (unsigned char)(p[n] % L);
        fwrite(buf,1,cnt,f); fclose(f); free(buf);
        fprintf(stderr,"l=%2d delta=%2d  %lld values -> %s\n", L, d, cnt, fn);
    }
    return 0;
}

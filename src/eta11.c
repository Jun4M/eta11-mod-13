/* Coefficients of prod_{k>=1}(1-q^k)^11 modulo 13.
   Exponents are m = 24*i + 11.  Writes a13.bin (one signed char per i).
   Checkpointed: rerun until it prints COMPLETE.
   Build:  gcc -O3 -march=native -fopenmp -o eta11 eta11.c
   Usage:  ./eta11 <MMAX> <passes_this_run>
   Memory: 2 bytes per coefficient, i.e. 2*(MMAX/24) bytes.
           MMAX=1e9  ->  0.08 GB ;  MMAX=1e10 -> 0.83 GB                */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
#define BSH 19
#define BS (1<<BSH)

int main(int argc, char **argv){
    if(argc<3){ fprintf(stderr,"usage: %s MMAX passes\n",argv[0]); return 1; }
    ll MMAX = atoll(argv[1]);
    int npass = atoi(argv[2]);
    ll NC = (MMAX-11)/24 + 1;
    signed char *f   = malloc(NC);
    signed char *out = malloc(NC);
    if(!f||!out){ fprintf(stderr,"alloc failed (need %.2f GB)\n", 2.0*NC/1e9); return 1; }

    ll cap = 4000000;
    ll *off = malloc(cap*sizeof(ll));
    signed char *sg = malloc(cap);
    ll nt = 0;
    for(ll k=1;;k++){
        ll g1 = k*(3*k-1)/2; if(g1>=NC) break;
        signed char s = (k&1)? -1 : 1;
        off[nt]=g1; sg[nt]=s; nt++;
        ll g2 = k*(3*k+1)/2; if(g2<NC){ off[nt]=g2; sg[nt]=s; nt++; }
    }
    fprintf(stderr,"NC=%lld  pentagonal terms=%lld  memory=%.2f GB\n",NC,nt,2.0*NC/1e9);

    int done = 0;
    FILE *st = fopen("state.bin","rb");
    if(st){
        if(fread(&done,sizeof(int),1,st)!=1 || fread(f,1,NC,st)!=(size_t)NC){
            fprintf(stderr,"state.bin does not match this MMAX; delete it and restart\n"); return 1; }
        fclose(st);
        fprintf(stderr,"resumed after pass %d\n",done);
    } else { memset(f,0,NC); f[0]=1; }

    for(int pass=0; pass<npass && done<11; pass++){
        ll nblocks = (NC+BS-1)/BS;
        #pragma omp parallel
        {
            int *acc = malloc(BS*sizeof(int));
            #pragma omp for schedule(dynamic)
            for(ll b=0; b<nblocks; b++){
                ll A = b*(ll)BS, hi = A+BS; if(hi>NC) hi=NC;
                ll len = hi-A;
                for(ll i=0;i<len;i++) acc[i]=f[A+i];
                for(ll i=0;i<nt;i++){
                    ll o=off[i]; if(o>=hi) break;
                    ll lo=(o>A)?o:A;
                    const signed char *src=f-o;
                    if(sg[i]>0) for(ll n=lo;n<hi;n++) acc[n-A]+=src[n];
                    else        for(ll n=lo;n<hi;n++) acc[n-A]-=src[n];
                }
                for(ll i=0;i<len;i++){ int t=acc[i]%13; if(t<0)t+=13; out[A+i]=(signed char)t; }
            }
            free(acc);
        }
        memcpy(f,out,NC); done++;
        fprintf(stderr,"pass %d/11 done\n",done);
    }
    st=fopen("state.bin","wb");
    fwrite(&done,sizeof(int),1,st); fwrite(f,1,NC,st); fclose(st);
    if(done==11){
        FILE *fp=fopen("a13.bin","wb"); fwrite(f,1,NC,fp); fclose(fp);
        fprintf(stderr,"COMPLETE: %lld coefficients in a13.bin\n",NC);
    } else fprintf(stderr,"checkpoint saved (%d/11). Run again to continue.\n",done);
    return 0;
}

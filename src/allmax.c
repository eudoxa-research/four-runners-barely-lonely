// Remark 9.12: all vertex maximizers of each subtorus and the slice-circle test.
// usage: ./allmax forms.txt ml.txt lo hi   (forms one per line, ml lines "p q X Y Den" from torusml or torusml2)
// For each 2-dim subtorus (line "s0 r0 s1 r1 s2 r2 s3 r3") with exact ML p/q (line "p q X Y Den" from torusml),
// enumerate ALL global maximizers that are vertices (>=2 independent active forms), and for each test whether
// some active pair (i,j) with matching signs has the maximizer on the identity component of ker(eps_i l_i - eps_j l_j).
// Output per subtorus: nmax ngood flags; prints subtori where some maximizer is bad.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static ll gcdll(ll a, ll b){ if(a<0)a=-a; if(b<0)b=-b; while(b){ll t=a%b;a=b;b=t;} return a; }
static ll md(ll a, ll m){ ll r=a%m; if(r<0) r+=m; return r; }
int S[4], R[4];
// point s = (X/N, Y/N) with N common denominator
typedef struct { ll X, Y, N; } pt;
static pt norm(pt a){ a.X=md(a.X,a.N); a.Y=md(a.Y,a.N); ll g=gcdll(gcdll(a.X,a.Y),a.N); if(g>1){a.X/=g;a.Y/=g;a.N/=g;} return a; }
// value of form i at s: (S X + R Y)/N mod 1 -> numerator mod N
static ll fnum(int i, pt s){ return md((ll)S[i]*s.X + (ll)R[i]*s.Y, s.N); }
int main(int argc, char **argv){
    FILE *fa=fopen(argv[1],"r"), *fb=fopen(argv[2],"r");
    double lo = argc>3 ? atof(argv[3]) : 1.0/3, hi = argc>4 ? atof(argv[4]) : 0.5;
    ll nsub=0, nbadsome=0, nbadall=0, nnovert=0, totmax=0, totbad=0;
    char line[512];
    static pt maxs[20000]; 
    while(fgets(line,sizeof line,fa)){
        if(sscanf(line,"%d %d %d %d %d %d %d %d",&S[0],&R[0],&S[1],&R[1],&S[2],&R[2],&S[3],&R[3])!=8) break;
        ll p,q,X0,Y0,D0; if(fscanf(fb,"%lld %lld %lld %lld %lld",&p,&q,&X0,&Y0,&D0)!=5) break;
        double dv=(double)p/q; if(!(dv>lo+1e-12 && dv<hi-1e-12)) continue;
        nsub++;
        int nm=0;
        for(int i=0;i<4;i++) for(int j=i+1;j<4;j++){
            ll det=(ll)S[i]*R[j]-(ll)S[j]*R[i]; if(det==0) continue;
            ll ad = det<0?-det:det;
            for(int ei=-1;ei<=1;ei+=2) for(int ej=-1;ej<=1;ej+=2){
                // solve S_i x + R_i y = ei p/q + ki, S_j x + R_j y = ej p/q + kj ; ki in [0,ad), kj in [0,ad) suffices (solutions mod 1)
                for(ll ki=0; ki<ad; ki++) for(ll kj=0; kj<ad; kj++){
                    // rhs scaled by q: bi = ei p + ki q, bj = ej p + kj q ; x = (R_j bi - R_i bj)/(det q), y = (S_i bj - S_j bi)/(det q)
                    ll bi = ei*p + ki*q, bj = ej*p + kj*q;
                    pt s; s.N = det*q; s.X = (ll)R[j]*bi - (ll)R[i]*bj; s.Y = (ll)S[i]*bj - (ll)S[j]*bi;
                    if(s.N<0){ s.N=-s.N; s.X=-s.X; s.Y=-s.Y; }
                    s = norm(s);
                    // G(s) == p/q ?
                    int ok=1;
                    for(int k=0;k<4;k++){ ll z=fnum(k,s); if(s.N-z<z) z=s.N-z; // ||.|| = z/N ; compare with p/q
                        if(z*q < p*s.N){ ok=0; break; } }
                    if(!ok) continue;
                    // dedupe
                    int dup=0; for(int t=0;t<nm;t++) if(maxs[t].X==s.X && maxs[t].Y==s.Y && maxs[t].N==s.N){dup=1;break;}
                    if(!dup && nm<20000) maxs[nm++]=s;
                }
            }
        }
        if(nm==0){ nnovert++; printf("NOVERTEX %d %d %d %d %d %d %d %d\n",S[0],R[0],S[1],R[1],S[2],R[2],S[3],R[3]); continue; }
        int nbad=0;
        for(int t=0;t<nm;t++){
            pt s=maxs[t]; int act[4], eps[4], na=0;
            for(int k=0;k<4;k++){ ll z=fnum(k,s); ll zz = (s.N-z<z)? s.N-z : z; act[k] = (zz*q == p*s.N); eps[k] = (z*2 < s.N) ? 1 : -1; if(act[k]) na++; }
            int good=0;
            for(int i=0;i<4 && !good;i++) for(int j=i+1;j<4 && !good;j++){
                if(!act[i]||!act[j]) continue;
                ll f0 = (ll)eps[i]*S[i] - (ll)eps[j]*S[j], f1 = (ll)eps[i]*R[i] - (ll)eps[j]*R[j];
                if(f0==0 && f1==0){ good=1; break; }
                ll g=gcdll(f0,f1); f0/=g; f1/=g;
                if(md(f0*s.X + f1*s.Y, s.N)==0) good=1;
            }
            if(!good) nbad++;
        }
        totmax+=nm; totbad+=nbad;
        if(nbad>0) nbadsome++;
        if(nbad==nm){ nbadall++; printf("ALLBAD %d %d %d %d %d %d %d %d  ML=%lld/%lld nmax=%d\n",S[0],R[0],S[1],R[1],S[2],R[2],S[3],R[3],p,q,nm); }
        else if(nbad>0 && nbadsome<=20) printf("SOMEBAD %d %d %d %d %d %d %d %d  ML=%lld/%lld nmax=%d nbad=%d\n",S[0],R[0],S[1],R[1],S[2],R[2],S[3],R[3],p,q,nm,nbad);
    }
    fprintf(stderr,"subtori in range: %lld ; total vertex maximizers %lld, bad %lld ; subtori with some bad maximizer %lld, all bad %lld, no vertex maximizer %lld\n",
        nsub, totmax, totbad, nbadsome, nbadall, nnovert);
    return 0;
}

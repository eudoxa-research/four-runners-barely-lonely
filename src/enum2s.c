// Enumerate all primitive rank-2 sublattices of Z^4 (row-HNF bases) with covolume^2 <= A2,
// not contained in a coordinate hyperplane (every column has a nonzero entry),
// and test whether ML(U) < THR (THR = thr_num/thr_den) using grid certificates; survivors are
// printed for exact processing.
// Usage: ./enum2 A thr_num thr_den [countonly]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef long long ll;
static ll gcdll(ll a, ll b){ if(a<0)a=-a; if(b<0)b=-b; while(b){ll t=a%b;a=b;b=t;} return a; }
static ll mod(ll a, ll m){ ll r=a%m; if(r<0) r+=m; return r; }

static ll A2; static ll thr_num, thr_den; static int countonly;
static ll total=0, survivors=0;

// quick certificate: is there a grid point (x,y) = (i/Q, j/Q) with min_k || s_k x + r_k y || >= THR ?
static int grid_cert(const int *s, const int *r, int Q){
    // value at (i,j): for each k, z = (s_k i + r_k j) mod Q, dist = min(z,Q-z)/Q >= thr  <=> dist*thr_den >= thr_num*Q
    for(int i=0;i<Q;i++) for(int j=0;j<Q;j++){
        int ok=1;
        for(int k=0;k<4;k++){
            ll z = mod((ll)s[k]*i + (ll)r[k]*j, Q);
            if(Q-z<z) z=Q-z;
            if(z*thr_den <= thr_num*Q){ ok=0; break; }
        }
        if(ok) return 1;
    }
    return 0;
}

static void process(const int *s, const int *r){
    total++;
    if(countonly) return;
    if(grid_cert(s,r,12)) return;
    if(grid_cert(s,r,60)) return;
    if(grid_cert(s,r,420)) return;
    if(grid_cert(s,r,16)) return;
    if(grid_cert(s,r,80)) return;
    if(grid_cert(s,r,208)) return;
    survivors++;
    printf("%d %d %d %d %d %d %d %d\n", s[0],r[0],s[1],r[1],s[2],r[2],s[3],r[3]);
}

int main(int argc, char **argv){
    ll A = atoll(argv[1]); A2 = A*A;
    thr_num = atoll(argv[2]); thr_den = atoll(argv[3]);
    countonly = argc>4;
    int s[4], r[4];
    for(int p=0;p<4;p++) for(int q=p+1;q<4;q++){
        for(int i=0;i<4;i++){ s[i]=0; r[i]=0; }
        for(ll a=1; a*a<=A2; a++){
            for(ll b=1; a*a*b*b<=A2; b++){
                ll rem0 = A2 - a*a*b*b;
                // columns: p (s=a, r=0), q (s=sq in [0,b), r=b), middle columns j in (p,q): s_j free, r_j=0;
                // columns j>q: s_j, r_j free.
                int mids[4], nm=0, highs[4], nh=0;
                for(int j=p+1;j<q;j++) mids[nm++]=j;
                for(int j=q+1;j<4;j++) highs[nh++]=j;
                for(ll sq=0; sq<b; sq++){
                    // middle columns: minor p_{jq} = s_j*b, p_{pj}=0 -> |s_j| <= sqrt(rem)/b
                    // enumerate middle s_j
                    ll smid[2]={0,0};
                    ll lim_mid = (ll)floor(sqrt((double)rem0)/b)+1;
                    // nested loops up to nm=2
                    ll m0lo = nm>0 ? -lim_mid : 0, m0hi = nm>0 ? lim_mid : 0;
                    for(ll m0=m0lo; m0<=m0hi; m0++){
                        ll rem1 = rem0 - (nm>0 ? m0*m0*b*b : 0); if(rem1<0) continue;
                        ll m1lo = nm>1 ? -lim_mid : 0, m1hi = nm>1 ? lim_mid : 0;
                        for(ll m1=m1lo; m1<=m1hi; m1++){
                            ll rem2 = rem1 - (nm>1 ? m1*m1*b*b : 0); if(rem2<0) continue;
                            smid[0]=m0; smid[1]=m1;
                            // high columns: for each j>q: p_{pj} = a*r_j, p_{qj} = sq*r_j - s_j*b, p_{mj} = s_m*r_j (for middle m)
                            // enumerate r_j and s_j for the high columns (nh <= 2)
                            ll lim_r = (ll)floor(sqrt((double)rem2)/a)+1;
                            ll h0rlo = nh>0 ? -lim_r : 0, h0rhi = nh>0 ? lim_r : 0;
                            for(ll r0=h0rlo; r0<=h0rhi; r0++){
                                ll c0 = (nh>0) ? (a*a*r0*r0 + (nm>0? smid[0]*smid[0]*r0*r0:0) + (nm>1? smid[1]*smid[1]*r0*r0:0)) : 0;
                                ll rem3 = rem2 - c0; if(rem3<0) continue;
                                // s_j0: p_{qj} = sq*r0 - s*b, need (sq*r0 - s*b)^2 <= rem3
                                ll s0lo=0, s0hi=0;
                                if(nh>0){ double sr=sqrt((double)rem3); s0lo=(ll)ceil((sq*r0 - sr)/b)-1; s0hi=(ll)floor((sq*r0 + sr)/b)+1; }
                                for(ll s0=s0lo; s0<=s0hi; s0++){
                                    ll pq0 = nh>0 ? (sq*r0 - s0*b) : 0;
                                    ll rem4 = rem3 - pq0*pq0; if(rem4<0) continue;
                                    ll lim_r1 = (ll)floor(sqrt((double)rem4)/a)+1;
                                    ll h1rlo = nh>1 ? -lim_r1 : 0, h1rhi = nh>1 ? lim_r1 : 0;
                                    for(ll r1=h1rlo; r1<=h1rhi; r1++){
                                        ll c1 = nh>1 ? (a*a*r1*r1 + (nm>0? smid[0]*smid[0]*r1*r1:0) + (nm>1? smid[1]*smid[1]*r1*r1:0)) : 0;
                                        ll rem5 = rem4 - c1; if(rem5<0) continue;
                                        ll s1lo=0, s1hi=0;
                                        if(nh>1){ double sr=sqrt((double)rem5); s1lo=(ll)ceil((sq*r1 - sr)/b)-1; s1hi=(ll)floor((sq*r1 + sr)/b)+1; }
                                        for(ll s1=s1lo; s1<=s1hi; s1++){
                                            ll pq1 = nh>1 ? (sq*r1 - s1*b) : 0;
                                            ll rem6 = rem5 - pq1*pq1; if(rem6<0) continue;
                                            // minor between the two high columns: s0*r1 - s1*r0
                                            ll phh = nh>1 ? (s0*r1 - s1*r0) : 0;
                                            if(phh*phh > rem6) continue;
                                            // assemble
                                            for(int i=0;i<4;i++){ s[i]=0; r[i]=0; }
                                            s[p]=(int)a; s[q]=(int)sq; r[q]=(int)b;
                                            if(nm>0) s[mids[0]]=(int)smid[0];
                                            if(nm>1) s[mids[1]]=(int)smid[1];
                                            if(nh>0){ s[highs[0]]=(int)s0; r[highs[0]]=(int)r0; }
                                            if(nh>1){ s[highs[1]]=(int)s1; r[highs[1]]=(int)r1; }
                                            // primitivity: gcd of all minors == 1 ; coordinate-free
                                            ll g=0; int cf=1;
                                            for(int i=0;i<4;i++){ if(s[i]==0 && r[i]==0) cf=0; }
                                            if(!cf) continue;
                                            for(int i=0;i<4;i++) for(int j=i+1;j<4;j++){ ll mnr=(ll)s[i]*r[j]-(ll)s[j]*r[i]; g=gcdll(g,mnr); }
                                            if(g!=1) continue;
                                            process(s,r);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    fprintf(stderr,"total lattices (covol<=%lld, primitive, coord-free): %lld ; survivors: %lld\n", A, total, survivors);
    return 0;
}

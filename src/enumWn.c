// Enumerate primitive rank-2 sublattices Lambda of Lambda_W for the two special 3-tori
//   W1 = {(2y, y, z, u)} : basis f1=(2,1,0,0), f2=(0,0,1,0), f3=(0,0,0,1)
//   W2 = {(a, b, a+b, d)}: basis f1=(1,0,1,0), f2=(0,1,1,0), f3=(0,0,0,1)
// Lambda = {c in Z^3 : n.c = 0} mapped by f, for primitive n in Z^3 (up to sign).
// covol^2 (in R^4 metric): W1: n1^2+5n2^2+5n3^2 ; W2: 2n1^2-2n1n2+2n2^2+3n3^2 (det(G) n^T G^{-1} n).
// Keep covol^2 < LIM2, coordinate-free, and print lattices failing the grid certificate for ML >= 5/16.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef long long ll;
static ll gcdll(ll a, ll b){ if(a<0)a=-a; if(b<0)b=-b; while(b){ll t=a%b;a=b;b=t;} return a; }
static ll mod(ll a, ll m){ ll r=a%m; if(r<0) r+=m; return r; }
static ll thr_num=5, thr_den=16;
static int grid_cert(const ll *s, const ll *r, int Q){
    for(int i=0;i<Q;i++) for(int j=0;j<Q;j++){
        int ok=1;
        for(int k=0;k<4;k++){
            ll z = mod(s[k]*i + r[k]*j, Q);
            if(Q-z<z) z=Q-z;
            if(z*thr_den <= thr_num*Q){ ok=0; break; }
        }
        if(ok) return 1;
    }
    return 0;
}
// kernel basis of n in Z^3: returns two integer vectors c1, c2 spanning {c : n.c=0} (primitive).
static void kernel3(ll n1, ll n2, ll n3, ll *c1, ll *c2){
    // Use extended Euclid: find basis of {c: n.c=0}. Approach: row-reduce [n^T | I3].
    ll E[3]={n1,n2,n3}; ll T[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    // reduce column vector E to (g,0,0) by unimodular row ops on (E,T)
    while(1){
        int nz=0, piv=-1; ll best=0;
        for(int i=0;i<3;i++) if(E[i]!=0){ nz++; if(piv<0 || llabs(E[i])<best){ piv=i; best=llabs(E[i]); } }
        if(nz<=1) break;
        for(int i=0;i<3;i++){ if(i==piv||E[i]==0) continue; ll q=E[i]/E[piv]; E[i]-=q*E[piv]; for(int k=0;k<3;k++) T[i][k]-=q*T[piv][k]; }
    }
    int idx=0; for(int i=0;i<3;i++){ if(E[i]==0){ ll *dst = idx==0? c1: c2; for(int k=0;k<3;k++) dst[k]=T[i][k]; idx++; } }
}
int main(int argc, char**argv){
    int which = atoi(argv[1]);        // 1 or 2
    ll LIM = atoll(argv[2]); ll LIM2 = LIM*LIM;
    ll f[3][4];
    if(which==1){ ll g[3][4]={{2,1,0,0},{0,0,1,0},{0,0,0,1}}; for(int i=0;i<3;i++)for(int k=0;k<4;k++) f[i][k]=g[i][k]; }
    else { ll g[3][4]={{1,0,1,0},{0,1,1,0},{0,0,0,1}}; for(int i=0;i<3;i++)for(int k=0;k<4;k++) f[i][k]=g[i][k]; }
    ll total=0, surv=0;
    ll N = LIM;
    for(ll n1=-N;n1<=N;n1++) for(ll n2=-N;n2<=N;n2++) for(ll n3=-N;n3<=N;n3++){
        if(n1==0&&n2==0&&n3==0) continue;
        // sign normalization: first nonzero positive
        ll first = n1?n1:(n2?n2:n3); if(first<0) continue;
        ll cv2 = (which==1)? (n1*n1+5*n2*n2+5*n3*n3) : (2*n1*n1-2*n1*n2+2*n2*n2+3*n3*n3);
        if(n1*n1+n2*n2+n3*n3 > LIM2) continue;
        if(gcdll(gcdll(n1,n2),n3)!=1) continue;
        ll c1[3],c2[3]; kernel3(n1,n2,n3,c1,c2);
        ll s[4],r[4];
        for(int k=0;k<4;k++){ s[k]=c1[0]*f[0][k]+c1[1]*f[1][k]+c1[2]*f[2][k]; r[k]=c2[0]*f[0][k]+c2[1]*f[1][k]+c2[2]*f[2][k]; }
        int cf=1; for(int k=0;k<4;k++) if(s[k]==0&&r[k]==0) cf=0;
        if(!cf) continue;
        // sanity: covolume^2 from Gram
        ll ss=0,rr=0,sr=0; for(int k=0;k<4;k++){ ss+=s[k]*s[k]; rr+=r[k]*r[k]; sr+=s[k]*r[k]; }
        if(ss*rr-sr*sr != cv2){ fprintf(stderr,"covol mismatch n=(%lld,%lld,%lld)\n",n1,n2,n3); return 1; }
        total++;
        if(grid_cert(s,r,12)) continue;
        if(grid_cert(s,r,60)) continue;
        if(grid_cert(s,r,420)) continue;
        if(grid_cert(s,r,16)) continue;
        if(grid_cert(s,r,80)) continue;
        if(grid_cert(s,r,208)) continue;
        surv++;
        printf("%lld %lld %lld %lld %lld %lld %lld %lld\n", s[0],r[0],s[1],r[1],s[2],r[2],s[3],r[3]);
    }
    fprintf(stderr,"W%d: lattices with covol<%lld: %lld, survivors %lld\n", which, LIM, total, surv);
    return 0;
}

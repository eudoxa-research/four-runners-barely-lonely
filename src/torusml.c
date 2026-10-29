// Exact ML of a 2-dim subtorus of (R/Z)^4 given 4 forms (s_i, r_i):
//   F(x,y) = min_i || s_i x + r_i y ||,  ML = max F over the 2-torus.
// Candidates: (a) three signed branches equal: s_i x + r_i y - sg_i d = k_i (3 eqs);
//             (b) two proportional forms at level d and a third form at a half-integer;
//             (c) two independent forms at half-integers (F=1/2 check).
// Input lines: s1 r1 s2 r2 s3 r3 s4 r4 ; output: p q X Y Den  (ML = p/q, witness (X/Den, Y/Den)).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static ll gcdll(ll a, ll b){ if(a<0)a=-a; if(b<0)b=-b; while(b){ll t=a%b;a=b;b=t;} return a; }
static ll mod(ll a, ll m){ ll r=a%m; if(r<0) r+=m; return r; }

static int ns[4], nr[4];
// evaluate F at (X/Den, Y/Den): returns numerator over Den (value = num/Den)
static ll Fnum(ll X, ll Y, ll Den){
    ll best = Den; // 1 (upper bound)
    for(int i=0;i<4;i++){
        ll z = mod((ll)ns[i]*X + (ll)nr[i]*Y, Den);
        if(Den - z < z) z = Den - z;
        if(z < best) best = z;
    }
    return best;
}
static ll det3(ll a11,ll a12,ll a13,ll a21,ll a22,ll a23,ll a31,ll a32,ll a33){
    return a11*(a22*a33-a23*a32) - a12*(a21*a33-a23*a31) + a13*(a21*a32-a22*a31);
}
int main(){
    while(scanf("%d %d %d %d %d %d %d %d",&ns[0],&nr[0],&ns[1],&nr[1],&ns[2],&nr[2],&ns[3],&nr[3])==8){
        ll bp=0, bq=1, bX=0, bY=0, bDen=1;
        int B[4]; for(int i=0;i<4;i++) B[i]=abs(ns[i])+abs(nr[i])+1;
        // (a) three active
        for(int i=0;i<4;i++) for(int j=i+1;j<4;j++) for(int m=j+1;m<4;m++)
        for(int sgi=-1;sgi<=1;sgi+=2) for(int sgj=-1;sgj<=1;sgj+=2) for(int sgm=-1;sgm<=1;sgm+=2){
            ll a11=ns[i],a12=nr[i],a13=-sgi, a21=ns[j],a22=nr[j],a23=-sgj, a31=ns[m],a32=nr[m],a33=-sgm;
            ll det = det3(a11,a12,a13,a21,a22,a23,a31,a32,a33);
            if(det==0) continue;
            // cofactors for Cramer: x = det(b|col1), etc. Precompute adjugate columns.
            // Solution: [x;y;d] = adj(A) * b / det, adj = transpose of cofactor matrix.
            ll C11=(a22*a33-a23*a32), C12=-(a21*a33-a23*a31), C13=(a21*a32-a22*a31);
            ll C21=-(a12*a33-a13*a32), C22=(a11*a33-a13*a31), C23=-(a11*a32-a12*a31);
            ll C31=(a12*a23-a13*a22), C32=-(a11*a23-a13*a21), C33=(a11*a22-a12*a21);
            // x*det = C11*b1 + C21*b2 + C31*b3 ; y*det = C12*b1 + C22*b2 + C32*b3 ; d*det = C13*b1+C23*b2+C33*b3
            ll adet = det<0? -det: det; ll sgn = det<0? -1: 1;
            for(ll ki=-B[i]; ki<=B[i]; ki++) for(ll kj=-B[j]; kj<=B[j]; kj++) for(ll km=-B[m]; km<=B[m]; km++){
                ll dn = sgn*(C13*ki + C23*kj + C33*km); // d = dn/adet
                if(dn<=0 || 2*dn>adet) continue;
                if(dn*bq <= bp*adet) continue; // cannot beat
                ll xn = sgn*(C11*ki + C21*kj + C31*km), yn = sgn*(C12*ki + C22*kj + C32*km);
                ll v = Fnum(mod(xn,adet), mod(yn,adet), adet);
                if(v*bq > bp*adet){ ll g=gcdll(v,adet); bp=v/g; bq=adet/g; bX=mod(xn,adet); bY=mod(yn,adet); bDen=adet; }
            }
        }
        // (b) proportional pair at level d, third at half integer: use doubled forms to keep integers:
        // equations: s_i x + r_i y - sg_i d = k_i ; s_j x + r_j y - sg_j d = k_j ; 2 s_m x + 2 r_m y = 2 k_m + 1
        for(int i=0;i<4;i++) for(int j=i+1;j<4;j++){
            if((ll)ns[i]*nr[j] != (ll)ns[j]*nr[i]) continue;
            for(int m=0;m<4;m++){ if(m==i||m==j) continue;
                for(int sgi=-1;sgi<=1;sgi+=2) for(int sgj=-1;sgj<=1;sgj+=2){
                    ll a11=ns[i],a12=nr[i],a13=-sgi, a21=ns[j],a22=nr[j],a23=-sgj, a31=2*ns[m],a32=2*nr[m],a33=0;
                    ll det = det3(a11,a12,a13,a21,a22,a23,a31,a32,a33);
                    if(det==0) continue;
                    ll C11=(a22*a33-a23*a32), C12=-(a21*a33-a23*a31), C13=(a21*a32-a22*a31);
                    ll C21=-(a12*a33-a13*a32), C22=(a11*a33-a13*a31), C23=-(a11*a32-a12*a31);
                    ll C31=(a12*a23-a13*a22), C32=-(a11*a23-a13*a21), C33=(a11*a22-a12*a21);
                    ll adet = det<0? -det: det; ll sgn = det<0? -1: 1;
                    for(ll ki=-B[i]; ki<=B[i]; ki++) for(ll kj=-B[j]; kj<=B[j]; kj++) for(ll km=-B[m]; km<=B[m]; km++){
                        ll b3 = 2*km+1;
                        ll dn = sgn*(C13*ki + C23*kj + C33*b3);
                        if(dn<=0 || 2*dn>adet) continue;
                        if(dn*bq <= bp*adet) continue;
                        ll xn = sgn*(C11*ki + C21*kj + C31*b3), yn = sgn*(C12*ki + C22*kj + C32*b3);
                        ll v = Fnum(mod(xn,adet), mod(yn,adet), adet);
                        if(v*bq > bp*adet){ ll g=gcdll(v,adet); bp=v/g; bq=adet/g; bX=mod(xn,adet); bY=mod(yn,adet); bDen=adet; }
                    }
                }
            }
        }
        // (c) two independent forms at half-integers
        for(int i=0;i<4;i++) for(int j=i+1;j<4;j++){
            ll det = (ll)ns[i]*nr[j] - (ll)ns[j]*nr[i];
            if(det==0) continue;
            ll adet = det<0?-det:det, sgn = det<0?-1:1;
            // 2 s_i x + 2 r_i y = 2ki+1 ; 2 s_j x + 2 r_j y = 2kj+1 ; solve with denominator 2*det
            for(ll ki=-B[i]; ki<=B[i]; ki++) for(ll kj=-B[j]; kj<=B[j]; kj++){
                ll b1=2*ki+1, b2=2*kj+1;
                ll xn = sgn*(b1*nr[j] - b2*nr[i]), yn = sgn*(ns[i]*b2 - ns[j]*b1); // over 2*adet
                ll Den = 2*adet;
                ll v = Fnum(mod(xn,Den), mod(yn,Den), Den);
                if(v*bq > bp*Den){ ll g=gcdll(v,Den); bp=v/g; bq=Den/g; bX=mod(xn,Den); bY=mod(yn,Den); bDen=Den; }
            }
        }
        printf("%lld %lld %lld %lld %lld\n", bp, bq, bX, bY, bDen);
    }
    return 0;
}

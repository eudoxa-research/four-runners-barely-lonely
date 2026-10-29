// g5enum2.c -- enumeration for Theorem 11.5 of "Four runners that are barely lonely" (five speeds). For a seed q in Z^5 (q[0]=1, positive entries: a tight quadruple {1,2,3,4} or {1,3,4,7}
// with one entry doubled) enumerate c=(0,c1..c4) with gcd(c)=1, first nonzero entry positive, and
// A2(c) = sum_{i<j} (q_i c_j - q_j c_i)^2 = covol(span(q,c))^2 <= MAXM. Keep W=span(q,c) only if every coordinate-free circle of W
// with a repeated absolute speed has exactly four distinct absolute speeds forming a multiple of {1,2,3,4}
// or {1,3,4,7}. Output: q, c, M for survivors; summary counts on stderr.
// Usage: ./g5enum q0 q1 q2 q3 q4 MAXM
#include <stdio.h>
#include <stdlib.h>
typedef long long ll;
static ll q[5], c[5]; static ll MAXM; static long long nvisit=0, nprim=0, nsurv=0;
static ll gcdll(ll a, ll b){ if(a<0)a=-a; if(b<0)b=-b; while(b){ll t=a%b;a=b;b=t;} return a; }
static ll absll(ll a){ return a<0?-a:a; }
static int tightset(ll *v){ // v[5] absolute values; need exactly 4 distinct, normalized {1,2,3,4} or {1,3,4,7}
    ll d[5]; int nd=0;
    for(int i=0;i<5;i++){ int f=0; for(int j=0;j<nd;j++) if(d[j]==v[i]) f=1; if(!f) d[nd++]=v[i]; }
    if(nd!=4) return 0;
    for(int i=0;i<4;i++) for(int j=i+1;j<4;j++) if(d[j]<d[i]){ ll t=d[i]; d[i]=d[j]; d[j]=t; }
    ll g=0; for(int i=0;i<4;i++) g=gcdll(g,d[i]);
    for(int i=0;i<4;i++) d[i]/=g;
    if(d[0]==1&&d[1]==2&&d[2]==3&&d[3]==4) return 1;
    if(d[0]==1&&d[1]==3&&d[2]==4&&d[3]==7) return 1;
    return 0;
}
static int filter(void){
    for(int i=0;i<5;i++) for(int j=i+1;j<5;j++) for(int e=-1;e<=1;e+=2){
        ll d1=q[i]-e*q[j], d2=c[i]-e*c[j];
        if(d1==0&&d2==0) continue;
        ll g=gcdll(d1,d2); ll u1=-d2/g, u2=d1/g;
        ll v[5]; int zero=0;
        for(int k=0;k<5;k++){ ll s=q[k]*u1+c[k]*u2; if(s==0){zero=1;break;} v[k]=absll(s); }
        if(zero) continue;
        if(!tightset(v)) return 0;
    }
    return 1;
}
static ll fpart(int k, ll x){ // sum_{i<k} (q_i x - q_k c_i)^2
    ll s=0; for(int i=0;i<k;i++){ ll d=q[i]*x-q[k]*c[i]; s+=d*d; } return s;
}
static void rec(int k, ll used){
    if(k==5){
        nvisit++;
        // sign normalization: first nonzero of c[1..4] positive ; gcd 1
        int fnz=0; for(int i=1;i<5;i++) if(c[i]!=0){ fnz=i; break; }
        if(fnz==0) return;
        if(c[fnz]<0) return;
        ll g=0; for(int i=1;i<5;i++) g=gcdll(g,c[i]); if(g!=1) return;
        nprim++;
        if(filter()){ nsurv++; printf("%lld %lld %lld %lld %lld  %lld %lld %lld %lld %lld  %lld\n",q[0],q[1],q[2],q[3],q[4],c[0],c[1],c[2],c[3],c[4],used); }
        return;
    }
    // f(x) = sum_{i<k} (q_i x - q_k c_i)^2 convex quadratic; integer minimizer near q_k*S1/S2
    ll S1=0,S2=0; for(int i=0;i<k;i++){ S1+=q[i]*c[i]; S2+=q[i]*q[i]; }
    ll t = (q[k]*S1) / S2; ll best=-1, bx=0;
    for(ll x=t-1;x<=t+1;x++){ ll f=fpart(k,x); if(best<0||f<best){best=f;bx=x;} }
    if(used+best>MAXM) return;
    for(ll x=bx; ; x++){ ll f=fpart(k,x); if(used+f>MAXM) break; c[k]=x; rec(k+1,used+f); }
    for(ll x=bx-1; ; x--){ ll f=fpart(k,x); if(used+f>MAXM) break; c[k]=x; rec(k+1,used+f); }
    c[k]=0;
}
int main(int argc, char **argv){
    if(argc<7){ fprintf(stderr,"usage\n"); return 1; }
    for(int i=0;i<5;i++) q[i]=atoll(argv[1+i]); MAXM=atoll(argv[6]);
    c[0]=0; rec(1,0);
    fprintf(stderr,"seed %lld %lld %lld %lld %lld MAXM %lld: visited %lld primitive %lld survivors %lld\n",q[0],q[1],q[2],q[3],q[4],MAXM,nvisit,nprim,nsurv);
    return 0;
}

#include <stdio.h>
#include <stdlib.h>
typedef long long ll;
static ll gcdll(ll a, ll b){ if(a<0)a=-a; if(b<0)b=-b; while(b){ll t=a%b;a=b;b=t;} return a; }
int main(int argc,char**argv){
  int A=atoi(argv[1]); ll A2=(ll)A*A; ll cnt=0, cntcf=0;
  for(int p12=-A;p12<=A;p12++) for(int p13=-A;p13<=A;p13++) for(int p14=-A;p14<=A;p14++)
  for(int p23=-A;p23<=A;p23++) for(int p24=-A;p24<=A;p24++) for(int p34=-A;p34<=A;p34++){
    ll n2=(ll)p12*p12+(ll)p13*p13+(ll)p14*p14+(ll)p23*p23+(ll)p24*p24+(ll)p34*p34;
    if(n2==0||n2>A2) continue;
    if((ll)p12*p34-(ll)p13*p24+(ll)p14*p23!=0) continue;
    ll g=0; g=gcdll(g,p12);g=gcdll(g,p13);g=gcdll(g,p14);g=gcdll(g,p23);g=gcdll(g,p24);g=gcdll(g,p34);
    if(g!=1) continue;
    // sign normalization: first nonzero positive
    int first = p12?p12:(p13?p13:(p14?p14:(p23?p23:(p24?p24:p34))));
    if(first<0) continue;
    cnt++;
    // coordinate-free: index i appears in some nonzero p
    int c1=(p12||p13||p14), c2=(p12||p23||p24), c3=(p13||p23||p34), c4=(p14||p24||p34);
    if(c1&&c2&&c3&&c4) cntcf++;
  }
  printf("A=%d: primitive decomposable Plucker vectors up to sign: %lld, coord-free: %lld\n",A,cnt,cntcf);
  return 0;
}

// mlc.c -- exact ML of a set of positive integers by the corner lemma (Lemma 2.1).
// Input lines "k v1 ... vk" (k <= 16); output lines "p q" with ML = p/q in lowest terms.
#include <stdio.h>
typedef long long ll;
static ll g(ll a,ll b){while(b){ll t=a%b;a=b;b=t;}return a;}
int main(){ int k; while(scanf("%d",&k)==1){ ll v[16]; for(int i=0;i<k;i++) if(scanf("%lld",&v[i])!=1) return 1;
  ll bn=0,bd=1;
  for(int i=0;i<k;i++) for(int j=i;j<k;j++){ ll D=v[i]+v[j]; for(ll m=1;2*m<=D;m++){ ll mn=D; for(int t=0;t<k;t++){ ll r=(v[t]*m)%D; if(D-r<r) r=D-r; if(r<mn) mn=r; if(mn*bd<=bn*D) break; } if(mn*bd>bn*D){bn=mn;bd=D;} } }
  ll gg=g(bn,bd); printf("%lld %lld\n",bn/gg,bd/gg); } return 0; }

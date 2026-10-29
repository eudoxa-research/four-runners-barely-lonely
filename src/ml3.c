#include <stdio.h>
#include <stdlib.h>
static long gcdl(long a,long b){while(b){long t=a%b;a=b;b=t;}return a;}
int main(int argc,char**argv){ int N=atoi(argv[1]); long cnt=0,low=0,norel=0;
 for(int a=1;a<=N;a++)for(int b=a+1;b<=N;b++)for(int c=b+1;c<=N;c++){
  if(gcdl(gcdl(a,b),c)!=1) continue; cnt++;
  int v[3]={a,b,c}; long bn=0,bd=1; int high=0;
  for(int i=0;i<3&&!high;i++)for(int j=i;j<3&&!high;j++){ long D=v[i]+v[j]; for(long m=1;2*m<=D;m++){ long mn=D; for(int k=0;k<3;k++){ long r=((long)v[k]*m)%D; if(D-r<r) r=D-r; if(r<mn) mn=r; } if(mn*bd>bn*D){bn=mn;bd=D;} if(3*bn>=bd){high=1;break;} } }
  if(high) continue; low++;
  int rel = (b==2*a)||(c==2*a)||(c==2*b)||(c==a+b);
  if(!rel){ norel++; printf("no special relation: %d %d %d ML=%ld/%ld\n",a,b,c,bn,bd); }
 }
 fprintf(stderr,"primitive triples max<=%d: %ld, ML<1/3: %ld, without special relation: %ld\n",N,cnt,low,norel); return 0; }

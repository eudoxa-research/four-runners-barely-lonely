// Exhaustive: all primitive 4-sets v1<v2<v3<v4 with v1^2+..+v4^2 < R^2; report those with ML < 1/4,
// with exact ML (max over corner times m/(u+w)). Quick certificate first: a list of times k/Q.
#include <stdio.h>
#include <stdlib.h>
typedef long long ll;
static ll gcdl(ll a, ll b){ while(b){ ll t=a%b; a=b; b=t;} return a; }
int main(int argc, char **argv){
    ll R = atoll(argv[1]); ll R2 = R*R;
    ll lo = argc>2 ? atoll(argv[2]) : 1, hi = argc>3 ? atoll(argv[3]) : R; // range for v1 (parallel split)
    // candidate times k/Q
    int Qs[] = {4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,72,80,96,100,120,128,144,160,192,240,256,320,384,480,512};
    int nQ = sizeof(Qs)/sizeof(int);
    static int ck[20000], cQ[20000]; int nc=0;
    for(int qi=0; qi<nQ; qi++){ int Q=Qs[qi]; for(int k=1; 2*k<=Q; k++){ if(gcdl(k,Q)!=1) continue; ck[nc]=k; cQ[nc]=Q; nc++; } }
    ll count_sets=0, count_low=0, count_exact=0;
    int v[4];
    for(v[0]=(int)lo; v[0]<=hi; v[0]++){
        if((ll)v[0]*v[0]*4 >= R2) break;
        for(v[1]=v[0]+1; ; v[1]++){
            ll s2=(ll)v[0]*v[0]+(ll)v[1]*v[1]; if(s2+(ll)(v[1]+1)*(v[1]+1)+(ll)(v[1]+2)*(v[1]+2) >= R2) break;
            for(v[2]=v[1]+1; ; v[2]++){
                ll s3=s2+(ll)v[2]*v[2]; if(s3+(ll)(v[2]+1)*(v[2]+1) >= R2) break;
                for(v[3]=v[2]+1; ; v[3]++){
                    ll s4=s3+(ll)v[3]*v[3]; if(s4>=R2) break;
                    if(gcdl(gcdl(v[0],v[1]),gcdl(v[2],v[3]))!=1) continue;
                    count_sets++;
                    // quick certificate
                    int cert=0;
                    for(int c=0;c<nc && !cert;c++){
                        int Q=cQ[c], k=ck[c]; int ok=1;
                        for(int i=0;i<4;i++){ ll z=((ll)v[i]*k)%Q; if(Q-z<z) z=Q-z; if(4*z<Q){ ok=0; break; } }
                        if(ok) cert=1;
                    }
                    if(cert) continue;
                    count_exact++;
                    // exact ML via corners
                    ll best_num=0, best_den=1; int found_high=0;
                    for(int a=0;a<4 && !found_high;a++) for(int b=a;b<4 && !found_high;b++){
                        ll D=v[a]+v[b];
                        for(ll m=1; 2*m<=D; m++){
                            ll mn=D;
                            for(int i=0;i<4;i++){ ll r=((ll)v[i]*m)%D; if(D-r<r) r=D-r; if(r<mn) mn=r; if(mn*best_den<=best_num*D) break; }
                            if(mn*best_den>best_num*D){ best_num=mn; best_den=D; }
                            if(4*best_num>=best_den){ found_high=1; break; }
                        }
                    }
                    if(!found_high){
                        ll g=gcdl(best_num,best_den); if(g==0) g=1;
                        printf("%d %d %d %d %lld %lld\n", v[0],v[1],v[2],v[3], best_num/g, best_den/g);
                        count_low++;
                    }
                }
            }
        }
    }
    fprintf(stderr,"v1 in [%lld,%lld]: primitive sets %lld, exact-checked %lld, ML<1/4: %lld\n", lo, hi, count_sets, count_exact, count_low);
    return 0;
}

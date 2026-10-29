# Independent pure-Python re-implementation of the G5 enumeration for a smaller covolume bound:
# enumerates c=(0,c1..c4) by a plain box search (no pruning) and compares the number of primitive
# normalized candidates and the survivors with the C program g5enum2 run at the same bound.
import sys, subprocess, itertools
from math import gcd, isqrt
B=int(sys.argv[1]) if len(sys.argv)>1 else 3000
SEEDS=[(1,1,2,3,4),(1,2,2,3,4),(1,2,3,3,4),(1,2,3,4,4),(1,1,3,4,7),(1,3,3,4,7),(1,3,4,4,7),(1,3,4,7,7)]
def gl(v):
    g=0
    for x in v: g=gcd(g,abs(x))
    return g
def tight(vals):
    d=sorted(set(vals))
    if len(d)!=4: return False
    g=gl(d); return [x//g for x in d] in ([1,2,3,4],[1,3,4,7])
def passes(q,c):
    for i in range(5):
        for j in range(i+1,5):
            for e in (1,-1):
                d1=q[i]-e*q[j]; d2=c[i]-e*c[j]
                if d1==0 and d2==0: continue
                g=gcd(abs(d1),abs(d2)); u1,u2=-d2//g,d1//g
                s=[abs(q[k]*u1+c[k]*u2) for k in range(5)]
                if 0 in s: continue
                if not tight(s): return False
    return True
ok=True
sel=[int(x) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else range(len(SEEDS))
for q in [SEEDS[i] for i in sel]:
    R=isqrt(B)  # |c_k| = |D_{0k}| <= covol
    prim=0; surv=[]
    import numpy as np
    rng=np.arange(-R,R+1,dtype=np.int64)
    C3,C4,C5=np.meshgrid(rng,rng,rng,indexing='ij')
    C3=C3.ravel(); C4=C4.ravel(); C5=C5.ravel()
    for c2 in range(-R,R+1):
        cols=[np.zeros_like(C3), np.full_like(C3,c2), C3, C4, C5]
        a2=np.zeros_like(C3)
        for i in range(5):
            for j in range(i+1,5):
                d=q[i]*cols[j]-q[j]*cols[i]; a2+=d*d
        for idx in np.nonzero(a2<=B)[0]:
            cc=(c2,int(C3[idx]),int(C4[idx]),int(C5[idx]))
            nz=[x for x in cc if x]
            if not nz or nz[0]<0 or gl(cc)!=1: continue
            prim+=1
            if passes(q,(0,)+cc): surv.append((0,)+cc)
    out=subprocess.run(['./g5enum2']+[str(x) for x in q]+[str(B)],capture_output=True,text=True)
    cs=sorted(tuple(int(x) for x in l.split()[5:10]) for l in out.stdout.strip().split('\n') if l.strip())
    cprim=int(out.stderr.split('primitive')[1].split()[0])
    agree = (cprim==prim) and (cs==sorted(surv))
    ok &= agree
    print(q, 'python primitive', prim, 'C primitive', cprim, 'survivors', len(surv), 'agree', agree, flush=True)
print('ALL AGREE' if ok else 'MISMATCH')

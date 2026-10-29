#!/usr/bin/env python3
"""check_obstr_c.py -- Section 10 of "Four runners that are barely lonely", checked with the C programs.
Usage: python3 src/check_obstr_c.py ./torusml2 ./mlc
 - Theorem 10.1(i): ML(U(A,B,C)) = 1/2 - 1/N for 255 random admissible (A,B,C) with N <= 126 (exact maximizer);
 - ML(1, N/2, N-1) = 1/2 - 1/N for 4 | N, 8 <= N <= 400;
 - Theorem 10.1(ii), (iii): the loneliness of every repeated-speed circle of U(A,B,C) for 10,958 admissible
   (A,B,C) with N <= 196 (a fixed residue subsample), and for the examples A = m, C = m+1, m <= 4, K <= 2;
 - Theorem 10.2 for 4 | N, 12 <= N <= 124, 1 <= l <= 6, and the comparison for (N, l) = (60, 61), (120, 121).
"""
import sys
TORUSML2 = sys.argv[1] if len(sys.argv) > 1 else './torusml2'
MLC = sys.argv[2] if len(sys.argv) > 2 else './mlc'
import itertools, subprocess, random
from fractions import Fraction as F
from math import gcd
def kernel3(n):
    E=list(n); T=[[1,0,0],[0,1,0],[0,0,1]]
    while sum(1 for e in E if e)>1:
        piv=min((i for i in range(3) if E[i]),key=lambda i:abs(E[i]))
        for i in range(3):
            if i==piv or E[i]==0: continue
            q=E[i]//E[piv]; E[i]-=q*E[piv]; T[i]=[x-q*y for x,y in zip(T[i],T[piv])]
    return [T[i] for i in range(3) if E[i]==0]
def forms_U(A,B,C):
    N=A+B+C; H=N//2
    c1,c2=kernel3((A,B,C))
    s=[c1[0],c1[1],c1[2],H*c1[1]]; r=[c2[0],c2[1],c2[2],H*c2[1]]
    return s,r
cases=[]
random.seed(1)
for N in range(12,130,4):
    for _ in range(12):
        A=random.randint(1,(N-2)//2-1); C=random.randint(1,(N-2)//2-A) if (N-2)//2-A>=1 else None
        if C is None: continue
        if N < 2*(A+C)+2: continue
        B=N-A-C
        if gcd(gcd(A,B),C)!=1: continue
        cases.append((A,B,C))
cases=sorted(set(cases))
lines=[]
for A,B,C in cases:
    s,r=forms_U(A,B,C); lines.append(' '.join(f"{s[i]} {r[i]}" for i in range(4)))
out=subprocess.run([TORUSML2],input='\n'.join(lines)+'\n',capture_output=True,text=True).stdout.split('\n')
bad=0
for (A,B,C),o in zip(cases,out):
    N=A+B+C; p,q=map(int,o.split()[:2])
    if F(p,q)!=F(1,2)-F(1,N): bad+=1; print("MISMATCH",A,B,C,F(p,q))
assert bad==0
print("PASS: Theorem 10.1(i) value on",len(cases),"random subtori U(A,B,C) with N <= 126")

import subprocess
def mlc(lists):
    inp='\n'.join(f"{len(v)} "+' '.join(map(str,v)) for v in lists)+'\n'
    out=subprocess.run([MLC],input=inp,capture_output=True,text=True).stdout.split('\n')
    return [F(int(o.split()[0]),int(o.split()[1])) for o in out if o.strip()]
# (11.3)
Ns=list(range(8,404,4)); vals=mlc([[1,N//2,N-1] for N in Ns])
assert all(v==F(1,2)-F(1,N) for N,v in zip(Ns,vals)); print("PASS: ML(1,N/2,N-1) = 1/2 - 1/N for 4 | N, 8 <= N <= 400")
# repeated-speed circles in U(A,B,C)
def rep_circles(A,B,C):
    s,r=forms_U(A,B,C); out=[]
    for i in range(4):
        for j in range(i+1,4):
            for sg in (1,-1):
                a=s[i]-sg*s[j]; b=r[i]-sg*r[j]
                if a==0 and b==0: continue
                gg=gcd(abs(a),abs(b)); bp,ap=b//gg,a//gg
                v=[s[k]*bp-r[k]*ap for k in range(4)]
                if any(x==0 for x in v): continue
                out.append(((i,j,sg),[abs(x) for x in v]))
    return out
def check_family(A,B,C):
    N=A+B+C; L=F(1,2)-F(1,N)
    circ=rep_circles(A,B,C); vals=mlc([v for _,v in circ])
    allg=all(gcd(x,N)>1 for x in (A,B,C))
    best=max(vals) if vals else None
    if allg: assert best<=L-F(1,N**3), (A,B,C,best)
    else: assert best==L, (A,B,C,best)
    return allg, best, L
cnt=[0,0]
for N in range(12,200,4):
    for A in range(1,N):
        for C in range(1,N):
            B=N-A-C
            if B<1 or N<2*(A+C)+2 or gcd(gcd(A,B),C)!=1: continue
            if (A*7+C*13+N)%5: continue   # subsample
            allg,best,L=check_family(A,B,C); cnt[allg]+=1
print("PASS: Theorem 10.1(ii),(iii) on",cnt[0]+cnt[1],"subtori (",cnt[1],"with all three gcds > 1 )")
# Corollary 11.3 examples
for m in (2,3,4):
    for K in (1,2):
        A,C=m,m+1; N=4*m*(m+1)*(2*m+1)*K; B=N-2*m-1
        allg,best,L=check_family(A,B,C); assert allg
        assert L-best>=F(1,N**3)
        print("PASS: example m=%d K=%d (N=%d): best repeated-speed circle %s, below L by at least N^-3" % (m,K,N,best))
# Theorem 13.3
vs=[]; keys=[]
for N in range(12,125,4):
    for l in range(1,7):
        vs.append([1,N//2,2*N*l-3,3*N*l-N//2-2]); keys.append((N,l))
vals=mlc(vs)
for (N,l),v in zip(keys,vals): assert v==F((N-2)*l-1,2*N*l-2),(N,l,v)
print("PASS: Theorem 10.2 formula for",len(keys),"pairs (N <= 124, l <= 6)")
# Corollary 13.4: N=60, l=61 and N=120, l=121
for N in (60,120):
    l=N+1; V=[1,N//2,2*N*l-3,3*N*l-N//2-2]; mv=mlc([V])[0]
    allg,best,L=check_family(2,N-5,3)
    assert mv>best and mv==F((N-2)*l-1,2*N*l-2)
    print("PASS: Theorem 10.2 comparison N=%d l=%d: ML(V)=%s exceeds the best repeated-speed circle %s of U(2,N-5,3)" % (N,l,mv,best))

# Section 10, Proposition 10.3 and Theorem 10.4 (the family U_H)
import random as _r
random=_r
def forms_H(A,B,C,H):
    c1,c2=kernel3((A,B,C))
    return [c1[0],c1[1],c1[2],H*c1[1]],[c2[0],c2[1],c2[2],H*c2[1]]
def tml2(items):
    inp='\n'.join(' '.join(f"{s[i]} {r[i]}" for i in range(4)) for s,r in items)+'\n'
    out=subprocess.run([TORUSML2],input=inp,capture_output=True,text=True).stdout.split('\n')
    return [F(int(o.split()[0]),int(o.split()[1])) for o in out if o.strip()]
def rep_best_H(A,B,C,H):
    s,r=forms_H(A,B,C,H); out=[]
    for i,j in itertools.combinations(range(4),2):
        for sg in (1,-1):
            a,b=s[i]-sg*s[j],r[i]-sg*r[j]
            if a==0 and b==0: continue
            g=gcd(abs(a),abs(b)); v=[s[k]*(b//g)-r[k]*(a//g) for k in range(4)]
            if all(v): out.append([abs(x) for x in v])
    return max(mlc(out)) if out else None
random.seed(7)

# Theorem 15.1 / Corollary 15.2
cases=[]
while len(cases)<300:
    A=random.randint(1,9); C=random.randint(1,9); B=random.randint(1,30); H=random.randint(1,25)
    S=A+C
    if gcd(gcd(A,B),C)!=1 or gcd(B,S)!=1: continue
    cases.append((A,B,C,H))
v1=tml2([forms_H(*c) for c in cases]); v2=mlc([[c[0]+c[2],c[1],c[3]*(c[0]+c[2])] for c in cases])
bad=[(c,a,b) for c,a,b in zip(cases,v1,v2) if a!=b]
assert not bad
print("PASS: Proposition 10.3, ML(U_H) = ML(S,B,HS) on",len(cases),"random subtori with gcd(B,S) = 1")
# dependence only on S (also when gcd(B,S) > 1)
grp={}
cases2=[]
for S in range(2,12):
    for B in range(1,25):
        for H in (1,2,3,5,8,13):
            for A in range(1,S):
                C=S-A
                if gcd(gcd(A,B),C)!=1: continue
                cases2.append((A,B,C,H))
v=tml2([forms_H(*c) for c in cases2])
for c,x in zip(cases2,v):
    key=(c[0]+c[2],c[1],c[3]); grp.setdefault(key,set()).add(x)
assert all(len(x)==1 for x in grp.values())
print("PASS: Proposition 10.3, independence of A and C on",len(grp),"parameter groups")
# Corollary 15.4 covolume
ok=True
for (A,B,C,H) in cases2[:400]:
    s,r=forms_H(A,B,C,H); c2=sum((s[i]*r[j]-s[j]*r[i])**2 for i,j in itertools.combinations(range(4),2))
    ok&= c2==B*B+(1+H*H)*(A*A+C*C)
assert ok
# Theorem 16.1 / 16.2
cases3=[]
for N in range(12,101,4):
    for A in range(1,N):
        for C in range(1,N):
            B=N-A-C; S=A+C
            if B<S or gcd(gcd(A,B),C)!=1: continue
            if (A*31+C*17+N)%7: continue
            for j in (0,1,2,5,9):
                cases3.append((A,B,C,N//2+j*N))
v=tml2([forms_H(*c) for c in cases3]); badv=0; tri=[]; tkeys=[]
for (A,B,C,H),x in zip(cases3,v):
    N=A+B+C; S=A+C; Ds=F(B,2*(H*S+B))
    if x!=F(1,2)-min(F(1,N),Ds): badv+=1
    g=gcd(B,S)
    if Ds<=F(1,N) and (N//2)%g==0: tri.append([B//g,S//g,H*S//g]); tkeys.append(x)
tv=mlc(tri)
assert badv==0 and all(a==b for a,b in zip(tv,tkeys))
print("PASS: Theorem 10.4 on",len(cases3),"subtori, second-phase triple values on",len(tri))
# table for (2,55,3)
for H in (30,90,270,330,390):
    x=tml2([forms_H(2,55,3,H)])[0]; b=rep_best_H(2,55,3,H)
    assert (b==x)==(H>319)
    print("PASS: (2,55,3), H=%d: ML=%s, attained by a repeated-speed circle: %s"%(H,x,b==x))

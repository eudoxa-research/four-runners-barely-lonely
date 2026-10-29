# check_five_c.py -- Section 9.1 of "Four runners that are barely lonely": the contact formula in the proof of
# Lemma 9.13 against the C maximizer for all subtori U_n of the five tori with |n|_inf <= N (default 20),
# and the claim that every value in [1/3, 2/5) satisfies 1 <= 2q - 5p <= 4.
# Usage: python3 src/check_five_c.py ./torusml2 [N]   (expected for N = 20: 142547 subtori, 0 mismatches)
import itertools, subprocess, sys
from fractions import Fraction as F
from math import gcd, floor, ceil
TYPES = {
 'A': [(1,0,0),(4,0,0),(0,1,0),(0,0,1)],   # (y,4y,z,u)
 'B': [(2,0,0),(3,0,0),(0,1,0),(0,0,1)],   # (2y,3y,z,u)
 'C': [(1,0,0),(0,1,0),(3,1,0),(0,0,1)],   # (a,b,3a+b,u)
 'D': [(1,0,0),(0,1,0),(2,2,0),(0,0,1)],   # (a,b,2a+2b,u)
 'E': [(1,0,0),(0,1,0),(0,0,1),(2,1,1)],   # (a,b,c,2a+b+c)
}
h=F(1,2); f25=F(2,5); f35=F(3,5)
# local polytope P(eps) = {c0 + c1*eps : vertices}, given as list of (base point, direction) vertex pairs (affine in eps)
def boxverts(ivals):
    # ivals: list of (lo0, lo1, hi0, hi1) meaning [lo0+lo1*e, hi0+hi1*e]
    out=[]
    for choice in itertools.product([0,1],repeat=len(ivals)):
        b=[]; d=[]
        for c,(l0,l1,h0,h1) in zip(choice,ivals):
            if c==0: b.append(l0); d.append(l1)
            else: b.append(h0); d.append(h1)
        out.append((b,d))
    return out
def tri(p, V, extra=None):
    out=[]
    for v in V:
        b=list(p); d=list(v)
        out.append((b,d))
    if extra is None: return out
    res=[]
    for (b,d) in out:
        for (lo0,lo1,hi0,hi1) in [extra]:
            res.append((b+[lo0],d+[lo1])); res.append((b+[hi0],d+[hi1]))
    return res
U=(f25,-1,f35,1)  # u in [2/5-e, 3/5+e]
LOCAL = {
 'A': boxverts([(f25,-1,f25,F(1,4)), U, U]),
 'B': boxverts([(F(1,5),-h,F(1,5),F(1,3)), U, U]),
 'C': tri([f25,f25],[(-1,-1),(F(2,3),-1),(-1,4)], U),
 'D': tri([f25,f25],[(-1,-1),(F(3,2),-1),(-1,F(3,2))], U),
 'E': tri([f25,f25,f25],[(-1,-1,-1),(F(3,2),-1,-1),(-1,4,-1),(-1,-1,4)]),
}
def local_eps(n, verts):
    vals=[(sum(F(ni)*bi for ni,bi in zip(n,b)), sum(F(ni)*di for ni,di in zip(n,d))) for b,d in verts]
    # range [min_j (v0_j + e*v1_j), max_j(...)]; for small e the extreme vertices are those maximizing v0 then v1
    lo0=min(v0 for v0,v1 in vals); hi0=max(v0 for v0,v1 in vals)
    a=max(-v1 for v0,v1 in vals if v0==lo0); b=max(v1 for v0,v1 in vals if v0==hi0)
    if floor(hi0)>=ceil(lo0): return F(0)
    cand=[]
    if b>0: cand.append((ceil(hi0)-hi0)/b)
    if a>0: cand.append((lo0-floor(lo0))/a)
    return min(cand)
def kernel3(n):
    E=list(n); T=[[1,0,0],[0,1,0],[0,0,1]]
    while sum(1 for e in E if e)>1:
        piv=min((i for i in range(3) if E[i]),key=lambda i:abs(E[i]))
        for i in range(3):
            if i==piv or E[i]==0: continue
            q=E[i]//E[piv]; E[i]-=q*E[piv]; T[i]=[x-q*y for x,y in zip(T[i],T[piv])]
    return [T[i] for i in range(3) if E[i]==0]
EXE=sys.argv[1] if len(sys.argv)>1 else './torusml2'
N=int(sys.argv[2]) if len(sys.argv)>2 else 20
for name,phi in TYPES.items():
    lines=[]; ns=[]
    for n in itertools.product(range(-N,N+1),repeat=3):
        if n==(0,0,0) or gcd(gcd(abs(n[0]),abs(n[1])),abs(n[2]))!=1: continue
        if n < tuple(-x for x in n): continue
        c1,c2=kernel3(n)
        s=[sum(p*c for p,c in zip(ph,c1)) for ph in phi]; r=[sum(p*c for p,c in zip(ph,c2)) for ph in phi]
        if any(s[i]==0 and r[i]==0 for i in range(4)): continue
        lines.append(' '.join(f"{s[i]} {r[i]}" for i in range(4))); ns.append(n)
    out=subprocess.run([EXE],input='\n'.join(lines)+'\n',capture_output=True,text=True).stdout.split('\n')
    bad=0; nlow=0; rcount={}
    for n,o in zip(ns,out):
        p,q=map(int,o.split()[:2]); v=F(p,q)
        e=local_eps(n,LOCAL[name])
        if e<=F(1,15):
            if v!=f25-e: bad+=1; print("MISMATCH",name,n,v,f25-e) if bad<5 else None
        else:
            if not v<F(1,3): bad+=1; print("NOT<1/3",name,n,v,f25-e) if bad<5 else None
        if F(1,3)<=v<f25:
            nlow+=1; r=2*v.denominator-5*v.numerator; rcount[r]=rcount.get(r,0)+1
    print(name,"subtori",len(ns),"mismatches",bad,"values in [1/3,2/5):",nlow,"r-distribution",dict(sorted(rcount.items())))

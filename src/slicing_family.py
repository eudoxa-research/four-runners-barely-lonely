# Remark 9.12: the family {lambda.y = 0, mu u = theta.y}; usage: python3 slicing_family.py 42 4 3 allg > fam.txt
# then ./torusml2 < fam.txt > fam_ml.txt ; ./allmax fam.txt fam_ml.txt 0.0 0.5 | grep -c ALLBAD   (expected 602 of 19000)
# Generate 2-dim subtori U = ker(m1) ∩ ker(m2), m1=(l1,l2,l3,0), m2=(-t1,-t2,-t3,mu), with lambda having all gcd(lambda_c, Lambda)>1
# or small, theta in a box, mu small; output forms (S_i,R_i) for torusml.
import sys, itertools
from math import gcd
def kernel_basis(M):
    # M: 2x4 integer; returns 4x2 integer matrix whose columns form a basis of Z^4 ∩ ker M
    rows=[list(r) for r in M]; n=4
    U=[[1 if i==j else 0 for j in range(n)] for i in range(n)]
    def colop(j,k,q):  # col_j -= q*col_k
        for r in rows: r[j]-=q*r[k]
        for r in U: r[j]-=q*r[k]
    def swap(j,k):
        for r in rows: r[j],r[k]=r[k],r[j]
        for r in U: r[j],r[k]=r[k],r[j]
    piv=0
    for i in range(2):
        # reduce row i on columns piv..n-1 to a single nonzero at column piv
        while True:
            nz=[j for j in range(piv,n) if rows[i][j]!=0]
            if not nz: break
            jm=min(nz,key=lambda j:abs(rows[i][j]))
            if jm!=piv: swap(piv,jm)
            done=True
            for j in range(piv+1,n):
                if rows[i][j]!=0:
                    q=rows[i][j]//rows[i][piv]; colop(j,piv,q)
                    if rows[i][j]!=0: done=False
            if done: break
        if any(rows[i][j]!=0 for j in range(piv,n)): piv+=1
    return [[U[r][c] for c in range(piv,n)] for r in range(n)]
seen=set(); out=[]
Lmax=int(sys.argv[1]); tbox=int(sys.argv[2]); mumax=int(sys.argv[3]); mode=sys.argv[4] if len(sys.argv)>4 else 'allg'
for l1 in range(1,Lmax):
    for l2 in range(l1,Lmax):
        for l3 in range(l2,Lmax):
            Lam=l1+l2+l3
            if Lam>Lmax or gcd(gcd(l1,l2),l3)!=1: continue
            gs=[gcd(l,Lam) for l in (l1,l2,l3)]
            if mode=='allg' and min(gs)==1: continue
            for mu in range(1,mumax+1):
                for t in itertools.product(range(-tbox,tbox+1),repeat=3):
                    M=[[l1,l2,l3,0],[-t[0],-t[1],-t[2],mu]]
                    # primitive row lattice: gcd of 2x2 minors == 1
                    g=0
                    for i,j in itertools.combinations(range(4),2): g=gcd(g,M[0][i]*M[1][j]-M[0][j]*M[1][i])
                    if g!=1: continue
                    K=kernel_basis(M)
                    if len(K[0])!=2: continue
                    S=[K[i][0] for i in range(4)]; R=[K[i][1] for i in range(4)]
                    if any(S[i]==0 and R[i]==0 for i in range(4)): continue
                    key=tuple(S+R)
                    if key in seen: continue
                    seen.add(key)
                    print(' '.join(f"{S[i]} {R[i]}" for i in range(4)))

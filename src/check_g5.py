# check_g5.py -- checks for Theorem 11.5 (five speeds), using the exact maximizer torusmlN.
# Usage (from the package root, after gcc -O2 -o torusmlN src/torusmlN.c): python3 src/check_g5.py
# 1. exact ML of all W = Zq + Zc with |c_i| <= 3 for the eight seeds q: every W with ML < 2/9 passes the
#    repeated-speed test of the enumeration and has ML = 1/5;
# 2. 20,000 random primitive planes in Z^5 with entries in [-3,3]: no value in (1/5, 2/9), and every
#    value 1/5 comes from a product subtorus.
import itertools, random, subprocess
from fractions import Fraction as F
from math import gcd
def gl(v):
    g=0
    for x in v: g=gcd(g,abs(x))
    return g
def ml_many(lattices):
    inp='\n'.join('5 '+' '.join(f'{b1[k]} {b2[k]}' for k in range(5)) for b1,b2 in lattices)+'\n'
    out=subprocess.run(['./torusmlN'],input=inp,capture_output=True,text=True).stdout.split('\n')
    return [F(int(o.split()[0]),int(o.split()[1])) for o in out if o.strip()]
def is_product(b1,b2):
    rows=[(b1[k],b2[k]) for k in range(5)]
    for k in range(5):
        others=[rows[i] for i in range(5) if i!=k]
        if all(a[0]*b[1]-a[1]*b[0]==0 for a in others for b in others):
            # other four rows proportional: speeds of that circle
            base=next(r for r in others)
            g=gl(base); u=(base[0]//g, base[1]//g)
            # ratio of each row to u
            vals=sorted({abs(r[0]//u[0]) if u[0] else abs(r[1]//u[1]) for r in others})
            g2=gl(vals); vals=[x//g2 for x in vals]
            return vals in ([1,2,3,4],[1,3,4,7])
    return False
SEEDS=[(1,1,2,3,4),(1,2,2,3,4),(1,2,3,3,4),(1,2,3,4,4),(1,1,3,4,7),(1,3,3,4,7),(1,3,4,4,7),(1,3,4,7,7)]
def tightset(v):
    d=sorted(set(v))
    if len(d)!=4: return False
    g=gl(d); d=[x//g for x in d]
    return d in ([1,2,3,4],[1,3,4,7])
def passes(q,c):
    for i in range(5):
        for j in range(i+1,5):
            for e in (1,-1):
                d1=q[i]-e*q[j]; d2=c[i]-e*c[j]
                if d1==0 and d2==0: continue
                g=gcd(abs(d1),abs(d2)); u1,u2=-d2//g,d1//g
                s=[q[k]*u1+c[k]*u2 for k in range(5)]
                if 0 in s: continue
                if not tightset([abs(x) for x in s]): return False
    return True
# 1. Soundness of the filter on a box: exact ML of every W=span(q,c), |c_k|<=3, c_0=0, gcd=1, sign-normalized
lat=[]; meta=[]
for q in SEEDS:
    for cc in itertools.product(range(-3,4),repeat=4):
        c=(0,)+cc
        if gl(cc)!=1: continue
        fnz=next(x for x in cc if x!=0)
        if fnz<0: continue
        lat.append((q,c)); meta.append(passes(q,c))
vals=ml_many(lat)
bad=0; below=0
for (q,c),ok,m in zip(lat,meta,vals):
    if m<F(2,9):
        below+=1
        if not ok: bad+=1; print('FILTER WOULD LOSE', q, c, m)
        if m!=F(1,5): bad+=1; print('G5 COUNTEREXAMPLE', q, c, m)
print(f'box check: {len(lat)} seeded subtori, {below} with ML<2/9, all pass the filter and have ML=1/5: {bad==0}')
# 2. Random primitive planes in Z^5 (no seed assumption): no value in (1/5,2/9); ML=1/5 only for products
random.seed(2026); lat=[]
while len(lat)<20000:
    b1=[random.randint(-3,3) for _ in range(5)]; b2=[random.randint(-3,3) for _ in range(5)]
    if any(b1[k]==0 and b2[k]==0 for k in range(5)): continue
    mins=[b1[i]*b2[j]-b1[j]*b2[i] for i in range(5) for j in range(i+1,5)]
    if gl(mins)!=1: continue
    lat.append((b1,b2))
vals=ml_many(lat); bad=0; crit=0; low=min(vals)
for (b1,b2),m in zip(lat,vals):
    if F(1,5)<m<F(2,9): bad+=1; print('GAP COUNTEREXAMPLE',b1,b2,m)
    if m<F(1,5): bad+=1; print('BELOW 1/5',b1,b2,m)
    if m==F(1,5):
        crit+=1
        if not is_product(b1,b2): bad+=1; print('NONPRODUCT CRITICAL',b1,b2)
nxt=min(v for v in vals if v>F(1,5))
print(f'random check: {len(lat)} planes, {crit} critical (all products), smallest value above 1/5 = {nxt}, failures {bad}')

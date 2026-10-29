#!/usr/bin/env python3
"""Checks Lemmas 9.2 and 9.4 (Section 9) against the exact C maximizer torusml on larger ranges
than the Python script: all primitive n = (n1, n2, n3) with |n1| <= 12, |n2|, |n3| <= 6 giving
coordinate-free subtori of W^1 and W^2 (3,356 and 3,354 subtori), and F(n1, n2) for |n_i| <= 12
(624 pairs).  Usage: python3 check_lemmas_c.py ./torusml   (torusml built from src/torusml.c)."""
import sys, subprocess
from fractions import Fraction as F
from math import gcd

def kernel3(n):
    E = list(n); T = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    while sum(1 for e in E if e) > 1:
        piv = min((i for i in range(3) if E[i]), key=lambda i: abs(E[i]))
        for i in range(3):
            if i == piv or E[i] == 0: continue
            q = E[i] // E[piv]; E[i] -= q * E[piv]
            T[i] = [a - q * b for a, b in zip(T[i], T[piv])]
    return [T[i] for i in range(3) if E[i] == 0]

def ml_exact(v):
    v = sorted(set(v)); best = F(0)
    for i, u in enumerate(v):
        for w in v[i:]:
            D = u + w
            for m in range(1, D // 2 + 1):
                mn = min(min(x * m % D, D - x * m % D) for x in v)
                if F(mn, D) > best: best = F(mn, D)
    return best

def subtori(f, A, B):
    out = []
    for n1 in range(-A, A + 1):
        for n2 in range(-B, B + 1):
            for n3 in range(-B, B + 1):
                nn = (n1, n2, n3)
                if nn == (0, 0, 0) or gcd(gcd(abs(n1), abs(n2)), abs(n3)) != 1: continue
                c1, c2 = kernel3(nn)
                s = tuple(sum(c1[i] * f[i][k] for i in range(3)) for k in range(4))
                r = tuple(sum(c2[i] * f[i][k] for i in range(3)) for k in range(4))
                if any(s[k] == 0 and r[k] == 0 for k in range(4)): continue
                out.append((nn, s, r))
    return out

def run(torusml, items):
    inp = "\n".join(" ".join(f"{s[k]} {r[k]}" for k in range(4)) for _, s, r in items) + "\n"
    res = subprocess.run([torusml], input=inp, capture_output=True, text=True, check=True)
    return [F(int(l.split()[0]), int(l.split()[1])) for l in res.stdout.strip().split("\n")]

def pred_F(n1, n2):
    if (n1 + n2) % 3: return F(1, 3)
    M = max(abs(n1 + n2), abs(2 * n1 - n2), abs(2 * n2 - n1)); return F(M, 3 * (M + 1))

torusml = sys.argv[1] if len(sys.argv) > 1 else "./torusml"
# W^1 = {(2y, y, z, u)}, Lemma 9.2
items = subtori([(2, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)], 12, 6); ex = run(torusml, items)
for (nn, s, r), e in zip(items, ex):
    n1, n2, n3 = nn
    if (n2, n3) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        pred = F(abs(n1) // 3, abs(n1) + 1) if abs(n1) % 3 == 0 else F(1, 3)
    else: pred = F(1, 3)
    assert e == pred, (nn, e, pred)
print(f"PASS: Lemma 9.2 on {len(items)} subtori of W^1")
# W^2 = {(a, b, a+b, u)}, Lemma 9.4
items = subtori([(1, 0, 1, 0), (0, 1, 1, 0), (0, 0, 0, 1)], 12, 6); ex = run(torusml, items)
for (nn, s, r), e in zip(items, ex):
    n1, n2, n3 = nn
    if abs(n3) >= 2: pred = F(1, 3)
    elif n3 == 0: pred = ml_exact([abs(n2), abs(n1), abs(n1 - n2)])
    else: pred = pred_F(-n3 * n1, -n3 * n2)
    assert e == pred, (nn, e, pred)
print(f"PASS: Lemma 9.4 on {len(items)} subtori of W^2")
# F(n1, n2), |n_i| <= 12
pairs = [(n1, n2) for n1 in range(-12, 13) for n2 in range(-12, 13) if (n1, n2) != (0, 0)]
items = [((n1, n2, 0), (1, 0, 1, n1), (0, 1, 1, n2)) for n1, n2 in pairs]; ex = run(torusml, items)
for (n1, n2), e in zip(pairs, ex):
    assert e == pred_F(n1, n2), (n1, n2, e)
print(f"PASS: formula for F(n1,n2) on {len(pairs)} pairs")

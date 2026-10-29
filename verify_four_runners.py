#!/usr/bin/env python3
"""Verification script for "Four runners that are barely lonely".

Standard library only.  Exact rational arithmetic throughout.

Checks performed (each prints PASS or raises AssertionError):
  1. The exact formula for ML(x, y, x+y, x+2y) against independent corner maximization.
  2. The formula for ML(a, b, 2b, 3b).
  3. The 21 classes of characters m with |m|^2 <= 19 and the certificates of Lemma 4.2.
  4. The certificates for all 123 surviving 2-dimensional subtori (file tori_survivors.json).
  5. Independent pure-Python enumeration of all primitive quadruples with max speed <= 60,
     compared with the classification of Theorem 1.1.
  6. The absence of m = 3, and the witnesses for m = 2, m = 5 and all m >= 4.
  7. Optional: re-check of the output of the exhaustive C search (file b512.txt) if present.
  8. The formulas of Section 9 (subtori of W^1 and W^2, and F(n1,n2)) against an exact
     two-dimensional maximizer written in Python.
  9. The by-product on triples: every primitive triple with max <= 80 and ML < 1/3 has a
     relation c = 2b or c = a + b (up to order).
 10. Lemmas 9.10 and 9.11 (generic case below 2/5): the arithmetic of the covolume bound, and the five types of
     characters with |m|_1 = 5 with their parametrizations.
 11. The subtorus of the remark on slicing: ML = 13/30, maximizers off every slice circle,
     and ML(1,8,29) = 13/30.
 12. Optional: the 18,960 survivors of Proposition 9.14 (file surv25.txt): witness points,
     values in R, reduction to orbits under signed permutations, and the exact value of every orbit.
 13. Lemma 9.13 (subtori of the five tori): the contact-interval formula against the exact
     two-dimensional maximizer, and 1 <= 2q-5p <= 4 for every value in [1/3, 2/5).
 14. Lemma 9.12 (superlevel sets of the five tori) on a rational grid.
 15. Theorem 9.8(ii): ML(2,3,5j) and ML(1,4,5j+3), and the two descriptions of R.
 16. The example at the level 3/7: a subtorus of W_(4,1,1,1) with ML = 43/103 = ML(1,102,25).
 17. Section 10: ML(1,N/2,N-1) = 1/2 - 1/N, ML(U(A,B,C)) for small N by the exact maximizer, the repeated-speed
     circles of U(A,B,C), the formula and comparison of Theorem 10.2, Proposition 10.3 and Theorem 10.4.
"""
import json, os, sys, itertools
from fractions import Fraction as F
from math import gcd

def dist(z):
    z = z - (z.numerator // z.denominator)
    return min(z, 1 - z)

def ml_exact(v):
    """Maximum loneliness of a list of positive integers (corner enumeration, exact)."""
    v = sorted(set(v)); best = F(0)
    dens = sorted({u + w for i, u in enumerate(v) for w in v[i:]})
    for D in dens:
        for m in range(1, D // 2 + 1):
            s = D
            for u in v:
                z = (u * m) % D; z = min(z, D - z)
                if z * best.denominator <= best.numerator * D: break
                s = min(s, z)
            else:
                if s * best.denominator > best.numerator * D: best = F(s, D)
    return best

# ---------- 1. the family V(x,y) = {x, y, x+y, x+2y} ----------
def ml_formula_U2(x, y):
    r = F((y - x) % 4, 4); r3 = F((y - 2 * x) % 4, 4)
    if r == 0 or r3 == 0: return F(1, 4)
    e = min(r / (x + y), (1 - r) / (x + 3 * y), r3 / (2 * x + y), (1 - r3) / (2 * x + 3 * y))
    return F(1, 4) - e

def m_of(x, y):
    r = (y - x) % 4; r3 = (y - 2 * x) % 4
    if r == 0 or r3 == 0: return None
    m1 = {1: (x + y - 1) // 2, 3: (x + 3 * y - 1) // 2, 2: (x + 3 * y - 2) // 4}[r]
    m3 = {1: (2 * x + y - 1) // 2, 3: (2 * x + 3 * y - 1) // 2, 2: (2 * x + 3 * y - 2) // 4}[r3]
    return max(m1, m3)

def check1(B=40):
    n = 0
    for x in range(1, B + 1):
        for y in range(1, B + 1):
            if x == y or gcd(x, y) != 1: continue
            v = [x, y, x + y, x + 2 * y]
            ex = ml_exact(v); fo = ml_formula_U2(x, y)
            assert ex == fo, (x, y, ex, fo)
            m = m_of(x, y)
            if m is None: assert ex == F(1, 4)
            else:
                assert ex == F(m, 4 * m + 2) and m >= 2 and m != 3, (x, y, ex, m)
            n += 1
    print(f"PASS 1: formula for V(x,y) agrees with exact ML on {n} coprime pairs with x,y <= {B}")

# ---------- 2. the family {a, b, 2b, 3b} ----------
def check2(A=60, Bb=10):
    n = 0
    for b in range(1, Bb + 1):
        for a in range(1, A + 1):
            if gcd(a, b) != 1 or a in (b, 2 * b, 3 * b): continue
            ex = ml_exact([a, b, 2 * b, 3 * b])
            pred = F(a, 4 * a + 4) if (b == 1 and a % 4 == 0) else None
            if pred is None: assert ex >= F(1, 4), (a, b, ex)
            else: assert ex == pred, (a, b, ex, pred)
            n += 1
    print(f"PASS 2: family {{a,b,2b,3b}}: {n} cases, below 1/4 exactly for b=1, 4|a, value a/(4a+4)")

# ---------- 3. characters of 3-dimensional subtori ----------
def canon_vec(m):
    best = None
    for perm in itertools.permutations(range(4)):
        for signs in itertools.product((1, -1), repeat=4):
            q = tuple(signs[i] * m[perm[i]] for i in range(4))
            for a in q:
                if a != 0:
                    if a < 0: q = tuple(-b for b in q)
                    break
            if best is None or q < best: best = q
    return best

def check3():
    classes = set()
    for m in itertools.product(range(-4, 5), repeat=4):
        if sum(a * a for a in m) > 19 or sum(1 for a in m if a) < 2: continue
        g = 0
        for a in m: g = gcd(g, abs(a))
        if g != 1: continue
        classes.add(canon_vec(m))
    assert len(classes) == 21, len(classes)
    special = []
    for m in sorted(classes):
        found = None
        for Q in (5, 10, 20, 15, 30):
            good = [a for a in range(Q) if 5 * min(a, Q - a) >= 2 * Q]
            for xs in itertools.product(good, repeat=4):
                if sum(mi * a for mi, a in zip(m, xs)) % Q == 0:
                    found = (Q, xs); break
            if found: break
        if found is None: special.append(m)
        else:
            Q, xs = found
            assert all(5 * min(a, Q - a) >= 2 * Q for a in xs)
    assert special == [(0, 0, 1, -2), (0, 1, -1, -1)], special
    print("PASS 3: 21 classes of characters with |m|^2<=19; all but x3=2x4 and x2=x3+x4 have a point with G>=2/5")

# ---------- 4. certificates for the surviving 2-dimensional subtori ----------
def check4(path):
    data = json.load(open(path))
    assert len(data) == 123
    low = [d for d in data if F(*d["ml"]) <= F(5, 16)]
    assert len(low) == 31
    tight = 0
    for d in data:
        s, r = d["s"], d["r"]; X, Y, D = d["point"]
        val = min(min((s[i] * X + r[i] * Y) % D, D - (s[i] * X + r[i] * Y) % D) for i in range(4))
        assert F(val, D) == F(*d["ml"]), d          # the point attains the claimed value
        assert all(not (s[i] == 0 and r[i] == 0) for i in range(4))   # coordinate-free
        cv2 = sum((s[i] * r[j] - s[j] * r[i]) ** 2 for i in range(4) for j in range(i + 1, 4))
        assert cv2 == d["covol2"]
        if F(*d["ml"]) == F(1, 4): tight += 1
        elif F(*d["ml"]) <= F(5, 16): assert F(val, D) >= F(2, 7)      # lower bound used in Section 6
        else: assert F(val, D) > F(5, 16)                                # certificate of exclusion
    assert tight == 2
    assert max(d["covol2"] for d in low) == 422
    vals = sorted({F(*d["ml"]) for d in low})
    assert vals == [F(1, 4), F(2, 7), F(3, 10), F(4, 13), F(5, 16)]
    print("PASS 4: 123 surviving subtori: 2 tight, 29 with certificate >= 2/7 and ML <= 5/16 (max covol^2 = 422), 92 excluded by a point > 5/16")

# ---------- 5. independent enumeration up to max speed 60 ----------
def in_U1(v):
    s = set(v)
    return any({b, 2 * b, 3 * b} <= s and len(s - {b, 2 * b, 3 * b}) == 1 for b in v)

def in_U2(v):
    sv = sorted(v)
    for x in range(1, max(v) + 1):
        for y in range(1, max(v) + 1):
            if sorted([x, y, x + y, x + 2 * y]) == sv: return True
    return False

def check5(N=60):
    found = []
    for v in itertools.combinations(range(1, N + 1), 4):
        if gcd(gcd(v[0], v[1]), gcd(v[2], v[3])) != 1: continue
        # quick certificate at t = 1/4 (speeds not divisible by 4) and t = 3/8
        if all(a % 4 for a in v): continue
        if all(min((3 * a) % 8, 8 - (3 * a) % 8) >= 2 for a in v): continue
        ml = ml_exact(list(v))
        if ml < F(1, 4): found.append((v, ml))
    for v, ml in found:
        if not in_U1(v) and not in_U2(v): assert v == (1, 3, 4, 14), v
        m = 2 * ml / (1 - 4 * ml)
        assert m.denominator == 1 and m >= 2 and m != 3, (v, ml)
    assert any(v == (1, 3, 4, 14) and ml == F(4, 17) for v, ml in found)
    print(f"PASS 5: all {len(found)} primitive quadruples with max speed <= {N} and ML < 1/4 are as classified")

# ---------- 6. the value set ----------
def check6():
    # m = 3 would require a pair with x + y <= 14 (Lemma 7.4); exhaustive check far beyond that
    for x in range(1, 61):
        for y in range(1, 61):
            if x == y or gcd(x, y) != 1: continue
            assert m_of(x, y) != 3
    assert ml_exact([1, 2, 3, 4]) == F(1, 5) and ml_exact([1, 3, 4, 7]) == F(1, 5)
    assert ml_exact([1, 7, 8, 15]) == F(5, 22)
    for s in range(0, 6):      # Fan-Sun family V(4s+3, 8): m = 2s+7
        assert ml_exact([8, 4 * s + 3, 4 * s + 11, 4 * s + 19]) == F(2 * s + 7, 8 * s + 30)
    for k in range(1, 8):      # {1,2,3,4k}: m = 2k
        assert ml_exact([1, 2, 3, 4 * k]) == F(k, 4 * k + 1)
    print("PASS 6: m=3 never occurs for x,y<=60; witnesses for m=2, m=5, odd m>=7 and even m verified")

# ---------- 7. optional re-check of the exhaustive search output ----------
def check7(path):
    rows = [tuple(map(int, l.split())) for l in open(path)]
    assert len(rows) == 14032
    for (a, b, c, d, p, q) in rows:
        v = (a, b, c, d)
        assert a * a + b * b + c * c + d * d < 512 ** 2 and gcd(gcd(a, b), gcd(c, d)) == 1
        ml = F(p, q); m = 2 * ml / (1 - 4 * ml)
        assert m.denominator == 1 and m >= 2 and m != 3
        if not in_U1(v) and not in_U2(v): assert v == (1, 3, 4, 14) and ml == F(4, 17)
    # recompute exact ML for a random sample
    import random; random.seed(0)
    for (a, b, c, d, p, q) in random.sample(rows, 200):
        assert ml_exact([a, b, c, d]) == F(p, q)
    print("PASS 7: exhaustive-search output: 14032 quadruples, all classified; 200 random exact re-checks")


# ---------- 8. exact loneliness of a 2-dimensional subtorus (small coefficients) ----------
def ml_torus(forms):
    """Exact max over T^2 of min_i ||s_i x + r_i y|| for forms = [(s_i, r_i)], by enumerating
    the candidate points where three signed branches are equal, where two proportional forms are at a
    common level and a third is at a half integer, or where two independent forms are at half integers."""
    forms = [(int(s), int(r)) for s, r in forms]
    def Fval(x, y): return min(dist(s * x + r * y) for s, r in forms)
    best = F(0)
    B = [abs(s) + abs(r) + 1 for s, r in forms]
    def det3(A): return (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]) + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))
    def solve3(A, b):
        d = det3(A)
        if d == 0: return None
        out = []
        for c in range(3):
            M = [list(row) for row in A]
            for i in range(3): M[i][c] = b[i]
            out.append(F(det3(M), d))
        return out
    for (i, j, m) in itertools.combinations(range(4), 3):
        for sg in itertools.product((1, -1), repeat=3):
            A = [[forms[i][0], forms[i][1], -sg[0]], [forms[j][0], forms[j][1], -sg[1]], [forms[m][0], forms[m][1], -sg[2]]]
            if det3(A) == 0: continue
            for ki in range(-B[i], B[i] + 1):
                for kj in range(-B[j], B[j] + 1):
                    for km in range(-B[m], B[m] + 1):
                        x, y, d = solve3(A, (ki, kj, km))
                        if d <= best or d > F(1, 2): continue
                        best = max(best, Fval(x, y))
    for (i, j) in itertools.combinations(range(4), 2):
        if forms[i][0] * forms[j][1] != forms[j][0] * forms[i][1]: continue
        for m in range(4):
            if m in (i, j): continue
            for sg in itertools.product((1, -1), repeat=2):
                A = [[forms[i][0], forms[i][1], -sg[0]], [forms[j][0], forms[j][1], -sg[1]], [forms[m][0], forms[m][1], 0]]
                if det3(A) == 0: continue
                for ki in range(-B[i], B[i] + 1):
                    for kj in range(-B[j], B[j] + 1):
                        for km in range(-B[m], B[m] + 1):
                            x, y, d = solve3(A, (ki, kj, F(2 * km + 1, 2)))
                            if d <= best or d > F(1, 2): continue
                            best = max(best, Fval(x, y))
    for (i, j) in itertools.combinations(range(4), 2):
        det = forms[i][0] * forms[j][1] - forms[j][0] * forms[i][1]
        if det == 0: continue
        for ki in range(-B[i], B[i] + 1):
            for kj in range(-B[j], B[j] + 1):
                bi, bj = F(2 * ki + 1, 2), F(2 * kj + 1, 2)
                x = F(bi * forms[j][1] - bj * forms[i][1], det); y = F(forms[i][0] * bj - forms[j][0] * bi, det)
                best = max(best, Fval(x, y))
    return best

def kernel3(n):
    """Basis of {c in Z^3 : n.c = 0} for primitive n, by integer row reduction."""
    E = list(n); T = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    while sum(1 for e in E if e) > 1:
        piv = min((i for i in range(3) if E[i]), key=lambda i: abs(E[i]))
        for i in range(3):
            if i == piv or E[i] == 0: continue
            q = E[i] // E[piv]; E[i] -= q * E[piv]
            T[i] = [a - q * b for a, b in zip(T[i], T[piv])]
    return [T[i] for i in range(3) if E[i] == 0]

def check8():
    # F(n1,n2) = {(a,b,a+b,n1 a+n2 b)}: forms (1,0),(0,1),(1,1),(n1,n2)
    n = 0
    for n1 in range(-4, 5):
        for n2 in range(-4, 5):
            if (n1, n2) == (0, 0): continue
            ex = ml_torus([(1, 0), (0, 1), (1, 1), (n1, n2)])
            if (n1 + n2) % 3: pred = F(1, 3)
            else:
                M = max(abs(n1 + n2), abs(2 * n1 - n2), abs(2 * n2 - n1)); pred = F(M, 3 * (M + 1))
            assert ex == pred, (n1, n2, ex, pred); n += 1
    # subtori of W^1 = {(2y,y,z,u)}: U_n for primitive n with small entries
    f = [(2, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    m = 0
    for n1 in range(-6, 7):
        for n2 in range(-3, 4):
            for n3 in range(-3, 4):
                nn = (n1, n2, n3)
                if nn == (0, 0, 0) or gcd(gcd(abs(n1), abs(n2)), abs(n3)) != 1: continue
                c1, c2 = kernel3(nn)
                s = tuple(sum(c1[i] * f[i][k] for i in range(3)) for k in range(4))
                r = tuple(sum(c2[i] * f[i][k] for i in range(3)) for k in range(4))
                if any(s[k] == 0 and r[k] == 0 for k in range(4)): continue
                ex = ml_torus(list(zip(s, r)))
                if (n2, n3) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    pred = F(abs(n1) // 3, abs(n1) + 1) if abs(n1) % 3 == 0 else F(1, 3)
                else: pred = F(1, 3)
                assert ex == pred, (nn, ex, pred); m += 1
    # subtori of W^2 = {(a,b,a+b,u)}: U_n for primitive n = (n1,n2,n3), all three cases of the lemma
    f = [(1, 0, 1, 0), (0, 1, 1, 0), (0, 0, 0, 1)]
    w = 0
    for n1 in range(-4, 5):
        for n2 in range(-4, 5):
            for n3 in range(-3, 4):
                nn = (n1, n2, n3)
                if nn == (0, 0, 0) or gcd(gcd(abs(n1), abs(n2)), abs(n3)) != 1: continue
                c1, c2 = kernel3(nn)
                s = tuple(sum(c1[i] * f[i][k] for i in range(3)) for k in range(4))
                r = tuple(sum(c2[i] * f[i][k] for i in range(3)) for k in range(4))
                if any(s[k] == 0 and r[k] == 0 for k in range(4)): continue
                ex = ml_torus(list(zip(s, r)))
                if abs(n3) >= 2: pred = F(1, 3)                                   # case (i)
                elif n3 == 0: pred = ml_exact([abs(n2), abs(n1), abs(n1 - n2)])    # case (ii): product with a circle
                else:                                                             # case (iii): F(-n3 n1, -n3 n2)
                    m1, m2 = -n3 * n1, -n3 * n2
                    if (m1 + m2) % 3: pred = F(1, 3)
                    else:
                        M = max(abs(m1 + m2), abs(2 * m1 - m2), abs(2 * m2 - m1)); pred = F(M, 3 * (M + 1))
                assert ex == pred, (nn, ex, pred); w += 1
    print(f"PASS 8: Lemma on F(n1,n2) verified on {n} parameter pairs, Lemma on subtori of W^1 on {m} subtori, Lemma on subtori of W^2 on {w} subtori")

def check9(N=80):
    low = 0
    for v in itertools.combinations(range(1, N + 1), 3):
        a, b, c = v
        if gcd(gcd(a, b), c) != 1: continue
        if ml_exact(list(v)) >= F(1, 3): continue
        low += 1
        assert b == 2 * a or c == 2 * a or c == 2 * b or c == a + b, v
    print(f"PASS 9: all {low} primitive triples with max <= {N} and ML < 1/3 have a relation c=2b or c=a+b")

def ml_above(v, thr):
    """True if some corner time gives min_i ||v_i t|| > thr (so ML(v) > thr)."""
    for i, u in enumerate(v):
        for w in v[i:]:
            D = u + w
            for m in range(1, D // 2 + 1):
                ok = True
                for x in v:
                    z = (x * m) % D; z = min(z, D - z)
                    if z * thr.denominator <= thr.numerator * D: ok = False; break
                if ok: return True
    return False

def check10():
    import math
    # covolume bound of the generic case: (35/2)^2 * 37 = 11331.25, and the other cases are smaller
    maxsq = {}
    for m in itertools.product(range(-15, 16), repeat=4):
        l1 = sum(abs(x) for x in m)
        if l1 in (7, 9, 11, 13, 15) and gcd(gcd(abs(m[0]), abs(m[1])), gcd(abs(m[2]), abs(m[3]))) == 1:
            maxsq[l1] = max(maxsq.get(l1, 0), sum(x * x for x in m))
    assert maxsq == {7: 37, 9: 65, 11: 101, 13: 145, 15: 197}, maxsq
    bound = F(35, 2) ** 2 * 37
    assert bound == F(45325, 4)                                   # 11331.25
    for s_ in (4, 5, 6):
        assert F(5 * (2 * s_ + 1) ** 2, 2 * (s_ - 2)) ** 2 < bound
    assert F(225, 2) ** 2 / 3 < bound                             # s = 7: A < (15/2)(15/sqrt 3)
    assert F(2500, 3) < bound                                     # even |m|_1: A < 50/sqrt 3
    assert all(3 * (s_ - 2) ** 2 >= 100 for s_ in range(8, 40))   # s >= 8 impossible
    # the five types: every primitive m with |m|_1 = 5 is a signed permutation of one of them,
    # and each parametrization is an injective homomorphism T^3 -> T^4 onto W_m
    types = {}
    for m in itertools.product(range(-5, 6), repeat=4):
        if sum(abs(x) for x in m) == 5 and gcd(gcd(abs(m[0]), abs(m[1])), gcd(abs(m[2]), abs(m[3]))) == 1:
            types.setdefault(tuple(sorted((abs(x) for x in m), reverse=True)), 0)
            types[tuple(sorted((abs(x) for x in m), reverse=True))] += 1
    assert set(types) == {(4, 1, 0, 0), (3, 2, 0, 0), (3, 1, 1, 0), (2, 2, 1, 0), (2, 1, 1, 1)}, types
    for mm, phi in FIVE.items():
        cols = [[phi[i][k] for i in range(4)] for k in range(3)]
        assert all(sum(mm[i] * phi[i][k] for i in range(4)) == 0 for k in range(3))
        minors = [abs(det3([phi[i] for i in rows])) for rows in itertools.combinations(range(4), 3)]
        g = 0
        for x in minors: g = gcd(g, x)
        assert g == 1, (mm, minors)
    print(f"PASS 10: generic covolume bound 11331.25; |m|_1 = 5 gives the five types ({sum(types.values())} vectors); parametrizations onto W_m")

def det3(A):
    return (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]) + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))

# the five tori: m and the coordinate forms phi_i (rows) of the parametrization T^3 -> W_m
FIVE = {
    (4, -1, 0, 0): [(1, 0, 0), (4, 0, 0), (0, 1, 0), (0, 0, 1)],    # A: (y,4y,z,u)
    (3, -2, 0, 0): [(2, 0, 0), (3, 0, 0), (0, 1, 0), (0, 0, 1)],    # B: (2y,3y,z,u)
    (3, 1, -1, 0): [(1, 0, 0), (0, 1, 0), (3, 1, 0), (0, 0, 1)],    # C: (a,b,3a+b,u)
    (2, 2, -1, 0): [(1, 0, 0), (0, 1, 0), (2, 2, 0), (0, 0, 1)],    # D: (a,b,2a+2b,u)
    (2, 1, 1, -1): [(1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 1, 1)],    # E: (a,b,c,2a+b+c)
}
H, Q25, Q35 = F(1, 2), F(2, 5), F(3, 5)
JINT = (Q25, -1, Q35, 1)                       # J(e) = [2/5 - e, 3/5 + e]
def _box(ivals):
    out = []
    for ch in itertools.product((0, 1), repeat=len(ivals)):
        b = []; d = []
        for c, (l0, l1, h0, h1) in zip(ch, ivals):
            b.append(l0 if c == 0 else h0); d.append(l1 if c == 0 else h1)
        out.append((b, d))
    return out
def _prism(p, V, J=None):
    out = [(list(p), list(v)) for v in V]
    if J is None: return out
    return [(b + [J[0]], d + [J[1]]) for b, d in out] + [(b + [J[2]], d + [J[3]]) for b, d in out]
CELLS = {   # vertices x0 + e*w of P_tau(e)
    'A': _box([(Q25, -1, Q25, F(1, 4)), JINT, JINT]),
    'B': _box([(F(1, 5), -H, F(1, 5), F(1, 3)), JINT, JINT]),
    'C': _prism((Q25, Q25), [(-1, -1), (F(2, 3), -1), (-1, 4)], JINT),
    'D': _prism((Q25, Q25), [(-1, -1), (F(3, 2), -1), (-1, F(3, 2))], JINT),
    'E': _prism((Q25, Q25, Q25), [(-1, -1, -1), (F(3, 2), -1, -1), (-1, 4, -1), (-1, -1, 4)]),
}
TYPEKEY = dict(zip('ABCDE', FIVE))

def contact_value(n, verts):
    """2/5 - e* from the contact interval n.P(e) = [alpha - e a, beta + e b], or None if e* > 1/15."""
    vals = [(sum(F(x) * y for x, y in zip(n, b)), sum(F(x) * y for x, y in zip(n, d))) for b, d in verts]
    alpha = min(v0 for v0, v1 in vals); beta = max(v0 for v0, v1 in vals)
    a = max(-v1 for v0, v1 in vals if v0 == alpha); b = max(v1 for v0, v1 in vals if v0 == beta)
    fl = alpha.numerator // alpha.denominator; ce = -((-beta.numerator) // beta.denominator)
    if fl + 1 <= beta or ce - 1 >= alpha or alpha.denominator == 1 or beta.denominator == 1:
        if any(alpha <= k <= beta for k in range(fl, ce + 1)): return Q25
    cand = []
    if b > 0: cand.append((ce - beta) / b)
    if a > 0: cand.append((alpha - fl) / a)
    e = min(cand)
    return Q25 - e if e <= F(1, 15) else None

def check11():
    forms = [(5, 1), (0, 1), (-4, -5), (-2, -10)]       # U = {4y1+21y2+5y3=0, u=2y1+3y2+3y3} in coordinates s
    # the parametrization: y1 = 5s1+s2, y2 = s2, y3 = -4s1-5s2, u = -2s1-10s2
    for (p, q) in forms[:3]: pass
    assert all(4 * f1 + 21 * f2 + 5 * f3 == 0 for f1, f2, f3 in [(5, 0, -4), (1, 1, -5)])
    assert (2 * 5 + 3 * 0 + 3 * (-4), 2 * 1 + 3 * 1 + 3 * (-5)) == (-2, -10)
    assert ml_torus(forms) == F(13, 30)
    xs = (F(3, 5), F(13, 30))
    vals = [(a * xs[0] + b * xs[1]) % 1 for a, b in forms]
    assert vals == [F(13, 30)] * 3 + [F(7, 15)]
    # slices y_a = y_b: identity components meet the diagonal only at r/30 with gcd(lambda_c,30) | r; 13 is not divisible by 2, 3, 5
    lam = (4, 21, 5)
    assert all(13 % gcd(l, 30) != 0 for l in lam)
    # direct check: for each pair, the primitive difference form does not vanish at x*
    for a, b in [(0, 1), (0, 2), (1, 2)]:
        f0, f1 = forms[a][0] - forms[b][0], forms[a][1] - forms[b][1]
        g = gcd(abs(f0), abs(f1)); assert ((f0 // g) * xs[0] + (f1 // g) * xs[1]) % 1 != 0
    assert ml_exact([1, 8, 29]) == F(13, 30)
    print("PASS 11: the slicing example has ML = 13/30, x* off all slice circles, ML(1,8,29) = 13/30")

def _hnf2(s, r):
    M = [list(s), list(r)]; piv = []; row = 0
    for col in range(4):
        if row == 2: break
        while True:
            nz = [i for i in range(row, 2) if M[i][col] != 0]
            if len(nz) <= 1: break
            i0 = min(nz, key=lambda i: abs(M[i][col]))
            for i in nz:
                if i != i0:
                    q = M[i][col] // M[i0][col]; M[i] = [x - q * y for x, y in zip(M[i], M[i0])]
        nz = [i for i in range(row, 2) if M[i][col] != 0]
        if not nz: continue
        if nz[0] != row: M[row], M[nz[0]] = M[nz[0]], M[row]
        if M[row][col] < 0: M[row] = [-x for x in M[row]]
        piv.append(col); row += 1
    q = piv[1]; t = M[0][q] // M[1][q]; M[0] = [x - t * y for x, y in zip(M[0], M[1])]
    return tuple(M[0] + M[1])

SIGNED_PERMS = [(pm, sg) for pm in itertools.permutations(range(4)) for sg in itertools.product((1, -1), repeat=4)]
def canon(s, r):
    return min(_hnf2([sg[i] * s[pm[i]] for i in range(4)], [sg[i] * r[pm[i]] for i in range(4)]) for pm, sg in SIGNED_PERMS)

def in_R(v):
    return F(1, 3) <= v < F(2, 5) and 1 <= 2 * v.denominator - 5 * v.numerator <= 4

def check12(path):
    """Survivors of the enumeration below 2/5: 's0 r0 s1 r1 s2 r2 s3 r3 p q X Y Den' per line."""
    n = 0; vals = {}; orbit_val = {}
    for line in open(path):
        t = list(map(int, line.split()))
        S, R = t[0:8:2], t[1:8:2]; p, q, X, Y, Den = t[8:]
        val = F(p, q)
        assert min(dist(F(S[i] * X + R[i] * Y, Den)) for i in range(4)) == val, t
        assert in_R(val), t
        assert sum((S[i] * R[j] - S[j] * R[i]) ** 2 for i, j in itertools.combinations(range(4), 2)) <= 11331
        c = canon(S, R)
        assert orbit_val.setdefault(c, val) == val
        vals[val] = vals.get(val, 0) + 1; n += 1
    for c, val in orbit_val.items():
        assert ml_torus(list(zip(c[:4], c[4:]))) == val, c
    print(f"PASS 12: {n} survivors, witnesses verified, {len(orbit_val)} orbits with exact values recomputed; values " + ", ".join(f"{v}: {k}" for v, k in sorted(vals.items())))

def check13(N=3):
    tot = 0; low = 0
    for tau, verts in CELLS.items():
        phi = FIVE[TYPEKEY[tau]]
        for n in itertools.product(range(-N, N + 1), repeat=3):
            if n == (0, 0, 0) or gcd(gcd(abs(n[0]), abs(n[1])), abs(n[2])) != 1 or n < tuple(-x for x in n): continue
            c1, c2 = kernel3(n)
            s = [sum(c1[k] * ph[k] for k in range(3)) for ph in phi]
            r = [sum(c2[k] * ph[k] for k in range(3)) for ph in phi]
            if any(s[i] == 0 and r[i] == 0 for i in range(4)): continue
            ex = ml_torus(list(zip(s, r))); pred = contact_value(n, verts)
            if pred is None: assert ex < F(1, 3), (tau, n, ex)
            else: assert ex == pred, (tau, n, ex, pred)
            if F(1, 3) <= ex < Q25: assert in_R(ex), (tau, n, ex); low += 1
            tot += 1
    print(f"PASS 13: contact formula of Lemma 9.13 equals the exact value on {tot} subtori of the five tori (|n|_inf <= {N}); {low} values in [1/3,2/5), all in R")

def check14(N=60):
    """Lemma 9.12 on the grid (1/N)Z^3: G_tau(y) >= 2/5 - e iff y lies in P_tau(e) or -P_tau(e) mod 1."""
    def inP(tau, y, e):
        if tau in 'AB':
            lo, hi = ((Q25 - e, Q25 + e / 4) if tau == 'A' else (F(1, 5) - e / 2, F(1, 5) + e / 3))
            return lo <= y[0] <= hi and all(Q25 - e <= t <= Q35 + e for t in y[1:])
        d = [(t - Q25) / e for t in (y[:2] if tau in 'CD' else y)]
        if tau in 'CD' and not (Q25 - e <= y[2] <= Q35 + e): return False
        if any(x < -1 for x in d): return False
        if tau == 'C': return 3 * d[0] + d[1] <= 1
        if tau == 'D': return d[0] + d[1] <= H
        return 2 * d[0] + d[1] + d[2] <= 1
    cnt = 0
    for tau in 'ABCDE':
        phi = FIVE[TYPEKEY[tau]]
        for e in (F(1, 15), F(1, 30), F(1, 40)):
            d = Q25 - e
            for y in itertools.product(range(N), repeat=3):
                yv = [F(t, N) for t in y]
                g = min(dist(sum(ph[k] * yv[k] for k in range(3))) for ph in phi)
                mine = [(1 - t) % 1 for t in yv]
                member = inP(tau, yv, e) or inP(tau, mine, e)
                assert (g >= d) == member, (tau, y, e)
                cnt += 1
    print(f"PASS 14: Lemma 9.12 confirmed at {cnt} grid points (N = {N}, three values of epsilon, five tori)")

def check15(J=40):
    for j in range(2, J + 1):
        assert ml_exact([2, 3, 5 * j]) == F(2 * j, 5 * j + 2), j
        assert ml_exact([1, 4, 5 * j + 3]) == F(2 * j + 1, 5 * j + 4), j
    R1 = {F(p, q) for q in range(1, 401) for p in range(1, q) if gcd(p, q) == 1 and in_R(F(p, q))}
    R2 = {F(2 * j, 5 * j + 2) for j in range(2, 400)} | {F(2 * j + 1, 5 * j + 4) for j in range(2, 400)}
    assert R1 == {x for x in R2 if x.denominator <= 400}
    print(f"PASS 15: ML(2,3,5j) = 2j/(5j+2) and ML(1,4,5j+3) = (2j+1)/(5j+4) for 2 <= j <= {J}; the two descriptions of R agree for q <= 400 ({len(R1)} fractions)")

def check16():
    forms = [(10, 1), (-49, -5), (9, 0), (0, 1)]
    m = (4, 1, 1, 1)
    assert all(sum(m[i] * forms[i][k] for i in range(4)) == 0 for k in range(2))
    assert ml_torus(forms) == F(43, 103) and ml_exact([1, 102, 25]) == F(43, 103)
    assert 3 * 103 - 7 * 43 == 8 and sum(abs(x) for x in m) == 7
    print("PASS 16: the subtorus {(10X+Y, -49X-5Y, 9X, Y)} of W_(4,1,1,1) has ML = 43/103 = ML(1,102,25)")

def _rep_circles(A, B, C):
    """Speed vectors of the identity components of U(A,B,C) cap {x_i = +-x_j} that are coordinate-free."""
    N = A + B + C; H = N // 2
    c1, c2 = kernel3((A, B, C))
    s = [c1[0], c1[1], c1[2], H * c1[1]]; r = [c2[0], c2[1], c2[2], H * c2[1]]
    out = []
    for i, j in itertools.combinations(range(4), 2):
        for sg in (1, -1):
            a, b = s[i] - sg * s[j], r[i] - sg * r[j]
            if a == 0 and b == 0: continue
            g = gcd(abs(a), abs(b)); v = [s[k] * (b // g) - r[k] * (a // g) for k in range(4)]
            if all(v): out.append([abs(x) for x in v])
    return s, r, out

def check17():
    for N in range(8, 201, 4):
        assert ml_exact([1, N // 2, N - 1]) == F(1, 2) - F(1, N), N
    n = 0
    for (A, B, C) in [(1, 9, 2), (2, 9, 1), (1, 13, 2), (1, 11, 4), (3, 11, 2)]:
        N = A + B + C; L = F(1, 2) - F(1, N)
        s, r, circ = _rep_circles(A, B, C)
        assert 4 * (N // 4) == N and N >= 2 * (A + C) + 2 and ml_torus(list(zip(s, r))) == L
        best = max(ml_exact(v) for v in circ)
        if all(gcd(x, N) > 1 for x in (A, B, C)): assert best <= L - F(1, N ** 3)
        else: assert best == L
        n += 1
    for N in range(12, 61, 4):
        for l in range(1, 4):
            V = [1, N // 2, 2 * N * l - 3, 3 * N * l - N // 2 - 2]
            assert ml_exact(V) == F((N - 2) * l - 1, 2 * N * l - 2), (N, l)
    N, l = 60, 61
    L = F(1, 2) - F(1, N)
    s, r, circ = _rep_circles(2, N - 5, 3)
    best = max(ml_exact(v) for v in circ)
    mv = ml_exact([1, N // 2, 2 * N * l - 3, 3 * N * l - N // 2 - 2])
    assert best <= L - F(1, N ** 3) < mv == F((N - 2) * l - 1, 2 * N * l - 2)
    # Proposition 10.3 (coalescence) on small cases, and Theorem 10.4 for N = 12 on both sides of the transition
    def forms_H(A, B, C, H):
        c1, c2 = kernel3((A, B, C))
        return list(zip([c1[0], c1[1], c1[2], H * c1[1]], [c2[0], c2[1], c2[2], H * c2[1]]))
    for (A, B, C, H) in [(1, 3, 1, 2), (1, 5, 2, 2), (2, 5, 1, 3), (1, 7, 3, 2)]:
        S = A + C
        assert gcd(B, S) == 1 and ml_torus(forms_H(A, B, C, H)) == ml_exact([S, B, H * S]), (A, B, C, H)
    assert ml_torus(forms_H(1, 7, 3, 2)) == ml_torus(forms_H(2, 7, 2, 2)) == ml_torus(forms_H(3, 7, 1, 2))
    # Theorem 10.4, N = 12: (1,10,1) with H = 6, 18 is in the first phase; (1,6,5) with H = 6 is in the second, g = 6
    for (A, B, C, H) in [(1, 10, 1, 6), (1, 10, 1, 18), (1, 6, 5, 6)]:
        N = A + B + C; S = A + C; Ds = F(B, 2 * (H * S + B))
        assert ml_torus(forms_H(A, B, C, H)) == F(1, 2) - min(F(1, N), Ds), (A, B, C, H)
    assert ml_exact([1, 1, 6]) == F(1, 2) - F(1, 14)
    assert ml_exact([11, 1, 330]) == F(15, 31) and ml_exact([11, 1, 390]) == F(195, 401)
    print(f"PASS 17: ML(1,N/2,N-1) for 4|N<=200, {n} subtori U(A,B,C) with their repeated-speed circles, Theorem 10.2 for N<=60, l<=3, the comparison in U(2,55,3) with l = 61, Proposition 10.3 and Theorem 10.4 on small cases")

if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    check1(); check2(); check3()
    check4(os.path.join(here, "tori_survivors.json"))
    check5(); check6(); check8(); check9(); check10(); check11()
    p12 = os.path.join(here, "surv25.txt")
    if os.path.exists(p12): check12(p12)
    else: print("SKIP 12: surv25.txt not present")
    check13(); check14(); check15(); check16(); check17()
    p = os.path.join(here, "b512.txt")
    if os.path.exists(p): check7(p)
    else: print("SKIP 7: b512.txt not present")
    print("ALL CHECKS PASSED")

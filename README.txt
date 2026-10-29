Four runners that are barely lonely
Source package, 28 October 2026 (v7: Section 11.1 on five speeds added, with g5enum2.c, torusmlN.c and two check scripts; v6: Section 11 added, an effective reduction of near-tight n-tuples to critical subtori for every n, proved without computation; v5: Section 10 on slicing obstructions of unbounded order added; Section 9.1 rewritten, two-dimensional subtori and three-runner values below 2/5,
shortest-character density bound, exact analysis of the five tori with |m|_1 = 5, the 9/25 enumeration replaced by a
smaller one; v4: two-dimensional subtori up to 9/25 and the value 7/50; v3: Section 9 on two-dimensional subtori below 1/3 added; v2: citations to Chen, Barajas-Serra and Cordella added; all computations rerun from a clean build, byte-identical)

FILES
  four_runners.tex / four_runners.pdf   the paper
  arxiv_abstract.txt                    the abstract for the arXiv form (1,808 characters)
  verify_four_runners.py                standard-library Python re-verification, checks 1 to 17 (about 16 min)
  tori_survivors.json                   the 123 surviving 2-dimensional subtori with certificate points
  b512.txt                              output of the exhaustive search: the 14032 primitive quadruples
                                        with |v| < 512 and ML < 1/4, each with its exact ML
  src/torusml.c      exact ML of a 2-dimensional subtorus (reads lines "s1 r1 s2 r2 s3 r3 s4 r4")
  src/enum2s.c       Hermite-normal-form enumeration of rank-2 sublattices of Z^4, covolume <= A,
                     with grid certificates; usage: ./enum2s 33 5 16 > surv33.txt
  src/enumWn.c       enumeration inside the special 3-tori; usage: ./enumWn 1 48 > survW1.txt ; ./enumWn 2 48 > survW2.txt
  src/pluckcount.c   independent count of primitive decomposable Pluecker vectors (completeness check)
  src/enum3s.c       enumeration of covolume <= 57 with the special-3-torus skip and threshold 1/3 (Section 9);
                     usage: ./enum3s 57 1 3 > surv57.txt   (expected: 63906892 lattices, 3267840 special, 0 survivors)
  src/ml4b.c         exhaustive search; usage: ./ml4b 512 [v1_lo v1_hi] > b512.txt
  src/ml3.c          by-product on triples (Section 9); usage: ./ml3 200   (expected: 1098601 primitive triples,
                     4639 with ML < 1/3, 0 without a relation c = 2b or c = a + b)
  src/check_lemmas_c.py   Lemmas 9.2 and 9.4 and the formula for F(n1,n2) against torusml on the larger ranges
                     stated in Section 12; usage: python3 src/check_lemmas_c.py ./torusml   (about 1 s)
  src/torusml2.c     the same exact maximizer as torusml.c, with the offsets of a base pair of forms reduced to a
                     fundamental domain (|det| classes); agrees with torusml.c on all 916,492 subtori of covolume <= 20,
                     on 19,300 others and on all 18,960 survivors of enum25
  src/enum25.c       enumeration of Proposition 9.14: squared covolume <= 11331, skip of the special 3-tori and of the
                     3-tori with |m|_1 = 5, grids Q in {2,10,20,30}, then the exact maximizer, threshold 2/5 (non-strict);
                     usage: ./enum25 11331 2 5 2 0 > v_0.txt ; ./enum25 11331 2 5 2 1 > v_1.txt ; cat v_0.txt v_1.txt | sort > surv25.txt
                     expected totals: 789460684 lattices, 21301056 special, 36832896 in |m|_1 = 5 tori, 44820720 discarded
                     by the exact maximizer, 18960 survivors (about 21 minutes per job)
  src/check_five_c.py   the contact formula in the proof of Lemma 9.13 against torusml2 for all U_n with |n|_inf <= 20 in
                     the five tori; usage: python3 src/check_five_c.py ./torusml2 20   (expected: 142547 subtori, 0 mismatches)
  src/allmax.c, src/slicing_family.py   the slicing statistics of Remark 9.16
  src/mlc.c          exact ML of a set of positive integers by the corner lemma (reads lines "k v1 ... vk")
  src/check_obstr_c.py   Section 10 against torusml2 and mlc; usage: python3 src/check_obstr_c.py ./torusml2 ./mlc
  src/torusmlN.c     the exact maximizer of torusml2.c for a two-dimensional subtorus of (R/Z)^n with n <= 8 forms;
                     input lines "n s1 r1 ... sn rn", output "p q X Y Den" as torusml2
  src/g5enum2.c      enumeration of Theorem 11.5 (five speeds): for a seed q, all normalized c with
                     covol(Zq+Zc)^2 <= 765275, keeping W = Zq+Zc only if every coordinate-free circle of W with a
                     repeated absolute speed is a tight quadruple circle; usage: ./g5enum2 1 1 2 3 4 765275 (and the
                     other seven seeds listed in g5_survivors.txt); expected: 16 survivors in total, see Section 12
  src/check_g5.py    exact ML (torusmlN) of all W = Zq+Zc with |c_i| <= 3 and of 20,000 random primitive planes
  src/check_g5_enum.py   pure-Python box search reproducing the candidate counts and survivors of g5enum2 at a
                     smaller bound; usage: python3 src/check_g5_enum.py 3000
  g5_survivors.txt   the 16 survivors of g5enum2 (seed q, vector c, squared covolume)
  surv25.txt         the 18,960 survivors of enum25 with exact ML and witness point (checked by verify_four_runners.py, check 12)

BUILD
  gcc -O2 -o torusml src/torusml.c ; gcc -O2 -o enum2s src/enum2s.c -lm ; gcc -O2 -o enumWn src/enumWn.c
  gcc -O2 -o pluckcount src/pluckcount.c ; gcc -O2 -o ml4b src/ml4b.c ; gcc -O2 -o enum3s src/enum3s.c -lm ; gcc -O2 -o ml3 src/ml3.c
  gcc -O2 -o torusml2 src/torusml2.c ; gcc -O2 -o enum25 src/enum25.c -lm ; gcc -O2 -o allmax src/allmax.c ; gcc -O2 -o mlc src/mlc.c
  gcc -O2 -o torusmlN src/torusmlN.c ; gcc -O2 -o g5enum2 src/g5enum2.c
  pdflatex four_runners.tex (twice)

RUN
  python3 verify_four_runners.py

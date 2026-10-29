# Four runners that are barely lonely

**Hyunjun Yi, Karam Gill, Bryan Bu**

October 2026. Preprint.

DOI: [10.5281/zenodo.23763104](https://doi.org/10.5281/zenodo.23763104)

## Abstract

In the lonely runner problem, ML(v) is the largest loneliness any time t achieves for speeds v, and Kravitz asked which values below 1/n occur. We settle the four-runner case: a primitive quadruple has ML < 1/4 exactly when it is {1, 2, 3, 4k}, or {x, y, x+y, x+2y} with gcd(x, y) = 1 and y ≢ x, 2x (mod 4), or {1, 3, 4, 14}. The values below 1/4 are therefore every m/(4m+2) with m ≥ 2 except 3/14. A density lemma bounds the length of a quadruple on neither critical subtorus by 512, and an exact search over 800,816,099 quadruples completes the proof. The same method shows that the accumulation points of the four-runner spectrum in (0, 2/5] are the three-runner values there, as Giri and Kravitz conjectured.

## Contents

- `four_runners.pdf`: the paper
- `four_runners.tex`: LaTeX source
- `verify_four_runners.py`: top-level verification script
- `src/`: C and Python enumeration and checking code
- `b512.txt, surv25.txt, g5_survivors.txt, tori_survivors.json`: search outputs used in the proof
- `README.txt`: how to rerun the computations

## Citation

```bibtex
@misc{fourrunners2026,
  title     = {Four runners that are barely lonely},
  author    = {Hyunjun Yi and Karam Gill and Bryan Bu},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23763104}
}
```

## License

The paper is released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

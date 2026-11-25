# The radius-two supported norm of the modular group acting on prime projective lines

`main.tex` and `main.pdf`: the paper. For every prime p it determines the largest norm of
A_p = pi(S) + pi(R) + pi(R^-1) on vectors supported in one radius-two ball of the Schreier graph of
PSL_2(Z) acting on P^1(F_p): 3 for p = 2, 3, 5, 7; sqrt(beta) = 2.90782... for p = 13; and
sqrt(alpha) = 2.890889... for every other prime, attained only at the centre 1.

## Checking the computations

Python 3.8 or later. From `code/`:

    python certificate.py      # standard library only, exact arithmetic; writes certificate_results.json
    python float_check.py 200  # NumPy; floating-point values of N_2(p) for p <= 200 and p = 257, 1009

`certificate.py` ends with `ALL CHECKS PASS` and takes under a minute. It checks every finite statement
used in the proof: the ball sizes, the obstruction values and the exceptional set E, the 22 commutation
classes and the 23 local models, the Gram matrix and characteristic polynomial of the maximal model, the
Sylvester certificates for the other models, the brackets for alpha, sqrt(alpha) and beta, and every centre
at every prime of E. It also compares the models with direct computation at p = 17, 19, 23, 101, 257, 1009.

`SHA256SUMS` lists the files of this repository.

## Licence

Text and code: CC BY 4.0 (see `LICENSE`).

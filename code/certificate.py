# Exact certificate for the radius-two supported norm of A_p = pi(S)+pi(R)+pi(R^-1) on P^1(F_p).
# Standard library only; every sign decision is made in integers or Fractions.
import json
import sys
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd

S = (0, -1, 1, 0)
R = (0, -1, 1, 1)
RI = (1, 1, -1, 0)
GENS = [S, R, RI]
I2 = (1, 0, 0, 1)
Q = [1, -19, 124, -331, 325, -60]
Q13 = [1, -18, 107, -242, 169, -24]
PRINTED_GRAM = [
    [3, 1, 0, 0, 1, 1, 1, 0], [1, 3, 2, 2, 0, 0, 1, 0], [0, 2, 5, 1, 1, 0, 0, 0], [0, 2, 1, 5, 1, 0, 0, 0],
    [1, 0, 1, 1, 3, 1, 1, 1], [1, 0, 0, 0, 1, 3, 0, 0], [1, 1, 0, 0, 1, 0, 3, 1], [0, 0, 0, 0, 1, 0, 1, 3]]
L_SQRT = F(2890889946779, 10 ** 12)
U_SQRT = F(2890889946780, 10 ** 12)


def norm(m):
    a, b, c, d = m
    first = next(v for v in m if v != 0)
    return m if first > 0 else (-a, -b, -c, -d)


def mul(m, n):
    a, b, c, d = m
    e, f, g, h = n
    return norm((a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h))


def inv(m):
    a, b, c, d = m
    return norm((d, -b, -c, a))


def balls(L):
    out = [{I2}]
    for _ in range(L):
        prev = out[-1]
        out.append(prev | {mul(s, g) for s in GENS for g in prev})
    return out


def delta(g, h):
    gh, hg = mul(g, h), mul(h, g)
    # use raw (unnormalised) products for the additive commutator
    a, b, c, d = g
    e, f, gg, hh = h
    p = (a * e + b * gg, a * f + b * hh, c * e + d * gg, c * f + d * hh)
    q = (e * a + f * c, e * b + f * d, gg * a + hh * c, gg * b + hh * d)
    m = [p[i] - q[i] for i in range(4)]
    return m[0] * m[3] - m[1] * m[2]


def commute(g, h):
    return mul(g, h) == mul(h, g)


def prime_factors(n):
    n, out, k = abs(n), set(), 2
    while k * k <= n:
        while n % k == 0:
            out.add(k)
            n //= k
        k += 1
    if n > 1:
        out.add(n)
    return out


def primes(n):
    return [k for k in range(2, n + 1) if all(k % j for j in range(2, int(k ** 0.5) + 1))]


def act(m, x, p):
    a, b, c, d = m
    if x == p:  # infinity
        return p if c % p == 0 else (a * pow(c, -1, p)) % p
    den = (c * x + d) % p
    return p if den == 0 else ((a * x + b) * pow(den, -1, p)) % p


def det(M):
    M = [[F(v) for v in row] for row in M]
    n, s = len(M), F(1)
    for i in range(n):
        piv = next((r for r in range(i, n) if M[r][i] != 0), None)
        if piv is None:
            return F(0)
        if piv != i:
            M[i], M[piv] = M[piv], M[i]
            s = -s
        s *= M[i][i]
        for r in range(i + 1, n):
            f = M[r][i] / M[i][i]
            if f:
                M[r] = [M[r][k] - f * M[i][k] for k in range(n)]
    return s


def charpoly(G):
    # Faddeev-LeVerrier over Fractions; returns monic integer coefficients, highest degree first
    n = len(G)
    A = [[F(v) for v in row] for row in G]
    coeffs = [F(1)]
    Mk = [[F(0)] * n for _ in range(n)]
    for k in range(1, n + 1):
        Mk = [[sum(A[i][t] * Mk[t][j] for t in range(n)) + (coeffs[-1] if i == j else 0) for j in range(n)]
              for i in range(n)]
        AM = [[sum(A[i][t] * Mk[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        coeffs.append(-sum(AM[i][i] for i in range(n)) / k)
    assert all(c.denominator == 1 for c in coeffs)
    return [int(c) for c in coeffs]


def polymul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def shifted_pd(G, t):
    # all leading principal minors of tI - G positive (t a Fraction) => every eigenvalue of G is < t
    n = len(G)
    T = [[(t if i == j else 0) - G[i][j] for j in range(n)] for i in range(n)]
    return all(det([row[:k] for row in T[:k]]) > 0 for k in range(1, n + 1))


def shifted_det(G, t):
    n = len(G)
    return det([[(t if i == j else 0) - G[i][j] for j in range(n)] for i in range(n)])


def gram(M):
    cols = len(M[0])
    return [[sum(row[i] * row[j] for row in M) for j in range(cols)] for i in range(cols)]


def model(B2, B3, C):
    # equivalence on B3: g ~ h iff h^-1 g is 1 or in C
    parent = {g: g for g in B3}

    def find(g):
        while parent[g] != g:
            parent[g] = parent[parent[g]]
            g = parent[g]
        return g
    for g, h in combinations(sorted(B3), 2):
        if mul(inv(h), g) in C:
            parent[find(g)] = find(h)
    rep = {}
    for g in sorted(B2):
        rep.setdefault(find(g), g)
    ins = sorted(rep.values())
    outs = sorted({find(mul(s, g)) for g in ins for s in GENS})
    oi = {y: i for i, y in enumerate(outs)}
    M = [[0] * len(ins) for _ in outs]
    for j, g in enumerate(ins):
        for s in GENS:
            M[oi[find(mul(s, g))]][j] += 1
    return M


def field_gram(p, x):
    nb = lambda u: [act(m, u, p) for m in GENS]
    ball, frontier = {x}, {x}
    for _ in range(2):
        frontier = {y for u in frontier for y in nb(u)} - ball
        ball |= frontier
    ball = sorted(ball)
    outs = sorted({y for u in ball for y in nb(u)})
    oi = {y: i for i, y in enumerate(outs)}
    M = [[0] * len(ball) for _ in outs]
    for j, u in enumerate(ball):
        for y in nb(u):
            M[oi[y]][j] += 1
    return gram(M), len(ball)


def main():
    res = {}
    B = balls(6)
    sizes = [len(b) for b in B]
    res["ball_sizes"] = sizes
    assert sizes == [1, 4, 8, 14, 22, 34, 50], sizes
    B6 = sorted(B[6])
    nonid = [g for g in B6 if g != I2]

    # Step 1: obstructions over distinct unordered pairs of B6
    pairs = list(combinations(B6, 2))
    ds = [delta(g, h) for g, h in pairs]
    sing = [(g, h) for (g, h), d in zip(pairs, ds) if d == 0]
    res["distinct_pairs"], res["singular_pairs"], res["max_abs_delta"] = len(pairs), len(sing), max(abs(d) for d in ds)
    res["pairs_with_repetition"] = len(pairs) + len(B6)
    res["singular_with_repetition"] = len(sing) + len(B6)
    res["distinct_pairs_without_identity"] = len(pairs) - len(nonid)
    res["singular_pairs_without_identity"] = sum(1 for g, h in sing if I2 not in (g, h))
    assert (res["distinct_pairs_without_identity"], res["singular_pairs_without_identity"]) == (1176, 50)
    assert all(commute(g, h) for g, h in sing), "singular non-commuting pair"
    cand = sorted(set().union(*[prime_factors(d) for d in ds if d != 0]))
    res["candidate_primes"] = cand
    assert max(cand) <= 100

    # E from actual common fixed points, every prime <= 101
    E = []
    for p in primes(101):
        bad = False
        for x in range(p + 1):
            Fx = [g for g in nonid if act(g, x, p) == x]
            if any(not commute(g, h) for g, h in combinations(Fx, 2)):
                bad = True
                break
        if bad:
            E.append(p)
    res["E"] = E
    assert E == [2, 3, 5, 7, 11, 13, 37, 41, 61, 73, 97], E
    assert set(E) <= set(cand) | {2, 3}

    # scalar reductions
    scal = sorted(set().union(*[prime_factors(gcd(gcd(b, c), a - d)) for a, b, c, d in nonid]))
    res["scalar_reduction_primes"] = scal
    assert set(scal) <= {2, 3}

    # Step 2: commutation classes by primitive (b, c, a-d)
    def key(m):
        a, b, c, d = m
        v = (b, c, a - d)
        g = gcd(gcd(abs(v[0]), abs(v[1])), abs(v[2]))
        v = tuple(t // g for t in v)
        return v if next(t for t in v if t) > 0 else tuple(-t for t in v)
    classes = {}
    for g in nonid:
        classes.setdefault(key(g), set()).add(g)
    for cl in classes.values():
        assert all(commute(g, h) for g, h in combinations(cl, 2))
    for k1, k2 in combinations(classes, 2):
        assert not any(commute(g, h) for g in classes[k1] for h in classes[k2])
    res["n_classes"] = len(classes)
    assert len(classes) == 22

    B2, B3 = sorted(B[2]), sorted(B[3])
    models = {"empty": model(B2, B3, set())}
    for k, cl in classes.items():
        models[str(k)] = model(B2, B3, cl)
    grams = {k: gram(M) for k, M in models.items()}
    polys = {k: charpoly(G) for k, G in grams.items()}
    res["model_dims"] = sorted(len(G) for G in grams.values())
    assert max(res["model_dims"]) <= 8

    # Step 3: the maximal model and the comparison
    cstar = {norm((0, 1, -1, 2)), norm((2, -1, 1, 0))}
    kstar = key(norm((0, 1, -1, 2)))
    assert classes[kstar] == cstar, classes[kstar]
    Gs = grams[str(kstar)]
    target = polymul(polymul([1, -4], [1, -5, 5]), Q)
    res["max_charpoly"] = polys[str(kstar)]
    assert polys[str(kstar)] == target == [1, -28, 320, -1942, 6784, -13740, 15285, -8000, 1200]
    n = len(PRINTED_GRAM)
    assert len(Gs) == n and any(all(Gs[pi[i]][pi[j]] == PRINTED_GRAM[i][j] for i in range(n) for j in range(n))
                                for pi in permutations(range(n))), "printed Gram not matched"
    smallest = None
    for k, G in grams.items():
        if k == str(kstar):
            continue
        T = [[(83 if i == j else 0) - 10 * G[i][j] for j in range(len(G))] for i in range(len(G))]
        minors = [det([r[:m] for r in T[:m]]) for m in range(1, len(G) + 1)]
        assert all(v > 0 for v in minors), k
        smallest = min(minors) if smallest is None else min(smallest, min(minors))
    res["smallest_competitor_minor_83I_10G"] = str(smallest)
    # alpha bracket: every eigenvalue < U^2, and det(L^2 I - G) < 0 gives an eigenvalue > L^2
    assert shifted_pd(Gs, U_SQRT ** 2) and shifted_det(Gs, L_SQRT ** 2) < 0
    assert shifted_pd(Gs, F(83572446844, 10 ** 10)) and shifted_det(Gs, F(83572446843, 10 ** 10)) < 0
    res["alpha_bracket"] = "8.3572446843 < alpha < 8.3572446844; 2.890889946779 < sqrt(alpha) < 2.890889946780"
    # x = 1 is fixed by (0 1;-1 2), infinity is not
    assert norm((0, 1, -1, 2)) in cstar and (1 - 1) ** 2 == 0

    # Step 4: primes in E, exhaustively over every centre
    model_polys = set(tuple(v) for v in polys.values())
    per_prime = {}
    for p in E:
        data = []
        for x in range(p + 1):
            G, size = field_gram(p, x)
            data.append((x, G, size))
        if p in (2, 3, 5, 7):
            full = [x for x, G, size in data if size == p + 1]
            assert full, p
            per_prime[p] = {"full_ball_centres": full}
        elif p == 13:
            pmax = polymul(polymul([1, -3], [1, -2]), Q13)
            tops = [x for x, G, s in data if charpoly(G) == pmax]
            assert tops == [2, 7], tops
            for x, G, s in data:
                if x not in tops:
                    assert shifted_pd(G, F(42, 5)), x
            G = data[2][1]
            assert shifted_pd(G, F(84554193489, 10 ** 10)) and shifted_det(G, F(84554193488, 10 ** 10)) < 0
            per_prime[p] = {"centres": tops}
        else:
            tops = [x for x, G, s in data if charpoly(G) == target]
            assert tops == [1], (p, tops)
            for x, G, s in data:
                if x != 1:
                    assert shifted_pd(G, F(83572446843, 10 ** 10)), (p, x)
            per_prime[p] = {"centres": tops}
    res["E_primes"] = per_prime

    # link test: at good primes every centre's Gram polynomial is one of the 23 model polynomials,
    # and the maximal one occurs only at x = 1
    link = {}
    for p in [17, 19, 23, 101, 257, 1009]:
        seen = set()
        for x in range(p + 1):
            G, _ = field_gram(p, x)
            cp = tuple(charpoly(G))
            assert cp in model_polys, (p, x)
            seen.add(cp)
            if cp == tuple(target):
                assert x == 1, (p, x)
        link[p] = len(seen)
    res["link_test_distinct_models_seen"] = link
    with open("certificate_results.json", "w", newline=chr(10)) as fh:
        json.dump(res, fh, indent=1, default=str)
        fh.write(chr(10))
    print("\nALL CHECKS PASS")


if __name__ == "__main__":
    main()

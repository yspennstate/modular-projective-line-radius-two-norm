# Floating-point check of the radius-two supported norm (NumPy eigenvalues).
# A_p = pi(S)+pi(R)+pi(R^-1) on P^1(F_p); N2(p) = max_x ||A_p P_{2,x}|| (input supported on the closed
# radius-two ball about x in the Schreier graph, output unrestricted). Claim: N2(p)^2 = alpha (top root of Q)
# for p not in {2,3,5,7,13}, unique center x=1; N2=3 for p in {2,3,5,7}; N2^2 = beta at p=13.
import sys
import numpy as np

INF = "inf"
MATS = {"S": (0, -1, 1, 0), "R": (0, -1, 1, 1), "Ri": (1, 1, -1, 0)}


def act(m, x, p):
    a, b, c, d = m
    if x == INF:
        return INF if c % p == 0 else (a * pow(c, -1, p)) % p
    num, den = (a * x + b) % p, (c * x + d) % p
    return INF if den == 0 else (num * pow(den, -1, p)) % p


def n2(p):
    pts = list(range(p)) + [INF]
    nb = {x: [act(m, x, p) for m in MATS.values()] for x in pts}
    best, arg = -1.0, []
    for x in pts:
        ball, frontier = {x}, {x}
        for _ in range(2):
            frontier = {y for u in frontier for y in nb[u]} - ball
            ball |= frontier
        ball = sorted(ball, key=lambda v: (v == INF, v if v != INF else 0))
        outs = sorted({y for u in ball for y in nb[u]}, key=lambda v: (v == INF, v if v != INF else 0))
        oi = {y: i for i, y in enumerate(outs)}
        M = np.zeros((len(outs), len(ball)))
        for j, u in enumerate(ball):
            for y in nb[u]:
                M[oi[y], j] += 1
        lam = np.linalg.eigvalsh(M.T @ M)[-1]
        if lam > best + 1e-9:
            best, arg = lam, [x]
        elif abs(lam - best) <= 1e-9:
            arg.append(x)
    return best, arg


def top_root(coeffs):
    return max(r.real for r in np.roots(coeffs) if abs(r.imag) < 1e-9)


def primes(n):
    s = [True] * (n + 1)
    s[0] = s[1] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = [False] * len(s[i * i::i])
    return [i for i in range(n + 1) if s[i]]


if __name__ == "__main__":
    hi = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    alpha = top_root([1, -19, 124, -331, 325, -60])
    beta = top_root([1, -18, 107, -242, 169, -24])
    print(f"alpha={alpha:.10f} sqrt={alpha**0.5:.12f} beta={beta:.10f}")
    bad = 0
    for p in primes(hi) + [257, 1009]:
        lam, arg = n2(p)
        if p in (2, 3, 5, 7):
            exp, ok_c = 9.0, True
        elif p == 13:
            exp, ok_c = beta, True
        else:
            exp, ok_c = alpha, arg == [1]
        ok = abs(lam - exp) < 1e-7 and ok_c
        bad += not ok
        if not ok or p in (2, 3, 5, 7, 11, 13, 97, 101, 257, 1009):
            print(f"p={p} N2^2={lam:.10f} expected={exp:.10f} centers={arg[:6]} {'OK' if ok else 'MISMATCH'}")
    print("mismatches:", bad)

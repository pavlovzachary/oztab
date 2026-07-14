"""Schur -> s~ coefficients, via Jacobi-Trudi on top of the h-side M-table.

The research object is

    s_lambda = sum_mu N_{lambda,mu} s~_mu,

but the tableau model counts the h-side (see coefficients.M_row).  Jacobi-Trudi

    s_lambda = det( h_{lambda_i - i + j} )_{1 <= i,j <= l}

writes s_lambda as a signed sum of products of h's, i.e. of h_nu's, so

    N_row(lambda) = sum over JT terms  sign * M_row(nu).

For two rows this is the familiar s_{(a,b)} = h_a h_b - h_{a+1} h_{b-1} that the
b = 2 and b = 3 analyses were built on; the determinant here is the general case.

Unlike M, N can be negative: nonnegativity of N is a hypothesis, not a theorem,
and `negative_terms` / scripts/sweep_nonneg.py exist to hunt for violations.
"""

from itertools import permutations

from .coefficients import M_row
from .multisets import integer_partitions


def _h_index(lam, i, j):
    """The Jacobi-Trudi matrix entry index lambda_i - i + j (0-based i, j)."""
    return lam[i] - (i + 1) + (j + 1)


def jacobi_trudi_terms(lam):
    """Expand det(h_{lambda_i - i + j}) as a list of (sign, nu) pairs, where nu
    is the partition of the h-product h_nu = prod_k h_{nu_k}.

    h_0 = 1 (index dropped from nu) and h_k = 0 for k < 0 (term dropped).
    """
    lam = tuple(p for p in lam if p > 0)
    size = len(lam)
    if size == 0:
        return [(1, ())]
    terms = []
    for perm in permutations(range(size)):
        sign = _permutation_sign(perm)
        indices = []
        for i in range(size):
            k = _h_index(lam, i, perm[i])
            if k < 0:
                break                      # h_negative = 0, term vanishes
            if k > 0:
                indices.append(k)          # h_0 = 1, contributes nothing
        else:
            terms.append((sign, tuple(sorted(indices, reverse=True))))
    return terms


def _permutation_sign(perm):
    """Sign of a permutation given as a tuple of images, by cycle count."""
    seen = [False] * len(perm)
    sign = 1
    for start in range(len(perm)):
        if seen[start]:
            continue
        length = 0
        node = start
        while not seen[node]:
            seen[node] = True
            node = perm[node]
            length += 1
        if length % 2 == 0:                # an even-length cycle is odd
            sign = -sign
    return sign


def N_row(lam):
    """Row of the N-matrix for s_lambda: {mu: coefficient of s~_mu}, nonzero mu.

    Coefficients may be negative.  Combines the JT terms' M-rows with signs.
    """
    row = {}
    for sign, nu in jacobi_trudi_terms(lam):
        for mu, count in M_row(nu).items():
            row[mu] = row.get(mu, 0) + sign * count
    return {mu: c for mu, c in row.items() if c != 0}


def N_table(n):
    """Full N-matrix for all s_lambda with |lambda| = n: {lambda: N_row}."""
    return {lam: N_row(lam) for lam in integer_partitions(n)}


def negative_terms(row):
    """The {mu: coefficient} entries of an N-row that are negative (the
    nonnegativity hypothesis says this is always empty)."""
    return {mu: c for mu, c in row.items() if c < 0}


def format_N_row(lam, row):
    """One-line 's_lambda = c1 s~_mu1 + ...' string (signs shown)."""
    def mu_str(mu):
        return "()" if not mu else "".join(str(p) for p in mu)

    ordered = sorted(row.items(), key=lambda kv: (sum(kv[0]), kv[0]))
    terms = [f"{c} s~_{mu_str(mu)}" for mu, c in ordered]
    lam_str = "".join(str(p) for p in lam) or "0"
    return f"s_{lam_str} = " + " + ".join(terms)

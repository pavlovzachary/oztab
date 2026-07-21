r"""Closed form for N_{(n,1),mu}: the s~_mu coefficient of the Schur function
s_{(n,1)}, for EVERY mu at once (not just the empty shape).

Why (n,1) is special.  Jacobi-Trudi collapses to two terms,

    s_{(n,1)} = h_n h_1 - h_{n+1},        so   N_{(n,1),mu} = M_{(n,1),mu} - M_{(n+1),mu},

and both M's are counted at the same m = |lambda| = n+1, so the body/first-row
split of tableaux.py applies to both with the same blank count.

The two contents involved are {{1^n, 2}} and {{1^(n+1)}}, i.e. at most ONE cell
of a tableau can contain the label 2.  That makes every ingredient a classical
generating function:

  * A body cell built from 1's only is {1^a}, and {1^a} < {1^b} iff a < b, so a
    pure body is exactly an SSYT of shape mu with positive integer entries (the
    entry = the multiplicity a).  Weighting by the entry sum gives the principal
    specialization

        F_mu(q) = s_mu(q, q^2, q^3, ...) = q^(|mu| + n(mu)) / prod_c (1 - q^h(c)),

    n(mu) = sum i*mu_i (0-based), h(c) = hook length.  [q^j] F_mu = # bodies of
    content {{1^j}}.

  * Every cell containing a 2 is LARGER than every pure cell ({1^a} < {1^b, 2}
    for all a, b, since a shorter tuple that is a prefix is smaller and 1 < 2).
    So in a body of content {{1^j, 2}} the unique 2-cell is greater than its
    row-left and column-top neighbours automatically, and cannot have anything
    to its right or below.  Hence it sits at an OUTER CORNER of mu, holds
    {1^a, 2} for any a >= 0, and the rest of the body is a pure filling of
    mu minus that corner:

        G_mu(q) = (1 / (1 - q)) * sum over outer corners c of F_{mu \ c}(q).

  * First rows are multiset partitions of the leftover.  For {{1^i}} that is
    p(i); for {{1^i, 2}} it is Q(i) = sum_{k<=i} p(k) (split off the block
    holding the 2).

Assembling, with P(q) = prod 1/(1-q^k) the partition generating function,

    sum_{n>=0} N_{(n,1),mu} q^n = F_mu(q) P(q) / (1-q)     [ pure body, 2 in row 1 ]
                                + G_mu(q) P(q)             [ 2 inside the body ]
                                - (F_mu(q) P(q))_+ / q     [ minus h_{n+1} ]

where X_+ drops the constant term before dividing by q.  `N_n1_series` is
literally this formula; `tests/test_n1_formula.py` pins it against the brute
force N_row for every mu up to size 6.

Specializations worth naming (all verified to n = 13):

    mu = ()      sum_n N q^n:  gives  N = sum_{k<=n} p(k) - p(n+1)
                 (the empty-shape theorem in proofs/null_coefficient_S_n-1_1.tex,
                 stated there for lambda = (n-1,1), i.e. shifted by one)

    mu = (k)     sum_n N q^n = P(q) q^k (2 - q^(k-1)) / ((1-q) prod_{i<=k}(1-q^i))
                 so N_{(n,1),(k)} = sum_j [2 p_k(j) - p_k(j-k+1)] * Q(n-k-j),
                 with p_k = partitions into parts <= k.  Nonnegative because
                 p_k is nondecreasing.  Stable top value N_{(n,1),(n)} = 2.

    mu = (1^k)   sum_n N q^n = P(q) q^(k(k-1)/2) (1 - q^(k-1) + q^k)
                                / ((1-q) prod_{i<=k}(1-q^i))
"""

from functools import lru_cache

from .multisets import integer_partitions

# --- truncated power series: plain lists, index = power of q -----------------


def _zero(deg):
    return [0] * deg


def _mul(a, b, deg):
    out = [0] * deg
    for i, x in enumerate(a):
        if not x:
            continue
        for j, y in enumerate(b):
            if i + j >= deg:
                break
            out[i + j] += x * y
    return out


def _add(a, b):
    return [x + y for x, y in zip(a, b)]


def _sub(a, b):
    return [x - y for x, y in zip(a, b)]


@lru_cache(maxsize=None)
def _inv_1mqk_cached(k, deg):
    out = [0] * deg
    for i in range(0, deg, k):
        out[i] = 1
    return tuple(out)


def _inv_1mqk(k, deg):
    """1 / (1 - q^k)."""
    return list(_inv_1mqk_cached(k, deg))


def _drop_and_divide_q(a, deg):
    """(a - a_0) / q, i.e. shift the series down by one after killing a_0."""
    return list(a[1:]) + [0]


@lru_cache(maxsize=None)
def _partition_series_cached(deg):
    out = [0] * deg
    out[0] = 1
    for k in range(1, deg):
        out = _mul(out, _inv_1mqk(k, deg), deg)
    return tuple(out)


def partition_series(deg):
    """P(q) = prod_k 1/(1-q^k) truncated below q^deg."""
    return list(_partition_series_cached(deg))


# --- shape combinatorics -----------------------------------------------------


def hook_lengths(mu):
    """Hook lengths of the cells of mu, in reading order."""
    mu = tuple(mu)
    if not mu:
        return []
    conj = [sum(1 for part in mu if part > j) for j in range(mu[0])]
    return [mu[i] - j + conj[j] - i - 1
            for i in range(len(mu)) for j in range(mu[i])]


def outer_corners(mu):
    """The shapes mu \\ c, one per outer (removable) corner c of mu."""
    mu = tuple(mu)
    out = []
    for i in range(len(mu)):
        if i + 1 == len(mu) or mu[i] > mu[i + 1]:
            nu = list(mu)
            nu[i] -= 1
            out.append(tuple(part for part in nu if part > 0))
    return out


def principal_specialization(mu, deg):
    """F_mu(q) = s_mu(q, q^2, q^3, ...) truncated below q^deg.

    [q^j] F_mu = number of SSYT of shape mu with positive integer entries
    summing to j = number of column-strict pure-1 multiset fillings of mu with
    content {{1^j}}.
    """
    mu = tuple(mu)
    out = [0] * deg
    exponent = sum(mu) + sum(i * part for i, part in enumerate(mu))
    if exponent >= deg:
        return out
    out[exponent] = 1
    for h in hook_lengths(mu):
        out = _mul(out, _inv_1mqk(h, deg), deg)
    return out


# --- the coefficient itself --------------------------------------------------


def N_n1_series(mu, nmax):
    """[N_{(n,1),mu} for n in 0..nmax], from the closed form in the module docstring.

    Costs no tableau enumeration, so it reaches n far beyond what N_row can.
    """
    mu = tuple(mu)
    # one spare degree: the -h_{n+1} term reads FP at q^(n+1), not q^n
    deg = nmax + 2
    P = partition_series(deg)
    inv_1mq = _inv_1mqk(1, deg)

    F = principal_specialization(mu, deg)
    G = _zero(deg)
    for nu in outer_corners(mu):
        G = _add(G, principal_specialization(nu, deg))
    G = _mul(G, inv_1mq, deg)

    FP = _mul(F, P, deg)
    series = _add(_mul(FP, inv_1mq, deg), _mul(G, P, deg))
    return _sub(series, _drop_and_divide_q(FP, deg))[:nmax + 1]


def N_n1(mu, n):
    """N_{(n,1),mu}, the coefficient of s~_mu in the Schur function s_{(n,1)}."""
    return N_n1_series(mu, n)[n]


def N_n1_row(n, nmax_mu=None):
    """{mu: N_{(n,1),mu}} over all mu with |mu| <= n+1, nonzero entries only.

    The closed-form counterpart of schur.N_row((n, 1)).
    """
    limit = n + 1 if nmax_mu is None else nmax_mu
    row = {}
    for size in range(limit + 1):
        for mu in integer_partitions(size):
            value = N_n1(mu, n)
            if value:
                row[mu] = value
    return row

"""The multiplicity-profile decomposition of M_{(a,b),()} -- pure Python.

This is the machinery from the Sage notebook (algebraic-combinatorics/), ported
so the injection programme runs without Sage.

Setting.  For the coefficient of s~_() in h_a h_b, the single-row multiset
tableau splits into sections indexed by the distinct parts k of a partition
lambda |- b: the cells whose b-content is k copies of the b-label.  The a-units
are distributed over a "pure" section (cells with no b-content) and one section
per distinct part.  Within the section for part k the cells are weakly ordered by
colex rank, and the multiplicity n_k = mult of k in lambda bounds the number of
cells, hence the number of parts of the partition of that section's a-units:

    C(a, lambda) = sum_{J_0 + sum_k J_k = a}  p(J_0) * prod_k p_{<= n_k}(J_k)

    M_{(a,b),()} = sum_{lambda |- b} C(a, lambda)

Profile invariance.  C(a, lambda) depends only on the multiset of multiplicities
{n_k}, not on the part values k -- immediately from the formula, since k appears
nowhere in it.  So partitions of b sharing a `profile` contribute equally and the
sum collapses onto distinct profiles:

    M_{(a,b),()} = sum_pi  count(b, pi) * C(a, pi)

`M_empty_via_profiles` is that collapsed sum; `cross_check` pins it against the
tableau enumerator in coefficients.M_row, which is the independent derivation.
"""

from collections import Counter
from functools import lru_cache
from itertools import product

from .multisets import integer_partitions


@lru_cache(maxsize=None)
def p(n):
    """Number of integer partitions of n."""
    if n < 0:
        return 0
    return p_leq(n, n) if n else 1


@lru_cache(maxsize=None)
def p_leq(n, m):
    """Number of partitions of n into at most m parts."""
    if n < 0 or m < 0:
        return 0
    if n == 0:
        return 1
    if m == 0:
        return 0
    # partitions of n into at most m parts = (at most m-1 parts) + (exactly m,
    # i.e. subtract 1 from each of the m parts)
    return p_leq(n, m - 1) + p_leq(n - m, m)


def P(a):
    """Cumulative partition count sum_{k=0}^{a} p(k)."""
    return sum(p(k) for k in range(a + 1))


def P_leq(m, n):
    """Cumulative sum_{j=0}^{m} p_{<= n}(j)."""
    return sum(p_leq(j, n) for j in range(m + 1))


def profile(lam):
    """Multiplicity profile of a partition: its multiplicities, sorted
    decreasing.  (2,1) and (3,1) both -> (1,1);  (1,1) and (2,2) both -> (2,)."""
    counts = Counter(p_ for p_ in lam if p_ > 0)
    return tuple(sorted(counts.values(), reverse=True))


def _compositions(total, parts):
    """All length-`parts` tuples of nonneg ints summing to `total`."""
    if parts == 0:
        return [()] if total == 0 else []
    if parts == 1:
        return [(total,)]
    out = []
    for first in range(total + 1):
        for rest in _compositions(total - first, parts - 1):
            out.append((first,) + rest)
    return out


@lru_cache(maxsize=None)
def contribution(a, lam):
    """C(a, lambda): the contribution of a single lambda |- b to M_{(a,b),()}.

    Depends only on profile(lam) -- see `contribution_by_profile`.
    """
    mults = profile(lam)
    if not mults:
        return p(a)
    total = 0
    for split in _compositions(a, len(mults) + 1):
        value = p(split[0])
        for i, n_k in enumerate(mults):
            value *= p_leq(split[i + 1], n_k)
            if value == 0:
                break
        total += value
    return total


@lru_cache(maxsize=None)
def contribution_by_profile(a, pi):
    """C(a, pi) for a profile pi directly (a profile is itself a valid multiset
    of multiplicities, so the same formula applies)."""
    if not pi:
        return p(a)
    total = 0
    for split in _compositions(a, len(pi) + 1):
        value = p(split[0])
        for i, n_k in enumerate(pi):
            value *= p_leq(split[i + 1], n_k)
            if value == 0:
                break
        total += value
    return total


def profile_counts(b):
    """{profile: number of lambda |- b with that profile}."""
    counts = Counter()
    for lam in integer_partitions(b):
        counts[profile(lam)] += 1
    return dict(counts)


def M_empty_via_profiles(a, b):
    """M_{(a,b),()} = sum_pi count(b,pi) * C(a,pi) -- the collapsed sum."""
    return sum(mult * contribution_by_profile(a, pi)
               for pi, mult in profile_counts(b).items())


def N_empty_via_profiles(a, b):
    """N_{(a,b),()} = M_{(a,b),()} - M_{(a+1,b-1),()}, the Jacobi-Trudi difference."""
    return M_empty_via_profiles(a, b) - M_empty_via_profiles(a + 1, b - 1)


def cross_check(max_n=8):
    """Check M_empty_via_profiles against the tableau enumerator (M_row) for all
    (a,b) with a >= b >= 0 and a + b <= max_n.  Returns list of failures."""
    from .coefficients import M_row

    failures = []
    for n in range(max_n + 1):
        for b in range(n // 2 + 1):
            a = n - b
            lam = (a, b) if b else (a,)
            enumerated = M_row(lam).get((), 0)
            formula = M_empty_via_profiles(a, b)
            if enumerated != formula:
                failures.append((a, b, enumerated, formula))
    return failures


def assignment_search(b, a_hi=45, show=6):
    """Search profile-level injection schemes for the Jacobi-Trudi difference.

    For each source lambda |- (b-1) we must injectively absorb its C(a+1, lambda)
    objects into a disjoint set of targets lambda' |- b.  A scheme is `feasible`
    from a0 if for every source, the summed target contributions dominate the
    source demand for all a in [a0, a_hi].  Returns schemes ranked by a0 (small
    a0 = the injection works further down, so less special-casing).

    This is the b=3 search from the notebook; it is exponential in p(b) and gets
    heavy past b = 5.
    """
    sources = [tuple(lam) for lam in integer_partitions(b - 1)]
    targets = [tuple(lam) for lam in integer_partitions(b)]
    n_src, n_tgt = len(sources), len(targets)

    @lru_cache(maxsize=None)
    def threshold(src_pi, tgt_pis):
        """Smallest a0 with sum_t C(a,t) >= C(a+1,src) throughout [a0, a_hi]."""
        ok_from = None
        for a in range(1, a_hi + 1):
            lhs = sum(contribution_by_profile(a, t) for t in tgt_pis)
            if lhs >= contribution_by_profile(a + 1, src_pi):
                if ok_from is None:
                    ok_from = a
            else:
                ok_from = None      # a later failure invalidates the earlier run
        return ok_from

    best = {}
    for bins in product(range(n_src + 1), repeat=n_tgt):
        groups = [[] for _ in range(n_src + 1)]     # last bin = unused targets
        for j, bin_index in enumerate(bins):
            groups[bin_index].append(j)
        worst = 1
        for i, src in enumerate(sources):
            if not groups[i]:
                worst = None
                break
            tgt_pis = tuple(sorted(profile(targets[j]) for j in groups[i]))
            t0 = threshold(profile(src), tgt_pis)
            if t0 is None:
                worst = None
                break
            worst = max(worst, t0)
        if worst is None:
            continue
        signature = tuple(sorted(
            (profile(src), tuple(sorted(profile(targets[j]) for j in groups[i])))
            for i, src in enumerate(sources)))
        if signature not in best or worst < best[signature][0]:
            scheme = [(sources[i], [targets[j] for j in groups[i]])
                      for i in range(n_src)]
            leftover = [targets[j] for j in groups[n_src]]
            best[signature] = (worst, scheme, leftover)

    ranked = sorted(best.values(), key=lambda entry: entry[0])
    return ranked[:show] if show else ranked

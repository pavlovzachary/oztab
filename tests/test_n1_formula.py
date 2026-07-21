"""The closed form for N_{(n,1),mu} must agree with the brute-force N_row.

N_row((n,1)) enumerates multiset tableaux; oztab.n1 evaluates a generating
function.  They are completely independent routes to the same integers, so this
is the test that would catch a wrong hook length, a missed outer corner, or an
off-by-one in the h_{n+1} subtraction.
"""

import pytest

from oztab import N_n1, N_n1_row, N_n1_series, principal_specialization
from oztab.multisets import integer_partitions_upto
from oztab.schur import N_row
from oztab.tableaux import count_column_strict_fillings

NMAX = 8


@pytest.mark.parametrize("n", range(1, NMAX + 1))
def test_closed_form_matches_enumeration(n):
    assert N_n1_row(n) == N_row((n, 1))


def test_empty_shape_is_the_proved_theorem():
    """N_{(n,1),()} = sum_{k<=n} p(k) - p(n+1), the result proved in
    proofs/null_coefficient_S_n-1_1.tex (stated there for lambda = (n-1,1))."""
    from oztab.n1 import partition_series

    P = partition_series(NMAX + 2)
    for n in range(0, NMAX + 1):
        assert N_n1((), n) == sum(P[:n + 1]) - P[n + 1]


@pytest.mark.parametrize("mu", integer_partitions_upto(5))
def test_principal_specialization_counts_pure_bodies(mu):
    """[q^j] s_mu(q,q^2,...) counts shape-mu fillings of content {{1^j}}."""
    deg = 10
    series = principal_specialization(mu, deg)
    for j in range(deg):
        expected = count_column_strict_fillings(tuple(mu), (1,) * j)
        assert series[j] == expected


def test_single_row_closed_form():
    """N_{(n,1),(k)}: coefficients of P(q) q^k (2 - q^{k-1}) / ((1-q) prod_{i<=k}(1-q^i))."""
    from oztab.n1 import _inv_1mqk, _mul, partition_series

    deg = 16
    for k in range(1, 6):
        series = _mul(partition_series(deg), _inv_1mqk(1, deg), deg)
        for i in range(1, k + 1):
            series = _mul(series, _inv_1mqk(i, deg), deg)
        numerator = [0] * deg
        numerator[k] += 2
        numerator[2 * k - 1] -= 1
        series = _mul(series, numerator, deg)
        assert series == N_n1_series((k,), deg - 1)


def test_stable_top_coefficients():
    """The unitriangularity fringe: 1 at mu = (n,1), 2 at mu = (n), 3 at (n-1,1)."""
    for n in range(3, NMAX + 1):
        assert N_n1((n, 1), n) == 1
        assert N_n1((n,), n) == 2
        assert N_n1((n - 1, 1), n) == 3

"""Schur-side (N) tests.

No Sage here, so the anchors are identities that must hold structurally:
  * s_(n) = h_n, so N_row((n,)) == M_row((n,));
  * two-row Jacobi-Trudi s_(a,b) = h_a h_b - h_{a+1} h_{b-1} must agree with the
    general determinant in schur.jacobi_trudi_terms;
plus the frozen N(a,3) column computed in Sage and recorded in the b=3 table of
algebraic-combinatorics/progress_summary.tex.
"""

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from oztab import M_row, N_row, jacobi_trudi_terms  # noqa: E402


def test_jt_terms_one_row():
    # s_(n) = h_n : a single term, no sign
    assert jacobi_trudi_terms((4,)) == [(1, (4,))]


def test_jt_terms_two_row():
    # s_(a,b) = h_a h_b - h_{a+1} h_{b-1}
    assert sorted(jacobi_trudi_terms((3, 2))) == sorted([(1, (3, 2)), (-1, (4, 1))])
    # b = 1: the negative term is h_4 h_0 = h_4
    assert sorted(jacobi_trudi_terms((3, 1))) == sorted([(1, (3, 1)), (-1, (4,))])


def test_one_row_schur_is_h():
    for n in range(1, 8):
        assert N_row((n,)) == M_row((n,))


def test_two_row_matches_explicit_jt():
    """The general determinant must reproduce the hand two-row identity."""
    for n in range(2, 9):
        for b in range(1, n // 2 + 1):
            a = n - b
            explicit = {}
            for mu, c in M_row((a, b)).items():
                explicit[mu] = explicit.get(mu, 0) + c
            for mu, c in M_row((a + 1, b - 1)).items():
                explicit[mu] = explicit.get(mu, 0) - c
            explicit = {mu: c for mu, c in explicit.items() if c}
            assert N_row((a, b)) == explicit, f"(a,b)=({a},{b})"


def test_N_a3_frozen_from_sage():
    """N(a,3) = coefficient of s~_() in s_(a,3), from the progress_summary table
    (computed in Sage as int(st(s[a,3]).coefficient([])))."""
    expected = {3: 2, 4: 10, 5: 20, 6: 44, 7: 76, 8: 134}
    for a, value in expected.items():
        assert N_row((a, 3)).get((), 0) == value, f"a={a}"


def test_null_coefficient_theorem():
    """proofs/null_coefficient_S_n-1_1.tex, Theorem 2:

        N_{(n-1,1), ()} = sum_{k=0}^{n-1} p(k)  -  p(n)

    together with its Proposition 3, mp({{1^(n-1), 2}}) = sum_{k=0}^{n-1} p(k).
    """
    from oztab.multisets import multiset_partitions
    from oztab.profiles import p

    for n in range(2, 12):
        cumulative = sum(p(k) for k in range(n))
        # Proposition 3
        content = tuple([1] * (n - 1) + [2])
        assert len(multiset_partitions(content)) == cumulative, f"Prop 3, n={n}"
        # Theorem 2
        assert N_row((n - 1, 1)).get((), 0) == cumulative - p(n), f"Thm 2, n={n}"


def test_N_a2_matches_solved_closed_form():
    """The solved case b = 2 (progress_summary, sec. 6):

        N_{(a,2),()} = C(a, (1,1)) - p(a+1)

    This is a genuine cross-check: the left side comes from the tableau
    enumerator via Jacobi-Trudi, the right side from the profile formula.
    """
    from oztab.profiles import contribution, p

    for a in range(2, 11):
        assert N_row((a, 2)).get((), 0) == contribution(a, (1, 1)) - p(a + 1)

"""Regression tests: pure Python (pytest), no Sage needed.

The target counts are taken from Orellana-Zabrocki, "Symmetric group characters
as symmetric functions": Example 10 (h_{2,1}, 20 tableaux) and Example 11
(h_{1^4}).  These freeze the enumerator against the paper.
"""

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from oztab import M_row, content_of, generate_tableaux, is_valid  # noqa: E402
from oztab.multisets import multiset_partitions  # noqa: E402


def test_multiset_partition_counts():
    assert len(multiset_partitions((1, 1, 2))) == 4
    assert len(multiset_partitions((1, 1, 1))) == 3
    assert len(multiset_partitions((1, 2, 3))) == 5      # Bell(3)
    assert len(multiset_partitions((1, 2, 3, 4))) == 15  # Bell(4)


def test_h_21_example10():
    assert M_row((2, 1)) == {(): 4, (1,): 7, (2,): 4, (1, 1): 3, (3,): 1, (2, 1): 1}
    assert sum(M_row((2, 1)).values()) == 20


def test_h_3():
    assert M_row((3,)) == {(): 3, (1,): 4, (2,): 2, (1, 1): 1, (3,): 1}


def test_h_1111_example11_visible_terms():
    row = M_row((1, 1, 1, 1))
    for mu, expected in [((), 15), ((1,), 37), ((1, 1), 31),
                         ((1, 1, 1), 10), ((1, 1, 1, 1), 1), ((2,), 31)]:
        assert row.get(mu, 0) == expected


def test_empty_shape_is_multiset_partitions():
    # coefficient of s~_() in h_lambda == number of multiset partitions of content
    for lam in [(2, 1), (3,), (1, 1, 1, 1), (3, 2)]:
        assert M_row(lam).get((), 0) == len(multiset_partitions(content_of(lam)))


def test_generated_tableaux_are_valid():
    for lam in [(2, 1), (3,), (2, 2), (3, 1)]:
        content = content_of(lam)
        for mu in [(), (1,), (2,), (1, 1), (2, 1)]:
            for tab in generate_tableaux(content, mu):
                assert is_valid(tab)
                assert tab.shape_index() == mu
                assert tab.content() == content

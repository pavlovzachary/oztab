"""The counting path must agree with the generator, exactly, always.

coefficients.M_row counts via tableaux.count_tableaux, which never builds the
tableaux (there are ~10x more of them at each n, so building them is what caps
the reachable n).  That speed is only safe if the count equals what the
generator actually produces -- these tests are the thing that guarantees it.
"""

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from oztab import content_of, generate_tableaux, count_tableaux  # noqa: E402
from oztab.multisets import integer_partitions, integer_partitions_upto  # noqa: E402
from oztab.tableaux import (column_strict_fillings,  # noqa: E402
                            count_column_strict_fillings)


def test_count_matches_generator_exhaustively():
    # every (lambda, mu) with |lambda| <= 6
    for n in range(7):
        for lam in integer_partitions(n):
            content = content_of(lam)
            for mu in integer_partitions_upto(n):
                built = len(generate_tableaux(content, mu))
                counted = count_tableaux(content, mu)
                assert built == counted, f"lambda={lam} mu={mu}: {built} vs {counted}"


def test_filling_count_matches_filling_list():
    for n in range(6):
        for lam in integer_partitions(n):
            content = content_of(lam)
            for mu in integer_partitions_upto(n):
                built = len(column_strict_fillings(mu, content))
                counted = count_column_strict_fillings(mu, content)
                assert built == counted, f"mu={mu} content={content}"

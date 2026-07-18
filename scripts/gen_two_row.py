"""Expansions of small two-row partitions lambda = (a, b), a >= b, a,b <= MAX.

For each two-row partition this prints two change-of-basis expansions into the
Orellana-Zabrocki s~ basis, both computed natively (no Sage required):

  h-side :  h_(a,b)  =  sum_mu  M_{(a,b),mu}  s~_mu     (M >= 0, counts tableaux)
  Schur  :  s_(a,b)  =  sum_mu  N_{(a,b),mu}  s~_mu     (N is signed via Jacobi-Trudi)

Usage:
    python gen_two_row.py          # a,b <= 5
    python gen_two_row.py 6        # a,b <= 6

Mirror of gen_n1.py: a thin driver over the public API, nothing new in src/.
To cross-check under a Sage kernel:  st(h[[a,b]]) should match M_row((a,b)).
"""

import sys

from oztab import M_row, format_row, N_row, format_N_row


def two_row_partitions(max_part):
    """(a, b) with max_part >= a >= b >= 1."""
    for a in range(1, max_part + 1):
        for b in range(1, a + 1):
            yield (a, b)


def main(max_part=5):
    for lam in two_row_partitions(max_part):
        print(f"\n{'=' * 60}")
        print(f"lambda = {lam}   |lambda| = {sum(lam)}")
        print(f"{'=' * 60}")
        print(f"  {format_row(lam, M_row(lam))}")
        print(f"  {format_N_row(lam, N_row(lam))}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

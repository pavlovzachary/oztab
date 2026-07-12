"""Coefficient tables built from the tableau generator.

M_lambda_mu is the coefficient of s~_mu in h_lambda; it equals the number of
Orellana-Zabrocki tableaux of content(lambda) with body shape mu.  Everything
here is derived from `generate_tableaux`, so the counts and the objects can
never disagree.
"""

from .multisets import integer_partitions, integer_partitions_upto
from .tableaux import generate_tableaux


def content_of(lam):
    """Content multiset of h_lambda: lambda_i copies of label i (sorted tuple).
    e.g. (2, 1) -> (1, 1, 2)."""
    elems = []
    for label, mult in enumerate(lam, start=1):
        elems.extend([label] * mult)
    return tuple(elems)


def M_row(lam):
    """Row of the M-matrix for h_lambda: {mu: coefficient} over nonzero mu."""
    content = content_of(lam)
    n = len(content)
    row = {}
    for mu in integer_partitions_upto(n):
        count = len(generate_tableaux(content, mu))
        if count:
            row[mu] = count
    return row


def M_table(n):
    """Full M-matrix for all h_lambda with |lambda| = n: {lambda: M_row}."""
    return {lam: M_row(lam) for lam in integer_partitions(n)}


def format_row(lam, row):
    """One-line human-readable 'h_lambda = c1 s~_mu1 + ...' string."""
    def mu_str(mu):
        return "()" if not mu else "".join(str(p) for p in mu)
    terms = [f"{c} s~_{mu_str(mu)}" for mu, c in sorted(row.items(), key=lambda kv: (sum(kv[0]), kv[0]))]
    return f"h_{''.join(str(p) for p in lam) or '0'} = " + " + ".join(terms)

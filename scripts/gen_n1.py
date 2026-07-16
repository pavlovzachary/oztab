"""Generate all Orellana-Zabrocki tableaux for lambda = (n, 1), n = 1..5.

For lambda = (n, 1) the content of h_lambda is the multiset {{1^n, 2}}.
The full set of tableaux is the union, over every body shape mu with
|mu| <= |lambda|, of generate_tableaux(content, mu).  Each tableau
contributes s~_{shape(body)}.
"""

from collections import Counter

from oztab import generate_tableaux
from oztab.multisets import integer_partitions_upto


def tableaux_for_lambda(lam):
    """All OZ tableaux for partition lam, grouped by body shape mu."""
    content = tuple(sorted(
        label
        for label, mult in enumerate(lam, start=1)
        for _ in range(mult)
    ))
    n = len(content)
    by_mu = {}
    for mu in integer_partitions_upto(n):
        ts = generate_tableaux(content, mu)
        if ts:
            by_mu[mu] = ts
    return content, by_mu


for n in range(1, 6):
    lam = (n, 1)
    content, by_mu = tableaux_for_lambda(lam)
    total = sum(len(v) for v in by_mu.values())
    print(f"\n{'='*60}")
    print(f"lambda = {lam}   content = {content}   ({total} tableaux total)")
    print(f"{'='*60}")
    for mu in sorted(by_mu, key=lambda m: (sum(m), m)):
        ts = by_mu[mu]
        print(f"\n  body shape mu = {mu or '()'}  ->  s~_{mu or '()'}   "
              f"({len(ts)} tableau{'x' if len(ts) != 1 else ''})")
        for t in ts:
            rendered = t.render(faithful=False).replace("\n", "\n      ")
            print(f"      {rendered}")

"""Test injections proving N_{(a,b),()} >= 0 by hand, for two-row lambda.

At mu = () the body is empty, so an OZ tableau IS just a multiset partition of
the whole content.  With s_(a,b) = h_(a,b) - h_(a+1,b-1) that makes the
empty-shape coefficient a difference of two counts of multiset partitions:

    N_{(a,b),()} = #MP{1^a 2^b} - #MP{1^(a+1) 2^(b-1)},

so proving nonnegativity means injecting  MP{1^(a+1) 2^(b-1)}  into
MP{1^a 2^b}.  The natural move is to turn one 1 into a 2 ("bump"); the question
is which cell to bump.  Since cells are sorted tuples, the lex order used
throughout oztab is Python's tuple order, and the lex-largest cell is the
rightmost one in the weakly-increasing display.

Reformulation worth knowing.  A cell {1^i 2^j} is just the pair (i,j), so a
multiset partition of {1^a 2^b} is a PARTITION OF THE BIPARTITE NUMBER (a,b)
(Auluck / Cheema-Motzkin style bipartite partitions), counted by

    sum_{a,b} #MP{1^a 2^b} x^a y^b = prod_{(i,j) != (0,0)} 1 / (1 - x^i y^j).

Writing c_k for the antidiagonal a+b = n, the whole mu = () two-row question is
exactly:  is (c_0, ..., c_n) symmetric and unimodal?  Symmetry is the cell
involution (i,j) -> (j,i); N_{(a,b),()} = c_a - c_{a+1} is one unimodality step.

Findings (run `python null_injection.py` to regenerate all of it):
  * "bump the rightmost 1" is injective for b <= 2 and FAILS first at (3,3).
    Every collision has the same shape -- see `report_collisions`.
  * No greedy single-cell rule tried here (either direction) gets past b = 2.
  * But a maximum bipartite matching SATURATES in every case tested, so a
    single-bump injection does exist for b >= 3; only a canonical rule for it
    is missing.  That is the gap to close.

Usage:  python null_injection.py [--upto N]     # default a+b <= 11
"""

import argparse
from collections import defaultdict

from oztab import N_row
from oztab.multisets import multiset_partitions

ones = lambda cell: sum(1 for x in cell if x == 1)
twos = lambda cell: sum(1 for x in cell if x == 2)


def content(a, b):
    return tuple([1] * a + [2] * b)


def partitions_of(a, b):
    """MP{1^a 2^b} = the partitions of the bipartite number (a, b)."""
    return multiset_partitions(content(a, b))


def show(cells):
    return " ".join("{" + ",".join(str(x) for x in c) + "}" for c in cells)


# --- the move ---------------------------------------------------------------


def bump(cell):
    """Turn one 1 into a 2.  Cells are sorted, so cell[0] is a 1 when present."""
    return tuple(sorted(cell[1:] + (2,)))


def apply_bump(P, pick):
    """Bump the cell chosen by `pick` among the cells of P that contain a 1."""
    candidates = [c for c in P if 1 in c]
    if not candidates:
        return None
    chosen = pick(candidates)
    rest = list(P)
    rest.remove(chosen)
    return tuple(sorted(rest + [bump(chosen)]))


def reachable(P):
    """Every partition obtainable from P by one bump, in any single cell."""
    out = set()
    for i, cell in enumerate(P):
        if 1 in cell:
            out.add(tuple(sorted(P[:i] + P[i + 1:] + (bump(cell),))))
    return out


RULES = {
    "rightmost (lex-max)": max,
    "leftmost (lex-min)": min,
    "pure 1-cells first": lambda cs: max([c for c in cs if twos(c) == 0] or cs),
    "most 1s, then fewest 2s": lambda cs: max(cs, key=lambda c: (ones(c), -twos(c))),
    "most 1s, then most 2s": lambda cs: max(cs, key=lambda c: (ones(c), twos(c))),
}


# --- the reverse move -------------------------------------------------------
# A canonical SURJECTION MP{1^a 2^b} ->> MP{1^(a+1) 2^(b-1)} proves the same
# inequality, so the un-bump direction is worth testing on its own.


def unbump(cell):
    """Turn one 2 into a 1.  Cells are sorted, so cell[-1] is a 2 when present."""
    return tuple(sorted(cell[:-1] + (1,)))


def apply_unbump(T, pick):
    candidates = [c for c in T if 2 in c]
    if not candidates:
        return None
    chosen = pick(candidates)
    rest = list(T)
    rest.remove(chosen)
    return tuple(sorted(rest + [unbump(chosen)]))


REVERSE_RULES = {
    "lex-max cell with a 2": max,
    "lex-min cell with a 2": min,
    "most 2s, then most 1s": lambda cs: max(cs, key=lambda c: (twos(c), ones(c))),
    "most 2s, then fewest 1s": lambda cs: max(cs, key=lambda c: (twos(c), -ones(c))),
    "fewest 1s, then most 2s": lambda cs: max(cs, key=lambda c: (-ones(c), twos(c))),
}


# --- experiments ------------------------------------------------------------


def two_row_cases(upto):
    return [(n - b, b) for n in range(2, upto + 1)
            for b in range(1, n // 2 + 1) if n - b >= b]


def check_rules(cases):
    """For each choice rule, the first (a,b) where the bump map is not injective."""
    print("Is 'bump one 1' injective, by choice-of-cell rule?\n")
    width = max(len(k) for k in RULES)
    print(f"  {'rule':{width}} | first failure | injective for")
    print("  " + "-" * (width + 34))
    for label, pick in RULES.items():
        first_bad, results = None, {}
        for (a, b) in cases:
            src = partitions_of(a + 1, b - 1)
            images = {apply_bump(P, pick) for P in src}
            results[(a, b)] = len(images) == len(src)
            if not results[(a, b)] and first_bad is None:
                first_bad = (a, b)
        good_b = [b for b in sorted({b for _, b in cases})
                  if all(v for (_, bb), v in results.items() if bb == b)]
        print(f"  {label:{width}} | {str(first_bad):13} | b in {good_b}")


def check_reverse_rules(cases):
    """For each un-bump rule, the first (a,b) where the map is not surjective."""
    print("\nIs 'un-bump one 2' surjective, by choice-of-cell rule?\n")
    width = max(len(k) for k in REVERSE_RULES)
    print(f"  {'rule':{width}} | first failure | surjective for")
    print("  " + "-" * (width + 34))
    for label, pick in REVERSE_RULES.items():
        first_bad, results = None, {}
        for (a, b) in cases:
            src = set(partitions_of(a + 1, b - 1))
            image = {apply_unbump(T, pick) for T in partitions_of(a, b)}
            assert image <= src
            results[(a, b)] = image == src
            if not results[(a, b)] and first_bad is None:
                first_bad = (a, b)
        good_b = [b for b in sorted({b for _, b in cases})
                  if all(v for (_, bb), v in results.items() if bb == b)]
        print(f"  {label:{width}} | {str(first_bad):13} | b in {good_b}")


def report_collisions(a, b, limit=4):
    """Show the colliding preimages of the rightmost-1 rule."""
    fibers = defaultdict(list)
    for P in partitions_of(a + 1, b - 1):
        fibers[apply_bump(P, max)].append(P)
    collisions = {k: v for k, v in fibers.items() if len(v) > 1}
    print(f"\nCollisions of the rightmost-1 rule at (a,b) = ({a},{b}): "
          f"{len(collisions)}")
    for image, preimages in sorted(collisions.items())[:limit]:
        print(f"    image  {show(image)}")
        for P in preimages:
            hit = max(c for c in P if 1 in c)
            print(f"      <-  {show(P):30s} (bumped {show([hit])} -> {show([bump(hit)])})")
    if len(collisions) > limit:
        print(f"    ... {len(collisions) - limit} more, all of the same shape")


def max_matching(left, adj):
    """Kuhn's algorithm; returns {right: left}."""
    match = {}

    def augment(u, seen):
        for v in adj[u]:
            if v in seen:
                continue
            seen.add(v)
            if v not in match or augment(match[v], seen):
                match[v] = u
                return True
        return False

    for u in left:
        augment(u, set())
    return match


def check_matching(cases):
    """Does ANY single-bump injection exist?  (Saturating matching = yes.)"""
    print("\nDoes any single-bump injection exist (maximum bipartite matching)?\n")
    print(f"  {'(a,b)':>8} {'|src|':>7} {'|tgt|':>7} {'matched':>8}  saturates?")
    for (a, b) in cases:
        src = partitions_of(a + 1, b - 1)
        adj = {P: reachable(P) for P in src}
        match = max_matching(src, adj)
        print(f"  {str((a, b)):>8} {len(src):>7} {len(partitions_of(a, b)):>7} "
              f"{len(match):>8}  {'YES' if len(match) == len(src) else 'NO':>9}")


def check_antidiagonals(upto):
    """The bipartite-partition reformulation: symmetry + unimodality."""
    print("\nAntidiagonals c_k = #partitions of the bipartite number (k, n-k):\n")
    for n in range(2, upto + 1):
        seq = [len(partitions_of(k, n - k)) for k in range(n + 1)]
        symmetric = seq == seq[::-1]
        unimodal = all(x <= y for x, y in zip(seq, seq[1:n // 2 + 1]))
        print(f"  n={n:2d}  {str(seq):52s} symmetric={symmetric} unimodal={unimodal}")

    print("\n  cross-check N_{(a,b),()} == c_a - c_(a+1):")
    bad = 0
    for (a, b) in two_row_cases(upto):
        lhs = N_row((a, b)).get((), 0)
        rhs = len(partitions_of(a, b)) - len(partitions_of(a + 1, b - 1))
        if lhs != rhs:
            bad += 1
            print(f"    MISMATCH ({a},{b}): {lhs} vs {rhs}")
    print(f"    {bad} mismatches")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--upto", type=int, default=11,
                        help="largest a+b to test (default 11)")
    args = parser.parse_args()
    cases = two_row_cases(args.upto)

    check_rules(cases)
    check_reverse_rules(cases)
    report_collisions(3, 3)
    check_matching(cases)
    check_antidiagonals(args.upto)


if __name__ == "__main__":
    main()

"""Multiset primitives and the lexicographic order used throughout.

A *multiset* is represented as a sorted-ascending tuple of positive ints,
e.g. {{1,1,2}} -> (1, 1, 2).  With that representation the lexicographic
order on multisets (Orellana-Zabrocki convention: 1<2, {1,2}<{2},
{1,1}<{2}, {1,1,2}<{1,2}) is exactly Python's native tuple comparison:

    (1,)     < (1, 1)      # "1"   < "11"
    (1, 1)   < (2,)        # "11"  < "2"
    (1, 2)   < (2,)        # "12"  < "2"
    (1, 1, 2)< (1, 2)      # "112" < "12"

so cells never need a custom comparator; `<` on the tuples is the order.
"""

from collections import Counter
from functools import lru_cache
from itertools import product


def as_multiset(elements):
    """Normalize an iterable of labels into the canonical sorted-tuple form."""
    return tuple(sorted(elements))


def _counter_to_tuple(counter):
    out = []
    for value in sorted(counter):
        out.extend([value] * counter[value])
    return tuple(out)


@lru_cache(maxsize=None)
def sub_multisets(multiset):
    """Every sub-multiset of `multiset` (including () and the whole), each as a
    sorted tuple, as a tuple of tuples.  Iterates multiplicity 0..k independently
    per label.

    Cached and materialized rather than a generator: the tableau backtracking
    asks for the sub-multisets of the same `remaining` content over and over
    across sibling branches, and rebuilding the list at every node was the
    dominant cost of the whole enumeration.
    """
    counts = Counter(multiset)
    labels = sorted(counts)
    out = []
    for choice in product(*[range(counts[v] + 1) for v in labels]):
        piece = []
        for value, take in zip(labels, choice):
            piece.extend([value] * take)
        out.append(tuple(piece))
    return tuple(out)


def multiset_difference(a, b):
    """Multiset a minus multiset b (b assumed contained in a for our uses)."""
    counts = Counter(a)
    counts.subtract(Counter(b))
    return _counter_to_tuple({v: n for v, n in counts.items() if n > 0})


@lru_cache(maxsize=None)
def multiset_partitions(multiset):
    """All multiset partitions of `multiset`, each as a sorted tuple of cells
    (a cell is itself a sorted tuple), with no duplicates.

    Method: remove one copy of the smallest label, recurse on the rest, then for
    each partition of the rest either drop the removed label into one of the
    existing blocks or start a new singleton block.  Canonicalizing every
    partition as a sorted tuple of cells and collecting into a set makes the
    indistinguishable-copy duplicates coincide instead of multiplying.
    """
    multiset = tuple(sorted(multiset))
    if not multiset:
        return [()]
    smallest = multiset[0]
    rest = multiset[1:]
    out = set()
    for partition in multiset_partitions(rest):
        # new singleton block
        out.add(tuple(sorted(partition + ((smallest,),))))
        # insert into each existing block
        for i in range(len(partition)):
            grown = tuple(sorted(partition[i] + (smallest,)))
            new_partition = partition[:i] + (grown,) + partition[i + 1:]
            out.add(tuple(sorted(new_partition)))
    return sorted(out)


def integer_partitions(n):
    """All integer partitions of n as non-increasing tuples (() for n == 0)."""
    if n == 0:
        return [()]
    out = []

    def rec(remaining, largest, current):
        if remaining == 0:
            out.append(tuple(current))
            return
        for part in range(min(remaining, largest), 0, -1):
            rec(remaining - part, part, current + [part])

    rec(n, n, [])
    return out


def integer_partitions_upto(n):
    """All integer partitions of sizes 0, 1, ..., n."""
    result = []
    for size in range(n + 1):
        result.extend(integer_partitions(size))
    return result

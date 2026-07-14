"""The multiplicity-profile decomposition of M_{(a,b),()}, and the injection search.

Reproduces (without Sage) the tables the b=2 / b=3 analysis was built on:

  * the C(a, lambda) contribution table for each lambda |- b,
  * M(a,b) vs M(a+1,b-1) and their difference N(a,b),
  * the profile collapse -- which partitions of b share a profile,
  * `assignment_search`: profile-level injection schemes ranked by the smallest
    a0 from which the domination inequalities hold persistently.

Usage:
    python scripts/profile_report.py --b 3 [--a-lo 3] [--a-hi 12] [--search]
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from oztab.multisets import integer_partitions  # noqa: E402
from oztab.profiles import (assignment_search, contribution,  # noqa: E402
                            M_empty_via_profiles, N_empty_via_profiles,
                            profile, profile_counts)


def part_str(lam):
    return "[" + ",".join(str(p) for p in lam) + "]" if lam else "[]"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--b", type=int, default=3)
    parser.add_argument("--a-lo", type=int, default=None)
    parser.add_argument("--a-hi", type=int, default=12)
    parser.add_argument("--search", action="store_true",
                        help="also run the injection assignment search")
    args = parser.parse_args()

    b = args.b
    a_lo = args.a_lo if args.a_lo is not None else b
    lams = [tuple(lam) for lam in integer_partitions(b)]

    print(f"=== profiles of the partitions of b = {b} ===")
    for pi, count in sorted(profile_counts(b).items()):
        sharers = [part_str(lam) for lam in lams if profile(lam) == pi]
        note = "  <-- shared" if count > 1 else ""
        print(f"  profile {str(pi):10} count {count}   {', '.join(sharers)}{note}")
    print(f"  p({b}) = {len(lams)} partitions collapse to "
          f"{len(profile_counts(b))} distinct profiles")

    print(f"\n=== contribution table, b = {b} ===")
    header = f"{'a':>4}" + "".join(f"{'C(a,' + part_str(l) + ')':>14}" for l in lams)
    header += f"{'M(a,b)':>10}{'M(a+1,b-1)':>12}{'N(a,b)':>9}"
    print(header)
    print("-" * len(header))
    for a in range(a_lo, args.a_hi + 1):
        line = f"{a:>4}"
        for lam in lams:
            line += f"{contribution(a, lam):>14}"
        line += f"{M_empty_via_profiles(a, b):>10}"
        line += f"{M_empty_via_profiles(a + 1, b - 1):>12}"
        line += f"{N_empty_via_profiles(a, b):>9}"
        print(line)

    if args.search:
        print(f"\n=== injection assignment search, b = {b} ===")
        print("each source lambda |- (b-1) must be absorbed into a disjoint set")
        print("of targets lambda' |- b; a0 = smallest a from which it holds.\n")
        for a0, scheme, leftover in assignment_search(b):
            print(f"  threshold a0 = {a0}")
            for src, targets in scheme:
                tgt = ", ".join(part_str(t) for t in targets)
                print(f"    {part_str(src):12} -> [{tgt}]")
            print(f"    leftover: [{', '.join(part_str(t) for t in leftover)}]")
            print()


if __name__ == "__main__":
    main()

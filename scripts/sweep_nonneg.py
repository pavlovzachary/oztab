"""Nonnegativity sweep of the Schur -> s~ coefficients N_{lambda,mu}.

The hypothesis under test: N_{lambda,mu} >= 0 for all lambda, mu.  Nothing in the
theory forces it -- N is a signed (Jacobi-Trudi) combination of the nonnegative
tableau counts M -- so a single negative entry would be the headline.

Usage:
    python scripts/sweep_nonneg.py --max-n 8 [--two-row] [--out data/]

Writes data/N_coefficients.csv with one row per nonzero coefficient, and prints
a per-n summary plus every negative entry found.  Cost grows ~7x per n; n <= 8
is seconds, n = 9 is minutes.
"""

import argparse
import csv
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from oztab import N_row, negative_terms  # noqa: E402
from oztab.multisets import integer_partitions  # noqa: E402


def part_str(lam):
    return "".join(str(p) for p in lam) if lam else "()"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=8,
                        help="sweep all |lambda| <= max_n (default 8)")
    parser.add_argument("--two-row", action="store_true",
                        help="restrict to two-row lambda (the research case)")
    parser.add_argument("--out", default="data",
                        help="output directory (default: data/)")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    path = os.path.join(args.out, "N_coefficients.csv")

    negatives = []
    rows = 0
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["n", "lambda", "mu", "N"])
        for n in range(args.max_n + 1):
            lams = [lam for lam in integer_partitions(n)
                    if not args.two_row or len(lam) <= 2]
            start = time.time()
            n_neg = 0
            for lam in lams:
                row = N_row(lam)
                for mu, coeff in sorted(row.items(), key=lambda kv: (sum(kv[0]), kv[0])):
                    writer.writerow([n, part_str(lam), part_str(mu), coeff])
                    rows += 1
                bad = negative_terms(row)
                for mu, coeff in bad.items():
                    negatives.append((lam, mu, coeff))
                    n_neg += 1
            flag = f"  <-- {n_neg} NEGATIVE" if n_neg else ""
            print(f"n={n:>2}: {len(lams):>3} partitions, "
                  f"{time.time() - start:6.2f}s{flag}", flush=True)

    print(f"\nwrote {rows} nonzero coefficients to {path}")
    if negatives:
        print(f"\n*** {len(negatives)} NEGATIVE COEFFICIENT(S) -- hypothesis violated ***")
        for lam, mu, coeff in negatives:
            print(f"    N_{{{part_str(lam)}, {part_str(mu)}}} = {coeff}")
        return 1
    print(f"\nNo negative coefficients for any lambda with |lambda| <= {args.max_n}.")
    print("Nonnegativity hypothesis holds on this range.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

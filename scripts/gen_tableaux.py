"""Print every Orellana-Zabrocki tableau for the requested partitions lambda.

The full set of tableaux for a lambda is the union, over every body shape mu
with |mu| <= |lambda|, of generate_tableaux(content(lambda), mu) -- there is no
single "all tableaux of lambda" call, since generate_tableaux fixes mu.  Each
tableau contributes one s~_{shape(body)} to the expansion of h_lambda.

Shape specs.  Each argument is a comma-separated partition; the letter `n`
stands for a free parameter swept over 1..UPTO (default 5):

    python gen_tableaux.py                  # n,1  -- the original (n,1) family
    python gen_tableaux.py 4,2 3,1,1        # just those two partitions
    python gen_tableaux.py n,2 --upto 4     # (2,2), (3,2), (4,2)
    python gen_tableaux.py n,n --upto 4     # (1,1), (2,2), (3,3), (4,4)
    python gen_tableaux.py all --upto 5     # every partition of size <= 5

Substitutions that are not weakly decreasing (n = 1 in `n,2`) are dropped.
Counts blow up fast -- |lambda| = 8 is already tens of thousands of tableaux --
so `--count-only` reports the totals without building anything.

Jacobi-Trudi pairs (`--pair`).  For a two-row lambda = (a, b),

    s_(a,b) = h_(a,b) - h_(a+1,b-1),      so   N_{(a,b),mu} = M_{(a,b),mu} - M_{(a+1,b-1),mu}

and the Schur coefficient is a difference of two tableau counts *shape by shape*.
`--pair` prints both sides of that difference interleaved by mu, which is the
view you want when hunting for the injection that proves nonnegativity:

    python gen_tableaux.py n,2 --upto 4 --pair
    python gen_tableaux.py 3,3 --pair

Both sides have the same |lambda|, hence the same blank count m, so the two
families of tableaux are directly comparable.  b = 1 is the already-solved case
(the partner is the single row (a+1)); b >= 2 is open.
"""

import argparse

from oztab import N_row, count_tableaux, generate_tableaux
from oztab.multisets import integer_partitions, integer_partitions_upto


def content_of_lambda(lam):
    """Content of h_lambda: {{1^lambda_1, 2^lambda_2, ...}} as a sorted tuple."""
    return tuple(sorted(
        label for label, mult in enumerate(lam, start=1) for _ in range(mult)
    ))


def expand_specs(specs, upto):
    """Turn shape specs into a de-duplicated list of partitions, in size order.

    A spec is `all` (every partition of size 1..upto) or a comma-separated list
    of parts in which the literal `n` is swept over 1..upto.
    """
    out = []
    for spec in specs:
        if spec == "all":
            for size in range(1, upto + 1):
                out.extend(integer_partitions(size))
            continue
        pieces = [p for p in spec.replace(" ", ",").split(",") if p]
        if not pieces:
            continue
        values = range(1, upto + 1) if "n" in pieces else [None]
        for n in values:
            lam = tuple(n if p == "n" else int(p) for p in pieces)
            if all(p > 0 for p in lam) and all(a >= b for a, b in zip(lam, lam[1:])):
                out.append(lam)
    seen, unique = set(), []
    for lam in out:
        if lam not in seen:
            seen.add(lam)
            unique.append(lam)
    return sorted(unique, key=lambda l: (sum(l), l))


def jt_partner(lam):
    """The subtracted term of s_(a,b) = h_(a,b) - h_(a+1,b-1).

    Returns None unless lam has exactly two rows.  b = 1 gives the single row
    (a+1), since h_0 = 1.
    """
    if len(lam) != 2:
        return None
    a, b = lam
    return (a + 1,) if b == 1 else (a + 1, b - 1)


def tableaux_for_lambda(lam):
    """All OZ tableaux for partition lam, grouped by body shape mu."""
    content = content_of_lambda(lam)
    by_mu = {}
    for mu in integer_partitions_upto(len(content)):
        ts = generate_tableaux(content, mu)
        if ts:
            by_mu[mu] = ts
    return content, by_mu


def total_for_lambda(lam):
    """Number of tableaux for lam, without building any (uses count_tableaux)."""
    content = content_of_lambda(lam)
    return sum(count_tableaux(content, mu)
               for mu in integer_partitions_upto(len(content)))


def _part_str(mu):
    return "(" + ",".join(str(p) for p in mu) + ")" if mu else "()"


def print_single(lam, faithful):
    content, by_mu = tableaux_for_lambda(lam)
    total = sum(len(v) for v in by_mu.values())
    print(f"\n{'=' * 66}")
    print(f"lambda = {lam}   content = {content}   ({total} tableaux total)")
    print(f"{'=' * 66}")
    for mu in sorted(by_mu, key=lambda m: (sum(m), m)):
        ts = by_mu[mu]
        print(f"\n  body shape mu = {_part_str(mu)}  ->  s~_{_part_str(mu)}   "
              f"({len(ts)} tableau{'x' if len(ts) != 1 else ''})")
        for t in ts:
            print("      " + t.render(faithful=faithful).replace("\n", "\n      "))


def print_pair(lam, partner, faithful):
    """The two sides of s_lam = h_lam - h_partner, interleaved by body shape."""
    left_content, left = tableaux_for_lambda(lam)
    right_content, right = tableaux_for_lambda(partner)
    lam_s, par_s = _part_str(lam), _part_str(partner)

    print(f"\n{'=' * 66}")
    print(f"s_{lam_s} = h_{lam_s} - h_{par_s}      (m = {len(left_content)} blanks on both sides)")
    print(f"  h_{lam_s}: content {left_content}   "
          f"{sum(len(v) for v in left.values())} tableaux")
    print(f"  h_{par_s}: content {right_content}   "
          f"{sum(len(v) for v in right.values())} tableaux")
    print(f"{'=' * 66}")

    reference = N_row(lam)
    running = {}
    for mu in sorted(set(left) | set(right), key=lambda m: (sum(m), m)):
        lts, rts = left.get(mu, []), right.get(mu, [])
        diff = len(lts) - len(rts)
        if diff:
            running[mu] = diff
        print(f"\n  mu = {_part_str(mu)}  ->  s~_{_part_str(mu)}      "
              f"{len(lts)} - {len(rts)} = {diff:+d}")
        for label, ts in ((f"h_{lam_s}", lts), (f"h_{par_s}", rts)):
            print(f"    {label}:" + ("" if ts else "  (none)"))
            for t in ts:
                print("      " + t.render(faithful=faithful).replace("\n", "\n      "))

    status = "OK" if running == reference else "MISMATCH"
    print(f"\n  band-by-band differences vs schur.N_row({lam}): {status}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("specs", nargs="*",
                        help="shape specs, e.g. 'n,1', '4,2', 'all' (default: n,1)")
    parser.add_argument("--upto", type=int, default=5,
                        help="range of the free parameter n (default 5)")
    parser.add_argument("--faithful", action="store_true",
                        help="render the implied blank cells of row 1")
    parser.add_argument("--count-only", action="store_true",
                        help="print totals per lambda, do not build the tableaux")
    parser.add_argument("--pair", action="store_true",
                        help="also show the Jacobi-Trudi partner (a+1,b-1) of each "
                             "two-row lambda, interleaved by body shape")
    args = parser.parse_args()

    lambdas = expand_specs(args.specs or ["n,1"], args.upto)
    if not lambdas:
        raise SystemExit("no partitions matched those specs")

    if args.pair:
        skipped = [lam for lam in lambdas if jt_partner(lam) is None]
        if skipped:
            print(f"--pair needs a two-row lambda; skipping {skipped}")
        lambdas = [lam for lam in lambdas if jt_partner(lam) is not None]
        if not lambdas:
            raise SystemExit("no two-row partitions to pair")

    if args.count_only:
        for lam in lambdas:
            if args.pair:
                partner = jt_partner(lam)
                print(f"s_{_part_str(lam)} = h_{_part_str(lam)} - h_{_part_str(partner)}   "
                      f"{total_for_lambda(lam)} - {total_for_lambda(partner)} tableaux")
            else:
                print(f"lambda = {lam}   |lambda| = {sum(lam)}   "
                      f"{total_for_lambda(lam)} tableaux")
        return

    for lam in lambdas:
        if args.pair:
            print_pair(lam, jt_partner(lam), args.faithful)
        else:
            print_single(lam, args.faithful)


if __name__ == "__main__":
    main()

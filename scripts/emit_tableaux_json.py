"""Enumerate OZ tableaux for the requested lambdas and splice them into the viewer.

Running this (1) enumerates every tableau of every requested partition lambda,
and (2) rewrites the embedded JSON data block of oz_tableaux.html *in place*, so
opening (or just refreshing) that file in a browser shows the current data.

Shape specs are the same as gen_tableaux.py: a comma-separated partition, with
the letter `n` swept over 1..UPTO, or the word `all`.

    python emit_tableaux_json.py                     # n,1  up to (5,1)
    python emit_tableaux_json.py n,1 --upto 6
    python emit_tableaux_json.py n,2 3,2,1 --upto 5  # mixed families
    python emit_tableaux_json.py all --upto 6        # every partition, |lambda| <= 6

With `--pair`, each two-row lambda = (a,b) is emitted together with its
Jacobi-Trudi partner (a+1,b-1) as a single tab, and the viewer renders the two
side by side within each body-shape band -- so the band difference on screen IS
N_{(a,b),mu} (see gen_tableaux.jt_partner):

    python emit_tableaux_json.py n,2 --upto 5 --pair

The JSON is inlined into the HTML, so the file grows with the data.  Past
roughly 20k tableaux the page gets unpleasant to load; `--limit` (default 20000)
refuses to write past that, and `--limit 0` disables the check.
"""

import argparse
import json
import re
from pathlib import Path

from gen_tableaux import (content_of_lambda, expand_specs, jt_partner,
                          total_for_lambda)
from oztab import generate_tableaux
from oztab.multisets import integer_partitions_upto

HTML_PATH = Path(__file__).with_name("oz_tableaux.html")
DATA_RE = re.compile(
    r'(<script type="application/json" id="data">)(.*?)(</script>)', re.S
)


def part_str(lam):
    return "(" + ",".join(str(p) for p in lam) + ")" if lam else "∅"


def side(lam):
    """(metadata, {mu: [tableau dicts]}) for one partition."""
    content = content_of_lambda(lam)
    m = len(content)
    by_mu = {}
    for mu in integer_partitions_upto(m):
        ts = generate_tableaux(content, mu)
        if not ts:
            continue
        by_mu[mu] = [
            {
                "first_row": [list(cell) for cell in t.first_row],
                "body": [[list(cell) for cell in row] for row in t.body],
            }
            for t in ts
        ]
    meta = {
        "lambda": list(lam),
        "label": part_str(lam),
        "content": list(content),
        "m": m,
        "total": sum(len(v) for v in by_mu.values()),
    }
    return meta, by_mu


def build(entries):
    """The {"entries": [...]} payload the viewer consumes.

    Each entry is one tab and holds one or two `sides`.  A two-sided entry is a
    Jacobi-Trudi pair: side 0 is added, side 1 subtracted, so a band's `diff` is
    the Schur coefficient N_{lambda,mu}.
    """
    data = {"entries": []}
    for lambdas in entries:
        metas, bodies = [], []
        for lam in lambdas:
            meta, by_mu = side(lam)
            metas.append(meta)
            bodies.append(by_mu)

        shapes = set()
        for by_mu in bodies:
            shapes |= set(by_mu)
        groups = []
        for mu in sorted(shapes, key=lambda m: (sum(m), m)):
            per_side = [by_mu.get(mu, []) for by_mu in bodies]
            counts = [len(ts) for ts in per_side]
            groups.append({
                "mu": list(mu),
                "counts": counts,
                "diff": counts[0] - (counts[1] if len(counts) > 1 else 0),
                "tableaux": per_side,
            })

        data["entries"].append({
            "paired": len(metas) == 2,
            "label": " − ".join(m["label"] for m in metas),
            "sides": metas,
            "total": sum(m["total"] for m in metas),
            "groups": groups,
        })
    return data


def splice(data):
    """Replace the <script id="data"> block of the HTML with `data`. Returns
    the JSON payload written."""
    payload = json.dumps(data)
    html = HTML_PATH.read_text()
    if not DATA_RE.search(html):
        raise SystemExit(f"data block not found in {HTML_PATH}")
    # function replacement -> payload is inserted literally (no backref parsing)
    new_html = DATA_RE.sub(
        lambda mo: mo.group(1) + "\n" + payload + "\n" + mo.group(3),
        html, count=1,
    )
    HTML_PATH.write_text(new_html)
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("specs", nargs="*",
                        help="shape specs, e.g. 'n,1', '4,2', 'all' (default: n,1)")
    parser.add_argument("--upto", type=int, default=5,
                        help="range of the free parameter n (default 5)")
    parser.add_argument("--limit", type=int, default=20000,
                        help="refuse to write more than this many tableaux (0 = no limit)")
    parser.add_argument("--pair", action="store_true",
                        help="pair each two-row lambda with its Jacobi-Trudi "
                             "partner (a+1,b-1) in one tab, shown side by side")
    args = parser.parse_args()

    lambdas = expand_specs(args.specs or ["n,1"], args.upto)
    if not lambdas:
        raise SystemExit("no partitions matched those specs")

    if args.pair:
        skipped = [lam for lam in lambdas if jt_partner(lam) is None]
        if skipped:
            print(f"--pair needs a two-row lambda; skipping {skipped}")
        entries = [(lam, jt_partner(lam))
                   for lam in lambdas if jt_partner(lam) is not None]
        if not entries:
            raise SystemExit("no two-row partitions to pair")
    else:
        entries = [(lam,) for lam in lambdas]

    # count first: counting is cheap, building is not
    counts = {lam: total_for_lambda(lam)
              for entry in entries for lam in entry}
    grand = sum(counts.values())
    if args.limit and grand > args.limit:
        listing = ", ".join(f"{lam}: {c:,}" for lam, c in counts.items())
        raise SystemExit(
            f"{grand:,} tableaux exceeds --limit {args.limit:,} ({listing}).\n"
            f"Narrow the specs, lower --upto, or pass --limit 0 to force it."
        )

    data = build(entries)
    payload = splice(data)
    totals = ", ".join(f"{e['label']}: {e['total']}" for e in data["entries"])
    print(f"Updated {HTML_PATH.name} - {len(entries)} tab(s)"
          f"{' (Jacobi-Trudi pairs)' if args.pair else ''}")
    print(f"  totals: {totals}")
    print(f"  data block: {len(payload):,} bytes")


if __name__ == "__main__":
    main()

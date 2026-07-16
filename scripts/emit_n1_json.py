"""Generate all OZ tableaux for lambda=(n,1) and splice them into the viewer.

Running this (1) enumerates every tableau for lambda=(1,1)..(NMAX,1), and
(2) rewrites the embedded JSON data block of n1_tableaux.html *in place*, so
opening (or just refreshing) that file in a browser shows the current data.

Usage:  python emit_n1_json.py [NMAX]      # default NMAX = 5
"""

import json
import re
import sys
from pathlib import Path

from oztab import generate_tableaux
from oztab.multisets import integer_partitions_upto

HTML_PATH = Path(__file__).with_name("n1_tableaux.html")
DATA_RE = re.compile(
    r'(<script type="application/json" id="data">)(.*?)(</script>)', re.S
)


def content_of_lambda(lam):
    return tuple(sorted(
        label for label, mult in enumerate(lam, start=1) for _ in range(mult)
    ))


def build(nmax):
    """The {"lambdas": [...]} payload the viewer consumes."""
    data = {"lambdas": []}
    for n in range(1, nmax + 1):
        lam = (n, 1)
        content = content_of_lambda(lam)
        m = len(content)
        groups = []
        total = 0
        for mu in integer_partitions_upto(m):
            ts = generate_tableaux(content, mu)
            if not ts:
                continue
            total += len(ts)
            groups.append({
                "mu": list(mu),
                "tableaux": [
                    {
                        "first_row": [list(cell) for cell in t.first_row],
                        "body": [[list(cell) for cell in row] for row in t.body],
                    }
                    for t in ts
                ],
            })
        # order groups by (|mu|, mu)
        groups.sort(key=lambda g: (sum(g["mu"]), g["mu"]))
        data["lambdas"].append({
            "lambda": list(lam),
            "content": list(content),
            "m": m,
            "total": total,
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
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    data = build(nmax)
    payload = splice(data)
    totals = ", ".join(
        f"({L['lambda'][0]},1): {L['total']}" for L in data["lambdas"]
    )
    print(f"Updated {HTML_PATH.name} - lambda=(1,1)..({nmax},1)")
    print(f"  totals: {totals}")
    print(f"  data block: {len(payload):,} bytes")


if __name__ == "__main__":
    main()

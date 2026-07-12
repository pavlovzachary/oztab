"""Correctness check: run under SageMath (`sage -python scripts/sage_check.py`
or inside a Sage session).

Compares our hand-rolled M-table against Sage's built-in irreducible symmetric
group character basis (the Orellana-Zabrocki s~ basis).  For every lambda with
|lambda| <= N it asserts that our {mu: count} equals the coefficients of
st(h[lambda]).  This is the independent oracle: our enumeration on one side,
Sage's algebra on the other.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from oztab import M_row, integer_partitions  # noqa: E402

try:
    from sage.all import SymmetricFunctions, QQ, Partition  # noqa: E402
except ImportError:
    sys.exit("This script must be run under SageMath (sage -python scripts/sage_check.py).")

N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

Sym = SymmetricFunctions(QQ)
h = Sym.h()
# The irreducible character basis; accessor name has varied across Sage versions.
try:
    st = Sym.st()
except AttributeError:
    st = Sym.irreducible_symmetric_group_character()


def sage_row(lam):
    """{mu (tuple): coeff} for st(h[lam])."""
    expansion = st(h[list(lam)] if lam else h.one())
    out = {}
    for part, coeff in expansion.monomial_coefficients().items():
        out[tuple(part)] = int(coeff)
    return {mu: c for mu, c in out.items() if c != 0}


failures = 0
checked = 0
for size in range(N + 1):
    for lam in integer_partitions(size):
        ours = M_row(lam)
        theirs = sage_row(lam)
        checked += 1
        if ours != theirs:
            failures += 1
            print(f"MISMATCH  h_{lam or '0'}")
            print(f"   ours  : {ours}")
            print(f"   sage  : {theirs}")

print(f"\nchecked {checked} partitions up to size {N}: "
      f"{'ALL MATCH' if failures == 0 else str(failures) + ' MISMATCH(ES)'}")
sys.exit(1 if failures else 0)

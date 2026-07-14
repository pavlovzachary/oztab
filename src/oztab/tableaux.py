"""Semistandard multiset tableaux for the Orellana-Zabrocki character basis.

Object.  To expand h_lambda in the irreducible character basis s~_mu, we count
column-strict tableaux whose first row starts with m blank cells (m fixed,
m >= |mu|; we take m = |lambda|), whose remaining cells hold nonempty multisets,
with total content {{1^lambda_1, 2^lambda_2, ...}}.  Each such tableau T
contributes s~_{shape(T) with its first row deleted}.  (Orellana-Zabrocki,
"Symmetric group characters as symmetric functions", Prop. 2 / Remark 8.)

Reduction we exploit.  At m = |lambda| the body (rows 2, 3, ...) has first-row
length mu_1 <= |mu| <= |lambda| = m, so the whole body sits underneath blank
cells.  A blank is smaller than every multiset, so no column constraint ever
links the body to the first row's filled cells.  Hence a tableau splits into two
INDEPENDENT pieces that only share the content:

    first row (filled part)  =  a multiset partition of some sub-content
                                (displayed weakly increasing, left to right)
    body (rows 2..)          =  a shape-mu filling with nonempty multiset cells,
                                rows weakly increasing, columns strictly
                                increasing, in the lex order on multisets.

`MultisetTableau` stores exactly these two pieces.  The m blank cells are implied
(and shown only in the faithful rendering); shape_index() returns mu.
"""

from functools import lru_cache

from .multisets import multiset_difference, multiset_partitions, sub_multisets


def _cell_str(cell):
    return "{" + ",".join(str(x) for x in cell) + "}"


class MultisetTableau:
    """A single tableau: `first_row` and `body`, both tuples of multiset cells.

    first_row : tuple of cells (each a sorted tuple), weakly increasing.
    body      : tuple of rows, each row a tuple of cells; its shape is mu.
    m         : number of implied blank cells at the start of row 1.
    """

    __slots__ = ("first_row", "body", "m")

    def __init__(self, first_row, body, m):
        self.first_row = tuple(first_row)
        self.body = tuple(tuple(row) for row in body)
        self.m = m

    def shape_index(self):
        """The partition mu indexing s~_mu (the shape of the body)."""
        return tuple(len(row) for row in self.body)

    def content(self):
        """Total content multiset over every filled cell (sorted tuple)."""
        elems = []
        for cell in self.first_row:
            elems.extend(cell)
        for row in self.body:
            for cell in row:
                elems.extend(cell)
        return tuple(sorted(elems))

    # -- rendering -----------------------------------------------------------
    def render(self, faithful=False):
        """Text rendering.  Compact by default; faithful shows the m blanks so
        the skew placement (body under blanks, first-row cells to the right) is
        visible exactly as in the Orellana-Zabrocki picture."""
        tail = " ".join(_cell_str(c) for c in self.first_row)
        if not faithful:
            lines = ["row1: " + (tail if tail else "(empty)")]
            if self.body:
                body_lines = [" ".join(_cell_str(c) for c in row) for row in self.body]
                lines.append("body: " + body_lines[0])
                lines.extend("      " + bl for bl in body_lines[1:])
            else:
                lines.append("body: (empty)")
            return "\n".join(lines)
        # faithful: m blank markers, then the tail; body rows aligned at column 0
        blanks = " ".join(["."] * self.m)
        first = (blanks + (" " + tail if tail else "")).rstrip()
        lines = [first if first else "."]
        for row in self.body:
            lines.append(" ".join(_cell_str(c) for c in row))
        return "\n".join(lines)

    def __str__(self):
        return self.render(faithful=False)

    def __repr__(self):
        return (f"MultisetTableau(shape={self.shape_index()}, "
                f"first_row={self.first_row}, body={self.body})")


def column_strict_fillings(shape, content):
    """All fillings of Young diagram `shape` (a partition tuple) with nonempty
    multiset cells using exactly `content`, rows weakly increasing and columns
    strictly increasing in the lex order.  Cells are filled in reading order
    (left to right, top to bottom) by backtracking."""
    shape = tuple(shape)
    coords = [(r, c) for r in range(len(shape)) for c in range(shape[r])]
    total = len(coords)
    if total == 0:
        return [()] if len(content) == 0 else []
    results = []
    grid = {}

    def backtrack(index, remaining):
        if index == total:
            if not remaining:
                rows = tuple(
                    tuple(grid[(r, c)] for c in range(shape[r]))
                    for r in range(len(shape))
                )
                results.append(rows)
            return
        r, c = coords[index]
        left = grid.get((r, c - 1))       # row: cell must be >= left  (weak)
        top = grid.get((r - 1, c))        # column: cell must be > top  (strict)
        for cell in sub_multisets(remaining):
            if not cell:
                continue
            if left is not None and cell < left:
                continue
            if top is not None and not (cell > top):
                continue
            grid[(r, c)] = cell
            backtrack(index + 1, multiset_difference(remaining, cell))
            del grid[(r, c)]

    backtrack(0, tuple(sorted(content)))
    return results


@lru_cache(maxsize=None)
def count_column_strict_fillings(shape, content):
    """Number of fillings counted by `column_strict_fillings`, without building
    them.  Same backtracking, but it accumulates a count, so the tableaux never
    have to be materialized -- the count is what the M-table needs, and there are
    ~10x more tableaux at each n.  Cached, since the same (shape, sub-content)
    pair recurs across many mu and many lambda.

    `tests/test_counting.py` pins this against len(column_strict_fillings(...)).
    """
    shape = tuple(shape)
    content = tuple(sorted(content))
    coords = [(r, c) for r in range(len(shape)) for c in range(shape[r])]
    total = len(coords)
    if total == 0:
        return 1 if len(content) == 0 else 0
    grid = {}

    def backtrack(index, remaining):
        if index == total:
            return 1 if not remaining else 0
        cells_left = total - index
        if len(remaining) < cells_left:
            return 0                    # every cell is nonempty: not enough left
        r, c = coords[index]
        left = grid.get((r, c - 1))
        top = grid.get((r - 1, c))
        found = 0
        for cell in sub_multisets(remaining):
            if not cell:
                continue
            if len(remaining) - len(cell) < cells_left - 1:
                continue                # would starve the remaining cells
            if left is not None and cell < left:
                continue
            if top is not None and not (cell > top):
                continue
            grid[(r, c)] = cell
            found += backtrack(index + 1, multiset_difference(remaining, cell))
            del grid[(r, c)]
        return found

    return backtrack(0, content)


def count_tableaux(content, mu):
    """Number of Orellana-Zabrocki tableaux of `content` with body shape `mu`,
    i.e. len(generate_tableaux(content, mu)) -- computed without building them.

    Uses the same body/first-row split: sum over the body's sub-content of
    (# shape-mu fillings) * (# multiset partitions of the leftover).
    """
    content = tuple(sorted(content))
    mu = tuple(mu)
    cells_needed = sum(mu)
    total = 0
    for body_content in sub_multisets(content):
        if len(body_content) < cells_needed:
            continue
        bodies = count_column_strict_fillings(mu, body_content)
        if not bodies:
            continue
        leftover = multiset_difference(content, body_content)
        total += bodies * len(multiset_partitions(leftover))
    return total


def generate_tableaux(content, mu, m=None):
    """All Orellana-Zabrocki tableaux of the given `content` (a multiset, the
    content of h_lambda) whose body shape is `mu`.  `m` defaults to |content|.

    Enumeration: choose the body's sub-content, fill the shape-mu body, and
    independently choose a multiset partition of the leftover for the first row.
    """
    content = tuple(sorted(content))
    mu = tuple(mu)
    if m is None:
        m = len(content)
    cells_needed = sum(mu)
    tableaux = []
    seen_bodies = set()
    for body_content in sub_multisets(content):
        if len(body_content) < cells_needed:
            continue
        if body_content in seen_bodies:
            continue
        seen_bodies.add(body_content)
        bodies = column_strict_fillings(mu, body_content)
        if not bodies:
            continue
        tails = multiset_partitions(multiset_difference(content, body_content))
        for body in bodies:
            for tail in tails:
                tableaux.append(MultisetTableau(tail, body, m))
    return tableaux


def is_valid(tableau):
    """Independently re-check a tableau's constraints (used in tests / to trust
    eyeballing): first row weakly increasing, body rows weak, body columns
    strict, all cells nonempty."""
    fr = tableau.first_row
    if any(len(c) == 0 for c in fr):
        return False
    if any(fr[i] > fr[i + 1] for i in range(len(fr) - 1)):
        return False
    body = tableau.body
    for row in body:
        if any(len(c) == 0 for c in row):
            return False
        if any(row[i] > row[i + 1] for i in range(len(row) - 1)):
            return False
    for r in range(len(body) - 1):
        for c in range(len(body[r + 1])):     # lower row is shorter or equal
            if not (body[r + 1][c] > body[r][c]):
                return False
    return True

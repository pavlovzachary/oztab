# oztab

Tools for the Orellana-Zabrocki **irreducible character basis** `s~_mu` of
symmetric functions, built to generate data and test nonnegativity hypotheses
for the Schur -> `s~` change of basis.

## The object

To expand `h_lambda` in the `s~` basis we count column-strict tableaux with
`m` blank cells at the start of the first row (`m = |lambda|` here), the
remaining cells holding nonempty multisets, total content
`{{1^lambda_1, 2^lambda_2, ...}}`. Each tableau `T` contributes `s~_mu` where
`mu` is `shape(T)` with its first row removed. (OZ, *Symmetric group characters
as symmetric functions*, Prop. 2 / Remark 8.)

At `m = |lambda|` the body sits entirely under blank cells, so a tableau splits
into two independent pieces sharing only the content: a **multiset partition**
(the first row's filled cells) and a **shape-`mu` multiset filling** (the body,
rows weakly increasing, columns strictly increasing in the lex order on
multisets: `1<2`, `{1,2}<{2}`, `{1,1,2}<{1,2}`). The enumeration is hand-rolled
on exactly that split.

## The two sides

The tableau model only counts the **`h`-side** (`M`). The research object is the
**Schur side** (`N`), reached by Jacobi-Trudi:

    s_lambda = det(h_{lambda_i - i + j})   ->   N_row = sum of signed M_rows

`M` entries are counts, so nonnegative by construction. `N` entries are a
*signed* combination of them, so nothing forces `N >= 0` — that it seems to hold
anyway is the **nonnegativity hypothesis**, and `scripts/sweep_nonneg.py` is the
hunt for a counterexample.

## Layout

    src/oztab/
      multisets.py     lex order, sub-multisets, multiset partitions
      tableaux.py      MultisetTableau + generator + validator + fast counter
      coefficients.py  M-row / M-table  (h-side; nonnegative counts)
      schur.py         N-row / N-table  (Schur side, via Jacobi-Trudi; signed)
      profiles.py      multiplicity-profile decomposition of the s~_() coefficient
    scripts/
      sweep_nonneg.py    sweep N for negative coefficients -> data/*.csv
      profile_report.py  contribution tables + the injection assignment search
      sage_check.py      compare M-table against Sage's s~ basis (run under Sage)
    notebooks/tour.ipynb   worked tour of all of the above
    tests/                 pytest, no Sage needed

## Use

    pip install -e .
    pytest                                   # 14 tests, no Sage needed

    from oztab import M_row, N_row, generate_tableaux, format_row
    print(format_row((2,1), M_row((2,1))))   # h_21 = 4 s~_() + 7 s~_1 + ...
    N_row((3,2))                             # {(): 6, (1,): 5, ...}
    for T in generate_tableaux((1,1,2), (1,1)):
        print(T); print()                    # compact view
        print(T.render(faithful=True))       # with the m blank cells shown

    python scripts/sweep_nonneg.py --max-n 8
    python scripts/profile_report.py --b 3 --search

Start with `notebooks/tour.ipynb`.

## Cost

The enumeration is exponential, roughly **7x per `n`**: `n <= 8` is seconds,
`n = 9` is ~2.5 minutes, `n = 10` is the better part of an hour. `M_row` counts
tableaux without building them (`count_tableaux`), which is what makes `n = 9`
reachable at all — there are 1.3M tableaux already at `n = 8`.

## Correctness

Four independent checks, none of which need Sage:

- **The paper.** `tests/test_regression.py` freezes OZ Example 10 (`h_{2,1}`,
  20 tableaux) and Example 11 (`h_{1^4}`).
- **Counter vs generator.** `M_row` counts tableaux without building them, so
  the counts are no longer *structurally* tied to the objects.
  `tests/test_counting.py` restores that guarantee by asserting
  `count_tableaux == len(generate_tableaux(...))` for every `(lambda, mu)` with
  `|lambda| <= 6`. **If you touch the enumeration, this is the test that matters.**
- **Enumerator vs profile formula.** `profiles.cross_check` derives
  `M_{(a,b),()}` a completely different way (the closed-form contribution sum)
  and asserts it agrees. Also pins the solved `b = 2` case,
  `N_{(a,2),()} = C(a,(1,1)) - p(a+1)`.
- **Sage.** `scripts/sage_check.py` is the external oracle: our counts vs
  `st(h[lambda])`. Needs a Sage install.

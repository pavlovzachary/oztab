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

## Layout

    src/oztab/
      multisets.py     lex order, sub-multisets, multiset partitions
      tableaux.py      MultisetTableau + generator (hand-rolled) + validator
      coefficients.py  M-row / M-table counted from the generator
    scripts/sage_check.py   compare M-table against Sage's s~ basis (run under Sage)
    tests/test_regression.py  frozen counts from the OZ paper (pytest, no Sage)

## Use

    pip install -e .
    pytest                                   # pure-Python regression tests

    from oztab import M_row, generate_tableaux, format_row
    print(format_row((2,1), M_row((2,1))))   # h_21 = 4 s~_() + 7 s~_1 + ...
    for T in generate_tableaux((1,1,2), (1,1)):
        print(T); print()                    # compact view
        print(T.render(faithful=True))       # with the m blank cells shown

## Correctness

`tests/` freezes the paper's Example 10 (`h_{2,1}`, 20 tableaux) and Example 11
(`h_{1^4}`). `scripts/sage_check.py` is the independent oracle: it asserts our
counts equal `st(h[lambda])` in Sage for all `lambda` up to a chosen size.

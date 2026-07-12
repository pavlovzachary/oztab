"""oztab: multiset tableaux for the Orellana-Zabrocki irreducible character basis.

Primitive object is the tableau (see tableaux.MultisetTableau); coefficient
tables are counted from it (see coefficients.M_row / M_table).
"""

from .tableaux import MultisetTableau, generate_tableaux, is_valid
from .coefficients import content_of, M_row, M_table, format_row
from .multisets import multiset_partitions, integer_partitions

__all__ = [
    "MultisetTableau", "generate_tableaux", "is_valid",
    "content_of", "M_row", "M_table", "format_row",
    "multiset_partitions", "integer_partitions",
]

"""oztab: multiset tableaux for the Orellana-Zabrocki irreducible character basis.

Primitive object is the tableau (see tableaux.MultisetTableau); coefficient
tables are counted from it (see coefficients.M_row / M_table).  The Schur side
(schur.N_row) sits on top of the h-side via Jacobi-Trudi, and profiles.py holds
the multiplicity-profile decomposition of the s~_() coefficient.
"""

from .tableaux import (MultisetTableau, generate_tableaux, is_valid,
                       count_tableaux)
from .coefficients import content_of, M_row, M_table, format_row
from .multisets import multiset_partitions, integer_partitions
from .schur import (N_row, N_table, jacobi_trudi_terms, negative_terms,
                    format_N_row)
from .profiles import (profile, contribution, contribution_by_profile,
                       profile_counts, M_empty_via_profiles,
                       N_empty_via_profiles, assignment_search, cross_check)
from .n1 import (N_n1, N_n1_row, N_n1_series, principal_specialization,
                 outer_corners, hook_lengths, partition_series)

__all__ = [
    # tableaux
    "MultisetTableau", "generate_tableaux", "is_valid", "count_tableaux",
    # h-side coefficients
    "content_of", "M_row", "M_table", "format_row",
    # Schur side
    "N_row", "N_table", "jacobi_trudi_terms", "negative_terms", "format_N_row",
    # profile decomposition
    "profile", "contribution", "contribution_by_profile", "profile_counts",
    "M_empty_via_profiles", "N_empty_via_profiles", "assignment_search",
    "cross_check",
    # closed form for lambda = (n, 1)
    "N_n1", "N_n1_row", "N_n1_series", "principal_specialization",
    "outer_corners", "hook_lengths", "partition_series",
    # multisets
    "multiset_partitions", "integer_partitions",
]

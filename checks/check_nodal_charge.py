#!/usr/bin/env python3
"""Exact linear-algebra checks for the nodal charge derivation.

This script checks only the finite matrix statements used in Section 3
of the paper.  The mixed matrix records the geometric X_9
relation after an integral change of charge basis.  Its E2 moment-map
interpretation is derived in Section 2.3 of the paper and checked
algebraically in checks/check_e2_gauging.py.
"""

from functools import reduce
from math import gcd

from sympy import Matrix, ilcm


def primitive_integer_vector(vector: Matrix) -> tuple[int, ...]:
    """Return the primitive integer representative of a rational vector."""
    common_denom = reduce(ilcm, (int(entry.q) for entry in vector), 1)
    integers = [int(entry * common_denom) for entry in vector]
    common_factor = reduce(gcd, (abs(entry) for entry in integers if entry), 0)
    integers = [entry // common_factor for entry in integers]
    first_nonzero = next(entry for entry in integers if entry)
    if first_nonzero < 0:
        integers = [-entry for entry in integers]
    return tuple(integers)


def main() -> None:
    # X_9's two nodal curve classes are independent.  After choosing them as
    # a basis of their span, their exact charge matrix is the identity.
    q_node = Matrix([[1, 0], [0, 1]])
    assert q_node.rank() == 2
    assert q_node.nullspace() == []

    # The branch-restricted geometric relation is (-1,-1,1), up to an
    # overall sign.  The E2 derivation identifies beta as the Cartan moment
    # map with embedding vector (1,1), so the effective mixed matrix is:
    q_mixed_effective = Matrix([[1, 0, 1], [0, 1, 1]])
    mixed_kernel = q_mixed_effective.nullspace()
    assert q_mixed_effective.rank() == 2
    assert len(mixed_kernel) == 1
    relation = primitive_integer_vector(mixed_kernel[0])
    assert relation in {(1, 1, -1), (-1, -1, 1)}

    print("Q_node rank:", q_node.rank())
    print("ker(Q_node): 0")
    print("Q_mixed_effective rank:", q_mixed_effective.rank())
    print("ker(Q_mixed_effective): span", relation)
    print("Arithmetic check passed; see Section 2.3 of the paper for the physics.")


if __name__ == "__main__":
    main()

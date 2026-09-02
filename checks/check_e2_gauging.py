#!/usr/bin/env python3
"""Exact algebra and lattice checks for the reduced E2 gauging derivation.

The conceptual inputs (the E2 chiral ring, its flavor interpretation, and the
M2/divisor charge formula) come from the sources cited in Section 2.3 of
the paper.  This script checks the polynomial relations,
the X_9 lattice reduction, and representative real-moment-map lifts.
"""

from __future__ import annotations

from cmath import exp, phase
from functools import reduce
from math import gcd, isclose, sqrt

from sympy import Matrix, groebner, symbols
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


def primitive_integer_vector(vector: Matrix) -> tuple[int, ...]:
    """Return the primitive integer representative of an integral null vector."""
    integers = [int(entry) for entry in vector]
    common_factor = reduce(gcd, (abs(entry) for entry in integers if entry), 0)
    integers = [entry // common_factor for entry in integers]
    first_nonzero = next(entry for entry in integers if entry)
    if first_nonzero < 0:
        integers = [-entry for entry in integers]
    return tuple(integers)


def check_higgs_ring() -> None:
    """Check the E2 relations and their Cartan-invariant normal forms."""
    instanton, meson, anti_instanton, fat = symbols("I M Itilde P")
    ideal = [
        meson**2 - instanton * anti_instanton,
        fat**2,
        fat * instanton,
        fat * meson,
        fat * anti_instanton,
    ]
    basis = groebner(
        ideal,
        instanton,
        anti_instanton,
        fat,
        meson,
        order="lex",
    )
    assert all(basis.reduce(relation)[1] == 0 for relation in ideal)

    # A monomial I^a M^b Itilde^c has Cartan weight 2(a-c).  At weight
    # zero, a=c and (I Itilde)^a M^b reduces exactly to M^(2a+b).
    for a in range(8):
        for b in range(8):
            monomial = instanton**a * meson**b * anti_instanton**a
            assert basis.reduce(monomial - meson ** (2 * a + b))[1] == 0

    # Any nonconstant monomial multiplied by P vanishes, as does P^2.
    for generator in (instanton, meson, anti_instanton, fat):
        assert basis.reduce(fat * generator)[1] == 0

    print("Cartan-invariant E2 ring: C[M,P]/(P^2, M P)")
    print("Altmann base identification: beta = M, alpha = P")


def check_real_moment_map_lifts() -> None:
    """Check explicit lifts beta=ab with |a|=|b| on representative values."""
    for beta in (0j, 1 + 0j, -1 + 0j, 2 + 3j, -4j):
        radius = sqrt(abs(beta))
        a = radius * exp(1j * phase(beta)) if beta else 0j
        b = radius + 0j
        assert isclose(abs(a * b - beta), 0.0, abs_tol=1e-12)
        assert isclose(abs(a) ** 2 - abs(b) ** 2, 0.0, abs_tol=1e-12)
    print("Representative E2 lifts satisfy beta=ab and mu_R=0")


def check_x9_matrix() -> None:
    """Reduce the printed ambient X_9 matrix to its two effective constraints."""
    q_ambient = Matrix(
        [
            [-1, 0, -1],
            [1, 0, 1],
            [-1, 1, 0],
            [1, -1, 0],
            [0, -1, -1],
            [0, 1, 1],
            [0, 0, 0],
            [0, 0, 0],
        ]
    )
    q_effective = Matrix([[1, 0, 1], [0, 1, 1]])

    assert q_ambient.rank() == q_effective.rank() == 2
    assert Matrix.vstack(q_ambient, q_effective).rank() == 2
    smith = smith_normal_form(q_ambient, domain=ZZ)
    nonzero_smith_entries = [
        int(smith[i, i])
        for i in range(min(smith.rows, smith.cols))
        if smith[i, i] != 0
    ]
    assert [abs(entry) for entry in nonzero_smith_entries] == [1, 1]

    kernel = q_ambient.nullspace()
    assert len(kernel) == 1
    relation = primitive_integer_vector(kernel[0])
    assert relation == (1, 1, -1)

    print("Q_X9 rank: 2")
    print("Smith invariants:", tuple(abs(x) for x in nonzero_smith_entries))
    print("ker(Q_X9): span", relation)
    print("Generic reduced quotient dimension: 3 - 2 = 1 quaternionic")


def main() -> None:
    check_higgs_ring()
    check_real_moment_map_lifts()
    check_x9_matrix()
    print("E2 algebra and X_9 lattice checks passed.")


if __name__ == "__main__":
    main()

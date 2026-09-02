#!/usr/bin/env python3
"""Exact checks for the E3 moment-map/deformation-ring derivation.

The physical inputs are the E3 chiral-ring relations and the M2/divisor
charge formula cited in Section 2.4 of the paper.  This script checks the
resulting invariant-ring presentation, the dP6 root lattice, and explicit
representative lifts of the complex Cartan moment maps to the zero set of all
real Cartan moment maps.
"""

from __future__ import annotations

from cmath import exp, phase
from math import isclose, sqrt

from sympy import Matrix, Rational, groebner, symbols


def close_zero(value: complex, tolerance: float = 1e-12) -> bool:
    """Return whether a real or complex floating-point value vanishes."""
    return isclose(abs(value), 0.0, abs_tol=tolerance)


def balanced_lift(products: tuple[complex, ...]) -> tuple[list[complex], list[complex]]:
    """Lift z_i=q_i*qt_i with |q_i|=|qt_i| term by term."""
    q: list[complex] = []
    qt: list[complex] = []
    for product in products:
        radius = sqrt(abs(product))
        q.append(radius * exp(1j * phase(product)) if product else 0j)
        qt.append(radius + 0j)
    return q, qt


def check_union_scheme() -> None:
    """Check that a plane and a transverse line have the claimed union ring."""
    elimination, h1, h2, ell = symbols("u H1 H2 L")

    # I(plane)=(L), I(line)=(H1,H2).  The standard elimination formula
    # I cap J = (u I + (1-u) J) cap C[H1,H2,L] computes their union ideal.
    basis = groebner(
        [elimination * ell, (1 - elimination) * h1, (1 - elimination) * h2],
        elimination,
        h1,
        h2,
        ell,
        order="lex",
    )
    eliminated = {
        polynomial.as_expr()
        for polynomial in basis.polys
        if elimination not in polynomial.as_expr().free_symbols
    }
    expected = {h1 * ell, h2 * ell}
    assert eliminated == expected

    target = groebner([h1 * ell, h2 * ell], h1, h2, ell, order="lex")
    assert target.reduce(h1 * ell)[1] == 0
    assert target.reduce(h2 * ell)[1] == 0
    assert target.reduce(h1 * h2)[1] != 0
    assert target.reduce(ell**2)[1] != 0

    # Degree d>0 consists of d+1 plane monomials and one line monomial.
    for degree in range(1, 10):
        plane_monomials = {
            h1**power * h2 ** (degree - power)
            for power in range(degree + 1)
        }
        line_monomial = {ell**degree}
        assert len(plane_monomials | line_monomial) == degree + 2

    print("E3 torus-invariant ring: C[H1,H2,L]/(H1 L,H2 L)")
    print("Hilbert function: 1 in degree 0 and d+2 in every degree d>0")


def check_torus_invariants() -> None:
    """Check the balanced-monomial description for SU(3) and SU(2)."""
    # If d_i is the q_i exponent minus the qt_i exponent, the rows are the
    # ADHM gauge weight and independent flavor-Cartan weights.  Their kernels
    # vanish, so invariance forces d_i=0 separately for every pair.
    su3_weight_matrix = Matrix([[1, 1, 1], [1, 0, -1], [0, 1, -1]])
    su2_weight_matrix = Matrix([[1, 1], [1, -1]])
    assert su3_weight_matrix.det() == 3
    assert su2_weight_matrix.det() == -2
    assert su3_weight_matrix.nullspace() == []
    assert su2_weight_matrix.nullspace() == []

    # The complex ADHM equations leave n-1 independent products z_i.
    assert Matrix([[1, 1, 1]]).nullspace() == [
        Matrix([-1, 1, 0]),
        Matrix([-1, 0, 1]),
    ]
    assert Matrix([[1, 1]]).nullspace() == [Matrix([-1, 1])]
    print("SU(3) Cartan quotient: C[z1,z2,z3]/(z1+z2+z3) = C[H1,H2]")
    print("Root-basis coordinates: H1=z1, H2=-z3")
    print("SU(2) Cartan quotient: C[w1,w2]/(w1+w2) = C[L], L=w1")


def check_charge_coordinate_pairing() -> None:
    """Check that root charges multiply root-basis moment coordinates directly."""
    h1, h2, d1, d2 = symbols("H1 H2 d1 d2")
    moment_diagonal = Matrix([h1, -h1 + h2, -h2])
    fundamental_coweights = Matrix(
        [
            [Rational(2, 3), Rational(1, 3)],
            [Rational(-1, 3), Rational(1, 3)],
            [Rational(-1, 3), Rational(-2, 3)],
        ]
    )
    compact_generator = fundamental_coweights * Matrix([d1, d2])
    coupling = (moment_diagonal.T * compact_generator)[0].expand()
    assert coupling == d1 * h1 + d2 * h2

    # The coweight generator evaluates to d_i on each simple root.
    assert (Matrix([[1, -1, 0]]) * compact_generator)[0].expand() == d1
    assert (Matrix([[0, 1, -1]]) * compact_generator)[0].expand() == d2

    ell, charge = symbols("L d")
    su2_moment = Matrix([ell, -ell])
    su2_generator = Rational(1, 2) * Matrix([charge, -charge])
    assert (su2_moment.T * su2_generator)[0].expand() == charge * ell
    print("Coweight embedding check: mu_D = d1 H1 + d2 H2 (and d L for A1)")


def check_real_moment_map_lifts() -> None:
    """Check explicit zero-real-moment lifts on both E3 components."""
    a2_samples = (
        (0j, 0j),
        (1 + 0j, 0j),
        (1 + 2j, -3 + 1j),
        (-2j, 4 - 3j),
    )
    for h1, h2 in a2_samples:
        # diag(z)=H1*diag(1,-1,0)+H2*diag(0,1,-1).
        z1 = h1
        z2 = -h1 + h2
        z3 = -h2
        products = (z1, z2, z3)
        q, qt = balanced_lift(products)

        assert close_zero(sum(products))
        assert all(close_zero(qi * qti - zi) for qi, qti, zi in zip(q, qt, products))
        real_diagonal = [abs(qi) ** 2 - abs(qti) ** 2 for qi, qti in zip(q, qt)]
        assert all(close_zero(entry) for entry in real_diagonal)
        assert close_zero(sum(real_diagonal))
        assert close_zero(products[0] - h1)
        assert close_zero(-products[2] - h2)

    a1_samples = (0j, 1 + 0j, -2j, 3 + 4j)
    for ell in a1_samples:
        products = (ell, -ell)
        q, qt = balanced_lift(products)
        assert close_zero(sum(products))
        assert all(close_zero(qi * qti - zi) for qi, qti, zi in zip(q, qt, products))
        real_diagonal = [abs(qi) ** 2 - abs(qti) ** 2 for qi, qti in zip(q, qt)]
        assert all(close_zero(entry) for entry in real_diagonal)
        assert close_zero(products[0] - ell)

    print("Representative A2-plane and A1-line lifts have all real moments zero")


def check_dp6_root_lattice() -> None:
    """Check K-perpendicularity and the A2 plus A1 intersection matrix."""
    intersection = Matrix.diag(1, -1, -1, -1)  # (ell,e1,e2,e3)
    canonical = Matrix([-3, 1, 1, 1])
    a2_root_1 = Matrix([0, 1, -1, 0])
    a2_root_2 = Matrix([0, 0, 1, -1])
    a1_root = Matrix([1, -1, -1, -1])
    roots = Matrix.hstack(a2_root_1, a2_root_2, a1_root)

    assert (canonical.T * intersection * roots) == Matrix([[0, 0, 0]])
    gram = roots.T * intersection * roots
    assert gram == Matrix([[-2, 1, 0], [1, -2, 0], [0, 0, -2]])
    assert abs(int(gram.det())) == 6

    print("K_perp(dP6) Gram matrix:")
    print(gram)
    print("Root lattice: A2 plus A1, discriminant 6")


def main() -> None:
    check_union_scheme()
    check_torus_invariants()
    check_charge_coordinate_pairing()
    check_real_moment_map_lifts()
    check_dp6_root_lattice()
    print("E3 algebra, lifting, and root-lattice checks passed.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact finite checks for the two magnetic-quiver components of the E3 Higgs branch.

The calculation uses the abelian monopole formula with x=t^2.  It checks the
affine A1 and affine A2 quivers printed in arXiv:2007.15600, Appendix B, and
compares their unrefined Hilbert functions with the minimal nilpotent orbit
closures a1 and a2.  It is independent of the compact smoothing matrices.
"""

from __future__ import annotations

import hashlib
import json


MAX_DEGREE = 12


def convolve(left: list[int], right: list[int], degree: int) -> list[int]:
    out = [0] * (degree + 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if i + j <= degree:
                out[i + j] += a * b
    return out


def a1_monopole_coefficients(degree: int) -> list[int]:
    """Coulomb HS of U(1)^2/U(1)_diag with two bifundamentals."""
    bare = [1] + [2] * degree
    dressing = [1] * (degree + 1)  # (1-x)^(-1)
    return convolve(bare, dressing, degree)


def a2_monopole_coefficients(degree: int) -> list[int]:
    """Coulomb HS of the rank-one affine-A2 triangle modulo its diagonal U(1)."""
    bare = [0] * (degree + 1)
    # Fix the diagonal shift by setting the third magnetic charge to zero.
    # The exponent of x=t^2 is
    # (|m1-m2|+|m2|+|m1|)/2.
    for m1 in range(-degree, degree + 1):
        for m2 in range(-degree, degree + 1):
            twice_exponent = abs(m1 - m2) + abs(m2) + abs(m1)
            assert twice_exponent % 2 == 0
            exponent = twice_exponent // 2
            if exponent <= degree:
                bare[exponent] += 1
    dressing = [n + 1 for n in range(degree + 1)]  # (1-x)^(-2)
    return convolve(bare, dressing, degree)


def mat_vec(matrix: tuple[tuple[int, ...], ...], vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in matrix)


def main() -> None:
    adjacency_a1 = ((0, 2), (2, 0))
    adjacency_a2 = ((0, 1, 1), (1, 0, 1), (1, 1, 0))
    ranks_a1 = (1, 1)
    ranks_a2 = (1, 1, 1)

    # Every node is balanced: 2 r_i = sum_j A_ij r_j.
    assert mat_vec(adjacency_a1, ranks_a1) == tuple(2 * r for r in ranks_a1)
    assert mat_vec(adjacency_a2, ranks_a2) == tuple(2 * r for r in ranks_a2)

    # Removing the decoupled diagonal U(1) gives the quaternionic dimensions.
    qdim_a1 = sum(ranks_a1) - 1
    qdim_a2 = sum(ranks_a2) - 1
    assert (qdim_a1, qdim_a2) == (1, 2)

    hs_a1 = a1_monopole_coefficients(MAX_DEGREE)
    hs_a2 = a2_monopole_coefficients(MAX_DEGREE)

    # Minimal nilpotent orbit closures:
    # a1: degree-n piece V_{2n} of SU(2), dimension 2n+1;
    # a2: degree-n piece V_{(n,n)} of SU(3), dimension (n+1)^3.
    expected_a1 = [2 * n + 1 for n in range(MAX_DEGREE + 1)]
    expected_a2 = [(n + 1) ** 3 for n in range(MAX_DEGREE + 1)]
    assert hs_a1 == expected_a1
    assert hs_a2 == expected_a2

    # The two cones meet only at the origin, so their union obeys the surgery formula.
    hs_e3 = [1] + [hs_a1[n] + hs_a2[n] for n in range(1, MAX_DEGREE + 1)]
    expected_e3 = [1] + [2 * n + 1 + (n + 1) ** 3 for n in range(1, MAX_DEGREE + 1)]
    assert hs_e3 == expected_e3

    payload = {
        "degree_variable": "x=t^2",
        "max_degree": MAX_DEGREE,
        "a1": {
            "adjacency": adjacency_a1,
            "ranks": ranks_a1,
            "quaternionic_dimension": qdim_a1,
            "hilbert_coefficients": hs_a1,
        },
        "a2": {
            "adjacency": adjacency_a2,
            "ranks": ranks_a2,
            "quaternionic_dimension": qdim_a2,
            "hilbert_coefficients": hs_a2,
        },
        "e3_union_hilbert_coefficients": hs_e3,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    print("E3 magnetic-quiver check SHA-256:", digest)
    print("a1: affine A1, qdim=1, HS coefficients:", hs_a1)
    print("a2: affine A2, qdim=2, HS coefficients:", hs_a2)
    print("E3=a2 union a1 coefficients:", hs_e3)
    print("All E3 magnetic-quiver checks passed.")


if __name__ == "__main__":
    main()

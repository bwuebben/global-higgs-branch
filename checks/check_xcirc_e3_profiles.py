#!/usr/bin/env python3
"""Exact compact-gauging checks for the four E3 profiles of X-circ.

The sparse 36 x 26 matrix below is the matrix independently constructed by
Paper 4's ``mixed_candidate.sage`` and ``global_kernel.py`` computations.  Its
rows are local deformation/moment-map coordinates and its columns are the
divisor basis complementary to rays (15, 22, 24, 26).  Thus the physical
charge map on a profile is the transpose of the corresponding row restriction.

This self-contained checker verifies ranks, Smith invariants, projected-kernel
dimensions, forced coordinates, inactive compact-vector combinations, and the
five primitive moment-map combinations printed in Appendix B of the paper.
"""

from __future__ import annotations

import hashlib
import json
from functools import reduce
from math import gcd

from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


DIVISOR_RAYS = (
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
    16, 17, 18, 19, 20, 21, 23, 25, 27, 28, 29, 30,
)

LIVE_MATRIX_SHA256 = "61835cdf70f3dbede5ee89fae8abf6e7c28be69562dbe4194b3fc23a801822ae"

SPARSE_ROWS: dict[str, tuple[tuple[int, int], ...]] = {
    "N_1": ((1, 1), (2, -1), (7, -1), (8, 1)),
    "N_2": ((1, 1), (2, -1), (13, -1), (14, 1)),
    "N_3": ((1, 1), (3, -1), (7, -1), (9, 1)),
    "N_4": ((2, 1), (4, -1), (8, -1), (10, 1)),
    "N_5": ((2, 1), (4, -1), (14, -1)),
    "N_6": ((3, 1), (5, -1), (9, -1), (11, 1)),
    "N_7": ((3, 1), (5, -1), (13, -1), (16, 1)),
    "N_8": ((3, 1), (5, -1), (23, -1), (25, 1)),
    "N_9": ((3, 1), (9, -1), (19, 1), (23, -1)),
    "N_10": ((4, 1), (6, -1), (10, -1), (12, 1)),
    "N_11": ((4, 1), (10, -1)),
    "N_12": ((5, 1), (6, -1), (11, -1), (12, 1)),
    "N_13": ((5, 1), (6, -1), (16, -1)),
    "N_14": ((5, 1), (6, -1), (25, -1)),
    "N_15": ((7, 1), (8, -1), (17, -1), (18, 1)),
    "N_16": ((7, 1), (9, -1), (17, -1), (19, 1)),
    "N_17": ((10, 1), (12, -1), (18, -1), (20, 1)),
    "N_18": ((10, 1), (12, -1)),
    "N_19": ((11, 1), (12, -1), (19, -1), (20, 1)),
    "N_20": ((11, 1), (12, -1), (25, -1)),
    "N_21": ((13, 1), (14, -1), (16, -1)),
    "N_22": ((13, 1), (16, -1), (23, -1), (25, 1)),
    "N_23": ((16, -1), (25, 1)),
    "N_24": ((17, 1), (18, -1), (19, -1), (20, 1)),
    "N_25": ((18, 1), (20, -1)),
    "N_26": ((19, 1), (20, -1), (25, -1)),
    "D_6,1:s1": ((1, 1), (2, -1), (5, -1), (6, 1)),
    "D_6,1:s2": ((2, 1), (3, -1), (4, -1), (5, 1)),
    "D_6,1:s3": ((1, 1), (2, -1), (3, -1), (4, 1), (5, 1), (6, -1)),
    "D_6,2:s1": ((7, 1), (8, -1), (11, -1), (12, 1)),
    "D_6,2:s2": ((8, 1), (9, -1), (10, -1), (11, 1)),
    "D_6,2:s3": ((7, 1), (8, -1), (9, -1), (10, 1), (11, 1), (12, -1)),
    "D_7,1:alpha": ((1, -2), (7, 2), (13, 1), (17, -1)),
    "D_7,1:beta": ((1, -1), (13, 1), (17, 1), (21, -1)),
    "D_7,2:alpha": ((2, 2), (8, -1), (14, -2)),
    "D_7,2:beta": ((8, 1), (14, -1), (18, -1)),
}

LOCAL_LABELS = tuple(SPARSE_ROWS)
DIVISOR_POSITION = {ray: index for index, ray in enumerate(DIVISOR_RAYS)}

PROFILE_FORBIDDEN = {
    "LL": {
        "D_6,1:s1", "D_6,1:s2", "D_6,2:s1", "D_6,2:s2",
        "D_7,1:alpha", "D_7,2:alpha",
    },
    "LP": {
        "D_6,1:s1", "D_6,1:s2", "D_6,2:s3",
        "D_7,1:alpha", "D_7,2:alpha",
    },
    "PL": {
        "D_6,1:s3", "D_6,2:s1", "D_6,2:s2",
        "D_7,1:alpha", "D_7,2:alpha",
    },
    "PP": {
        "D_6,1:s3", "D_6,2:s3",
        "D_7,1:alpha", "D_7,2:alpha",
    },
}

EXPECTED = {
    "LL": {
        "allowed": 30,
        "rank": 21,
        "kernel": 9,
        "forced": ("N_9", "N_11", "D_7,1:beta", "D_7,2:beta"),
    },
    "LP": {
        "allowed": 31,
        "rank": 21,
        "kernel": 10,
        "forced": (
            "N_9", "N_11", "D_6,1:s3", "D_7,1:beta", "D_7,2:beta",
        ),
    },
    "PL": {
        "allowed": 31,
        "rank": 21,
        "kernel": 10,
        "forced": (
            "N_9", "N_11", "D_6,2:s3", "D_7,1:beta", "D_7,2:beta",
        ),
    },
    "PP": {
        "allowed": 32,
        "rank": 20,
        "kernel": 12,
        "forced": ("N_9", "N_11", "D_7,1:beta", "D_7,2:beta"),
    },
}


def dense_matrix() -> Matrix:
    """Return B with local rows and divisor columns."""
    rows = []
    for label in LOCAL_LABELS:
        row = [0] * len(DIVISOR_RAYS)
        for ray, coefficient in SPARSE_ROWS[label]:
            row[DIVISOR_POSITION[ray]] = coefficient
        rows.append(row)
    return Matrix(rows)


def divisor_vector(coefficients: dict[int, int]) -> Matrix:
    """Return a column vector in the printed divisor basis."""
    result = [0] * len(DIVISOR_RAYS)
    for ray, coefficient in coefficients.items():
        result[DIVISOR_POSITION[ray]] = coefficient
    return Matrix(result)


def nonzero_moment(moment: Matrix) -> tuple[tuple[str, int], ...]:
    """Return the nonzero labeled entries of a local moment-map vector."""
    return tuple(
        (LOCAL_LABELS[index], int(entry))
        for index, entry in enumerate(moment)
        if entry != 0
    )


def smith_invariants(matrix: Matrix) -> tuple[int, ...]:
    """Return the absolute nonzero Smith invariants."""
    smith = smith_normal_form(matrix, domain=ZZ)
    return tuple(
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if smith[index, index] != 0
    )


def primitive(vector: Matrix) -> bool:
    """Return whether an integral vector is primitive."""
    entries = [abs(int(entry)) for entry in vector if entry]
    return bool(entries) and reduce(gcd, entries) == 1


B = dense_matrix()

VECTORS = {
    "V_9": divisor_vector(
        {ray: 1 for ray in tuple(range(1, 13)) + (17, 18, 19, 20)}
    ),
    "V_11": divisor_vector(
        {ray: -1 for ray in tuple(range(7, 13)) + (17, 18, 19, 20, 21)}
    ),
    "V_beta1": divisor_vector({21: -1}),
    "V_beta2": divisor_vector(
        {**{ray: 1 for ray in range(1, 13)}, 21: -1}
    ),
    "V_A1": divisor_vector(
        {
            1: 1, 3: 1, 5: 1, 7: 1, 9: 1, 11: 1, 13: 1,
            16: 1, 17: 1, 19: 1, 21: 1, 23: 1, 25: 1,
        }
    ),
}

EXPECTED_MOMENTS = {
    "V_9": (("N_9", 1), ("D_7,1:alpha", -1), ("D_7,2:alpha", 1)),
    "V_11": (("N_11", 1), ("D_7,1:alpha", -1), ("D_7,2:alpha", 1)),
    "V_beta1": (("D_7,1:beta", 1),),
    "V_beta2": (("D_7,2:alpha", 1), ("D_7,2:beta", 1)),
    "V_A1": (("D_6,1:s3", 1), ("D_6,2:s3", 1)),
}


def check_full_matrix() -> None:
    """Check the unrestricted matrix and primitive isolating vectors."""
    assert B.shape == (36, 26)
    payload = {
        "labels": LOCAL_LABELS,
        "rows": [[int(entry) for entry in B.row(index)] for index in range(B.rows)],
    }
    digest = hashlib.sha256(
        json.dumps(payload, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()
    assert digest == LIVE_MATRIX_SHA256
    assert B.rank() == 21
    assert smith_invariants(B) == (1,) * 21

    for name, vector in VECTORS.items():
        assert primitive(vector)
        assert nonzero_moment(B * vector) == EXPECTED_MOMENTS[name]

    print("X-circ matrix SHA-256:", digest)
    print("X-circ matrix: 36 x 26, rank 21, Smith invariants 1^21")
    for name in VECTORS:
        print(f"  {name}: {EXPECTED_MOMENTS[name]}")


def check_profiles() -> None:
    """Check all four branch-restricted charge maps."""
    for profile, forbidden in PROFILE_FORBIDDEN.items():
        allowed_indices = tuple(
            index
            for index, label in enumerate(LOCAL_LABELS)
            if label not in forbidden
        )
        labels = tuple(LOCAL_LABELS[index] for index in allowed_indices)
        profile_matrix = B[list(allowed_indices), :]
        charge_map = profile_matrix.T
        rank = charge_map.rank()
        kernel_basis = charge_map.nullspace()
        forced = tuple(
            label
            for coordinate, label in enumerate(labels)
            if all(vector[coordinate] == 0 for vector in kernel_basis)
        )

        expected = EXPECTED[profile]
        assert len(labels) == expected["allowed"]
        assert rank == expected["rank"]
        assert len(kernel_basis) == expected["kernel"]
        assert forced == expected["forced"]
        assert smith_invariants(profile_matrix) == (1,) * rank
        assert len(profile_matrix.nullspace()) == 26 - rank

        # On every reduced profile, the first four vectors isolate the four
        # universally forced coordinates because alpha_1=alpha_2=0.
        reduced_moments = {
            name: tuple(
                (label, coefficient)
                for label, coefficient in EXPECTED_MOMENTS[name]
                if label in labels
            )
            for name in ("V_9", "V_11", "V_beta1", "V_beta2")
        }
        assert reduced_moments == {
            "V_9": (("N_9", 1),),
            "V_11": (("N_11", 1),),
            "V_beta1": (("D_7,1:beta", 1),),
            "V_beta2": (("D_7,2:beta", 1),),
        }

        active_a1_moment = tuple(
            item for item in EXPECTED_MOMENTS["V_A1"] if item[0] in labels
        )
        expected_a1 = {
            "LL": (("D_6,1:s3", 1), ("D_6,2:s3", 1)),
            "LP": (("D_6,1:s3", 1),),
            "PL": (("D_6,2:s3", 1),),
            "PP": (),
        }
        assert active_a1_moment == expected_a1[profile]

        print(
            f"{profile}: local={len(labels)}, effective rank={rank}, "
            f"projected kernel={len(kernel_basis)}, inactive vectors={26-rank}"
        )
        print("  forced zero:", forced)
        print("  V_A1 moment:", active_a1_moment or "inactive")


def check_inactive_vector_lattices() -> None:
    """Check the five universal inactive directions and the extra PP one."""
    inactive = [divisor_vector({ray: 1}) for ray in (27, 28, 29, 30)]
    anticanonical = divisor_vector(
        {
            1: 3, 2: 2, 3: 3, 4: 1, 5: 2, 6: 1,
            7: 3, 8: 2, 9: 3, 10: 1, 11: 2, 12: 1,
            13: 2, 14: 1, 16: 1, 17: 2, 18: 1, 19: 2,
            20: 1, 21: 1, 23: 2, 25: 1,
        }
    )
    inactive.append(anticanonical)
    inactive_matrix = Matrix.hstack(*inactive)
    assert inactive_matrix.rank() == 5
    assert B * inactive_matrix == Matrix.zeros(36, 5)

    pp_allowed = tuple(
        index
        for index, label in enumerate(LOCAL_LABELS)
        if label not in PROFILE_FORBIDDEN["PP"]
    )
    pp_matrix = B[list(pp_allowed), :]
    pp_inactive = Matrix.hstack(inactive_matrix, VECTORS["V_A1"])
    assert pp_inactive.rank() == 6
    assert pp_matrix * pp_inactive == Matrix.zeros(len(pp_allowed), 6)

    print("Universal inactive vector lattice: rank 5")
    print("PP inactive vector lattice: rank 6, with primitive extra generator V_A1")


def main() -> None:
    check_full_matrix()
    check_profiles()
    check_inactive_vector_lattices()
    print("All X-circ E3 profile-gauging checks passed.")


if __name__ == "__main__":
    main()

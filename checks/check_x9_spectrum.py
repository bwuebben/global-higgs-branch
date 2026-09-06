#!/usr/bin/env python3
"""Exact lattice and topology checks for the compact X_9 transition.

This check is self-contained.  It verifies:

1. the Batyrev Hodge numbers (6,122) of the MPCP resolution;
2. a saturated six-divisor basis of the toric Picard lattice, which is a
   complete basis of H^2 of the resolved hypersurface over Q;
3. the two primitive compact flavor gaugings of (u_1,u_2,beta);
4. the one-dimensional reduced Higgs kernel C(-1,-1,1);
5. the Mayer--Vietoris calculation giving (h11,h21)=(3,123) on the
   smooth side; and
6. the resulting five-dimensional massless multiplet counts.

Only SymPy is required.  Run from global_higgs_branch with

    ./venv/bin/python checks/check_x9_spectrum.py
"""

from __future__ import annotations

from collections import Counter
from functools import reduce
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form


V = (
    (1, 0, 0, 0),
    (0, 1, 0, 0),
    (-1, -1, 0, 0),
    (0, 0, 1, 0),
    (0, 0, 0, 1),
    (-4, -2, -1, 0),
    (-4, -2, 0, -1),
    (3, 1, 1, 1),
    (2, 1, 1, 1),
)

# The unique interior point of the dP_7 pentagon.  Its star subdivision
# produces the exceptional divisor E = D_p.
PENTAGON_INTERIOR = (-2, -1, 0, 0)

# This is the sole facet-interior lattice point.  The corresponding toric
# divisor misses a generic anticanonical hypersurface and is removed before
# forming H^2 of the Calabi--Yau.
FACET_INTERIOR = (-1, 0, 0, 0)

# Divisor rays that meet the resolved Calabi--Yau, in the paper's order.
RAYS = V + (PENTAGON_INTERIOR,)


def dot(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    return sum(x * y for x, y in zip(a, b))


def primitive(v: list[int]) -> tuple[int, ...]:
    g = reduce(gcd, (abs(x) for x in v if x), 0)
    assert g
    return tuple(x // g for x in v)


def facet_inequalities(vertices: tuple[tuple[int, ...], ...]):
    """Return the primitive inequalities <n,x> <= 1 of a reflexive 4-polytope."""
    facets = {}
    for inds in combinations(range(len(vertices)), 4):
        base = Matrix(vertices[inds[0]])
        diffs = Matrix.hstack(
            *(Matrix(vertices[i]) - base for i in inds[1:])
        ).T
        if diffs.rank() != 3:
            continue
        ns = diffs.nullspace()
        assert len(ns) == 1
        q = ns[0]
        den = reduce(lcm, (int(x.q) for x in q), 1)
        n = primitive([int(x * den) for x in q])
        c = dot(n, vertices[inds[0]])
        vals = [dot(n, v) for v in vertices]
        if all(x >= c for x in vals):
            n = tuple(-x for x in n)
            c = -c
            vals = [-x for x in vals]
        if not all(x <= c for x in vals):
            continue
        if c < 0:
            n = tuple(-x for x in n)
            c = -c
            vals = [-x for x in vals]
        if c != 1:
            continue
        support = tuple(i for i, x in enumerate(vals) if x == c)
        if Matrix(
            [[vertices[i][j] - vertices[support[0]][j] for j in range(4)]
             for i in support[1:]]
        ).rank() != 3:
            continue
        facets[n] = support
    return tuple(sorted((n, support) for n, support in facets.items()))


def lattice_points(vertices, facets):
    lo = [min(v[j] for v in vertices) for j in range(4)]
    hi = [max(v[j] for v in vertices) for j in range(4)]
    return tuple(
        p
        for p in product(*(range(lo[j], hi[j] + 1) for j in range(4)))
        if all(dot(n, p) <= 1 for n, _ in facets)
    )


def tight_set(point, facets):
    return tuple(i for i, (n, _) in enumerate(facets) if dot(n, point) == 1)


def edge_length(a, b):
    return reduce(gcd, (abs(x - y) for x, y in zip(a, b)), 0)


def batyrev_side(vertices):
    """Return l(Delta), sum facet interiors, and the codimension-two correction."""
    facets = facet_inequalities(vertices)
    points = lattice_points(vertices, facets)
    counts = Counter(tight_set(p, facets) for p in points)
    facet_interiors = sum(n for tight, n in counts.items() if len(tight) == 1)
    correction = 0
    for tight, npoints in counts.items():
        if len(tight) != 2:
            continue
        n1, n2 = facets[tight[0]][0], facets[tight[1]][0]
        correction += npoints * (edge_length(n1, n2) - 1)
    return facets, points, facet_interiors, correction


def hodge_numbers():
    facets, points, sfac, corr = batyrev_side(V)
    dual_vertices = tuple(n for n, _ in facets)
    _, dual_points, sfac_dual, corr_dual = batyrev_side(dual_vertices)
    h11 = len(points) - 5 - sfac + corr
    h21 = len(dual_points) - 5 - sfac_dual + corr_dual
    return {
        "facets": facets,
        "points": points,
        "dual_points": dual_points,
        "sfac": sfac,
        "corr": corr,
        "sfac_dual": sfac_dual,
        "corr_dual": corr_dual,
        "h11": h11,
        "h21": h21,
    }


H = hodge_numbers()
assert len(H["facets"]) == 12
assert len(H["points"]) == 12
assert len(H["dual_points"]) == 162
assert H["sfac"] == 1 and H["corr"] == 0
assert H["sfac_dual"] == 35 and H["corr_dual"] == 0
assert (H["h11"], H["h21"]) == (6, 122)

boundary_nonfacet = {
    p for p in H["points"]
    if p != (0, 0, 0, 0) and len(tight_set(p, H["facets"])) >= 2
}
facet_interior_points = {
    p for p in H["points"] if len(tight_set(p, H["facets"])) == 1
}
assert boundary_nonfacet == set(RAYS)
assert facet_interior_points == {FACET_INTERIOR}


# Principal-divisor map M -> Z^{10}.  The four columns v_0,v_1,v_3,v_4
# form the identity, so eliminating D_0,D_1,D_3,D_4 is integral and leaves
# the saturated toric Picard basis (D_2,D_5,D_6,D_7,D_8,D_p).  The Batyrev
# correction vanishes and its rank equals h11, so it is a complete H^2-basis
# over Q. Equality with the full free integral H^2 lattice is established
# separately by check_x9_integral_screening.sage, not by this count.
ray_matrix = Matrix.hstack(*(Matrix(v) for v in RAYS))
pivot = (0, 1, 3, 4)
free = (2, 5, 6, 7, 8, 9)
assert ray_matrix[:, pivot].det() == 1
assert ray_matrix.rank() == 4
PICARD_LABELS = ("D2", "D5", "D6", "D7", "D8", "E=Dp")


# Rows are (u_1,u_2,beta,alpha); columns are the ten divisor rays above.
# The first two are the node circuit relations.  beta is the dP_7 root
# moment map and alpha is its nilpotent complementary tangent direction.
B_FULL = Matrix([
    [0, 0,  1, 1,  0,  0, -1, -1,  0, 0],
    [0, 0, -1, 0, -1,  1,  0,  1,  0, 0],
    [0, 0,  0, 1, -1,  1, -1,  0,  0, 0],
    [0, 0,  0, 0,  3, -2,  1,  0, -2, 0],
])
assert B_FULL * ray_matrix.T == Matrix.zeros(4, 4)

B_REDUCED_PICARD = B_FULL[:3, list(free)]
Q_PICARD = B_REDUCED_PICARD.T
assert Q_PICARD == Matrix([
    [ 1, -1,  0],
    [ 0,  1,  1],
    [-1,  0, -1],
    [-1,  1,  0],
    [ 0,  0,  0],
    [ 0,  0,  0],
])

# Columns are the new integral basis in old Picard coordinates:
#   V1=-D6, V2=D5,
#   S1=D2+D7, S2=D2+D5+D6, S3=D8, Phi=E.
U = Matrix.hstack(
    Matrix([0, 0, -1, 0, 0, 0]),
    Matrix([0, 1,  0, 0, 0, 0]),
    Matrix([1, 0,  0, 1, 0, 0]),
    Matrix([1, 1,  1, 0, 0, 0]),
    Matrix([0, 0,  0, 0, 1, 0]),
    Matrix([0, 0,  0, 0, 0, 1]),
)
assert abs(U.det()) == 1
Q_ISOLATED = U.T * Q_PICARD
assert Q_ISOLATED == Matrix([
    [1, 0, 1],
    [0, 1, 1],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
])

snf = smith_normal_form(Q_PICARD, domain=ZZ)
nonzero_snf = tuple(abs(int(snf[i, i])) for i in range(min(snf.shape))
                    if snf[i, i] != 0)
assert nonzero_snf == (1, 1)

Q_EFFECTIVE = Q_ISOLATED[:2, :]
kernel = Q_EFFECTIVE.nullspace()
assert len(kernel) == 1
kernel_generator = tuple(int(x) for x in kernel[0])
assert kernel_generator == (-1, -1, 1)

# The full scheme has extra alpha couplings, but alpha vanishes on every
# ordinary reduced Higgs trajectory.  In the isolated basis the two effective
# vectors have alpha coefficients -1 and -2; Phi=E remains completely neutral.
Q_FULL_PICARD = B_FULL[:, list(free)].T
Q_FULL_ISOLATED = U.T * Q_FULL_PICARD
assert tuple(Q_FULL_ISOLATED[0, :]) == (1, 0, 1, -1)
assert tuple(Q_FULL_ISOLATED[1, :]) == (0, 1, 1, -2)
assert tuple(Q_FULL_ISOLATED[5, :]) == (0, 0, 0, 0)


# Rational Mayer--Vietoris bookkeeping.  L is the sum of the three links, R
# the resolution pieces, and M the Milnor fibres.  The branch kernel has
# dimension one by the charge calculation above.
b2_Y = H["h11"]
b2_L = 1 + 1 + 2                 # two nodes and the dP_7 link
b2_R = 1 + 1 + 3                 # two P^1s and dP_7
b2_M = 0 + 0 + 1                 # node Milnor fibres and dP_7 Milnor fibre
branch_kernel_dim = len(kernel)

# L -> R is injective, so b2(Y)=b2(W)+b2(R)-b2(L).
b2_W = b2_Y - b2_R + b2_L

# ker[L -> W+M] is the branch-restricted global relation kernel.
rank_L_to_WM = b2_L - branch_kernel_dim
b2_Z = b2_W + b2_M - rank_L_to_WM
assert (b2_L, b2_R, b2_M, b2_W, rank_L_to_WM, b2_Z) == (4, 5, 1, 5, 3, 3)

# Local Euler replacement: each node contributes 2-0, while dP_7 contributes
# chi(dP_7)-chi(M_7)=5-1.  The dP_7 Milnor fibre has (b2,b3)=(1,1).
chi_Y = 2 * (H["h11"] - H["h21"])
chi_local_change = 2 * (2 - 0) + (5 - 1)
chi_Z = chi_Y - chi_local_change
h21_Z = b2_Z - chi_Z // 2
assert (chi_Y, chi_local_change, chi_Z, h21_Z) == (-232, 8, -240, 123)

# Five-dimensional M-theory counts on a smooth compact Calabi--Yau threefold.
nV_Y, nH_Y = H["h11"] - 1, H["h21"] + 1
nV_Z, nH_Z = b2_Z - 1, h21_Z + 1
assert (nV_Y, nH_Y, nV_Z, nH_Z) == (5, 123, 2, 124)
assert (nV_Z - nV_Y, nH_Z - nH_Y) == (-3, 1)


payload = repr((RAYS, tuple(map(tuple, B_FULL.tolist())), tuple(map(tuple, U.tolist()))))
digest = sha256(payload.encode("ascii")).hexdigest()

print(f"X_9 lattice check SHA-256: {digest}")
print("MPCP resolution: (h11,h21)=(6,122), chi=-232")
print("Picard basis:", PICARD_LABELS)
print("Isolated reduced charge matrix:")
print(Q_ISOLATED)
print("Effective vectors: V1=-D6, V2=D5")
print("Spectators on the reduced branch: D2+D7, D2+D5+D6, D8, E=Dp")
print("Smith invariants:", nonzero_snf)
print("Reduced Higgs kernel:", kernel_generator)
print("Mayer-Vietoris ranks: b2(L,R,M,W,Z)=", (b2_L, b2_R, b2_M, b2_W, b2_Z))
print("Smooth fibre: (h11,h21)=(3,123), chi=-240")
print("5d endpoints: (nV,nH)=(5,123) -> (2,124)")
print("All X_9 compact-spectrum checks passed.")

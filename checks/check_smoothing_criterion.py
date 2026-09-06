#!/usr/bin/env python3
"""Exact linear and polynomial checks for the mixed smoothing translation.

Standard library only. The counterexample checks use the three printed
divisor identities, not a replacement for the full toric construction in
Smoothing Calabi--Yau Threefolds with Nodes and del Pezzo Cone Points, §5.
No finite computation here proves analytic integrability or sheaf stability.
"""

from fractions import Fraction


def rank(rows):
    rows = [[Fraction(x) for x in row] for row in rows]
    if not rows:
        return 0
    pivot_row = 0
    for col in range(len(rows[0])):
        pivot = next((i for i in range(pivot_row, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        value = rows[pivot_row][col]
        rows[pivot_row] = [x / value for x in rows[pivot_row]]
        for i in range(pivot_row + 1, len(rows)):
            value = rows[i][col]
            rows[i] = [x - value * y for x, y in zip(rows[i], rows[pivot_row])]
        pivot_row += 1
    return pivot_row


def smoothing_test(matrix, factors):
    """No discriminant functional may vanish identically on ker(matrix)."""
    base_rank = rank(matrix)
    return all(rank(matrix + [factor]) > base_rank for factor in factors)


def product_of_linear_forms(forms):
    polynomial = {(0, 0): 1}
    for a, b in forms:
        out = {}
        for (i, j), coefficient in polynomial.items():
            out[i + 1, j] = out.get((i + 1, j), 0) + coefficient * a
            out[i, j + 1] = out.get((i, j + 1), 0) + coefficient * b
        polynomial = {power: coefficient for power, coefficient in out.items() if coefficient}
    return polynomial


def main():
    checks = 0

    def check(condition):
        nonlocal checks
        assert condition
        checks += 1

    # H1 R1 + H2 R2 = H1 e1 + (H2-H1)e2 - H2 e3.
    e3_factors = [(1, 0), (-1, 1), (0, -1)]
    check(tuple(map(sum, zip(*e3_factors))) == (0, 0))
    check(product_of_linear_forms(e3_factors) == {(2, 1): 1, (1, 2): -1})
    check(smoothing_test([[0, 0]], e3_factors))
    for line in e3_factors:
        check(not smoothing_test([line], e3_factors))
    check(smoothing_test([[2, -1]], e3_factors))
    # An inactive cone is excluded; a nonzero image alone is insufficient.
    check(not smoothing_test([[1, 0], [0, 1]], e3_factors))
    check(rank([[1, -1]]) == 1)
    check(not smoothing_test([[1, -1]], e3_factors))

    q9 = [[1, 0, 1], [0, 1, 1]]
    coordinate_factors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    check(rank(q9) == 2)
    check(smoothing_test(q9, coordinate_factors))
    check(all(sum(x * y for x, y in zip(row, [-1, -1, 1])) == 0 for row in q9))
    check(not smoothing_test(q9 + [[0, 0, 1]], coordinate_factors))

    # Counterexample coordinates: n13,n16,lambda0,a0,b0,lambda1,a1,b1.
    identities = [
        [1, 0, 0, -1, 1, 0, -1, 1],
        [0, 0, 0, 0, 1, 0, 0, 1],
        [0, 1, 1, 0, 0, 1, 0, 0],
    ]
    profiles = {
        'LL': ([0, 1, 2, 5], 0),
        'LP': ([0, 1, 2, 6, 7], 7),
        'PL': ([0, 1, 3, 4, 5], 4),
        'PP': ([0, 1, 3, 4, 6, 7], 1),
    }
    for profile, (columns, forced) in profiles.items():
        restricted = [[row[j] for j in columns] for row in identities]
        factor = [int(j == forced) for j in columns]
        check(not smoothing_test(restricted, [factor]))
        if profile in ('LP', 'PL'):
            # Both cone blocks and both displayed nodal coordinates can be nonzero.
            witness = [1, -3, 3, 1, 0] if profile == 'LP' else [1, -3, 1, 0, 3]
            check(all(sum(x * y for x, y in zip(row, witness)) == 0 for row in restricted))
            check(all(witness[:2]))
            line_value = witness[2] if profile == 'LP' else witness[4]
            plane_value = witness[3:5] if profile == 'LP' else witness[2:4]
            check(bool(line_value) and any(plane_value))
    check(product_of_linear_forms([(-1, 0), (0, 1), (1, -1)]) ==
          {(2, 1): -1, (1, 2): 1})
    # In the paper's root convention H1=-a, H2=b-a.
    root_change = [(-1, 0), (-1, 1)]
    composed = [tuple(sum(f[k] * root_change[k][j] for k in range(2))
                      for j in range(2)) for f in e3_factors]
    check(composed == [(-1, 0), (0, 1), (1, -1)])
    print(f'PASS: {checks} exact smoothing-criterion checks.')


if __name__ == '__main__':
    main()

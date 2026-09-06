#!/usr/bin/env python3
"""Exact conifold-threshold and special-geometry identities for paper 2.

Run with Python and SymPy 1.14.0. The four marked charges are inputs,
not proved by this script. Checks cover derivatives, unit monodromy,
primitive versus monodromy lattices, neutral thresholds, and the
homogeneous special-geometry reconstruction. No BPS stability or
finite-Planck hypermultiplet metric is inferred.
"""

from sympy import I, Matrix, Rational, ZZ, diff, exp, expand_func, eye, limit, log, pi, polylog, simplify, symbols, zeros
from sympy.matrices.normalforms import smith_normal_form


def main():
    count = 0

    def check(condition):
        nonlocal count
        assert condition
        count += 1

    x = Matrix(symbols("x1 x2 x3", positive=True))
    zstar = symbols("zstar", positive=True)
    q = Matrix([[1, 0, 1, 0], [0, 1, 1, 0], [0, 0, 1, 1]])
    z = q.T * x
    check(z == Matrix([x[0], x[1], sum(x), x[2]]))
    check(q * Matrix([-1, -1, 1, -1]) == zeros(3, 1))
    check(q.rank() == 3)
    check([abs(d) for d in smith_normal_form(q, domain=ZZ).diagonal()] == [1, 1, 1])
    beta = q * q.T
    check(beta == Matrix([[2, 1, 1], [1, 2, 1], [1, 1, 2]]))
    check(beta.det() == 4)
    check([abs(d) for d in smith_normal_form(beta, domain=ZZ).diagonal()] == [1, 1, 4])
    check(beta.eigenvals() == {1: 2, 4: 1})

    fsing = sum(v**2 * (log(v / zstar) - Rational(3, 2)) for v in z) / (4 * pi * I)
    expected = sum((q[:, a] * q[:, a].T * log(z[a] / zstar) for a in range(4)), zeros(3)) / (2 * pi * I)
    hessian = Matrix(3, 3, lambda i, j: diff(fsing, x[i], x[j]))
    check((hessian - expected).applyfunc(simplify) == zeros(3))
    check(simplify(diff(fsing, zstar) + sum(v**2 for v in z) / (4 * pi * I * zstar)) == 0)
    # The manuscript uses -Li_3/(2*pi*i)^3 for a unit conifold index.
    v = symbols("v")
    li3_term = -polylog(3, exp(2 * pi * I * v)) / (2 * pi * I)**3
    check(simplify(expand_func(diff(li3_term, v, 2)) - log(1 - exp(2 * pi * I * v)) / (2 * pi * I)) == 0)
    check(limit((1 - exp(2 * pi * I * v)) / v, v, 0) == -2 * pi * I)

    # A positive continuation adds 2*pi*i to one logarithm.
    omega = zeros(3).row_join(eye(3)).col_join((-eye(3)).row_join(zeros(3)))
    monodromies = []
    for a in range(4):
        delta_f = z[a]**2 / 2
        check(Matrix([diff(delta_f, v) for v in x]) == q[:, a] * z[a])
        shift = q[:, a] * q[:, a].T
        check(Matrix(3, 3, lambda i, j: diff(delta_f, x[i], x[j])) == shift)
        monodromy = eye(3).row_join(zeros(3)).col_join(shift.row_join(eye(3)))
        check(monodromy.T * omega * monodromy == omega)
        check((monodromy - eye(6))**2 == zeros(6))
        monodromies.append(monodromy)
    check(all(a * b == b * a for a in monodromies for b in monodromies))
    combined = eye(6)
    for monodromy in monodromies:
        combined *= monodromy
    check(combined[3:, :3] == beta)
    check((combined - eye(6)).rank() == 3)

    # The surviving coordinates t2,t3,t4 annihilate all four charges.
    qfull = Matrix([[0, 1, 1, 0], [0, 0, 0, 0], [0, 0, 0, 0],
                    [0, 0, 0, 0], [1, 0, 1, 0], [0, 0, 1, 1]])
    neutral = eye(6)[:, 1:4]
    check(qfull.T * neutral == zeros(4, 3))
    for a in range(4):
        check(neutral.T * qfull[:, a] * qfull[:, a].T == zeros(3, 6))

    # Homogenization explains why an affine constant changes the D6 period.
    X0, X1, X2, X3 = symbols("X0 X1 X2 X3", nonzero=True)
    X = Matrix([X0, X1, X2, X3])
    c = symbols("c")
    delta_hom = c * X0**2
    check(diff(delta_hom, X0) == 2 * c * X0)
    check(diff(delta_hom, X0, X0) == 2 * c)
    check(all(diff(delta_hom, v) == 0 for v in X[1:]))

    # Algebraic identity N_AB X^B = F_AB X^B for symmetric F_AB=U+iV.
    # Verify for fully symbolic symmetric real U,V, not a numerical sample.
    def symmetric(prefix):
        entries = symbols(" ".join(f"{prefix}{i}{j}" for i in range(4) for j in range(i, 4)), real=True)
        data = dict(zip(((i, j) for i in range(4) for j in range(i, 4)), entries))
        return Matrix(4, 4, lambda i, j: data[min(i, j), max(i, j)])
    real, imag = symmetric("u"), symmetric("v")
    denom = (X.T * imag * X)[0]
    period_matrix = real - I * imag + 2 * I * (imag * X) * (X.T * imag) / denom
    check(period_matrix == period_matrix.T)
    check(((period_matrix - real - I * imag) * X).applyfunc(simplify) == zeros(4, 1))

    # Common real quadratic shifts leave K unchanged but shift Re N.
    shift = symmetric("a")
    check(period_matrix + shift == (real + shift) - I * imag + 2 * I * (imag * X) * (X.T * imag) / denom)
    Xbar = Matrix(symbols("xb0 xb1 xb2 xb3"))
    check(simplify((Xbar.T * shift * X - X.T * shift * Xbar)[0]) == 0)
    cim = symbols("cim", real=True)
    check(simplify(I * (Xbar[0] * 2 * I * cim * X[0] - X[0] * (-2 * I * cim) * Xbar[0]) + 4 * cim * X[0] * Xbar[0]) == 0)
    print(f"{count} exact threshold and special-geometry checks passed.")


if __name__ == "__main__":
    main()

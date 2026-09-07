#!/usr/bin/env python3
"""Exact local lattice and moment-map identities for the E1/E4 extension.

Run with SymPy 1.14.0, e.g. sage -python checks/check_e1_e4_dictionary.py.
The geometric local-base and Higgs-orbit descriptions are cited inputs.
These checks do not prove global smoothability or determine a compact gauge group.
"""

from sympy import Abs, I, Matrix, Poly, conjugate, diag, eye, prod, simplify, sqrt, symbols, zeros
from sympy.matrices.normalforms import smith_normal_form
from sympy import ZZ


def main():
    count = 0

    def check(statement):
        nonlocal count
        assert statement
        count += 1

    # Mark S5 by h,e1,e2,e3,e4; columns below are its four simple roots.
    form = diag(1, -1, -1, -1, -1)
    canonical = Matrix([-3, 1, 1, 1, 1])
    roots = Matrix([[0, 0, 0, 1], [1, 0, 0, -1], [-1, 1, 0, -1],
                    [0, -1, 1, -1], [0, 0, -1, 0]])
    pencils = Matrix([[1, 1, 1, 1, 2], [-1, 0, 0, 0, -1],
                      [0, -1, 0, 0, -1], [0, 0, -1, 0, -1],
                      [0, 0, 0, -1, -1]])
    incidence = Matrix(5, 4, lambda i, j: int(i == j) - int(i == j + 1))
    check(canonical.T * form * roots == zeros(1, 4))
    check(pencils * Matrix([1] * 5) == -2 * canonical)
    check(pencils.T * form * roots == incidence)
    check(roots.T * form * roots == -incidence.T * incidence)
    check((roots.T * form * roots).det() == 5)
    check([abs(d) for d in smith_normal_form(incidence, domain=ZZ).diagonal()] == [1] * 4)
    check(incidence[:4, :].det() == 1)
    c = Matrix(symbols('c1:5'))
    kappa = roots * c
    z = incidence * c
    check(z == Matrix([c[0], c[1]-c[0], c[2]-c[1], c[3]-c[2], -c[3]]))
    inverse = Matrix([-z[4]] + [z[i] + z[4] for i in range(4)])
    check(inverse == kappa)
    check(simplify((kappa.T * form * kappa + z.T * z)[0]) == 0)
    check(simplify(prod(z) + c[0]*c[3]*(c[1]-c[0])*(c[2]-c[1])*(c[3]-c[2])) == 0)

    # Trace pairing equals divisor pairing; trace-zero generator may be fractional.
    divisor = Matrix(symbols('d0:5'))
    m = roots.T * form * divisor
    tau = Matrix(list(incidence.T.col_join(Matrix([[1]*5])).inv() * m.col_join(Matrix([0]))))
    check(incidence.T * tau == m)
    check(sum(tau) == 0)
    check(simplify((tau.T * z - divisor.T * form * kappa)[0]) == 0)
    check(all(all(coefficient.is_Integer for coefficient in
                  Poly(simplify(5*v), *divisor).coeffs()) for v in tau))
    central = symbols('central')
    check(simplify(((tau + Matrix([central]*5)).T * z - tau.T*z)[0]) == 0)

    # E1: the primitive ruling difference is the negative A1 root.
    ruling_form = Matrix([[0, 1], [1, 0]])
    root = Matrix([1, -1])
    check((root.T*ruling_form*root)[0] == -2)
    check((Matrix([-2, -2]).T*ruling_form*root)[0] == 0)
    u, a, b = symbols('u a b')
    check(simplify((Matrix([a,b]).T*ruling_form*(u*root))[0] - u*(b-a)) == 0)

    # Combined ADHM and maximal-flavor-torus weights force equal exponents pairwise.
    for n in (2, 3, 5):
        w = Matrix([[1]*n] + [[int(i == j)-int(i == n-1) for i in range(n)] for j in range(n-1)])
        check(abs(w.det()) == n)
        check(w.rank() == n)
        check(w.row_join(-w) * eye(n).col_join(eye(n)) == zeros(n))

    # Exact balanced representatives, including zero coordinates and complex phases.
    samples = ([1, -1], [0, 0], [1, 1, -2], [1, 1, 1, 1, -4],
               [1+I, -2, 3*I, 0, 1-4*I])
    for values in samples:
        check(sum(values) == 0)
        pairs = [(v/sqrt(Abs(v)), sqrt(Abs(v))) if v != 0 else (0, 0) for v in values]
        check(all(simplify(q*t-v) == 0 for (q,t),v in zip(pairs,values)))
        check(all(simplify(q*conjugate(q)-t*conjugate(t)) == 0 for q,t in pairs))
    check(prod(samples[3]) != 0)  # equal diagonal values need not lie on the discriminant
    check(prod(samples[4]) == 0)
    print(f'{count} exact E1/E4 lattice, charge and moment-map checks passed.')


if __name__ == '__main__':
    main()

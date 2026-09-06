#!/usr/bin/env python3
"""The smooth side X_t of the X_9 transition: surviving classes, full cubic form, c_2.

Input: the triple intersections kappa_{IJK} of the resolution Xhat_9 on the basis
(D2, D5, D6, D7, D8, E) and c_2 . D_I (checks/check_x9_prepotential.sage), and the
local geometry of the pentagon cone: its Milnor fibre is V_7 minus a hyperplane
section S_7, with V_7 = Bl_p P^3 embedded by |2H - Q| (H^3 = 1, Q^3 = 1, H.Q = 0,
c_2(V_7).H = 6, c_2(V_7).Q = 0, c_1(V_7) = 4H - 2Q), so that H^2(V_7) -> H^2(S_7)
has image R^perp = span(f_1 + f_2, e), the classes orthogonal to the A_1 root.

What this script establishes, exactly:
  (1) the classes of Xhat_9 neutral under the two nodal charges and the root are
      span(D0, D1, D8, E); modulo the Coulomb class E they are 3-dimensional,
      which is h^{1,1}(X_t);
  (2) products involving at least one inherited class (D0 or D1) are unchanged
      by the transition, since inherited classes are compactly supported away
      from the singular points; so all of the cubic form and c_2 of X_t except
      nu^3 and c_2.nu are read off Xhat_9 with nu = D8;
  (3) the local correction: Z = F_7 cup (-N_E) has cubic form Phi^3 on V_7 and
      p_1(Z).zeta = p_1(V_7).Phi for the lift Phi of D|_E in R^perp; for
      D8|_E = e the lift is Q, so nu^3 = D8^3 + 1 = 9 and c_2.nu = c_2.D8 - 2 = -6;
      the result is independent of the representative D8 + tE;
  (4) these are the invariants of a smooth P^2 with normal bundle O(-3):
      S^3 = K_S^2 = 9, c_2(X).S = c_2(S) - K_S^2 = -6;
  (5) Riemann-Roch integrality holds on the lattice spanned by D0, D1, nu.
"""
from fractions import Fraction as Fr
from itertools import combinations_with_replacement, product
from collections import Counter
from math import factorial

names = ['D2', 'D5', 'D6', 'D7', 'D8', 'E']
kap = {('D2','D2','D2'):8,('D2','D2','D5'):-2,('D2','D2','D6'):-3,('D2','D2','D7'):-2,('D2','D5','D6'):1,
       ('D2','D6','D6'):1,('D2','D6','D7'):1,('D5','D6','E'):1,('D5','E','E'):-2,('D6','D6','D6'):-1,
       ('D6','D6','D7'):-1,('D6','D7','D7'):-1,('D6','E','E'):-2,('D7','D7','D7'):-3,('D7','D7','D8'):4,
       ('D7','D8','D8'):-6,('D8','D8','D8'):8,('D8','D8','E'):-1,('D8','E','E'):-1,('E','E','E'):7}
K = {}
for (x, y, z), v in kap.items():
    for p in {(x,y,z),(x,z,y),(y,x,z),(y,z,x),(z,x,y),(z,y,x)}:
        K[p] = v
def tri(u, v, w):
    return sum(u[i]*v[j]*w[k]*K.get((names[i], names[j], names[k]), 0)
               for i in range(6) for j in range(6) for k in range(6))
def cls(**kw):
    return [Fr(kw.get(n, 0)) for n in names]
c2 = [-4, 24, 26, 18, -4, -2]
c2dot = lambda D: sum(a*b for a, b in zip(c2, D))

D2, D5, D6, D7, D8, E = (cls(D2=1), cls(D5=1), cls(D6=1), cls(D7=1), cls(D8=1), cls(E=1))
D3, D4 = cls(D5=1, D7=-1, D8=-1), cls(D6=1, D7=-1, D8=-1)
D0 = cls(D2=1, D5=4, D6=4, D7=-3, D8=-2, E=2)
D1 = cls(D2=1, D5=2, D6=2, D7=-1, D8=-1, E=1)

# (1) charges: C1 = D6.D7, C2 = D2.D4, root R = R_{v4} - R_{v3} = E.D4 - E.D3
def charges(D):
    return (tri(D, D6, D7), tri(D, D2, D4), tri(D, E, D4) - tri(D, E, D3))
assert charges(D2) == (1, -1, 0) and charges(D5) == (0, 1, 1) and charges(D6) == (-1, 0, -1)
assert charges(D7) == (-1, 1, 0) and charges(D8) == (0, 0, 0) and charges(E) == (0, 0, 0)
assert charges(D0) == (0, 0, 0) and charges(D1) == (0, 0, 0)
# the neutral subspace is 4-dimensional: rank of the 6x3 charge matrix is 2
rows = [charges(cls(**{n: 1})) for n in names]
import itertools
def rank(M):
    M = [list(map(Fr, r)) for r in M]; r = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]; M[i] = [a - f*b for a, b in zip(M[i], M[r])]
        r += 1
    return r
assert rank(rows) == 2
assert rank([D0, D1, D8, E]) == 4
print("(1) neutral classes = span(D0, D1, D8, E); modulo E: h^{1,1}(X_t) = 3")

# (2) products with an inherited class, and E-shift independence
print("(2) cubic form of X_t on (D0, D1, nu), nu = extension of D8:")
cub = {}
B = {'D0': D0, 'D1': D1, 'nu': D8}
for a, b, c in combinations_with_replacement(['D0', 'D1', 'nu'], 3):
    cub[(a, b, c)] = tri(B[a], B[b], B[c])
assert all(tri(D0, E, X) == 0 and tri(D1, E, X) == 0 for X in [D0, D1, D8, E])   # plane . E . anything = 0
# (3) local correction: D8|_E = e (E.D8.D8 = e^2 = -1, E.E.D8 = K.e = -1), lift Phi = Q
assert tri(E, D8, D8) == -1 and tri(E, E, D8) == -1
V7 = {'H3': 1, 'Q3': 1, 'c2H': 6, 'c2Q': 0}
def local_cubic(a, c):      # Phi = aH + cQ
    return a**3 * V7['H3'] + c**3 * V7['Q3']
def local_p1(a, c):         # p_1(V_7).Phi = c_1^2.Phi - 2 c_2.Phi, c_1 = 4H - 2Q
    return a * (16 - 2*V7['c2H']) + c * (4 - 2*V7['c2Q'])
# lift of e + tK: K = -2(f1+f2) + e  ->  Phi = -2t H + (1+t) Q
for t in range(-3, 4):
    Dt = [D8[i] + t*E[i] for i in range(6)]
    assert tri(Dt, Dt, Dt) + local_cubic(-2*t, 1 + t) == 9
    assert c2dot(Dt) - Fr(1, 2) * local_p1(-2*t, 1 + t) == -6
cub[('nu', 'nu', 'nu')] = tri(D8, D8, D8) + local_cubic(0, 1)
c2t = {'D0': c2dot(D0), 'D1': c2dot(D1), 'nu': c2dot(D8) - Fr(1, 2) * local_p1(0, 1)}
expected = {('D0','D0','D0'): 83, ('D0','D0','D1'): 49, ('D0','D1','D1'): 27, ('D1','D1','D1'): 14,
            ('D0','D0','nu'): 0, ('D0','D1','nu'): 0, ('D1','D1','nu'): 1,
            ('D0','nu','nu'): 0, ('D1','nu','nu'): -3, ('nu','nu','nu'): 9}
assert cub == expected, cub
assert c2t == {'D0': 146, 'D1': 80, 'nu': -6}, c2t
for k, v in cub.items():
    print("    ", " ".join(k), "=", v)
print("     c2 . (D0, D1, nu) =", c2t['D0'], c2t['D1'], c2t['nu'])
print("(3) local correction on V_7 = Bl_p P^3: nu^3 = 8 + 1 = 9, c2.nu = -4 - 2 = -6, independent of D8 -> D8 + tE")
# (4) P^2 with O(-3)
assert cub[('nu','nu','nu')] == 9 and c2t['nu'] == 3 - 9
print("(4) nu has the invariants of a rigid P^2: S^3 = 9, c2.S = 3 - 9 = -6")
# (5) Riemann-Roch integrality: 2 D^3 + c2.D = 0 mod 12 on the lattice
def cubic(x, y, z):
    v = {'D0': x, 'D1': y, 'nu': z}; tot = 0
    for (a, b, c), val in cub.items():
        mult = 6
        for m in Counter((a, b, c)).values():
            mult //= factorial(m)
        tot += mult * val * v[a] * v[b] * v[c]
    return tot
assert all((2*cubic(x, y, z) + c2t['D0']*x + c2t['D1']*y + c2t['nu']*z) % 12 == 0
           for x, y, z in product(range(-4, 5), repeat=3))
print("(5) chi(O(D)) = D^3/6 + c2.D/12 is an integer on the lattice span(D0, D1, nu)")
print("ALL CHECKS PASSED")

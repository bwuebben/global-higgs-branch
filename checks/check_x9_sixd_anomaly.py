#!/usr/bin/env python3
"""Six-dimensional anomaly bookkeeping for the genus-one fibration of X_9.

Verified geometric input (checks/x9_fibrations.sage, checks/x9_weierstrass.sage):
  base P^2; two (4,6,12) points q_3, q_4 on the line L_0; an I_2 locus on L_0
  with three type-III and one type-IV point; 76 isolated I_2 fibres off L_0;
  a rank-one Mordell-Weil group; h^{2,1} = 122.
The six-dimensional model lives on B' = Bl_{q_3,q_4} P^2 (T = 2), with su(2) on
the proper transform L~ = H - E_3 - E_4 of L_0 and a u(1) from the section.
Conventions: Kumar-Morrison-Taylor and Park-Taylor; for su(2): lambda = 1,
A_fund = 1, C_fund = 1/2, A_adj = 4, C_adj = 8; a hyper in a representation R
of su(2) is counted dim R times in the gravitational and abelian sums.

What this script establishes, exactly:
  (1) su(2) on L~ carries ten fundamental hypermultiplets (a.b and b.b);
  (2) H - V + 29T = 273 forces 96 charged hypermultiplets: 20 from the
      doublets and 76 su(2)-singlets, matching the 76 isolated I_2 fibres;
  (3) the three abelian conditions have a UNIQUE solution up to the overall
      normalization of the u(1): doublets of charge 1, singlets of charge 2,
      height-pairing class b = 22 H - 6 E_3 - 6 E_4 (equivalently charges
      1/2 and 1 with b = (11/2) H - (3/2)(E_3 + E_4)).  Hence the gauge group
      is U(2) = (SU(2) x U(1))/Z_2, and b.H = 4 * (11/2), where 11/2 is the
      height pairing of the section E against D_2 on the resolved threefold
      including the fibral correction from the I_2 locus over L_0 (E meets
      D_8, D_2 meets D_7), computed here from the triple intersections.
"""
from fractions import Fraction as Fr
from itertools import product, combinations_with_replacement

def dot(u, v):          # classes h*H - e3*E3 - e4*E4 encoded as (h, e3, e4)
    return u[0]*v[0] - u[1]*v[1] - u[2]*v[2]

a = (-3, -1, -1)        # K_{B'} = -3H + E_3 + E_4
L = ( 1,  1,  1)        # proper transform of L_0
T, V, h21 = 2, 3 + 1, 122
assert dot(a, a) == 9 - T and dot(a, L) == -1 and dot(L, L) == -1

# (1) su(2) on a (-1)-curve
A_fund, C_fund, A_adj, C_adj = 1, Fr(1, 2), 4, 8
n_f = A_adj - 6 * dot(a, L)
assert n_f == 10 and Fr(1, 3) * (n_f * C_fund - C_adj) == dot(L, L)
print(f"(1) su(2) on L~: {n_f} fundamental hypermultiplets")

# (2) gravitational anomaly
H = 273 - 29 * T + V
H_neutral = h21 + 1
H_charged = H - H_neutral
assert (H, H_charged) == (219, 96) and 2 * n_f + 76 == H_charged
print(f"(2) H = {H} = {H_neutral} neutral + 2*{n_f} (doublets) + 76 (singlets at the isolated I_2 fibres)")

# (3) abelian conditions, solved exactly for every charge assignment.
#     Given the charges, L.b, -a.b and b.b are fixed; with b = (h, c3, c4):
#     h - c3 - c4 = L.b,  3h - c3 - c4 = -a.b,  h^2 - c3^2 - c4^2 = b.b.
charges = [Fr(k, 2) for k in range(0, 5)]            # 0, 1/2, ..., 2
sols = set()
for d in product(range(n_f + 1), repeat=len(charges)):     # doublet charge multiplicities
    if sum(d) != n_f:
        continue
    Lb = sum(n * q * q for n, q in zip(d, charges))          # b_su2 . b = sum_R x_R A_R q_R^2
    dq2 = 2 * Lb
    dq4 = 2 * sum(n * q**4 for n, q in zip(d, charges))
    for qa, qb in combinations_with_replacement(charges[1:], 2):  # two singlet charge values
        for s in range(0, 77):                                    # s singlets of charge qa, 76-s of qb
            sq2 = s * qa**2 + (76 - s) * qb**2
            sq4 = s * qa**4 + (76 - s) * qb**4
            m = Fr(sq2 + dq2, 6)                                  # = -a.b
            bb = Fr(sq4 + dq4, 3)                                 # = b.b
            h = (m - Lb) / 2
            S = h - Lb                                            # c3 + c4
            P = h * h - bb                                        # c3^2 + c4^2
            disc = 2 * P - S * S                                  # (c3 - c4)^2
            if disc < 0:
                continue
            from math import isqrt
            num, den = disc.numerator, disc.denominator
            rn, rd = isqrt(num), isqrt(den)
            if rn * rn != num or rd * rd != den:
                continue
            root = Fr(rn, rd)
            for c3 in {(S + root) / 2, (S - root) / 2}:
                c4 = S - c3
                if 2 * h % 1 or 2 * c3 % 1 or 2 * c4 % 1:      # classes in (1/2) Z
                    continue
                key = (h, min(c3, c4), max(c3, c4), tuple(d), (qa, qb, s))
                sols.add(key)
# canonical form: (b, doublet multiplicities by charge, singlet multiplicities by charge)
from collections import Counter
canon = set()
for (h, c3, c4, d, (qa, qb, s)) in sols:
    cnt = Counter(); cnt[qa] += s; cnt[qb] += 76 - s
    singlets = tuple((q, n) for q, n in sorted(cnt.items()) if n > 0)
    canon.add((h, c3, c4, d, singlets))
print(f"(3) abelian anomaly solutions with charges in (1/2)Z, at most 2, classes in (1/2)Z: {len(canon)}")
for sol in sorted(canon):
    print("    ", sol)
expected = {
    (Fr(11, 2), Fr(3, 2), Fr(3, 2), (0, 10, 0, 0, 0), ((Fr(1), 76),)),
    (Fr(22),    Fr(6),    Fr(6),    (0, 0, 10, 0, 0), ((Fr(2), 76),)),
}
assert canon == expected, "solution set differs from the unique normalized solution"
print("    unique up to normalization: doublets charge 1, singlets charge 2, b = 22H - 6E_3 - 6E_4")
print("    (-1, e^{i pi}) acts trivially: gauge group U(2)")

# Shioda-corrected height pairing on the resolved threefold, from the triple intersections
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
w  = cls(D5=1, D6=1, D7=-1, D8=-1)      # w = D_4 + D_5, the genus-one fibration class
E, D2, D7, D8 = cls(E=1), cls(D2=1), cls(D7=1), cls(D8=1)
assert tri(w, w, w) == 0
assert (tri(E, D7, w), tri(E, D8, w), tri(D2, D7, w), tri(D2, D8, w)) == (0, 1, 1, 0)
naive = [E[i] - D2[i] for i in range(6)]
assert tri(naive, naive, w) == -6
corr = [naive[i] - Fr(1, 2) * D7[i] for i in range(6)]     # orthogonal to the fibral component D_7
assert tri(corr, D7, w) == 0
height = -tri(corr, corr, w)
assert height == Fr(11, 2)
print(f"(4) Shioda-corrected height pairing on the resolved threefold: {height}; 4 * {height} = {4*height} = b.H")
print("ALL CHECKS PASSED")

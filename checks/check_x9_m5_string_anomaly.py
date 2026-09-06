#!/usr/bin/env python3
"""Wrapped-M5 anomaly across the X_9 F_1 -> P^2 surface transition.

Standard-library exact arithmetic. Run from the project root:
    python3 checks/check_x9_m5_string_anomaly.py

The independent surface input is Pic(F_1) = <h,e>, intersection diag(1,-1),
K=-3h+e, and Pic(P^2)=<h>, h^2=1, K=-3h. The restriction of (D0,D1,D8)
is (0,h,K). Compare this worldvolume lattice with the bulk cubic already
computed by the smooth-side check. Normalizations follow eq. (2.21) of
Katz--Kim--Tarazi--Vafa, arXiv:2004.14401. The SU(2) is transverse rotation,
not an assumed infrared R-symmetry of a tensionless local string.
"""
import contextlib
import io
from fractions import Fraction as Fr
with contextlib.redirect_stdout(io.StringIO()):
    import check_x9_smooth_side as bulk


def gram(vectors, signature):
    return [[sum(a*b*s for a,b,s in zip(u,v,signature))
             for v in vectors] for u in vectors]


restr_before = [(0,0),(1,0),(-3,1)]
restr_after = [(0,),(1,),(-3,)]
k_before = gram(restr_before, (1,-1))
k_after = gram(restr_after, (1,))
assert k_before == [[0,0,0],[0,1,-3],[0,-3,8]]
assert k_after == [[0,0,0],[0,1,-3],[0,-3,9]]
divisors = [bulk.D0, bulk.D1, bulk.D8]
assert k_before == [[bulk.tri(bulk.D8,u,v) for v in divisors] for u in divisors]
names = ['D0','D1','nu']
assert k_after == [[bulk.cub[tuple(sorted(('nu',a,b)))] for b in names] for a in names]
assert [bulk.tri(bulk.D8,bulk.E,d) for d in divisors] == [0,0,-1]
assert bulk.tri(bulk.D8,bulk.E,bulk.E) == -1
assert bulk.tri(bulk.D8,bulk.D1,bulk.D1) == 1
delta_k = [[b-a for a,b in zip(ra,rb)] for ra,rb in zip(k_before,k_after)]
assert delta_k == [[0,0,0],[0,0,0],[0,0,1]]


def anomaly_coefficients(canonical_square, euler):
    # Adjunction in a CY3: c2(X).S = chi_top(S) - K_S^2.
    c2dot = euler-canonical_square
    # Coefficients of (c2(N_SU2), p1(TSigma)); the gauge part is -k_AB F^A F^B/2.
    return (-Fr(2*canonical_square+c2dot,12), Fr(c2dot,48))


before = anomaly_coefficients(8,4)
after = anomaly_coefficients(9,3)
assert bulk.c2dot(bulk.D8) == 4-8 == -4
assert bulk.c2t['nu'] == 3-9 == -6
assert before == (-1, -Fr(1,12))
assert after == (-1, -Fr(1,8))
delta = tuple(b-a for a,b in zip(before,after))
assert delta == (0,-Fr(1,24))

# Tensor zero modes on rational rigid surfaces: b2+=1, b2-=rho-1;
# three translation bosons, four real right-moving fermions, no left fermions.
# These count semiclassical zero modes, NOT a claimed interacting IR central charge.
def zero_mode_central_charges(b2_minus):
    return (3+b2_minus, 3+1+4*Fr(1,2))


assert zero_mode_central_charges(1) == (4,6)
assert zero_mode_central_charges(0) == (3,6)
for bminus, coeffs in [(1,before),(0,after)]:
    cl,cr=zero_mode_central_charges(bminus)
    assert Fr(cl-cr,24)==coeffs[1]

# Lost negative lattice direction has degree -1 under nu and zero under D0,D1.
# Its left-moving anomaly is +(F^nu)^2/2 + p1(TSigma)/24 in these conventions.
lost_charge=(0,0,-1)
lost_gauge=[[Fr(a*b,2) for b in lost_charge] for a in lost_charge]
assert [[-Fr(x,2) for x in row] for row in delta_k] == [
    [-x for x in row] for row in lost_gauge]
assert delta[1] == -Fr(1,24)
print('M5 string: surface lattice F1 diag(1,-1) -> P2 (1)')
print('k_AB before:', k_before)
print('k_AB after: ', k_after)
print('I4 coefficients (c2(N_SU2), p1(TSigma)):', before, '->', after)
print('Delta I4 = -(F^nu)^2/2 - p1(TSigma)/24; no transverse SU(2) jump.')
print('Exactly minus the anomaly of one lost unit-charge left-moving lattice boson.')
print('ALL CHECKS PASSED (no full BPS index or tensionless IR spectrum asserted).')

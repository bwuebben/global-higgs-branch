"""Integral compact charges, E2 center lift, and E0 screening for X9.

The resolved fan is reconstructed by the prepotential check. Integer
matrices then verify the printed divisor/curve duality and local gluing.
The smooth cubic is also checked against an independently implemented
surface correction. No BPS degeneracy or finite-Planck Higgs metric is
computed here.

Run: sage -c "load('checks/check_x9_integral_screening.sage')"
"""

import contextlib
import io
from itertools import combinations_with_replacement, product

with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_prepotential.sage')

# 1. Actual integral compact curve classes pairing unimodularly with the
# six divisors. This proves equality with H^2/tors, not just toric saturation.
witness_curves = [(6,7),(2,4),(4,5),(8,9),(3,7),(0,5)]
pairing = matrix(ZZ, [curve(i,j) for i,j in witness_curves])
assert pairing == matrix(ZZ, [[1,0,-1,-1,0,0],[-1,1,0,1,0,0],
    [1,0,0,0,0,1],[0,0,0,0,-1,-1],[0,0,1,-1,2,0],[2,0,3,0,0,0]])
assert abs(pairing.det()) == 1
assert pairing.inverse().base_ring() == ZZ or all(x in ZZ for x in pairing.inverse().list())
all_pairings = matrix(ZZ, [cv for cv in eff.values()])
assert all_pairings.smith_form()[0].diagonal() == [1]*6

# 2. Restriction to E in the divisor basis (h,e1,e2), intersection diag(1,-1,-1).
res = matrix(ZZ, [[0,1,1,0,1,-3],[0,-1,0,0,-1,1],[0,0,-1,0,-1,1]])
surface_form = diagonal_matrix(ZZ, [1,-1,-1])
assert res.transpose()*surface_form*res == L(CLS['E'])
assert res.smith_form()[0].diagonal() == [1,1,1]
assert res*CLS['D4'] == vector(ZZ,[0,1,0])
assert res*CLS['D3'] == vector(ZZ,[0,0,1])
assert res*CLS['D8'] == vector(ZZ,[1,-1,-1])

# On gamma=d*h+m1*e1+m2*e2 the integral Coulomb, SU2-weight and auxiliary
# Abelian charges are (3d+m1+m2,m1-m2,d). Their index is two, with the
# unique parity condition q_g+w+q_a even. We claim a faithful action on
# this M2 charge lattice, not a determination of every UV flavor observable.
local_charges = matrix(ZZ, [[3,1,1],[0,1,-1],[1,0,0]])
assert local_charges.smith_form()[0].diagonal() == [1,1,2]
assert abs(local_charges.det()) == 2
assert all(sum(local_charges.column(i)) % 2 == 0 for i in range(3))
cochars = (res.transpose()*surface_form*local_charges.inverse()).transpose()
assert cochars == matrix(QQ, [[0,1/2,1/2,0,1,-1],
                             [0,1/2,-1/2,0,0,0],
                             [0,-1/2,-1/2,0,-2,0]])
assert cochars*(-CLS['D6']) == vector(QQ,[-1/2,1/2,1/2])
assert cochars*CLS['D5'] == vector(QQ,[1/2,1/2,-1/2])
assert 2*cochars.row(1) == R
for column in cochars.columns():
    # A loop is integral, or has simultaneous half-turns in all three factors.
    assert len({q-floor(q) for q in column}) == 1
assert cochars.transpose()*local_charges == res.transpose()*surface_form

# The geometric nilpotent column differs from the U(1)_a current coefficient
# by -3/2 of the reduced root column. The local invariant ring permits the
# corresponding shear beta_geom=J0+3U/2, alpha_geom=U.
alpha_geom = vector(QQ,[0,-2,1,0,-2,0])
assert alpha_geom == cochars.row(2)-3/2*R
Pr = PolynomialRing(QQ,'J0,U'); Jzero, Unil = Pr.gens()
Inil = Pr.ideal([Unil^2,Unil*Jzero])
assert Pr.ideal([Unil^2,Unil*(Jzero+3/2*Unil)]) == Inil

# 3. Integral Mayer--Vietoris data. H|E=2h-e1-e2 and Q|E=h-e1-e2
# form an integral basis of R-perp (not an index-two sublattice).
v7_restriction = matrix(ZZ, [[2,1],[-1,-1],[-1,-1]])
assert matrix(ZZ, [[2,1],[-1,-1]]).det() == -1
Ksurf = vector(ZZ,[-3,1,1])
assert v7_restriction*vector(ZZ,[-2,1]) == Ksurf
assert gcd(Ksurf) == gcd(vector(ZZ,[2,-1])) == 1

Aclass = CLS['D2']+CLS['D7']
Bclass = CLS['D5']+CLS['D6']-CLS['D7']
neutral = matrix(ZZ,[Aclass,Bclass,CLS['D8'],CLS['E']]).transpose()
extended = matrix(ZZ,[-CLS['D6'],CLS['D5'],Aclass,Bclass,CLS['D8'],CLS['E']]).transpose()
assert abs(extended.det()) == 1
assert Q*neutral == zero_matrix(QQ,3,4)
assert CLS['D0'] == Aclass+4*Bclass-2*CLS['D8']+2*CLS['E']
assert CLS['D1'] == Aclass+2*Bclass-CLS['D8']+CLS['E']
old_in_new = matrix(ZZ, [[1,1,0],[4,2,0],[-2,-1,1]])
assert abs(old_in_new.det()) == 2

# 4. Primitive-basis cubic: add the V7 correction to the resolved tensor.
# A|E=0, B|E=H, nu|E=Q, so Delta cubic=B^3+nu^3 and Delta c2=(0,-2,-2).
pb = [Aclass,Bclass,CLS['D8']]
expected = {(0,0,0):-1,(0,0,1):1,(0,0,2):4,(0,1,1):-1,(0,1,2):-4,
            (0,2,2):-6,(1,1,1):3,(1,1,2):4,(1,2,2):6,(2,2,2):9}
smooth_tensor = {}
for ijk in combinations_with_replacement(range(3),3):
    correction = 1 if ijk in [(1,1,1),(2,2,2)] else 0
    smooth_tensor[ijk] = trip(*(pb[i] for i in ijk))+correction
assert smooth_tensor == expected
c2_primitive = vector(ZZ,[sum(c2D[j]*pb[i][j] for j in range(6))
                          for i in range(3)])+vector(ZZ,[0,-2,-2])
assert c2_primitive == vector(ZZ,[14,30,-6])

def primitive_trip(u,v,w):
    return sum(u[i]*v[j]*w[k]*smooth_tensor[tuple(sorted((i,j,k)))]
               for i,j,k in product(range(3),repeat=3))
for a,b,c in product(range(12),repeat=3):
    vec=vector(ZZ,[a,b,c])
    assert (2*primitive_trip(vec,vec,vec)+c2_primitive*vec) % 12 == 0
assert [primitive_trip(old_in_new.column(i),old_in_new.column(j),old_in_new.column(k))
        for i,j,k in [(0,0,0),(0,0,1),(0,1,1),(1,1,1),(1,1,2),(1,2,2),(2,2,2)]] == [83,49,27,14,1,-3,9]

# 5. The local P2 line has nu charge -3, whereas compact D1^2 has nu charge 1.
# A divisor D restricts to k*h on P2, so D*nu^2=-3*k.
line_charges = vector(ZZ,[-smooth_tensor[(i,2,2)]/3 for i in range(3)])
assert line_charges == vector(ZZ,[2,-2,-3])
d1 = old_in_new.column(1)
screening_charge = vector(ZZ,[primitive_trip(vector(ZZ,[1 if j==i else 0 for j in range(3)]),d1,d1)
                              for i in range(3)])
assert screening_charge == vector(ZZ,[1,7,1])
assert gcd(abs(line_charges[2]),abs(screening_charge[2])) == 1

print('Resolved integral divisor/curve pairing det =',pairing.det())
print(pairing)
print('E2 local charge lattice index two: q_g+w+q_a=0 mod 2.')
print('V1,V2 cocharacters:',list(cochars*(-CLS['D6'])),list(cochars*CLS['D5']))
print('Primitive smooth basis A=D2+D7, B=D5+D6-D7, nu; old basis index two.')
print('Primitive cubic entries:',smooth_tensor,'; c2:',list(c2_primitive))
print('Local E0 line charge:',list(line_charges),'; compact D1^2 charge:',list(screening_charge))
print('Local Z3 electric line charge is screened in the compact model.')
print('ALL INTEGRAL AND SCREENING CHECKS PASSED; no BPS multiplicities asserted.')

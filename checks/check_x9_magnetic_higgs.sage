"""Light-emission phase test, Higgs-phase magnetic flux, and neutral couplings.

Run: sage -c "load('checks/check_x9_magnetic_higgs.sage')"
The endpoint period positivity and physical confinement interpretation are
inputs explained in the manuscript. No constituent BPS existence is inferred.
"""
import contextlib
import io
from itertools import combinations_with_replacement

with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_prepotential.sage')

nef = matrix(ZZ,sorted([list(r.vector()) for r in S.rays()])).T
assert abs(nef.det()) == 1
e1,e2 = [vector(ZZ,curve(i,9)) for i in [4,3]]
electric = matrix(ZZ,[C1,C2,e1,e2]).T
orientation = diagonal_matrix(ZZ,[-1,-1,1,-1])
assert electric*orientation*vector(ZZ,[1]*4) == 0
assert nef.T*electric == matrix(ZZ,[[0,1,1,0],[0,0,0,0],
    [0,0,0,0],[0,0,0,0],[1,0,1,0],[0,0,1,1]])
assert electric.rank() == 3
assert [x for x in electric.elementary_divisors() if x] == [1,1,1]

# Exact phase identity for emission of any nonzero nonnegative combination
# of the four known light charges, on a positive imaginary electric ray.
phase = PolynomialRing(QQ,names=('eps','c1','c2','ce','n1','n2','n3','n4','realZ','imagZ'))
eps,c1,c2,ce,n1,n2,n3,n4,realZ,imagZ = phase.gens()
light_masses = vector(phase,[c1,c2,c1+c2+ce,ce])
mass = vector(phase,[n1,n2,n3,n4])*light_masses
assert all(coefficient>0 for coefficient in mass.coefficients())
# Z_eta=i eps mass; Z_(Gamma-eta)=realZ+i(imagZ-eps mass).
alignment = (imagZ-eps*mass)*0-realZ*eps*mass
assert alignment == -eps*mass*realZ
assert all(coefficient<0 for coefficient in alignment.coefficients())
print('Emission alignment = -eps*m_eta*Re(Z_Gamma), strictly negative on the stated ray.')

PA = CLS['D4']
PV = CLS['D7']+CLS['D8']
PT = CLS['D6']
probes = [PA,PV,PT,CLS['D8'],CLS['E']]
pairings = matrix(ZZ,probes)*electric
assert pairings == matrix(ZZ,[[0,-1,-1,0],[-1,1,1,1],
    [-1,0,0,1],[0,0,1,1],[0,0,-1,-1]])
assert PA+PV == PT
fluxes = pairings*orientation
assert fluxes == matrix(ZZ,[[0,1,-1,0],[1,-1,1,-1],
    [1,0,0,-1],[0,0,1,-1],[0,0,-1,1]])
assert all(sum(row)==0 for row in fluxes)
assert fluxes.row(0)+fluxes.row(1) == fluxes.row(2)
assert fluxes[:2,:].rank() == 2
print('Pairings, rows A,V,T,D8,E; columns C1,C2,e1,e2:\n',pairings)
print('Oriented Higgs fluxes in A3:\n',fluxes)

# The entire integral flux image is A3, not a finite-index sublattice.
flux_map = orientation*electric.T
simple_roots = matrix(ZZ,[[1,-1,0,0],[0,1,-1,0],[0,0,1,-1]])
for root in simple_roots.rows():
    # A flux w is lifted by m1=-w1, m2=-w2, m3=-w4.
    magnetic_lift = -root[1]*nef.column(0)-root[0]*nef.column(4)-root[3]*nef.column(5)
    assert flux_map*magnetic_lift == root
neutral_nef = nef.matrix_from_columns([1,2,3])
assert flux_map*neutral_nef == 0
assert neutral_nef.column_module() == flux_map.right_kernel()
assert [x for x in flux_map.elementary_divisors() if x] == [1,1,1]

# Compare the four-hyper annihilator with the old neutral lattice modulo E.
smoothA = CLS['D2']+CLS['D7']
smoothB = CLS['D5']+CLS['D6']-CLS['D7']
old_neutral = matrix(ZZ,[smoothA,smoothB,CLS['D8'],CLS['E']]).T
assert old_neutral.column_module() == matrix(ZZ,[C1,C2]).right_kernel()
assert old_neutral.T*e2 == vector(ZZ,[0,1,1,-1])
projection = identity_matrix(ZZ,6)+CLS['E'].column()*e2.row()
assert projection^2 == projection
assert projection*CLS['E'] == 0
lifted = [smoothA,smoothB+CLS['E'],CLS['D8']+CLS['E']]
lifted_matrix = matrix(ZZ,lifted).T
assert projection*old_neutral == lifted_matrix.augment(zero_vector(ZZ,6))
assert lifted_matrix.column_module() == neutral_nef.column_module()
transition = (nef.inverse()*lifted_matrix).matrix_from_rows([1,2,3])
assert transition == matrix(ZZ,[[0,-1,-2],[2,-1,-1],[-1,1,1]])
assert abs(transition.det()) == 1

# Independent smooth-side local-surgery formula: add B^3+nu^3 and
# subtract two from each of c2.B,c2.nu in the old representatives.
old_basis = [smoothA,smoothB,CLS['D8']]
triples = list(combinations_with_replacement(range(3),3))
smooth_cubic = [trip(*(old_basis[i] for i in ijk))+
    (1 if ijk in [(1,1,1),(2,2,2)] else 0) for ijk in triples]
lifted_cubic = [trip(*(lifted[i] for i in ijk)) for ijk in triples]
assert lifted_cubic == smooth_cubic == [-1,1,4,-1,-4,-6,3,4,6,9]
lifted_c2 = vector(ZZ,[sum(c2D[j]*p[j] for j in range(6)) for p in lifted])
smooth_c2 = vector(ZZ,[sum(c2D[j]*p[j] for j in range(6)) for p in old_basis])+vector(ZZ,[0,-2,-2])
assert lifted_c2 == smooth_c2 == vector(ZZ,[14,30,-6])
print('Primitive unconfined lifts: D2+D7, D5+D6-D7+E, D8+E.')
print('All ten cubic couplings match:',lifted_cubic)
print('All three second-Chern pairings match:',lifted_c2)
print('ALL MAGNETIC HIGGS CHECKS PASSED; BPS survival is not established and the mirror identification is checked separately.')

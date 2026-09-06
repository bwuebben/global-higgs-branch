"""The middle Newton slice: a candidate K3 for heterotic stable degeneration.

Run: sage checks/check_x9_heterotic_k3_candidate.sage
Checks the integral two-piece subdivision and its reflexive 30-point middle
slice. Prints a characteristic-zero witness with fibers I2 + 22 I1.
This does NOT construct the semistable total space, the two heterotic bundles,
or their massless U(1); those remain separate tasks.
"""
from collections import defaultdict

vertices = [(1,0,0,0),(0,1,0,0),(-1,-1,0,0),(0,0,1,0),(0,0,0,1),
            (-4,-2,-1,0),(-4,-2,0,-1),(3,1,1,1),(2,1,1,1)]
polar = LatticePolytope(vertices).polar()
assert polar.npoints() == 162

# The exponents at v4 and v6 agree precisely on this affine-linear slice.
# The untwisted guess m4=0 misses two fiber monomials and is not reflexive.
delta = lambda m: m[3]+2*m[0]+m[1]
middle = [m for m in polar.points() if delta(m) == 0]
assert len(middle) == 30
assert all(m[3]+1 == 1-4*m[0]-2*m[1]-m[3] for m in middle)
by_fiber = defaultdict(list)
for m in middle:
    by_fiber[tuple(m[:2])].append(int(m[2])+1)
by_fiber = {key:sorted(exps) for key,exps in by_fiber.items()}
expected_support = {
    (1,-1):[0], (0,1):[0], (0,0):list(range(3)), (0,-1):list(range(5)),
    (-1,2):list(range(1,3)), (-1,1):list(range(1,5)),
    (-1,0):list(range(1,7)), (-1,-1):list(range(1,9))}
assert by_fiber == expected_support
middle_poly = LatticePolytope([tuple(m[:3]) for m in middle])
assert middle_poly.is_reflexive()
dual = middle_poly.polar()
assert set(map(tuple,dual.vertices())) == {
    (-4,-2,-1),(0,1,0),(1,0,1),(1,0,0),(-1,-1,0),(0,0,1)}
rho_toric = dual.npoints()-4-sum(len(fc.interior_points()) for fc in dual.facets())
rho_toric += sum(len(ed.interior_points())*len(ed.dual().interior_points()) for ed in dual.edges())
assert rho_toric == 4
newton = Polyhedron(vertices=[tuple(m) for m in polar.points()])
for sign,nv in [(1,12),(-1,11)]:
    half = newton.intersection(Polyhedron(ieqs=[[0,2*sign,sign,0,sign]]))
    assert half.n_vertices() == nv
    assert all(all(x in ZZ for x in v) for v in half.vertices())
print('middle slice: m4+2m1+m2=0; 30 lattice points; reflexive; toric Picard count 4')
print('both half-polytopes integral; 12 and 11 vertices')

# A fixed rational witness, including every supported coefficient.
ring = PolynomialRing(QQ,'s')
s = ring.gen()
coefficients = {
    'a0':4, 'b2':2, 'b1':3+s+5*s**2,
    'b0':7+4*s+4*s**2+9*s**3+5*s**4,
    'c3':7*s+4*s**2,
    'c2':s+8*s**2+5*s**3+6*s**4,
    'c1':2*s+3*s**2+s**3+5*s**4+9*s**5+4*s**6,
    'c0':5*s+9*s**2+8*s**3+6*s**4+5*s**5+s**6+3*s**7+8*s**8}
roles = {(1,-1):'a0',(0,1):'b2',(0,0):'b1',(0,-1):'b0',
         (-1,2):'c3',(-1,1):'c2',(-1,0):'c1',(-1,-1):'c0'}
for key,role in roles.items():
    coefficients[role] = ring(coefficients[role])
    assert sorted(coefficients[role].exponents()) == by_fiber[key]

# Complete the square in a0*w^2+(b0*u^2+b1*u*v+b2*v^2)*w
# +u*(c0*u^3+c1*u^2*v+c2*u*v^2+c3*v^3)=0.
a0,b0,b1,b2,c0,c1,c2,c3 = [coefficients[k] for k in ('a0','b0','b1','b2','c0','c1','c2','c3')]
A = b0**2-4*a0*c0
B = 2*b0*b1-4*a0*c1
C = b1**2+2*b0*b2-4*a0*c2
D = 2*b1*b2-4*a0*c3
E = b2**2
I = 12*A*E-3*B*D+C**2
J = 72*A*C*E+9*B*C*D-27*A*D**2-27*E*B**2-2*C**3
disc = 4*I**3-J**2
assert (I.degree(),J.degree(),disc.degree()) == (8,12,24)
assert (I.valuation(),J.valuation(),disc.valuation()) == (0,0,2)
residual = disc//s**2
assert residual.degree() == 22 and residual(0) != 0
assert residual.gcd(residual.derivative()) == 1
assert residual.gcd(I) == 1 and residual.gcd(J) == 1
# Leading degree 24 guarantees no discriminant zero at infinity.
assert disc[24] != 0
print('rational witness (polynomials in the affine base coordinate s):')
for key in ('a0','b2','b1','b0','c3','c2','c1','c0'):
    print('  %s = %s' % (key,coefficients[key]))
print('degrees (f,g,Delta) = (8,12,24); orders at 0 = (0,0,2)')
print('Delta/s^2: degree 22, squarefree, coprime to f and g; infinity smooth')
print('candidate K3 fibers: I2 + 22 I1, over characteristic zero')
print('ALL CHECKS PASSED; global semistable degeneration and bundle construction not asserted')

"""Integral marking of the compact four-node mirror family.

Run: sage -c "load('checks/check_x9_integral_vanishing_marking.sage')"
The proof uses uniform resummation of the seven electric periods,
Picard--Lefschetz invariance, and the normalized node residue.
Finite coefficient checks below support, but do not replace, that proof.
No regular magnetic periods or bound-state stability are computed.
"""
import contextlib
import io
import os
os.environ.setdefault('MIRROR_DEGREE','2')
with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_mirror_periods.sage')

def factorial_arguments(i,L,A,B,j,k):
    return [B,A,-i-B+j,L+j-k,L-i,i,-j+k,i+2*A-j,-L-2*A+k,-L+A+B-k]

def local_bounds(L,A,B):
    return [max(L,A-L),max(L+2*A,A+B-L),max(2*L+2*A,A+B-L)]

def reciprocal_factorial(n):
    return QQ(0) if n<0 else 1/factorial(n)

# Exact polyhedral certificates for the all-degree bounds: at most one
# denominator argument is negative. A violation of either maximum bound
# means both homogeneous linear inequalities are strict; scaling permits
# replacing their positive right-hand sides by 1.
coordinate_rows=identity_matrix(QQ,6).rows()
arg_rows=glsm.T.rows()
irow,Lrow,Arow,Brow,jrow,krow=coordinate_rows
spectator_degree=Lrow+Arow+Brow
violations=[(irow-Lrow,irow-Arow+Lrow),
            (jrow-Lrow-2*Arow,jrow-Arow-Brow+Lrow),
            (krow-2*Lrow-2*Arow,krow-Arow-Brow+Lrow)]
for exception in range(10):
    inequalities=[[0]+list(row) for row in coordinate_rows]
    inequalities += [[0]+list(row) for index,row in enumerate(arg_rows) if index!=exception]
    for first,second in violations:
        assert Polyhedron(ieqs=inequalities+[[-1]+list(first),[-1]+list(second)]).is_empty()
    cone=Polyhedron(ieqs=inequalities)
    assert not cone.lines()
    for ray in cone.rays():
        degree=spectator_degree*vector(ray)
        assert (irow+jrow+krow)*vector(ray)<=5*degree
        assert all(abs(row*vector(ray))<=3*degree for row in arg_rows)
    # Even after allowing signed local degrees, a negative local index
    # is incompatible with having only this one possible negative argument.
    signed_inequalities=[[0]+list(row) for row in [Lrow,Arow,Brow]]
    signed_inequalities += [[0]+list(row) for index,row in enumerate(arg_rows) if index!=exception]
    for row in [irow,jrow,krow]:
        assert Polyhedron(ieqs=signed_inequalities+[[-1]+list(-row)]).is_empty()
print('Exact polyhedral certificates verify the degree bounds and signed-index extension for all degrees.')

# At fixed (L,A,B), only finitely many local degrees contribute to the
# scalar or first jets. Independently test the bounds in a wider box,
# including negative indices that are absent from the original Mori sum.
for d in range(4):
  for L,A,B in IntegerVectors(d,3):
    bounds=local_bounds(L,A,B)
    assert bounds[0]<=d and bounds[1]<=2*d and bounds[2]<=2*d
    for i in range(-2,3*d+3):
      for j in range(-2,3*d+3):
        for k in range(-2,3*d+3):
          args=factorial_arguments(i,L,A,B,j,k)
          if sum(entry<0 for entry in args)<=1:
            assert all(0<=n<=bound for n,bound in zip([i,j,k],bounds))
            assert max(abs(n) for n in args)<=3*d if d else args==[0]*10
print('First-jet local-degree bounds and bilateral extension checked independently.')

def binomial_jet(N,k):
    """B_N(k), B_N'(k), where B_N(t)=1/(Gamma(1+t) Gamma(1+N-t))."""
    N,k=ZZ(N),ZZ(k)
    if k>=0 and N-k>=0:
        c=1/(factorial(k)*factorial(N-k))
        return c,c*(harmonic(N-k)-harmonic(k))
    if k<0 and N-k<0:
        return QQ(0),QQ(0)
    if k<0:
        return QQ(0),(-1)^(-k-1)*factorial(-k-1)/factorial(N-k)
    return QQ(0),-(-1)^(k-N-1)*factorial(k-N-1)/factorial(k)

# The elementary convolution proof works even when a pair total is negative.
for N in range(7):
  for k in range(-15,16):
    c,d=binomial_jet(N,k)
    assert c==QQ(binomial(N,k))/factorial(N), ('binomial scalar',N,k,c)
    expected_derivative=sum(QQ(binomial(N,j))/factorial(N)*binomial_jet(0,k-j)[1] for j in range(N+1))
    assert d==expected_derivative, ('binomial convolution',N,k,d,expected_derivative)
for m in range(1,7):
  for k in range(-15,16):
    c,d=binomial_jet(-m,k)
    assert c==0
    assert d==(-1)^(m+k)*prod(k+j for j in range(1,m)), ('negative total',m,k,d)
  for M in range(m,m+4):
    for k in range(-5,6):
      assert sum(binomial(M,j)*binomial_jet(-m,k-j)[1] for j in range(M+1))==0
print('Binomial derivative convolution and negative-total finite differences checked.')

# Exact coefficient resummation on Sigma: z1=z5=z6=1.
# Test all spectator degrees L+A+B <= 5, keeping every local degree,
# not merely a total six-variable degree truncation.
tested=0
for d in range(6):
  for L,A,B in IntegerVectors(d,3):
    scalar=QQ(0)
    gradients=zero_vector(QQ,6)
    limits=local_bounds(L,A,B)
    for i in range(limits[0]+1):
      for j in range(limits[1]+1):
        for k in range(limits[2]+1):
          args=factorial_arguments(i,L,A,B,j,k)
          assert args==list(vector(ZZ,[i,L,A,B,j,k])*glsm)
          if sum(entry<0 for entry in args)>1:
            continue
          coefficient,gradient,hessian=gamma_coefficients(glsm,(i,L,A,B,j,k))
          scalar+=coefficient
          gradients+=gradient
    closed=factorial(2*A+B)*prod(reciprocal_factorial(n) for n in
             [B,L,L,L,A-L,2*A-B,B-A-2*L])
    assert scalar==closed
    assert all(gradients[j]==0 for j in [0,4,5])
    if d==0:
        assert scalar==1 and gradients==0
    tested+=1
assert tested==56
print('All 56 spectator coefficients through degree five: C1=C5=C6=0; W0 matches the closed formula.')

# Four factorial pairs and the shift cancellation underlying the all-orders
# proof. The constrained convolution has total index A-L and total order A.
RP=PolynomialRing(QQ,names=('i','j','k','L','A','B','rho1','rho5','rho6'))
i,j,k,L,A,B,rho1,rho5,rho6=RP.gens()
totals=[L,2*A-B,L,B-A-2*L]
indices=[i,j-i-B,k-j,A+B-L-k]
shifts=[rho1,rho5-rho1,rho6-rho5,-rho6]
assert sum(totals)==A and sum(indices)==A-L and sum(shifts)==0
assert sum([totals[t]-indices[t] for t in range(4)])==L

# The local conifold integral has unit multiplicity:
# normalized Res(du dv ds dt/(uv+st-mu)) integrates to +/-mu/(2*pi*i).
# A T^2 fibration over q=uv in [0,mu] gives (2*pi*i)^2*mu before
# the common (2*pi*i)^(-3) normalization. Thus there is no undetermined
# integer multiple when comparing residue derivatives to flat logarithms.

# Make the limiting residue rows regular, avoiding subtraction of radicals.
KR=PolynomialRing(QQ,names=('v','b')).fraction_field()
v,b=KR.gens()
# At r=0, d=sqrt(1+4v). At r=1/b, S=sqrt((1+b)^2+4vb^2).
RR=PolynomialRing(QQ,names=('d','S','b')).fraction_field()
d,S,b=RR.gens()
local_e2=vector(RR,[(1/d-1)/2,(1/d-1)/2,1/d])
local_e1=vector(RR,[(1+1/d)/2,(1+1/d)/2,1/d])
remote_C1=vector(RR,[(1-(1+b)/S)/2,(1+(1-b)/S)/2,-b/S])
remote_C2=vector(RR,[(1+(1+b)/S)/2,(1-(1-b)/S)/2,b/S])
assert local_e1-local_e2==remote_C1+remote_C2==vector(RR,[1,1,0])
limits=[vector(QQ,[entry.subs({d:1,S:1,b:0}) for entry in row])
        for row in [remote_C1,remote_C2,local_e1,local_e2]]
assert limits==[vector(QQ,[0,1,0]),vector(QQ,[1,0,0]),
                vector(QQ,[1,1,1]),vector(QQ,[0,0,1])]
print('Unit limiting charge rows in (beta1,beta5,beta6), ordered C1,C2,e1,e2:',limits)

e1,e2=[vector(ZZ,curve(i,9)) for i in [4,3]]
charges=matrix(ZZ,[C1,C2,e1,e2]).T
assert nef.T*charges==matrix(ZZ,[[0,1,1,0],[0,0,0,0],[0,0,0,0],
                                [0,0,0,0],[1,0,1,0],[0,0,1,1]])
assert charges.rank()==3
assert charges*vector(ZZ,[-1,-1,1,-1])==0
assert [n for n in charges.elementary_divisors() if n]==[1,1,1]

# The full 14-dimensional symplectic lattice includes D0 and its D6 dual.
# All four charges have zero D0 and magnetic entries, not just matching
# projections to a local flavor lattice.
I7=identity_matrix(ZZ,7)
Z7=zero_matrix(ZZ,7)
symplectic=block_matrix([[Z7,I7],[-I7,Z7]])
matrices=[]
for column in (nef.T*charges).columns():
    electric=vector(ZZ,[0]+list(column))
    matrix_M=block_matrix([[I7,Z7],[electric.column()*electric.row(),I7]])
    assert matrix_M.T*symplectic*matrix_M==symplectic
    matrices.append(matrix_M)
assert all(left*right==right*left for left in matrices for right in matrices)
assert (prod(matrices)-identity_matrix(ZZ,14)).rank()==3
print('Full integral charges in the original divisor-dual basis:',charges)
print('Saturated rank-three span; zero D0/D4/D6 components; commuting rank-three combined monodromy.')

# Batyrev--Kreuzer, Corollaries 1.9 and 3.9: vertex generation of
# the entire fan lattice suffices to kill pi_1 and the Brauer group.
# Check both sides, so no mirror H_3 torsion is invisible to periods.
fan_vertices=[(1,0,0,0),(0,1,0,0),(-1,-1,0,0),(0,0,1,0),(0,0,0,1),
              (-4,-2,-1,0),(-4,-2,0,-1),(3,1,1,1),(2,1,1,1)]
fan_polytope=Polyhedron(vertices=fan_vertices)
assert all(h.b()==1 for h in fan_polytope.inequalities())
polar_vertices=[vector(ZZ,h.A()) for h in fan_polytope.inequalities()]
polar_basis=matrix(ZZ,[(-1,2,1,-1),(1,-1,-1,-1),(0,-1,-1,1),(-1,-1,-1,4)])
assert all(row in polar_vertices for row in polar_basis.rows())
assert abs(polar_basis.det())==1
assert matrix(ZZ,fan_vertices).elementary_divisors()[:4]==[1,1,1,1]
assert matrix(ZZ,polar_vertices).elementary_divisors()[:4]==[1,1,1,1]
print('Original and polar vertices generate their full lattices; no mirror H3 torsion ambiguity.')
print('ALL INTEGRAL VANISHING-MARKING CHECKS PASSED; regular magnetic periods and stability are not computed by this check.')

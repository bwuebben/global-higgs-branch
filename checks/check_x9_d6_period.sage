"""Integral reduced period frame and D6/prepotential matching data.

Run: sage -c "load('checks/check_x9_d6_period.sage')"
The analytic proof uses polarized conifold descent, the seven-period
identification, and reality of the marked Abel continuation. This check
verifies their lattice/normalization inputs, not analytic continuation by
finite truncation. A separate zero-spectator calculation checks the four
trilogarithm constants without truncating their multiple covers.
"""
import contextlib
import io
with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_neutral_period_descent.sage')

def integral_frame(kappas,chern):
    """Demirtas et al. (2.20)--(2.24), with geometric integral generators.

    m0=O_X, ma=O_X(Ha)-O_X, e0=-O_point, and ea a pure D2 class.
    Pairing is minus the Euler pairing. These are Chern-character
    coordinates; c2/24 is added only when forming the Mukai charge.
    """
    r=len(kappas)
    ambient=2*r+2
    zero=zero_matrix(QQ,r)
    # ch ordering: (rank, ch1[1:r], ch2[1:r], ch3).
    # -chi(E,F) = -r*s' + s*r' + D.q' - q.D'
    #              - (r*D' - r'*D).c2/12.
    pairing=zero_matrix(QQ,ambient)
    pairing[0,-1]=-1
    pairing[-1,0]=1
    for a in range(r):
        pairing[1+a,1+r+a]=1
        pairing[1+r+a,1+a]=-1
        pairing[0,1+a]=-chern[a]/12
        pairing[1+a,0]=chern[a]/12
    cols=[]
    m0=zero_vector(QQ,ambient); m0[0]=1; cols.append(m0)
    for a in range(r):
        ma=zero_vector(QQ,ambient)
        ma[1+a]=1
        for b in range(r):
            ma[1+r+b]=kappas[b][a,a]/2
        ma[-1]=kappas[a][a,a]/6
        cols.append(ma)
    e0=zero_vector(QQ,ambient); e0[-1]=-1; cols.append(e0)
    for a in range(r):
        ea=zero_vector(QQ,ambient); ea[1+r+a]=1; cols.append(ea)
    charges=matrix(QQ,cols).T
    indices=vector(QQ,[chern[a]/12+kappas[a][a,a]/6 for a in range(r)])
    h=matrix(QQ,r,r,lambda a,b:(kappas[a][b,b]-kappas[a][a,b])/2)
    assert all(x in ZZ for x in indices) and all(x in ZZ for x in h.list())
    lam=zero_matrix(QQ,r+1)
    for a in range(r):
        lam[0,a+1]=-indices[a]
        lam[a+1,0]=indices[a]
        for b in range(r):
            lam[a+1,b+1]=h[a,b]
    eye=identity_matrix(ZZ,r+1)
    z=zero_matrix(ZZ,r+1)
    assert charges.T*pairing*charges==block_matrix([[lam,eye],[-eye,z]])
    # Column convention: corrected magnetic generators acquire electric
    # components. The D6 generator is unchanged.
    correction=zero_matrix(ZZ,r+1)
    for a in range(r):
        correction[a+1,0]=indices[a]
        for b in range(a+1,r):
            correction[a+1,b+1]=h[a,b]
    transform=block_matrix([[eye,z],[correction.T,eye]])
    assert transform.det()==1
    symp=block_matrix([[z,eye],[-eye,z]])
    new_charges=charges*transform
    assert new_charges.T*pairing*new_charges==symp
    acoef=matrix(QQ,r,r,lambda a,b:
                 kappas[a][a,b]/2 if a>=b else kappas[a][b,b]/2)
    assert acoef==acoef.T
    for a in range(r):
        # -integral exp(-t) ch(E) sqrt(Td) gives the stated Fa polynomial.
        q=new_charges[1+r:1+2*r,a+1]
        assert vector(q)==vector(acoef.row(a))
        assert -new_charges[-1,a+1]-chern[a]/24==chern[a]/24
    return indices,h,acoef,transform,new_charges

old_c2=vector(QQ,vector(c2D)*nef)
new_c2=vector(QQ,sc2)
old_frame=integral_frame(nef_kappa,old_c2)
new_frame=integral_frame(smooth_kappa,new_c2)
assert old_frame[2].matrix_from_rows_and_columns([1,2,3],[1,2,3])==new_frame[2]
assert new_c2==old_c2[1:4]
print('Old RR indices:',old_frame[0])
print('Smooth RR indices:',new_frame[0])
print('Common surviving quadratic matrix a_ij:\n',new_frame[2])

# Divisor intersections need not generate the integral curve lattice.
# The old index is one, but the smooth index is ten. The integral electric
# basis uses all topological D2 charges, not just complete intersections.
curve_indices=[]
for qmat,kappas in [(glsm,nef_kappa),(Qsm,smooth_kappa)]:
    curve_columns=matrix(ZZ,[[u*k*v for k in kappas]
                            for u,v in combinations_with_replacement(qmat.columns(),2)]).T
    curve_indices.append(prod(x for x in curve_columns.elementary_divisors() if x))
assert curve_indices==[1,10]
# Vertex generation on both sides kills pi1 and the Brauer group by
# Batyrev--Kreuzer. Thus integral cohomology is torsion-free. The untwisted
# K-theory AHSS degenerates (its differentials are rationally zero), and
# supplies a lift of every integral H4 class. For rank=c1=0, RR makes ch3
# integral, removable by a point class. No BPS existence claim is made.
sm_polytope=Polyhedron(vertices=smooth_rays)
sm_polar_vertices=[vector(ZZ,h.A()) for h in sm_polytope.inequalities()]
assert all(h.b()==1 for h in sm_polytope.inequalities())
for vertices in [smooth_rays,sm_polar_vertices]:
    assert [x for x in matrix(ZZ,vertices).elementary_divisors() if x]==[1]*4
print('Smooth-side vertex lattices exclude torsion; complete-intersection D2 index is 10, not 1.')

# The integral reduced frame is a coordinate subquotient, after the
# curvature/flux corrections above, not in unshifted Mukai coordinates.
assert surviving.T*J14*surviving==J8
seven=identity_matrix(ZZ,8).matrix_from_columns([1,2,3,4,5,6,7])
assert (seven.T*J8).right_kernel().basis_matrix()==matrix(ZZ,[[0,0,0,0,1,0,0,0]])
print('An eighth-period difference pairing trivially with the other seven is a scalar-period multiple.')

# Euler characteristics from the charge matrices and cubic tensors,
# independently of the Hodge-number subtraction.
def euler_from_glsm(qmat,kappas):
    total=sum(qmat.columns())
    cubic=lambda x:sum(x[a]*(x*kappas[a]*x) for a in range(len(kappas)))
    return -(cubic(total)-sum(cubic(q) for q in qmat.columns()))/3
chi_old=euler_from_glsm(glsm,nef_kappa)
chi_sm=euler_from_glsm(Qsm,smooth_kappa)
assert (chi_old,chi_sm)==(-232,-240)
instanton_shift=(chi_sm-chi_old)/2
assert instanton_shift==-4
assert chi_old+2*instanton_shift==chi_sm

# Independent all-local-degree boundary calculation. With spectator
# degree zero, at most two negative arguments occur only on four active
# rays and two zero-intersection rays. Three-negative sectors contract
# to zero by the cubic tensor. Polyhedra certify this for every degree.
local_map=matrix(ZZ,[[1,0,0],[0,0,0],[0,0,0],
                     [0,0,0],[0,1,0],[0,0,1]])
local_arguments=(local_map.T*glsm).T.rows()
for negative in combinations(range(10),3):
    if not trip(*[vector(CLS[names[i]]) for i in negative]):
        continue
    inequalities=[[0]+list(row) for row in identity_matrix(ZZ,3).rows()]
    inequalities += [[-1]+list(-local_arguments[i]) if i in negative
                     else [0]+list(local_arguments[i]) for i in range(10)]
    assert Polyhedron(ieqs=inequalities).is_empty()

# If exactly two reciprocal Gamma zeros are present, the contracted
# third Taylor coefficient equals prefactor*(Di.Dj).gradient_log.
# At zero spectator degree each active ray has prefactor 1/h^2 and
# this last contraction -2/h. The harmonic pieces cancel identically.
ray_data=[((2,4),(1,0,0)),((6,7),(0,1,0)),
          ((4,9),(1,1,1)),((3,9),(0,0,1))]
for pair,ray in ray_data:
    slopes=vector(ZZ,ray)*local_map.T*glsm
    assert [i for i,s in enumerate(slopes) if s<0]==list(pair)
    assert sorted(slopes)==[-1,-1,0,0,0,0,0,0,1,1]
    charge=vector(QQ,[glsm.column(pair[0])*k*glsm.column(pair[1]) for k in nef_kappa])
    positive_sum=sum((glsm.column(i) for i,s in enumerate(slopes) if s>0),zero_vector(QQ,6))
    negative_sum=sum((glsm.column(i) for i in pair),zero_vector(QQ,6))
    assert charge*positive_sum==2
    assert charge*(positive_sum+negative_sum)==0
    # gradient=-positive_sum*H_h-negative_sum*H_(h-1)
    # H_h-H_(h-1)=1/h, giving exactly -2/h^3 on every ray.
assert -2*len(ray_data)==2*instanton_shift==-8
# Direct multivariate Taylor multiplication checks the contraction's
# combinatorial normalization independently of the harmonic formula.
third_ring=PowerSeriesRing(QQ,names=tuple('r'+str(i) for i in range(6)),default_prec=4)
rho_vars=vector(third_ring,third_ring.gens())
def direct_rational_cubic(n):
    scalar=third_ring(1).add_bigoh(4)
    total=sum(glsm.columns())
    total_index=vector(ZZ,n)*total
    for j in range(1,total_index+1):
        scalar*=j+rho_vars*total
    for column in glsm.columns():
        index=vector(ZZ,n)*column
        linear=rho_vars*column
        if index>=0:
            for j in range(1,index+1):
                scalar/=j+linear
        else:
            for j in range(-index):
                scalar*=linear-j
    value=QQ(0)
    for a,b,c in combinations_with_replacement(range(6),3):
        exponent=tuple(ZZ(i==a)+ZZ(i==b)+ZZ(i==c) for i in range(6))
        value+=nef_kappa[a][b,c]*scalar[exponent]
    return value
for pair,ray in ray_data:
    for h in range(1,5):
        assert direct_rational_cubic(h*local_map*vector(ZZ,ray))==-2/QQ(h)^3
for local in IntegerVectors(3,3):
    full=local_map*vector(ZZ,local)
    expected=QQ(0)
    for pair,ray in ray_data:
        if vector(local)==3*vector(ray):
            expected=-2/QQ(3)^3
        if tuple(local)==ray:
            expected=-2
    assert direct_rational_cubic(full)==expected
print('Each zero-spectator cubic Gamma tail is -2/h^3; total -8*zeta(3).')
print('chi_old=-232, chi_sm=-240; F_inst shift is -4*zeta(3)/(2*pi*i)^3.')
print('ALL D6 FRAME/CONSTANT DATA CHECKS PASSED; the analytic equality uses the printed pairing and reality proof.')

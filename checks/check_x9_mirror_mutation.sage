"""Exact birational identification with the smooth-side Batyrev mirror.

Run: sage -c "load('checks/check_x9_mirror_mutation.sage')"
All arithmetic is exact. The printed geometric argument, not this script,
identifies the companion smoothing's ambient fan and invokes the flop theorem.
No chosen crepant chamber, BPS survival or five-dimensional limit is asserted.
"""

rays = [vector(ZZ,r) for r in [(-1,-3,-1,0),(-1,-3,0,-1),
    (-1,-1,0,0),(-1,0,0,0),(0,1,0,0),(1,-1,0,0),(1,3,1,1)]]
cones = [[int(c) for c in row] for row in
    ['1456','1346','0456','0346','0145','0134','01256','1236','0236','0123']]
normals = [vector(ZZ,m) for m in [(-2,-1,-2,6),(1,-1,-2,3),
    (-2,-1,6,-2),(1,-1,3,-2),(-2,-1,6,6),(1,-1,3,3),
    (0,1,-2,-2),(1,0,-2,0),(1,0,0,-2),(1,0,0,0)]]
P = Polyhedron(vertices=rays)
dual = Polyhedron(ieqs=[[1]+list(r) for r in rays])
assert P.dim()==dual.dim()==4 and dual.is_compact()
assert set(map(tuple,dual.vertices_list()))==set(map(tuple,normals))
for cone,m in zip(cones,normals):
    assert all(m*rays[i]==-1 if i in cone else m*rays[i]>-1 for i in range(7))
    if len(cone)==4:
        assert abs(matrix(ZZ,[rays[i] for i in cone]).det())==1
assert rays[0]+rays[1]+rays[6]==2*rays[2]+rays[5]
extra = vector(ZZ,[0,-1,0,0])
origin = vector(ZZ,[0,0,0,0])
points = list(P.integral_points())
dpoints = list(dual.integral_points())
assert set(map(tuple,points))==set(map(tuple,rays+[extra,origin]))
assert len(dpoints)==162

def interior_count(face, lattice_points):
    return sum(face.relative_interior_contains(p) for p in lattice_points)

def h21(poly, polar, lattice_points, polar_points):
    facet_sum = sum(interior_count(f.as_polyhedron(),lattice_points)
                    for f in poly.faces(3))
    correction = 0
    for face in poly.faces(2):
        fpoly = face.as_polyhedron()
        lp = interior_count(fpoly,lattice_points)
        if lp:
            df = Polyhedron(vertices=[m for m in polar.vertices_list()
                if all(vector(m)*vector(v)==-1 for v in fpoly.vertices_list())])
            correction += lp*interior_count(df,polar_points)
    return len(lattice_points)-5-facet_sum+correction, facet_sum, correction

assert h21(P,dual,points,dpoints)==(3,1,0)
assert h21(dual,P,dpoints,points)==(123,34,0)
print('Reflexive mirror Newton polytope: 9 points, dual 162; (h11,h21)=(123,3).')

# Original ten monomials, independent of the four-node check.
old_rays = [vector(ZZ,r) for r in [(1,0,0,0),(0,1,0,0),(-1,-1,0,0),
    (0,0,1,0),(0,0,0,1),(-4,-2,-1,0),(-4,-2,0,-1),
    (3,1,1,1),(2,1,1,1),(-2,-1,0,0)]]
basis = matrix(ZZ,[(0,1,0,0),(1,0,0,0),(0,1,1,0),(0,1,0,1)]).T
assert abs(basis.det())==1
K = PolynomialRing(QQ,names=('T','u','p','q','a','b','v','h','lam')).fraction_field()
T,u,p,q,a,b,v,h,lam = K.gens()
def monomial(exponent, variables):
    return prod(z^n for z,n in zip(variables,exponent))

old_coeff = [b,a,1,v,1,1,v,1,1,1]
old = sum(c*monomial(basis.inverse()*r,[h,u,p,q])
          for c,r in zip(old_coeff,old_rays))-1
G = 1+u^4*p*q
assert old==a*h+b*u-1+G/h*(u^-1+u^-4*p^-1+v*u^-4*q^-1+u^-2)
new_coeff = [1,v,1,1,b,a,a]
new = sum(c*monomial(r,[T,u,p,q]) for c,r in zip(new_coeff,rays))-1
assert old.subs(h=G*T/u)==new
assert new.subs(T=u*h/G)==old

# The mutation is a dense chart of the existing blowup of (s,1+xy).
x,y,r,s = u^2*p,u^2*q,u,G*T
f,g = 1+x*y,1+1/x+v/y+r
chart = f*g+s*(r*(b*r-1)+a*s)
strict = g+T*r*(b*r-1)+a*T^2*f
assert chart==G*strict and strict/(u*T)==new
# The original Laurent equation is chart/(r*s), with h=s/r.
assert old.subs(h=s/r)==chart/(r*s)
new_variables = [T,u,p,q]
old_variables = [G*T/u,u,p,q]
log_jacobian = matrix(K,[[z*w.derivative(z)/w for z in new_variables]
                        for w in old_variables])
assert log_jacobian.det()==1
print('Full polynomial, inverse map, small-resolution chart and log-volume preservation: PASS.')

# Remove the unique facet-interior coefficient by its root automorphism.
normal = vector(ZZ,[0,1,-2,-2])
assert normal in normals and normal*extra==-1
assert all(m*extra>=0 for m in normals if m!=normal)
assert P.interior_contains(origin)
assert sum(m*extra==-1 for m in normals)==1
root_transformed = sum(c*monomial(r,new_variables)*(1+lam/u)^(1+normal*r)
                       for c,r in zip(new_coeff,rays))-(1+lam/u)
assert root_transformed-new==2*b*lam+(b*lam^2-lam)/u+lam/(T*u)
# A pre-existing c/u is unchanged. Its transformed coefficient has
# derivative -1 in lam at our slice: the omitted coefficient is a gauge.
assert (b*lam^2-lam).derivative(lam).subs(lam=0)==-1
gauge = matrix(ZZ,[rays[0],rays[2],rays[3],rays[6]-rays[5]])
assert gauge.det()==1
print('One root gauge and four unimodular torus gauges leave exactly (v,a,b).')

# Exact nondegenerate member, checked on every face including the full torus.
# Passing to rational face coordinates is a finite torus cover and cannot
# remove critical points in characteristic zero.
sample_coeff = [1,QQ(1)/2,1,1,QQ(1)/10,QQ(1)/100,QQ(1)/100,-1]
assert all(sample_coeff[i]!=0 for i in range(7))  # No zero-dimensional face meets the hypersurface.
support = rays+[origin]
checked = 0
for dim in range(1,5):
    faces = [P] if dim==4 else [f.as_polyhedron() for f in P.faces(dim)]
    for face in faces:
        indices = [i for i,r in enumerate(support) if face.contains(r)]
        differences = matrix(QQ,[support[i]-support[indices[0]] for i in indices])
        span = differences.row_space().basis_matrix()
        coords = [span.solve_left(row) for row in differences.rows()]
        denom = lcm(c.denominator() for row in coords for c in row)
        exponents = [vector(ZZ,row*denom) for row in coords]
        minimum = vector(ZZ,[min(row[j] for row in exponents) for j in range(dim)])
        RF = PolynomialRing(QQ,names=tuple('x%d'%j for j in range(dim))+('inv',))
        xx = RF.gens()
        polynomial = sum(sample_coeff[i]*monomial(e-minimum,xx[:dim])
                         for i,e in zip(indices,exponents))
        ideal = RF.ideal([polynomial]+[polynomial.derivative(z) for z in xx[:dim]]
                         +[xx[-1]*prod(xx[:dim])-1])
        assert ideal==RF.ideal(1), (dim,indices)
        checked += 1
print('Exact nondegeneracy sample: all %d positive-dimensional faces including the full torus PASS.'%checked)

# Integral coefficient coordinates and all-orders period index map.
relations = matrix(ZZ,[[1,1,0,0,4,1,1],[0,0,1,0,2,1,0],[0,0,0,1,1,1,0]])
assert relations*matrix(ZZ,rays)==0
assert relations.row_module()==matrix(ZZ,rays).left_kernel()
assert [prod(K(c)^n for c,n in zip(new_coeff,row)) for row in relations]==[v*a^2*b^4,a*b^2,a*b]
parameter_map = matrix(ZZ,[[1,2,4],[0,1,2],[0,1,1]])
assert parameter_map.det()==-1
index_map = parameter_map.T
assert index_map.inverse()==matrix(ZZ,[[1,0,0],[-2,-1,1],[0,2,-1]])
count = 0
for alpha in range(6):
    for beta in range(6-alpha):
        for gamma in range(6-alpha-beta):
            ell,m,n = index_map*vector(ZZ,[alpha,beta,gamma])
            counts = vector(ZZ,[alpha,beta,gamma])*relations
            constant_term = factorial(sum(counts))/prod(factorial(k) for k in counts)
            restricted = factorial(2*m+n)/(factorial(n)*factorial(ell)^3*
                factorial(m-ell)*factorial(2*m-n)*factorial(n-m-2*ell))
            assert constant_term==restricted
            count += 1
assert count==56
print('Unimodular parameter/index maps and %d exact period coefficients: PASS.'%count)
print('All mirror-mutation checks PASS; no chosen-chamber or 5d-limit claim.')

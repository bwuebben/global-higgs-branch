"""Exact four-node specialization of the compact X9 mirror.

Run: sage -c "load('checks/check_x9_finite_bulk_nodes.sage')"
Independent of the truncated Frobenius expansion. This checks actual
threefold singularities, not just zeros of candidate charge periods.
It does not compute global integral period transport or magnetic stability.
"""

vertices = [vector(ZZ,p) for p in
    [(1,0,0,0),(0,1,0,0),(-1,-1,0,0),(0,0,1,0),(0,0,0,1),
     (-4,-2,-1,0),(-4,-2,0,-1),(3,1,1,1),(2,1,1,1),(-2,-1,0,0)]]
poly = Polyhedron(vertices=vertices)
glsm_bulk = matrix(ZZ,[[0,0,-1,0,-1,1,0,1,0,0],
    [0,0,0,1,1,0,0,0,-1,-1],[0,1,0,0,0,0,0,2,-2,1],
    [1,0,-1,0,0,0,0,0,0,1],[0,0,1,1,0,0,-1,-1,0,0],
    [0,0,0,-1,0,0,1,0,1,-1]])

# A unimodular monomial change, with x,y invertible and A,T regular.
monomials = matrix(ZZ,[(2,1,1,0),(2,1,0,1),(1,0,0,0),(1,1,0,0)]).T
assert abs(monomials.det()) == 1
normal_A = vector(ZZ,[1,-1,-1,-1])
normal_T = vector(ZZ,[0,1,-1,-1])
assert monomials.inverse()[2,:] == matrix(ZZ,[normal_A])
assert monomials.inverse()[3,:] == matrix(ZZ,[normal_T])
assert all(1+normal_A*p>=0 and 1+normal_T*p>=0 for p in vertices)
assert gcd(normal_A-normal_T)==1
assert gcd(matrix(ZZ,[normal_A,normal_T]).minors(2))==1
assert [i for i,p in enumerate(vertices) if normal_A*p==normal_T*p==-1] == [3,4,5,6,8,9]

PB = PolynomialRing(QQ,names=('z1','z5','w','v','a','b','x','y','A','T'))
z1,z5,w,v,a,b,x,y,A,T = PB.gens()
KB = PB.fraction_field()
coefficients = list(map(KB,[b,a*w^2*z5^2,1,v,1,z1*z5*w,v*w,1/(w*z5),1,1]))
assert [prod(coefficients[j]^glsm_bulk[i,j] for j in range(10))
        for i in range(6)] == [z1,v,a,b,z5,w]
chart = sum(coefficients[i]*prod(KB(t)^n for t,n in zip([x,y,A,T],
                monomials.inverse()*(p-vertices[9]))) for i,p in enumerate(vertices))-A*T
expected = 1+v*x+y+z1*z5*w/x+v*w/y+x*y+A*(1+x*y/(w*z5))-A*T+b*A^2*T+a*w^2*z5^2*T^2
assert chart==expected
special = chart.subs({z1:1,z5:1,w:1})
f = 1+x*y
g = 1+1/x+v/y+A
h = A*(b*A-1)+a*T
assert special==f*g+T*h
print('Full compact chart: F =',chart)
print('On z1=z5=w=1: F=(1+xy)(1+1/x+v/y+A)+T(A(bA-1)+aT).')

# The critical points have T=0, xy=-1, A=0 or 1/b,
# and v*x^2-(1+A)*x-1=0. Use x,A,b free, eliminating v.
node_sub = {z1:1,z5:1,w:1,y:-1/x,T:0,v:(1+(1+A)*x)/x^2}
assert chart.subs(node_sub)==0
for variable in [x,y,A]:
    assert chart.derivative(variable).subs(node_sub)==0
assert chart.derivative(T).subs(node_sub)==A*(b*A-1)
Hessian = matrix(KB,[[chart.derivative(i).derivative(j) for j in [x,y,A,T]] for i in [x,y,A,T]])
det_node = Hessian.det().subs(node_sub)
disc_node = (1+A+2/x)^2
for location in [KB(0),1/b]:
    assert det_node.subs(A=location)==disc_node.subs(A=location)
print('Node Hessian determinant = (1+A)^2+4v, at A=0 or 1/b.')

# Four smoothing covectors; rank exactly three when the nodes are distinct.
covector = vector(KB,[chart.derivative(parameter).subs(node_sub)
                     for parameter in [z1,z5,w]])
assert covector==vector(KB,[1/x,1/x+A,-1])
for parameter in [v,a,b]:
    assert chart.derivative(parameter).subs(node_sub)==0
PC = PolynomialRing(QQ,names=('r','s','t','B'))
r,s,t,B = PC.gens()
minor = matrix(PC,[[r,r,-1],[s,s,-1],[t,t+B,-1]]).det()
assert minor==B*(r-s)
print('Smoothing rank: a 3x3 minor is (r-s)/b, nonzero.')

# Residue normalization and the unit-coefficient relation. At each A,
# D=1/x+vx=+/-sqrt((1+A)^2+4v), 1/x=(D-1-A)/2.
PD = PolynomialRing(QQ,names=('D','AA')).fraction_field()
D,AA = PD.gens()
def residue_row(D,AA):
    r=(D-1-AA)/2
    return vector(PD,[r,r+AA,-1])/D
assert residue_row(D,AA)+residue_row(-D,AA)==vector(PD,[1,1,0])
assert special.derivative(A).subs(T=0)==f
assert h.derivative(A)==2*b*A-1
print('Opposite A-roots have h_A=-1,+1: the four oriented residue rows sum to zero.')

# Square charts really are toric charts; the remote pair approaches the
# C1 and C2 square critical points as b -> 0, rather than disappearing.
for index in [6,5]:
    square = matrix(ZZ,[vertices[index]-vertices[2],vertices[7]-vertices[2],
                       (-1,0,0,0),(2,1,0,0)]).T
    assert abs(square.det())==1
    expected_normal = [-1,2,-1,1] if index==6 else [-1,2,1,-1]
    assert list(square.inverse().row(2))==expected_normal
    assert square.inverse().row(3)==normal_T
    assert all(n>=0 for p in vertices for n in (square.inverse()*(p-vertices[2]))[2:])
    for power in [0,1]:
        character=square.column(power)
        expr=prod(KB(t)^n for t,n in zip([x,y,A,T],monomials.inverse()*character))
        assert expr == ([1/(A*y),x*y] if index==6 else [1/(A*x),x*y])[power]
    square_polynomial=sum(coefficients[i]*prod(KB(t)^n for t,n in zip([x,y,A,T],
        square.inverse()*(p-vertices[2]))) for i,p in enumerate(vertices))
    square_polynomial-=prod(KB(t)^n for t,n in zip([x,y,A,T],
        square.inverse()*(-vertices[2])))
    assert square_polynomial.derivative(T).subs({A:0,T:0})==b

# On the compact rational curve f=g=T=0, the two missing x-chart
# endpoints lie in the square charts. The normal derivative there is b.
# The four interior zeros of h are a quartic with nonzero end coefficients.
intersection_polynomial=PB(x^2*h.subs({A:-1-1/x+v*x,T:0}))
assert intersection_polynomial==(v*x^2-x-1)*(b*(v*x^2-x-1)-x)
assert intersection_polynomial.coefficient({x:0})==b
assert intersection_polynomial.coefficient({x:4})==b*v^2

# Exact sample: v=1/2, a=1/100, b=1/10. The saturated critical
# scheme of this whole chart consists of precisely the four stated nodes.
PS = PolynomialRing(QQ,names=('xx','yy','AA','TT','inv'),order='degrevlex')
xx,yy,AA,TT,inv = PS.gens()
Fs=(1+xx*yy)*(xx*yy+yy+QQ(1)/2*xx+AA*xx*yy)+xx*yy*TT*(AA*(AA/10-1)+TT/100)
critical_ideal=PS.ideal([Fs]+[Fs.derivative(t) for t in [xx,yy,AA,TT]]+[inv*xx*yy-1])
predicted_ideal=PS.ideal([TT,xx*yy+1,AA*(AA-10),QQ(1)/2*xx^2-(1+AA)*xx-1,inv+1])
assert critical_ideal==predicted_ideal
assert critical_ideal.vector_space_dimension()==4
print('Exact rational sample: saturated critical scheme has length four and is reduced.')

# Check every other toric stratum at the same sample. Eleven facets have
# no face-critical points. The remaining facet is T=0, already treated;
# the only nonsimplicial two-faces are the pentagon and the two squares.
sample_coefficients=[QQ(1)/10,QQ(1)/100,1,QQ(1)/2,1,1,QQ(1)/2,1,1,1]
facet_dimensions=[]
for face in poly.faces(3):
    on=[i for i,p in enumerate(vertices) if face.as_polyhedron().contains(p)]
    differences=matrix(QQ,[vertices[i]-vertices[on[0]] for i in on])
    basis=differences.row_space().basis_matrix()
    coordinates=[basis.solve_left(row) for row in differences.rows()]
    denominator=lcm([c.denominator() for row in coordinates for c in row])
    exponents=[vector(ZZ,row*denominator) for row in coordinates]
    minimum=vector(ZZ,[min(row[k] for row in exponents) for k in range(3)])
    RF=PolynomialRing(QQ,names=('r','s','t','inverse'))
    rf=RF.gens()
    polynomial=sum(sample_coefficients[i]*prod(rf[k]^(e[k]-minimum[k]) for k in range(3))
                   for i,e in zip(on,exponents))
    ideal=RF.ideal([polynomial]+[polynomial.derivative(t) for t in rf[:3]]+
                  [rf[3]*prod(rf[:3])-1])
    critical_dimension=ideal.dimension()
    expected_dimension=1 if on==[2,3,4,5,6,7,8,9] else -1
    assert critical_dimension==expected_dimension
    facet_dimensions.append((on,critical_dimension))
print('Facet critical dimensions (-1 means empty):',facet_dimensions)
nonsimplicial=[]
for dimension in [1,2]:
    for face in poly.faces(dimension):
        on=tuple(i for i,p in enumerate(vertices) if face.as_polyhedron().contains(p))
        if len(on)>dimension+1:
            nonsimplicial.append(on)
assert sorted(nonsimplicial)==[(2,3,6,7),(2,4,5,7),(3,4,5,6,8,9)]
print('No other singularities at the sample: facets, squares and lower strata checked.')
print('ALL FINITE-BULK NODE CHECKS PASSED; magnetic stability and global period marking are not computed.')

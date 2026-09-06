"""Local I2 vanishing charges, exact endpoint and compact charge transport.

Run from the project root:
sage -c "load('checks/check_x9_local_vanishing_charges.sage')"
The all-orders endpoint proof is the printed Chu--Vandermonde identity,
not this finite coefficient check. No finite-bulk stability is asserted.
"""
import contextlib
import io
with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_mirror_continuation.sage')

PRv = PolynomialRing(QQ,names=('u','v','w','x','y','d','c'))
u,v,w,x,y,d,c = PRv.gens()
F = x*(1+x)*y^2+(v*x^2+x+u*w)*y+v*w*x
Qv = F.discriminant(y)
Dv = Qv.discriminant(x)
residual = 16*v^2*w+27*v*w^2-18*v*w-v+1
assert Dv.subs(u=1) == -256*v^3*w^4*(w-1)^2*residual
assert F.subs({u:1,w:1}) == (x*y+1)*(v*x+x*y+y)
assert Qv.subs({u:1,w:1}) == (v*x^2-x-1)^2
assert (27*w^2-18*w-1)^2-64*w == (w-1)*(9*w-1)^3
assert Dv.coefficient({v:3}) == -256*u^2*w^4*(w-1)*(u*w-1)
assert Dv.subs({u:1+d,v:1+c*d,w:1+c*d}).coefficient({d:2}) == -6400*(c^2+c-1)
print('Root-mass splitting: w=1+(-1 +/- sqrt(5))*(u-1)/2+O((u-1)^2).')

# Match the conventional circle E2 curve, Closset--Magureanu (7.1),
# on v=w: U=-1/w, lambda=u, M1=1, with Weierstrass scale 2w.
field_v = PRv.fraction_field()
a4,b4,c4,d4,e4 = [Qv.coefficient({x:i}).subs(v=w) for i in [4,3,2,1,0]]
Iu = 12*a4*e4-3*b4*d4+c4^2
Ju = 72*a4*c4*e4+9*b4*c4*d4-27*a4*d4^2-27*b4^2*e4-2*c4^3
U = -1/field_v(w)
standard_g2 = (U^4-8*(1+u)*U^2-24*u*U+16*(1-u+u^2))/12
standard_g3 = -(U^6-12*(1+u)*U^4-36*u*U^3+144*u*(1+u)*U+
                   24*(2+u+2*u^2)*U^2-8*(8-12*u-39*u^2+8*u^3))/216
assert 4*Iu/3==(2*w)^4*standard_g2
assert 4*Ju/27==(2*w)^6*standard_g3
print('Exact match to the standard circle E2 curve: U=-1/w, lambda=u, M1=1.')

# Both nodes are transverse; their smoothing parameters are independent.
K5 = QuadraticField(5,'r5')
r5 = K5.gen()
nodes = [(1+r5)/2,(1-r5)/2]
rows = []
for root in nodes:
    ynode = -1/root
    evaluate = PRv.hom([K5(1),K5(1),K5(1),root,ynode,K5(0),K5(0)],K5)
    assert evaluate(F)==evaluate(F.derivative(x))==evaluate(F.derivative(y))==0
    jac = matrix(K5,[[evaluate((x*y+1).derivative(variable)) for variable in [x,y]],
                    [evaluate((x*y+x+y).derivative(variable)) for variable in [x,y]]])
    assert jac.det()!=0
    rows.append([evaluate(F.derivative(w)+F.derivative(v)),evaluate(F.derivative(u))])
assert matrix(K5,rows).det()!=0
print('Two distinct ordinary nodes; independent first-order smoothing parameters:',matrix(K5,rows))

# General local logarithmic correction, expanded in v. All u,w powers
# at each v degree are included. Summing m uses ordinary Vandermonde;
# the remaining polynomial is terminating 2F1(-ell,2ell;ell+1;w).
RW = PolynomialRing(QQ,'ww')
ww = RW.gen()
for ell in range(1,21):
    direct = RW(0)
    mass_derivative = RW(0)
    for k in range(ell,2*ell+1):
        coefficient_sum = QQ(0)
        for m in range(k-ell,ell+1):
            coefficient = (-1)^(ell+k-1)*factorial(ell+k-1)/(
                factorial(m)*factorial(ell+m-k)*factorial(ell-m)*
                factorial(k-m)*factorial(k-ell))
            reflected_m = k-m
            reflected = (-1)^(ell+k-1)*factorial(ell+k-1)/(
                factorial(reflected_m)*factorial(ell+reflected_m-k)*
                factorial(ell-reflected_m)*factorial(k-reflected_m)*factorial(k-ell))
            assert coefficient==reflected
            coefficient_sum += coefficient
            mass_derivative += m*coefficient*ww^k
        assert coefficient_sum == (-1)^(ell+k-1)*factorial(ell+k-1)*binomial(2*ell,k)/(
            factorial(ell)^2*factorial(k-ell))
        direct += coefficient_sum*ww^k
    prefactor = -factorial(2*ell-1)*factorial(2*ell)/factorial(ell)^4
    hypergeom = sum(rising_factorial(-ell,j)*rising_factorial(2*ell,j)/(
        rising_factorial(ell+1,j)*factorial(j))*ww^j for j in range(ell+1))
    assert direct == prefactor*ww^ell*hypergeom
    assert hypergeom(1) == rising_factorial(1-ell,ell)/rising_factorial(ell+1,ell) == 0
    assert direct(1)==0
    assert 2*mass_derivative == ww*direct.derivative()
print('Endpoint coefficient cancellation and Weyl reflection checked through v^20.')

# The positive-real triangle 0<v<=w<1 has no further discriminant.
# For w>=1/9 the residual quadratic has nonpositive discriminant.
# For w<=1/9 it decreases over 0<=v<=w and its endpoint exceeds 2/3.
small_w_derivative_bound = 59*QQ(1)/81-1
assert small_w_derivative_bound<0
assert QQ(1)-QQ(1)/9-18*QQ(1)/81 == QQ(2)/3
assert residual.subs(v=w) == 43*w^3-18*w^2-w+1
assert residual.subs(w=1) == (4*v+1)^2

# Differential normalization at the split point. The mass reflection
# h(u,v,w)=h(1/u,v,u*w) gives 2 h_u=h_w at u=w=1.
omega_endpoint = 1/r5
logq_mass_derivative = -(1-omega_endpoint)/2
slopes = [(-1-r5)/2,(-1+r5)/2]
assert omega_endpoint*slopes[0]+logq_mass_derivative+1==0  # e1
assert omega_endpoint*slopes[1]+logq_mass_derivative==0    # e2

# Primitive compact lifts in the established six-dimensional charge basis.
e1,e2 = [vector(ZZ,curve(i,9)) for i in [4,3]]
assert e1-e2==R==C1+C2
assert nef.T*curve(8,9)==vector(ZZ,[0,1,0,0,0,0])
assert e1 == vector(ZZ,[0,1,0,0,1,-1])
assert e2 == vector(ZZ,[0,0,1,0,1,-1])
assert gcd(e1)==gcd(e2)==1
charges = matrix(ZZ,[C1,C2,e1,e2]).T
assert charges.rank()==3
assert charges*vector(ZZ,[-1,-1,1,-1])==0
assert [entry for entry in charges.elementary_divisors() if entry] == [1,1,1]
print('Compact charges, columns C1,C2,e1,e2:',charges)
pairings_local = matrix(ZZ,[[p*charge for charge in [e1,e2]]
                          for p in [PA,PV,PT,CLS['E']]])
assert pairings_local == matrix(ZZ,[[-1,0],[1,1],[0,1],[-1,-1]])
print('Local-curve pairings, rows A,V,T,E:',pairings_local)

Me1,Me2 = [transvection(charge) for charge in [e1,e2]]
assert Me1*Me2==Me2*Me1
assert (Me1*Me2-1).rank()==2
assert (M1*M2*Me1*Me2-1).rank()==3
assert (Me1*Me2*seed-seed)[6:]==e2
local_electric_projection = matrix(ZZ,[[1,1]])
assert (local_electric_projection*local_electric_projection.T)[0,0]==2

# The rank-three effective compact charge matrix is integrally equivalent
# to the four-hypermultiplet model with relation (-1,-1,1,-1).
basis3 = matrix(ZZ,[C1,C2,e2]).T
Qeffective = matrix(ZZ,[[1,0,1,0],[0,1,1,0],[0,0,1,1]])
assert basis3*Qeffective==charges
assert [entry for entry in basis3.elementary_divisors() if entry] == [1,1,1]
print('Effective compact charge matrix:',Qeffective)
print('ALL LOCAL VANISHING-CHARGE CHECKS PASSED; finite-bulk magnetic stability remains open.')

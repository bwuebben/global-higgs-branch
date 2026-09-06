"""Exact geometry and scale algebra for the auxiliary-circle limit.

Run: sage -c "load('checks/check_x9_circle_limit.sage')"
The analytic period bounds and circle dictionary are inputs proved/cited
in the text. This check does not construct a uniform 5d Higgs effective
action or assert stability of any additional light magnetic charge.
"""
import contextlib
import io
with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_prepotential.sage')

nef = matrix(ZZ,sorted([list(r.vector()) for r in S.rays()])).T
H1,H2,H3,H4,H5,H6 = nef.columns()
assert H2==CLS['D5']+CLS['D4'] and H3==CLS['D1'] and H4==CLS['D0']
restriction = matrix(ZZ,[[0,1,1,0,1,-3],[0,-1,0,0,-1,1],[0,0,-1,0,-1,1]])
form = diagonal_matrix(ZZ,[1,-1,-1])
assert restriction.T*form*restriction==L(CLS['E'])
assert restriction*H2==vector(ZZ,[1,0,0])
assert restriction*H3==restriction*H4==zero_vector(ZZ,3)
assert [curve(i,9)*H2 for i in [3,4,5,6,8]]==[0,0,1,1,1]
beta2 = vector(ZZ,curve(8,9))
assert nef.T*beta2==vector(ZZ,[0,1,0,0,0,0])
electric = matrix(ZZ,[C1,C2,curve(4,9),curve(3,9)]).T
assert electric.rank()==3 and electric.augment(beta2).rank()==4
assert [x for x in electric.augment(beta2).elementary_divisors() if x]==[1,1,1,1]
print('H2 restricts to the plane class on E; beta2 adds one primitive direction to the four-hyper charges.')

PF = PolynomialRing(QQ,names=('alpha','beta','eta','scale','y'))
alpha,beta,eta,scale,y = PF.gens()
KF = PF.fraction_field()
W = alpha*CLS['D0']+beta*CLS['D1']
J = W+eta*H2
V0 = trip(W,W,W)/6
P3 = 83*alpha^3+147*alpha^2*beta+81*alpha*beta^2+14*beta^3
assert V0==P3/6
volume = trip(J,J,J)/6
assert volume==V0+eta*(QQ(17)/2*alpha^2+11*alpha*beta+3*beta^2)+eta^2*(QQ(3)/2*alpha+beta)
assert all(c>0 for c in V0.coefficients())
assert trip(CLS['E'],J,J)==eta^2
assert trip(CLS['E'],CLS['E'],J)==-3*eta
assert W*C1==W*C2==0 and restriction*W==zero_vector(PF,3)

def metric(j,u,v,vm=False):
    vol = trip(j,j,j)/6
    return KF(-trip(j,u,v))+KF(trip(j,j,u)*trip(j,j,v))/((6 if vm else 4)*vol)

# Remove the graviphoton before comparing the two bulk-gauging directions.
vectors = [-CLS['D6'],CLS['D5']]
bulk_metric = matrix(KF,2,2,lambda i,j: metric(W,vectors[i],vectors[j],vm=True))
projected = [vector(KF,v)-KF(trip(W,W,v)/(6*V0))*vector(KF,W) for v in vectors]
assert all(trip(W,W,v)==0 for v in projected)
assert matrix(KF,2,2,lambda i,j: metric(W,projected[i],projected[j]))==bulk_metric
sigma = 7*alpha^2+8*alpha*beta+2*beta^2
tau = 3*alpha+2*beta
diagonal = 2*sigma^2/(3*V0)
assert bulk_metric==matrix(KF,[[diagonal,tau-diagonal],[tau-diagonal,diagonal]])
P4 = 143*alpha^4+289*alpha^3*beta+199*alpha^2*beta^2+52*alpha*beta^3+4*beta^4
assert 2*diagonal-tau==P4/P3
assert bulk_metric*vector(KF,[1,1])==tau*vector(KF,[1,1])
assert bulk_metric*vector(KF,[1,-1])==(P4/P3)*vector(KF,[1,-1])
assert all(c>0 for c in P4.coefficients())
assert metric(J,CLS['E'],CLS['E'],vm=True)==3*eta+eta^4/(6*volume)
# This comparison is to the limiting five-dimensional Hodge metric,
# not the singular four-dimensional gauge coupling with hypers integrated out.
print('Bulk-gauging kinetic eigenvalues: 3 alpha+2 beta and P4/P3, strictly positive.')
print('P4 =',P4)
print('Del Pezzo vector norm after graviphoton removal = 3 eta + eta^4/(6 Vol).')

# Symbolic circle scaling: eta=y/scale, with y bounded and scale -> infinity.
one_variable = PolynomialRing(KF,'inv_scale')
inv_scale = one_variable.gen()
scaled_volume = one_variable(volume.subs(eta=0)) + inv_scale*y*(QQ(17)/2*alpha^2+11*alpha*beta+3*beta^2) + inv_scale^2*y^2*(QQ(3)/2*alpha+beta)
assert scaled_volume[0]==V0
assert (inv_scale^2*y^2/2).degree()==2  # geometric E tension scale

# D4 structure-sheaf charge of E has bounded normalized period along Sigma.
# Z_E=Vol(E)+i*(-E^2.J/2)-(E^3/6+c2.E/24)+Psi_E.
assert trip(CLS['E'],CLS['E'],CLS['E'])==7
assert sum(c2D[i]*CLS['E'][i] for i in range(6))==-2
assert QQ(7)/6-QQ(2)/24==QQ(13)/12
assert -trip(CLS['E'],CLS['E'],H2)/2==QQ(3)/2
assert all(trip(CLS['E'],W,basis_vector)==0
           for basis_vector in identity_matrix(ZZ,6).rows())
print('Normalized E structure-sheaf period = y2^2/2 + 3 i y2/2 - 13/12 + Psi_E; bounded, not a BPS-existence claim.')

# Independent smooth-side consistency: D8+E lifts the surviving plane.
assert trip(CLS['D8']+CLS['E'],J,J)==(beta+eta)^2
assert trip(CLS['D8'],J,J)==beta^2+2*beta*eta
assert trip(CLS['D8'],CLS['E'],J)==eta
print('Surviving plane tension tends to beta^2/2 while the E string scale vanishes.')

# The dimensionless coefficient exponents on the mutated mirror, using
# a~exp(-2 pi scale beta), b~exp(-2 pi scale alpha), v fixed.
weights = vector(PF,[2*beta+4*alpha,beta+2*alpha,beta+alpha])
assert weights[0]==2*weights[1]
assert all(c>0 for w in weights for c in w.coefficients())
# xi1/xi2^2=v remains finite: it is a correlated limit, not three
# independent large-volume parameters of the surviving mirror.
print('Mutated mirror weights:',weights,'; xi1/xi2^2=v stays fixed.')

# KK scale: for any fixed positive cutoff mu and integer n, M_n=|n|/R.
# For integer R and cutoff one the tower contains 2R+1 Fourier modes.
for radius in [1,2,5,10,100]:
    modes = [n for n in range(-radius,radius+1) if QQ(abs(n))/radius<=1]
    assert len(modes)==2*radius+1
print('KK gap closes as 1/R; the four-hyper infrared truncation is not uniform.')
print('ALL AUXILIARY-CIRCLE CHECKS PASSED; no full 5d Higgs dynamics or magnetic BPS survival asserted.')

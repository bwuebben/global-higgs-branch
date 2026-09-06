"""Exact nodal monodromy and the del Pezzo face discriminant of X9.

Run: sage -c "load('checks/check_x9_mirror_continuation.sage')"
This does not continue the full compact magnetic period pair to E2.
"""
import contextlib
import io

with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_mirror_periods.sage')

# The only nonsimplicial two-faces: two squares and the E pentagon.
face_sets = []
for face in P.faces(2):
    face_poly = face.as_polyhedron()
    on = tuple(i for i,p in enumerate(pts[:-1]) if face_poly.contains(p))
    if len(on)>3:
        face_sets.append(on)
assert sorted(face_sets) == [(2,3,6,7),(2,4,5,7),(3,4,5,6,8,9)]
ex = pts[3]-pts[9]
ey = pts[4]-pts[9]
assert pts[5]-pts[9] == -ex
assert pts[6]-pts[9] == -ey
assert pts[8]-pts[9] == ex+ey
assert gcd(matrix(ZZ,[ex,ey]).minors(2)) == 1

# Elliptic face curve, with u=z1*z5, v=z2, w=z6.
# After torus rescaling:
# 1+v*x+y+u*w/x+v*w/y+x*y=0.
# Clearing x*y makes this quadratic in y.
PR = PolynomialRing(QQ,names=('u','v','w','x','y'))
u,v,w,x,y = PR.gens()
Fface = x*(1+x)*y^2+(v*x^2+x+u*w)*y+v*w*x
quartic = (v*x^2+x+u*w)^2-4*v*w*x^2*(1+x)
assert Fface.discriminant(y) == quartic
disc = quartic.discriminant(x)
print('Elliptic quartic:',quartic)
print('Quartic discriminant factorization:',disc.factor())

# Independent critical-point (Horn) parametrization for the face.
KR = PolynomialRing(QQ,names=('alpha','beta','gamma')).fraction_field()
alpha,beta,gamma = KR.gens()
Qface = matrix(ZZ,[[1,-1,1,-1,0,0],
                  [1,1,0,0,-1,-1],
                  [-1,0,0,1,1,-1]])
face_points = matrix(ZZ,[[1,0,1],[0,1,1],[-1,0,1],
                       [0,-1,1],[1,1,1],[0,0,1]])
assert Qface*face_points == 0
critical = vector(KR,[alpha,beta,gamma])*Qface
uh,vh,wh = [prod(critical[j]^Qface[i,j] for j in range(6)) for i in range(3)]
evaluate_face = PR.hom([uh,vh,wh,KR(0),KR(0)],KR)
assert evaluate_face(disc) == 0
print('Face Horn coordinates (u,v,w):',uh,vh,wh)
print('Equal-area specialization v=w:',disc.subs(v=w).factor())
print('Trivial root holonomy u=1:',disc.subs(u=1).factor())

# Full principal torus-critical locus, without a large implicit elimination.
LR = PolynomialRing(QQ,names=tuple('ell'+str(i+1) for i in range(6))).fraction_field()
ells = vector(LR,LR.gens())
critical_full = ells*glsm
central = sum(critical_full)
assert central == 2*ells[2]+ells[3]
assert critical_full*matrix(ZZ,pts[:-1]) == 0
full_horn = [prod((critical_full[j]/central)^glsm[i,j]
                  for j in range(10)) for i in range(6)]
print('Full torus-critical parametrization: a_i=(ell*Q)_i, a_*=2ell3+ell4.')

# Electric vanishing charges have zero D0 component. The ordered magnetic
# and electric coordinates are the dual integral six-dimensional bases.
# A positive loop in t.C gives Delta Psi_P=(P.C)*(t.C).
I6 = identity_matrix(ZZ,6)
Z6 = zero_matrix(ZZ,6)
omega = block_matrix([[Z6,I6],[-I6,Z6]])
def transvection(c):
    c = vector(ZZ,c)
    return block_matrix([[I6,Z6],[c.column()*c.row(),I6]])

M1,M2,MR = [transvection(c) for c in [C1,C2,R]]
for monodromy in [M1,M2,MR]:
    assert monodromy.T*omega*monodromy == omega
    assert (monodromy-1)^2 == 0
assert M1*M2 == M2*M1
assert (M1*M2-1).rank() == 2
assert (MR-1).rank() == 1
assert M1*M2 != MR
assert (M1*M2-MR)[6:,0:6] == -(C1.column()*C2.row()+C2.column()*C1.row())

PA = CLS['D4']
PV = CLS['D7']+CLS['D8']
PT = CLS['D6']
assert PA+PV == PT
pairings = matrix(ZZ,[[p*c for c in [C1,C2]] for p in [PA,PV,PT]])
assert pairings == matrix(ZZ,[[0,-1],[-1,1],[-1,0]])
seed = vector(QQ,list(PT)+[0]*6)
assert (M1*M2*seed-seed)[6:] == -C1
assert (MR*seed-seed)[6:] == -R
for relative in [C2,R]:
    charged = seed+vector(QQ,[0]*6+list(relative))
    for monodromy in [M1,M2,M1*M2]:
        assert monodromy*charged-monodromy*seed == charged-seed
print('Nodal magnetic pairings, rows A,V,T and columns C1,C2:\n',pairings)
print('Combined nodal monodromy has rank 2; a hypothetical R transvection has rank 1.')

# Exact normalized electric periods on the previous flat-coordinate path.
SRpath = PolynomialRing(QQ,names=('s','epsilon','lam','a','b','kk','ii'))
s,epsilon,lam,a,b,kk,ii = SRpath.gens()
tpath = -s*vector(SRpath,CLS['D5'])+ii*lam*(a*vector(SRpath,CLS['D0'])+
              b*vector(SRpath,CLS['D1'])+epsilon*vector(SRpath,J0))
assert tpath*C1 == ii*lam*epsilon
assert tpath*C2 == -s+ii*lam*epsilon
assert tpath*R == -s+2*ii*lam*epsilon
assert tpath*curve(4,9) == -s+3*ii*lam*epsilon
assert tpath*curve(3,9) == ii*lam*epsilon
assert tpath*curve(8,9) == ii*lam*epsilon
print('C2 tower: Z_k=k-s+i*lambda*epsilon; a massless endpoint requires s in ZZ.')

# Sign check by direct continuation of the dilogarithm derivative.
# The logarithm of (1-exp(2pi*i*t))/(-2pi*i*t) is single valued on this
# small disk; the explicit i*theta term follows the logarithm's sheet.
import mpmath as mp
with mp.workdps(50):
    tbase = mp.j/50
    def continued_dilog_derivative(theta):
        t = tbase*mp.exp(mp.j*theta)
        regular = (1-mp.exp(2*mp.pi*mp.j*t))/(-2*mp.pi*mp.j*t)
        continued_log = mp.log(-2*mp.pi*mp.j*tbase)+mp.j*theta+mp.log(regular)
        return 2*mp.pi*t*continued_log
    change = mp.quad(continued_dilog_derivative,[0,mp.pi,2*mp.pi])/(4*mp.pi^2)
    assert abs(change-tbase)<mp.mpf('1e-45')
    print('Positive-loop Li2 continuation: Delta Psi/(P.C)=',mp.nstr(change,20))
print('ALL CONTINUATION-DATA CHECKS PASSED; no compact E2 period endpoint asserted.')

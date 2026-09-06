"""Exact polarization and B-field slope tests for the isolated D6 probes.

Run: sage -c "load('checks/check_x9_probe_transport.sage')"
No claim of exact quantum stability, a complete DT jump, or a particle
index is made. The reduced-pencil genericity hypothesis is not used.
"""
import contextlib
import io

with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_prepotential.sage')

T = CLS['D6']
A, V = CLS['D4'], CLS['D7']+CLS['D8']
parts = [A,CLS['D7'],CLS['D8']]
ell = CLS['D5']
F = L(T)*ell
assert F == vector(QQ,[1,0,0,0,0,1])
assert L(ell)*ell == 0 and ell*C1 == 0 and ell*C2 == ell*(C1+C2) == 1
assert L(ell)*A == F and L(ell)*V == 0

PR = PolynomialRing(QQ,names=('a','b','eps'))
a,b,eps = PR.gens()
J = a*vector(PR,CLS['D0'])+b*vector(PR,CLS['D1'])+eps*vector(PR,J0)
assert J*C1 == J*C2 == eps
assert J*F == 3*a+2*b+13*eps
areaT = trip(J,J,T)/2
assert areaT == (14*a^2+16*a*b+132*a*eps+4*b^2+72*b*eps+297*eps^2)/2
assert areaT.subs(eps=0) == 7*a^2+8*a*b+2*b^2

subsets = [ids for size in [1,2] for ids in combinations(range(3),size)]
monomials = [a^3,a^2*b,a*b^2,b^3,a^2*eps,a*b*eps,b^2*eps,a*eps^2,b*eps^2,eps^3]
# Coefficientwise common lower bounds for the two sheaves, printed in Appendix D.
lower_rows = [
    [42,76,44,8,553,662,190,2385,1420,3388],
    [42,90,60,12,639,882,279,3141,2073,4974],
    [0,42,48,12,70,476,237,660,1255,1488],
    [0,42,48,12,70,476,233,660,1239,1476],
    [42,90,60,12,595,842,275,2781,1937,4282],
    [42,76,44,8,597,702,198,2745,1572,4092],
]
denominator = 2*(J*F)*trip(J,J,V)
assert all(c>0 for c in denominator.coefficients())
assert denominator(0,0,1) == 4576

for label,weights,wall_at_J0 in [('F2',[1,0,0],QQ(77)/104),
                                ('F_R',[1,1,0],QQ(165)/208)]:
    records = []
    for ids,row in zip(subsets,lower_rows):
        U = sum(parts[i] for i in ids)
        complement = T-U
        base = trip(J,U,U+2*complement)*trip(J,J,T)-trip(J,T,T)*trip(J,J,U)
        delta = PR(base+2*eps*(sum(weights)*trip(J,J,U)
                              -sum(weights[i] for i in ids)*trip(J,J,T)))
        lower = sum(c*m for c,m in zip(row,monomials))
        assert all(c>=0 for c in (delta-lower).coefficients())
        assert all(c>=0 for c in delta.coefficients())
        assert delta.monomial_coefficient(eps^3)>0
        # B=-s*ell: the slope comparison is delta+s*k.
        k = PR(2*((J*F)*trip(J,J,U)-trip(J,ell,U)*trip(J,J,T)))
        records.append((delta,k))
    d1,k1 = records[0]
    assert k1 == -denominator
    assert d1(0,0,1)/denominator(0,0,1) == wall_at_J0
    assert all(c<0 for c in (d1+k1).coefficients())
    assert (d1+k1).monomial_coefficient(eps^3)<0
    for delta,k in records[1:]:
        residual = PR(denominator*delta+d1*k)
        assert all(c>=0 for c in residual.coefficients())
        assert residual.monomial_coefficient(eps^6)>0
    face_num = 7*a^2+8*a*b+2*b^2
    face_den = 11*a^2+10*a*b+2*b^2
    assert d1.subs(eps=0)*face_den == denominator.subs(eps=0)*face_num
    print(label, 'stable for all a,b>=0, eps>0 at B=0.')
    print('  first fixed-charge B=-s D5 slope wall at J*: s=',wall_at_J0)
    print('  D4 comparison at s=1:',d1+k1)
assert QQ(7+8+2)/QQ(11+10+2) == QQ(17)/23
print('Common limiting slope wall: (7a^2+8ab+2b^2)/(11a^2+10ab+2b^2).')

def flow(p,q,q0,line):
    return p,q+L(p)*line,q0+line*q+trip(p,line,line)/2

seed = (T,-C1/2,QQ(11)/12)
charged = [(T,-C1/2+C2,-QQ(1)/12),(T,C1/2+C2,-QQ(1)/12)]
flowed_seed = flow(*seed,ell)
assert flowed_seed == (T,seed[1]+F,seed[2])
for charge,gamma in zip(charged,[C2,C1+C2]):
    flowed = flow(*charge,ell)
    assert flowed[1]-flowed_seed[1] == gamma
    assert flowed[2]-flowed_seed[2] == 0
    # Gauge-invariant relative degree-six component at B'=ell is still -1.
    assert flowed[2]-flowed_seed[2]-ell*gamma == -1
print('Spectral flow by D5 sets the relative raw D0 to zero and adds F to the core.')

def divisor_sheaf(p,line):
    c2p = sum(c2D[k]*p[k] for k in range(6))
    q = L(p)*line-L(p)*p/2
    q0 = trip(p,p,p)/6-trip(p,p,line)/2+trip(p,line,line)/2+c2p/24
    return p,q,q0

G = divisor_sheaf(A,CLS['D2']+A)
H2 = divisor_sheaf(V,zero_vector(QQ,6))
HR = divisor_sheaf(V,T)
assert G[2] == -QQ(1)/2 and H2[2] == HR[2] == QQ(5)/12
for H,charge in zip([H2,HR],charged):
    assert tuple(G[i]+H[i] for i in range(3)) == charge
    assert G[0]*H[1]-H[0]*G[1] == -1
assert G[2]+sum(c2D[k]*A[k] for k in range(6))/24 == 0
print('Magnetic extensions: G on D4, H on D7+D8; each pairing is -1.')
print('ALL EXACT CHECKS PASSED. Quantum endpoint and full index remain undetermined.')

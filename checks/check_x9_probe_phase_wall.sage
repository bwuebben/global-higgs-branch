"""Finite-volume magnetic phase alignment, not numerical quantum periods.

Run: sage -c "load('checks/check_x9_probe_phase_wall.sage')"
The analytic large-volume persistence argument and the Ext/cohomology
arguments are in the manuscript. This checks their intersection arithmetic,
the exact polynomial central charges, and the constituent slope tests.
"""
import contextlib
import io

with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_probe_transport.sage')

vA, vV = trip(J,J,A)/2, trip(J,J,V)/2
g, f = J*G[1], J*F
assert vA == (3*a^2+6*a*b+2*b^2+42*a*eps+34*b*eps+121*eps^2)/2
assert vV == (11*a^2+10*a*b+2*b^2+90*a*eps+38*b*eps+176*eps^2)/2
assert g == -(3*a+2*b+11*eps)/2
assert f == 3*a+2*b+13*eps
assert ell*G[1] == 1
assert ell*H2[1] == ell*HR[1] == 0
assert G[2] == -QQ(1)/2 and H2[2] == HR[2] == QQ(5)/12

QS = PolynomialRing(QQ,names=('a','b','eps','lam','s'))
aa,bb,ee,lam,s = QS.gens()
inc = PR.hom([aa,bb,ee],QS)
av, vv, gg, ff = map(inc,[vA,vV,g,f])
field = QS.fraction_field()

def central_parts(charge):
    """Expand -P.t^2/2+q.t-q0 directly from the cubic and Mukai charge."""
    p,q,q0 = charge
    re = lam^2*inc(PR(trip(J,J,p)/2))-s^2*trip(ell,ell,p)/2-s*(ell*q)-q0
    im = lam*(inc(PR(J*q))+s*inc(PR(trip(J,ell,p))))
    return re,im

for label,H,h_expected,expected_n,expected_d in [
    ('F2',H2,(3*a+2*b+12*eps)/2,20328*lam^2+17,27456*lam^2+14),
    ('F_R',HR,(3*a+2*b+14*eps)/2,21780*lam^2+29,27456*lam^2+38),
]:
    h = J*H[1]
    assert h == h_expected
    hh = inc(h)
    # Real/imaginary parts in Z=-P.t^2/2+q.t-q0, B=-s D5.
    reG, imG = lam^2*av-s+QQ(1)/2, lam*(gg+s*ff)
    reH, imH = lam^2*vv-QQ(5)/12, lam*hh
    assert (reG,imG) == central_parts(G)
    assert (reH,imH) == central_parts(H)
    determinant = imG*reH-reG*imH
    K, C = ff*vv, av*hh-vv*gg
    slope_den = lam^2*K+hh-5*ff/12
    numerator = lam^2*C+5*gg/12+hh/2
    assert determinant == lam*(s*slope_den-numerator)
    assert determinant.degree(s) == 1  # all quadratic B terms cancel
    root = field(numerator)/field(slope_den)
    assert determinant(s=root) == 0
    assert 4*C == inc(PR(
        trip(J,A,A+2*V)*trip(J,J,T)-trip(J,T,T)*trip(J,J,A)
        +2*((J*(G[1]+H[1]+C1/2))*trip(J,J,A)
            -eps*trip(J,J,T))))
    assert numerator(0,0,1,lam,0)*expected_d == slope_den(0,0,1,lam,0)*expected_n

    # At J* and lambda>=1: the root lies in (0,1), is transverse,
    # and positive real parts exclude anti-alignment.
    n0 = numerator(0,0,1,lam,0)
    d0 = slope_den(0,0,1,lam,0)
    assert all(c>0 for c in n0.coefficients())
    assert all(c>0 for c in d0.coefficients())
    assert (d0-n0)(lam=1)>0
    assert (d0-n0).monomial_coefficient(lam^2)>0
    assert vA(0,0,1)-QQ(1)/2 == 60
    assert vV(0,0,1)-QQ(5)/12 == QQ(1051)/12

    # Strict slope stability of H at the first F wall. A line bundle
    # on V has maximal subsheaf O_U(L-(V-U)) for U=D7 or D8.
    at_J0 = [422,282] if label == 'F2' else [419,285]
    for U,expected_gap in zip([CLS['D7'],CLS['D8']],at_J0):
        W = V-U
        line = zero_vector(QQ,6) if label == 'F2' else T
        subq = L(U)*(line-W)-L(U)*U/2
        # D5.V=0, so the H comparisons do not depend on s.
        gapH = PR(h*trip(J,J,U)-(J*subq)*trip(J,J,V))
        assert all(c0>=0 for c0 in gapH.coefficients())
        assert gapH.monomial_coefficient(eps^3)>0
        assert gapH(0,0,1) == expected_gap
        print(label, 'quotient subunion',list(U),'strict slope gap:',gapH)
    print(label, 'polynomial phase wall at J*: s =',expected_n,'/',expected_d)

for charge,multiplicity in zip(charged,[1,2]):
    re,im = central_parts(charge)
    re_seed,im_seed = central_parts(seed)
    assert re-re_seed == 1-s
    assert im-im_seed == multiplicity*lam*ee
    assert charge[0] == seed[0]  # the same magnetic period cancels

# Charge, cohomology, and Euler inputs to the rigidity argument.
assert trip(A,A,A) == 0
assert SURF['D4'][:3] == (0,12,1)
assert G[0]*H2[1]-H2[0]*G[1] == -1
assert G[0]*HR[1]-HR[0]*G[1] == -1
print('ALL EXACT PHASE-WALL CHECKS PASSED.')
print('Quantum periods are not evaluated; their smallness is an asymptotic input.')

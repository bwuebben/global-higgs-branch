"""Exact data for two isolated charged sheaves on T=D4+D7+D8.

Run: sage -c "load('checks/check_x9_glued_probe.sage')"
The geometric proof (duality, stability, spherical twist and rigidity)
is in papers/03-magnetic-sheaves/sections/02_sheaves.tex. This checks its algebraic inputs;
it does not classify the full moduli or certify a physical chamber.
The additional reduced-pencil hypothesis is NOT used here.
"""
import contextlib
import io

with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_prepotential.sage')

T = CLS['D6']
A, B, C = [CLS[k] for k in ['D4', 'D7', 'D8']]
D = CLS['D2']
parts = [A, B, C]
assert sum(parts) == T
assert L(D)*A == C2 and L(T)*B == C1 and L(T)*T == C1
assert [part*C2 for part in parts] == [-1, 1, 0]
assert T*C2 == 0 and D*C2 == A*C2 == -1
assert J0*C1 == J0*C2 == 1
assert trip(T,A,J0) == trip(T,C,J0) == 0
assert [trip(J0,J0,part) for part in parts] == [121,173,3]
assert trip(J0,J0,T) == 297 and trip(J0,T,T) == 1
print('C2 complete-intersection normal degrees (-1,-1); component incidences (-1,1,0).')

# Two-term locally free resolution of I_(C2,T):
# 0 -> O(-D2-D4) + O(-T) -> O(-D2) + O(-D4) -> I -> 0.
# With a=z4, d=z2, b=z7*z8, the map is [(-a,0),(d,b)].
# Its transpose presents Ext^1_X(I,O_X)=F_R; determinant cuts out T.
Rloc = PolynomialRing(QQ, names=('a','d','b'))
a,d,b = Rloc.gens()
resolution_map = matrix(Rloc, [[-a,0],[d,b]])
assert resolution_map.det() == -a*b
assert resolution_map.subs({a:0,d:0,b:0}).rank() == 0
assert resolution_map.subs({a:0,d:0,b:1}).rank() == 1
print('Local resolution determinant -a*b; dual is non-locally-free at the double-curve point.')

def ch_line_sum(positive, negative):
    """Rank-zero Chern character of a difference of line bundles."""
    assert len(positive) == len(negative)
    p = sum(positive)-sum(negative)
    q = (sum(L(v)*v for v in positive)-sum(L(v)*v for v in negative))/2
    s = (sum(trip(v,v,v) for v in positive)-sum(trip(v,v,v) for v in negative))/6
    return p,q,s

pI,qI,sI = ch_line_sum([-D,-A],[-D-A,-T])
pR,qR,sR = ch_line_sum([D+A,T],[D,A])
assert pI == pR == T and qI == -C1/2-C2 and qR == C1/2+C2
assert qR == -qI and sR == sI
c2T = sum(c2D[k]*T[k] for k in range(6))
q0R = sR+c2T/24
assert q0R == -QQ(1)/12
q2 = qR-L(T)*T
q02 = q0R-T*qR+trip(T,T,T)/2
assert q2 == -C1/2+C2 and q02 == q0R
seed_q = -C1/2
seed_q0 = trip(T,T,T)/6+c2T/24
assert seed_q0 == QQ(11)/12 and q0R-seed_q0 == -1
assert qR-seed_q == C1+C2 and q2-seed_q == C2
print('F2: Q2=-C1/2+C2; F_R: Q2=C1/2+C2; both Q0=-1/12.')
print('Relative to the D6 seed: (C2,-1) and (R,-1), respectively.')

expected = {'I': [4092,4530,1476,1494,4974,3630],
            'F2': [3388,5222,1488,1482,4282,4334],
            'F_R': [3630,4974,1494,1476,4530,4092]}
for name, flux_degrees in [('I',[-1,0,0]),('F2',[1,0,0]),('F_R',[1,1,0])]:
    gaps = []
    for size in [1,2]:
        for indices in combinations(range(3),size):
            U = sum(parts[i] for i in indices)
            V = T-U
            base = trip(J0,U,U+2*V)*trip(J0,J0,T)-trip(J0,T,T)*trip(J0,J0,U)
            gap = base+2*(sum(flux_degrees)*trip(J0,J0,U)
                          -sum(flux_degrees[i] for i in indices)*trip(J0,J0,T))
            gaps.append(gap)
    assert gaps == expected[name] and min(gaps)>0
    print(name, 'stability numerators:', gaps)

# P1 cohomology input: normal bundle has no sections; zero-D0 O(-1)
# is acyclic, while the quotient O(-2) has Euler characteristic -1.
h0_P1 = lambda n: max(n+1,0)
h1_P1 = lambda n: max(-n-1,0)
assert 2*h0_P1(-1) == 0
assert h0_P1(-1) == h1_P1(-1) == 0
assert h0_P1(-2)-h1_P1(-2) == -1
assert 2-h0_P1(0) == 1  # kernel of the surjective pencil restriction
print('Rigidity inputs: h0(N_C2/X)=0; one section of I_C2(T); chi(O_C2(-2))=-1.')
print('ALL ALGEBRAIC CHECKS PASSED. Each proved isolated stable point has DT weight +1.')

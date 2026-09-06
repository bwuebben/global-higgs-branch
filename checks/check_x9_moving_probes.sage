"""Exact arithmetic for moving D5/D6 magnetic probes of X9.

Run: sage -c "load('checks/check_x9_moving_probes.sage')"
The seed-component proof is in the manuscript and assumes reduced pencils
whose only reducible members are the displayed toric ones. This check does
NOT certify that genericity hypothesis or a physical stability condition.
"""
import contextlib
import io

with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_prepotential.sage')

P5, P6 = CLS['D5'], CLS['D6']
parts5 = [CLS[k] for k in ['D3', 'D7', 'D8']]
parts6 = [CLS[k] for k in ['D4', 'D7', 'D8']]
assert sum(parts5) == P5 and sum(parts6) == P6
assert L(P5)*P5 == 0 and L(P6)*P6 == C1
assert trip(P6, CLS['D4'], J0) == trip(P6, CLS['D8'], J0) == 0
assert trip(P6, CLS['D7'], J0) == C1*J0 == 1
assert trip(P6, P6, CLS['D7']) == trip(P6, CLS['D7'], CLS['D7']) == -1

# Each ambient divisor has exactly the two displayed toric sections.
for idx, direction in [(5, (0,0,1,0)), (6, (0,0,0,1))]:
    poly = Polyhedron(ieqs=[[int(j == idx)] + list(p) for j, p in enumerate(pts)])
    assert sorted(map(tuple, poly.integral_points())) == sorted([(0,0,0,0), direction])
    monomials = [tuple(vector(ZZ,m)*p + int(j == idx) for j,p in enumerate(pts))
                 for m in poly.integral_points()]
    print(names[idx], 'section exponents:', monomials)

# J0 is the sum of all primitive ambient nef generators, hence ample.
assert J0 == vector(QQ, [6,18,16,-11,-8,7])
expected = {'D5': [2952,3690,1230,1230,3690,2952],
            'D6': [3740,4876,1482,1488,4628,3982]}
for name, total, pieces in [('D5',P5,parts5), ('D6',P6,parts6)]:
    gaps = []
    for size in [1,2]:
        for indices in combinations(range(3),size):
            A = sum(pieces[i] for i in indices)
            B = total-A
            assert trip(J0,J0,A) > 0 and trip(J0,A,B) > 0
            gap = trip(J0,A,A+2*B)*trip(J0,J0,total)-trip(J0,total,total)*trip(J0,J0,A)
            gaps.append(gap)
    assert gaps == expected[name] and min(gaps) > 0
    print(name, 'proper-subunion stability numerators:', gaps)

# Mukai charges of the zero/one-point ideal sheaves on each pencil member.
for name, divisor, expected_q0 in [('D5',P5,1), ('D6',P6,QQ(11)/12)]:
    cubic = trip(divisor,divisor,divisor)
    c2value = sum(c2D[k]*divisor[k] for k in range(6))
    q2 = -L(divisor)*divisor/2
    q0 = cubic/6 + c2value/24
    assert q0 == expected_q0
    assert q2 == (zero_vector(QQ,6) if name == 'D5' else -C1/2)
    print(name, 'Q2:',q2,'Q0 at n=0,1:',q0,q0-1)

# Euler and Behrend-sign arithmetic for the PROVED conditional components.
euler_X = 2*(6-122)
euler_blowup = euler_X + (2-1)*2  # exceptional P1-bundle over C1=P1
assert euler_X == -232 and euler_blowup == -230
assert [-2,-euler_X] == [-2,232]
assert [-2,-euler_blowup] == [-2,230]

# Charge restrictions are unconditional intersection statements.
A6 = CLS['D7']-P6
assert A6 == -CLS['D4']-CLS['D8'] and L(P6)*A6 == 0
assert P5*C2 == P5*(C1+C2) == 1
assert A6*C1 == 0 and A6*C2 == A6*(C1+C2) == 1
assert A6*(-C1/2) == 0
print('D5-horizontal charge: Q_D5(C2)=Q_D5(R)=1.')
print('D6 boundary selector A6=D7-D6=-(D4+D8): (C1,C2,R) -> (0,1,1).')
print('Conditional seed components: P1,X and P1,Bl_C1(X); indices (-2,232),(-2,230).')
print('ALL ARITHMETIC CHECKS PASSED; reduced-pencil hypothesis and physical chamber not certified.')

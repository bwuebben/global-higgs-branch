"""Exact protected-ring checks for the formal X9 connected-torus quotient.

Physical input: the reduced E2 current ring Jp*Jm = beta^2, with
current charges +/-R and R=C1+C2. This is NOT a computation of a
finite-Planck index or a construction of a rigid limit of X9.

Run: sage -c "load('checks/check_x9_dressed_instantons.sage')"
"""

from collections import defaultdict
from itertools import product

# 1. The primitive effective action follows from the compact integral charges.
C = matrix(ZZ, [[1,-1], [0,1], [-1,0], [-1,1], [0,0], [0,0]])
assert C.rank() == 2
assert C.smith_form()[0].diagonal() == [1,1]
root = C.column(0) + C.column(1)
assert root == vector(ZZ, [0,1,-1,0,0,0])
# Raising current maps the e2 state to the e1 state, so its charge is
# their difference, R, not 2R. The double-cover variables are not operators.
charges = matrix(ZZ, [[1,-1,0,0,1,-1,0], [0,0,1,-1,1,-1,0]])
degrees = (1,1,1,1,2,2,2)  # h1,k1,h2,k2,Jp,Jm,beta
x_exp = vector(ZZ, [1,0,1,0,0,1,0])
y_exp = vector(ZZ, [0,1,0,1,1,0,0])
assert charges*x_exp == charges*y_exp == vector(ZZ, [0,0])
assert vector(ZZ, degrees)*x_exp == vector(ZZ, degrees)*y_exp == 4
# In the root-normalized action, b has charge (-1/2,-1/2).
assert vector(QQ, [1,1]) - vector(QQ, [1/2,1/2]) != 0

# 2. Eliminate genuine operators, retaining their quadratic E2 relation.
A = PolynomialRing(QQ, 'h1,k1,h2,k2,Jp,Jm,beta,X,Y', order='lex')
h1,k1,h2,k2,Jp,Jm,beta,X,Y = A.gens()
relations = [Jp*Jm-beta^2, h1*k1+beta, h2*k2+beta,
             X-h1*h2*Jm, Y-k1*k2*Jp]
eliminated = A.ideal(relations).elimination_ideal(A.gens()[:6])
assert eliminated == A.ideal(X*Y-beta^4)
# Removing the local instanton currents leaves no positive-dimensional
# reduced branch; the nonreduced remnant beta^2=0 is not a finite branch.
finite_coupling = A.ideal(relations + [Jp,Jm]).elimination_ideal(A.gens()[:6])
assert finite_coupling == A.ideal([X,Y,beta^2])
physical_basis = A.ideal(relations).groebner_basis()

# Completeness: after removing h_i*k_i and Jp*Jm pairs, charge neutrality
# forces a monomial to be a power of X or Y, times a power of beta.
for exponents in product(range(4), repeat=6):
    a,b,c,d,p,m = exponents
    if a-b+p-m != 0 or c-d+p-m != 0:
        continue
    n = p-m
    if n >= 0:
        beta_power = 2*min(p,m) + min(a,b) + min(c,d)
        reduced = (-1)^(min(a,b)+min(c,d))*Y^n*beta^beta_power
    else:
        beta_power = 2*min(p,m) + min(a,b) + min(c,d)
        reduced = (-1)^(min(a,b)+min(c,d))*X^(-n)*beta^beta_power
    monomial = h1^a*k1^b*h2^c*k2^d*Jp^p*Jm^m
    assert (monomial-reduced).reduce(physical_basis) == 0

# 3. Independent integral auxiliary hyperkahler quotient: four hypers / U(1)^3.
# It is only a presentation of the ring, not a claimed 5d UV Lagrangian.
Qaux = matrix(ZZ, [[1,0,0,1], [0,1,0,1], [0,0,1,1]])
assert Qaux.smith_form()[0].diagonal() == [1,1,1]
assert Qaux*vector(ZZ, [-1,-1,-1,1]) == 0
B = PolynomialRing(QQ, 'q1,k1,q2,k2,q3,k3,q4,k4,beta,X,Y', order='lex')
q1,k1,q2,k2,q3,k3,q4,k4,beta,X,Y = B.gens()
aux_relations = [q1*k1+q4*k4, q2*k2+q4*k4, q3*k3+q4*k4,
                 beta-q4*k4, X+q1*q2*q3*k4, Y-k1*k2*k3*q4]
assert B.ideal(aux_relations).elimination_ideal(B.gens()[:8]) == B.ideal(X*Y-beta^4)

# 4. Independent finite Molien expansion in the PHYSICAL seven generators.
# Each key is (2R-degree, gauge charge 1, gauge charge 2, instanton grading).
# Numerator (1-t^4)(1-t^2)^2 removes the E2 and two moment-map relations.
cutoff = 16
series = {(0,0,0,0): ZZ(1)}
for i, degree in enumerate(degrees):
    q1,q2 = charges.column(i)
    n = 1 if i == 4 else (-1 if i == 5 else 0)
    expanded = defaultdict(ZZ)
    for (d,c1,c2,j), coefficient in series.items():
        for k in range((cutoff-d)//degree+1):
            expanded[(d+k*degree,c1+k*q1,c2+k*q2,j+k*n)] += coefficient
    series = expanded
constant_term = defaultdict(ZZ)
for (d,c1,c2,n), coefficient in series.items():
    if c1 == c2 == 0:
        constant_term[(d,n)] += coefficient
for shift in (4,2,2):
    updated = defaultdict(ZZ, constant_term)
    for (d,n), coefficient in constant_term.items():
        if d+shift <= cutoff:
            updated[(d+shift,n)] -= coefficient
    constant_term = updated
constant_term = {key: value for key,value in constant_term.items() if value}
normal_forms = {(4*abs(n)+2*j,n): ZZ(1)
                for n in range(-cutoff//4,cutoff//4+1)
                for j in range((cutoff-4*abs(n))//2+1)}
assert constant_term == normal_forms

# Exact rational equality: A3 complete intersection = dressed-instanton sum.
F = FractionField(PolynomialRing(QQ, 't,z'))
t,z = F.gens()
hilbert = (1-t^8)/((1-t^2)*(1-z*t^4)*(1-t^4/z))
sector_sum = (1 + z*t^4/(1-z*t^4) + (t^4/z)/(1-t^4/z))/(1-t^2)
assert hilbert == sector_sum

print('Effective compact charges primitive; J+ has charge R=C1+C2.')
print('Physical-current elimination: XY=beta^4 (A3, C^2/Z4).')
print('Independent integral four-hyper quotient: same ring; Smith invariants (1,1,1).')
print('Refined Molien coefficients agree through degree', cutoff, 'with all instanton sectors.')
print('H(t,z)=(1-t^8)/[(1-t^2)(1-z*t^4)(1-z^(-1)*t^4)].')
print('Unrefined coefficients at degrees 0,2,...,16:',
      [sum(v for (d,n),v in constant_term.items() if d == degree)
       for degree in range(0,cutoff+1,2)])
print('Nonzero beta requires both nodal bilinears and both instanton currents.')
print('All formal X9 dressed-instanton checks passed.')

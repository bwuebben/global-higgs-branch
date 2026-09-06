from random import Random
V = [(1,0,0,0),(0,1,0,0),(-1,-1,0,0),(0,0,1,0),(0,0,0,1),(-4,-2,-1,0),(-4,-2,0,-1),(3,1,1,1),(2,1,1,1)]
p = (-2,-1,0,0); f = (-1,0,0,0)
rays = [vector(ZZ,x) for x in V + [p, f]]
P = LatticePolytope(V); polar_pts = [vector(ZZ,m) for m in P.polar().points()]
fiber_idx = [0,1,2,9,10]; base_idx = [3,4,5,6,7,8]      # x0,x1,x2,e,w ; z3,z4,z5,z6,z7,z8
R = PolynomialRing(QQ, 'z3,z4,z5,z6,z7,z8'); z = R.gens()
rng = Random(int(20260902))
sections = {}
for m in polar_pts:
    fib = tuple(int(m*rays[i])+1 for i in fiber_idx); bas = tuple(int(m*rays[i])+1 for i in base_idx)
    coeff = QQ(rng.randint(-9, 9)) or QQ(1)
    sections[fib] = sections.get(fib, R(0)) + coeff * prod(z[j]**bas[j] for j in range(6))
# fiber monomial exponents (x0,x1,x2,e,w) -> role
role = {(2,0,1,0,0):'a0', (1,2,0,0,1):'b2', (1,1,1,1,1):'b1', (1,0,2,2,1):'b0',
        (0,3,0,1,2):'c3', (0,2,1,2,2):'c2', (0,1,2,3,2):'c1', (0,0,3,4,2):'c0'}
S = {role[k]: v for k, v in sections.items()}
assert len(S) == 8
def analyse(chart_name, subs, vars2):
    Rc = PolynomialRing(QQ, vars2); a, b = Rc.gens()
    sub = {}
    for name, val in subs.items():
        sub[R(name)] = val if val in (0,1) else (a if val == 'a' else b)
    s = {k: Rc(v.subs({R(n): (Rc(1) if vv == 1 else (a if vv == 'a' else b)) for n, vv in subs.items()})) for k, v in S.items()}
    A0, B0, B1, B2, C0, C1, C2, C3 = (s[k] for k in ('a0','b0','b1','b2','c0','c1','c2','c3'))
    # quartic Q(u,v) = q2^2 - 4 a0 u c(u,v):  coefficients of u^4..v^4
    A = B0**2 - 4*A0*C0; Bq = 2*B0*B1 - 4*A0*C1; C = B1**2 + 2*B0*B2 - 4*A0*C2; D = 2*B1*B2 - 4*A0*C3; E = B2**2
    I = 12*A*E - 3*Bq*D + C**2
    J = 72*A*C*E + 9*Bq*C*D - 27*A*D**2 - 27*E*Bq**2 - 2*C**3
    Dl = 4*I**3 - J**2
    def ord_along(poly, var_index):
        if poly == 0: return 'inf'
        return min(e[var_index] for e in poly.exponents())
    def mult0(poly):
        if poly == 0: return 'inf'
        return min(sum(e) for e in poly.exponents())
    print(f"chart {chart_name} (coords {vars2}): degrees I,J,Delta = {I.degree()},{J.degree()},{Dl.degree()}")
    for vi, vname in enumerate(vars2.split(',')):
        print(f"   along {vname}=0: ord(f,g,Delta) = ({ord_along(I,vi)}, {ord_along(J,vi)}, {ord_along(Dl,vi)})")
    print(f"   at the origin: mult(f,g,Delta) = ({mult0(I)}, {mult0(J)}, {mult0(Dl)})")
    return I, J, Dl
# chart around q4 = L0 ∩ L5 : z3=z4=z6=z8=1, coordinates (z5, z7)
analyse("q4", {'z3':1,'z4':1,'z6':1,'z8':1,'z5':'a','z7':'b'}, 'z5,z7')
# chart around q3 = L0 ∩ L6 : z3=z4=z5=z8=1, coordinates (z6, z7)
analyse("q3", {'z3':1,'z4':1,'z5':1,'z8':1,'z6':'a','z7':'b'}, 'z6,z7')
# chart around L5 ∩ L6 : z3=z4=z7=z8=1, coordinates (z5, z6)
analyse("L5∩L6", {'z3':1,'z4':1,'z7':1,'z8':1,'z5':'a','z6':'b'}, 'z5,z6')

print("\n=== matter counting ===")
Fp = GF(32003)
# chart z7 != 0 (P^2 minus L0): singular points of Delta = cusps (f=g=0) + I2 nodes
Rc = PolynomialRing(Fp, 'a,b'); a, b = Rc.gens()
def specialize(chartsubs):
    s = {k: Rc(v.subs({R(n): (R(1) if vv == 1 else (R('z5') if vv=='a' else R('z6'))) for n, vv in chartsubs.items()}).change_ring(Fp) if False else Rc(str(v.subs({R(n): (R(1) if vv == 1 else R(n)) for n, vv in chartsubs.items()})).replace('z5','a').replace('z6','b').replace('z7','b'))) for k, v in S.items()}
    return s
# simpler: rebuild with z7=1, z3=z4=z8=1, a=z5, b=z6
s = {k: Rc(str(v.subs({R('z3'):R(1), R('z4'):R(1), R('z8'):R(1), R('z7'):R(1)})).replace('z5','a').replace('z6','b')) for k, v in S.items()}
A0, B0, B1, B2, C0, C1, C2, C3 = (s[k] for k in ('a0','b0','b1','b2','c0','c1','c2','c3'))
A = B0**2 - 4*A0*C0; Bq = 2*B0*B1 - 4*A0*C1; C = B1**2 + 2*B0*B2 - 4*A0*C2; D = 2*B1*B2 - 4*A0*C3; E = B2**2
I = 12*A*E - 3*Bq*D + C**2; J = 72*A*C*E + 9*Bq*C*D - 27*A*D**2 - 27*E*Bq**2 - 2*C**3; Dl = 4*I**3 - J**2
sing = Rc.ideal([Dl, Dl.derivative(a), Dl.derivative(b)])
cusp = Rc.ideal([I, J])
print("chart z7!=0: deg Delta =", Dl.degree(), "; length of singular scheme of Delta:", sing.vector_space_dimension(), "; length of (f,g):", cusp.vector_space_dimension())
print("  => I2 nodes off L0 = length(sing) - 2*length(f,g) =", sing.vector_space_dimension() - 2*cusp.vector_space_dimension())
# on L0 (z7 = 0), chart z6 != 0: coordinates z5 (along L0); use q4-chart polynomials restricted to z7=0
Rq = PolynomialRing(Fp, 'a,b'); a, b = Rq.gens()
s4 = {k: Rq(str(v.subs({R('z3'):R(1), R('z4'):R(1), R('z8'):R(1), R('z6'):R(1)})).replace('z5','a').replace('z7','b')) for k, v in S.items()}
A0, B0, B1, B2, C0, C1, C2, C3 = (s4[k] for k in ('a0','b0','b1','b2','c0','c1','c2','c3'))
A = B0**2 - 4*A0*C0; Bq = 2*B0*B1 - 4*A0*C1; C = B1**2 + 2*B0*B2 - 4*A0*C2; D = 2*B1*B2 - 4*A0*C3; E = B2**2
I = 12*A*E - 3*Bq*D + C**2; J = 72*A*C*E + 9*Bq*C*D - 27*A*D**2 - 27*E*Bq**2 - 2*C**3; Dl = 4*I**3 - J**2
Dl2 = Rq(Dl // b**2) if Dl % b**2 == 0 else None
assert Dl2 is not None and Dl2 % b != 0
U = PolynomialRing(Fp, 'a'); t = U.gen()
d_on = U(str(Dl2.subs({b: 0})).replace('a','t')) if False else Dl2.subs({b: Rq(0)}).univariate_polynomial()
f_on = I.subs({b: Rq(0)}).univariate_polynomial(); g_on = J.subs({b: Rq(0)}).univariate_polynomial()
print("on L0 (chart z6!=0, coordinate z5; q4 at z5=0): deg Delta'|L0 =", d_on.degree(), " factorization multiplicities:", sorted([m for _, m in d_on.factor()], reverse=True)[:8], " total roots (with mult) in this chart:", sum(m*fa.degree() for fa, m in d_on.factor()))
print("   ord_{z5} at q4 of (f,g,Delta') on L0:", f_on.valuation(), g_on.valuation(), d_on.valuation())
print("   gcd(f|L0, g|L0) degree (III/IV points on L0):", f_on.gcd(g_on).degree(), " gcd(Delta'|L0, f|L0) degree:", d_on.gcd(f_on).degree())

#!/usr/bin/env sage
"""Fibration structure of Xhat_9 used in the heterotic/F-theory comparison.

Asserts, with exact rational arithmetic:
  * the K3 fibration Xhat_9 -> P^1 (m = -e_3^*): fiber polytope reflexive with
    8 points, generic-fiber Picard rank 3, toric fiber lattice U + <-4>;
  * the genus-one fibration phi_w, w = D4 + D5 (w^3 = 0, w^2 = F): the toric
    projection along span(v0, v1) is NOT a fan, the morphism target is P^2;
    E and D2 are sections, D0 a 3-section, D1 a 2-section; D3 and D4 lie over
    points (non-flat fibers) and are rational with K^2 = 1 and 0; D7, D8 lie
    over one line with a two-component fiber; D5, D6 lie over lines;
  * the incidences: C2 and R_{v4} are disjoint (-1)-curves of D4|X, C1 meets
    D3 once, R_{v3} is a (-1)-curve of D3|X, and R = C1 + C2 in H_2.
Run: sage -c "load('checks/x9_fibrations.sage')"
"""
from itertools import combinations, combinations_with_replacement, product
from collections import Counter

V = [(1,0,0,0),(0,1,0,0),(-1,-1,0,0),(0,0,1,0),(0,0,0,1),(-4,-2,-1,0),(-4,-2,0,-1),(3,1,1,1),(2,1,1,1)]
pts = [vector(ZZ, x) for x in V + [(-2,-1,0,0), (-1,0,0,0)]]
NR = 11
P = Polyhedron(vertices=V)
facets = []
for ineq in P.inequalities():
    n = -vector(ineq.A()); assert ineq.b() == 1
    facets.append((n, [i for i, x in enumerate(pts) if n * x == 1]))
H = {0: 10, 1: 11, 2: 3, 3: 9, 4: 2, 5: 7, 6: 1, 7: 4, 8: 8, 9: -5, 10: -3}
H = {i: QQ(h) + QQ(i + 1) / QQ(1009) for i, h in H.items()}
simplices = set()
for n, on in facets:
    if len(on) == 4:
        simplices.add(tuple(sorted(on))); continue
    j = next(k for k in range(4) if n[k] != 0); keep = [k for k in range(4) if k != j]
    lifted = [tuple([pts[i][k] for k in keep] + [H[i]]) for i in on]
    Q = Polyhedron(vertices=lifted)
    for Fc in Q.faces(3):
        ieq = [e for e in Fc.ambient_Hrepresentation() if e.is_inequality()]
        if -vector(ieq[0].A())[-1] < 0:
            simplices.add(tuple(sorted(on[lifted.index(tuple(v))] for v in Fc.vertices())))
simplices = sorted(simplices)
assert len(simplices) == 29
fan = Fan(cones=[list(s) for s in simplices], rays=pts, check=False)
X = ToricVariety(fan)
HH = X.cohomology_ring(); D = [HH(X.divisor(i)) for i in range(NR)]; aK = sum(D)
kap = {}
for i, j, k in combinations_with_replacement(range(NR), 3):
    kap[(i, j, k)] = X.integrate(D[i] * D[j] * D[k] * aK)
def kappa(i, j, k): return kap[tuple(sorted((i, j, k)))]
BASIS = [2, 5, 6, 7, 8, 9]
CL = {'D0': [1,4,4,-3,-2,2], 'D1': [1,2,2,-1,-1,1], 'D2': [1,0,0,0,0,0], 'D3': [0,1,0,-1,-1,0],
      'D4': [0,0,1,-1,-1,0], 'D5': [0,1,0,0,0,0], 'D6': [0,0,1,0,0,0], 'D7': [0,0,0,1,0,0],
      'D8': [0,0,0,0,1,0], 'E': [0,0,0,0,0,1]}
CL = {k: vector(QQ, v) for k, v in CL.items()}
def trip(a, b, c):
    return sum(kappa(BASIS[I], BASIS[J], BASIS[K]) * a[I] * b[J] * c[K]
               for I in range(6) for J in range(6) for K in range(6))
def curve(a, b):
    return vector(QQ, [trip(a, b, vector(QQ, [1 if k == K else 0 for k in range(6)])) for K in range(6)])
unit = lambda i: vector(QQ, [kappa(i, i, b) for b in BASIS])   # not used

# ---- K3 fibration ---------------------------------------------------------
m = vector(ZZ, (0, 0, -1, 0))
assert not any(any(m * pts[i] > 0 for i in s) and any(m * pts[i] < 0 for i in s) for s in simplices)
slice_idx = [i for i, x in enumerate(pts) if x[2] == 0]
LP = LatticePolytope([(pts[i][0], pts[i][1], pts[i][3]) for i in slice_idx])
assert LP.is_reflexive() and LP.npoints() == 8 and LP.polar().npoints() == 34
facet_int = sum(len(Fc.interior_points()) for Fc in LP.facets())
edge_corr = sum(len(e.interior_points()) * len(e.dual().interior_points()) for e in LP.edges())
assert (facet_int, edge_corr) == (1, 0)
assert LP.npoints() - 4 - facet_int + edge_corr == 3          # Picard rank of the generic K3 fiber
G = matrix(QQ, [[trip(CL[a], CL[b], CL['D5']) for b in ('D2', 'D6', 'E')] for a in ('D2', 'D6', 'E')])
assert G == matrix(QQ, [[-2, 1, 0], [1, 0, 1], [0, 1, -2]]) and G.det() == 4
s, f, sp = (vector(QQ, e) for e in ([1,0,0], [0,1,0], [0,0,1]))
v = sp - s - 2 * f
assert v * G * v == -4 and v * G * f == 0 and v * G * s == 0     # <s,f,s'> = U + <-4>
assert (s + f) * G * (s + f) == 0 and f * G * (s + f) == 1
Gall = matrix(QQ, [[trip(CL[a], CL[b], CL['D5']) for b in CL] for a in CL])
assert Gall.rank() == 3

# ---- genus-one fibration ----------------------------------------------------
w = CL['D4'] + CL['D5']
assert w == CL['D3'] + CL['D6']
assert trip(w, w, w) == 0
F = curve(w, w)
assert F == curve(CL['D4'], CL['D5']) == vector(QQ, [1, 0, 0, 0, 0, 1])
deg = {k: F * c for k, c in CL.items()}
assert deg == {'D0': 3, 'D1': 2, 'D2': 1, 'E': 1, 'D3': 0, 'D4': 0, 'D5': 0, 'D6': 0, 'D7': 0, 'D8': 0}
# the toric projection along span(v0, v1) does not define a fan
Lsat = (ZZ**4).submodule([pts[0], pts[1]]).saturation()
Qm = (ZZ**4).quotient(Lsat)
cones2 = set()
for sg in simplices:
    pr = [vector(ZZ, Qm(pts[k])) for k in sg if pts[k] not in Lsat]
    cones2.add(tuple(sorted(tuple(r) for r in Cone(pr, lattice=ZZ**2).rays())))
assert ((0, -1), (1, 1)) in cones2 and ((0, -1), (1, 0)) in cones2      # overlapping cones
try:
    Fan([Cone([vector(ZZ, r) for r in c], lattice=ZZ**2) for c in cones2], check=True)
    raise AssertionError("projected cones unexpectedly form a fan")
except ValueError:
    pass
# fiber polygon: 4 vertices, 5 boundary points, dual edge lengths (1,1,3,2)
Bm = matrix(ZZ, Lsat.basis())
poly2 = LatticePolytope([tuple(int(x) for x in Bm.solve_left(vector(QQ, pts[k]))) for k in range(NR) if pts[k] in Lsat])
assert poly2.is_reflexive() and poly2.npoints() == 6 and poly2.nvertices() == 4
pol = poly2.polar(); assert pol.npoints() == 8
lengths = sorted(len([mm for mm in pol.points() if vector(mm) * vector(vt) == -1]) - 1 for vt in poly2.vertices())
assert lengths == [1, 1, 2, 3]
# divisors over points versus over curves
assert curve(w, CL['D3']) == 0 and curve(w, CL['D4']) == 0            # D3, D4 over points
for k in ('D5', 'D6', 'D7', 'D8'):
    assert curve(w, CL[k]) != 0 and trip(w, w, CL[k]) == 0             # over curves
assert curve(w, CL['D5']) == F and curve(w, CL['D6']) == F              # D5, D6 over lines, flat there
assert curve(w, CL['D7']) * CL['D2'] == 1 and curve(w, CL['D7']) * CL['E'] == 0
assert curve(w, CL['D8']) * CL['E'] == 1 and curve(w, CL['D8']) * CL['D2'] == 0
assert curve(CL['D7'], CL['D8']) != 0 and kappa(7, 7, 8) == 4 and kappa(7, 8, 8) == -6   # I2-type double curve
assert kappa(7, 8, 8) + kappa(7, 7, 8) == -2                              # genus 0
# surfaces over points: D3 rational K^2 = 1, D4 rational elliptic K^2 = 0
c2 = HH(X.Chern_class(2))
chi = lambda i: X.integrate(c2 * D[i] * aK) + kappa(i, i, i)
assert (kappa(3, 3, 3), chi(3)) == (1, 11) and (kappa(4, 4, 4), chi(4)) == (0, 12)
assert kappa(4, 4, 5) == 0 and kappa(4, 5, 5) == 0                        # K3 slices of D4 are its elliptic fibers
assert all(kappa(3, 4, b) == 0 for b in BASIS)                            # D3 and D4 disjoint
# incidences of the transition curves
C1 = vector(QQ, [kappa(b, 6, 7) for b in BASIS]); C2 = vector(QQ, [kappa(b, 2, 4) for b in BASIS])
Rv3 = vector(QQ, [kappa(b, 3, 9) for b in BASIS]); Rv4 = vector(QQ, [kappa(b, 4, 9) for b in BASIS])
assert Rv4 - Rv3 == C1 + C2                                               # R = C1 + C2
assert C2 * CL['D5'] == 1 and C1 * CL['D5'] == 0                          # C2 section, C1 vertical for the K3 fibration
assert C1 * CL['D3'] == 1 and C1 * CL['D4'] == 0
assert kappa(2, 2, 4) == -1 and kappa(4, 9, 9) == -1 and kappa(2, 4, 9) == 0   # C2, R_v4 disjoint (-1)-curves in D4|X
assert kappa(3, 3, 9) == -1                                               # R_v3 a (-1)-curve of D3|X
assert (CL['E'] - CL['D2']) * F == 0 and trip(CL['E'] - CL['D2'], CL['E'] - CL['D2'], w) == -6
print("Xhat_9 fibration facts: all assertions passed.")
print("K3 fiber lattice: U + <-4>, Picard rank 3;  genus-one fibration over P^2 with sections E, D2,")
print("non-flat fibers D3|X (K^2=1) and D4|X (K^2=0) over two points on the I2 line pi(D7) = pi(D8).")

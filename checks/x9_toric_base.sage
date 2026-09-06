# The six-dimensional base of X_9 is toric, and the flat model is the flop of C_1.
#
# The projection N -> Z^2 onto the last two coordinates sends the rays of the
# resolved fan of Xhat_9 to the fan of S_7 = Bl_2 P^2 (rays (1,0),(0,1),(-1,0),
# (0,-1),(1,1)) or to 0.  The MPCP fan used in the paper (diagonal v6v7 in the
# square {v2,v3,v6,v7}) has exactly four maximal cones that do not map into a
# cone of S_7, all containing the curve C_1 = D6.D7; replacing that diagonal by
# v2v3 (the flop of C_1) gives a smooth complete fan mapping cone-to-cone, i.e.
# a toric morphism to S_7 whose generic fiber is the toric surface with rays
# (1,0),(0,1),(-1,-1),(-2,-1),(-1,0) (a dP_7), so that the anticanonical
# hypersurface is a flat elliptic fibration over S_7 with sections E and D_2.
import sys, io
_stdout = sys.stdout; sys.stdout = io.StringIO()
load('checks/check_x9_prepotential.sage')      # defines pts (11 rays) and simplices (29 cones)
sys.stdout = _stdout

proj = lambda v: (v[2], v[3])
S7_rays  = [(1,0),(1,1),(0,1),(-1,0),(0,-1)]                 # cyclic order
S7_cones = [(S7_rays[i], S7_rays[(i+1) % 5]) for i in range(5)]
S7 = Fan(cones=[[i, (i+1) % 5] for i in range(5)], rays=S7_rays, check=True)
assert S7.is_smooth() and S7.is_complete()
# self-intersections: three (-1)-curves in a chain and two 0-curves
selfint = []
for i in range(5):
    a, b, c = vector(S7_rays[i-1]), vector(S7_rays[i]), vector(S7_rays[(i+1) % 5])
    k = [m for m in range(-3, 4) if a + c == m * b]
    selfint.append(-k[0])
assert selfint == [-1, -1, -1, 0, 0], selfint
# linear equivalences: R3 ~ R1 + R5, R4 ~ R2 + R5 (characters e1*, e2*)
R = {r: i for i, r in enumerate(S7_rays)}
assert [vector(r)[0] for r in S7_rays] == [1, 1, 0, -1, 0]    # e1*: R(1,0)+R(1,1)-R(-1,0) = 0
assert [vector(r)[1] for r in S7_rays] == [0, 1, 1, 0, -1]    # e2*: R(1,1)+R(0,1)-R(0,-1) = 0

def in_one_cone(vs):
    vs = [v for v in vs if v != (0, 0)]
    return (not vs) or any(all(Cone(list(c)).contains(vector(v)) for v in vs) for c in S7_cones)
def bad_cones(simps):
    return sorted(s for s in simps if not in_one_cone([proj(pts[i]) for i in s]))

simps = set(tuple(sorted(s)) for s in simplices)
assert len(simps) == 29
bad = bad_cones(simps)
assert len(bad) == 4 and all({6, 7} <= set(s) for s in bad), bad
print("paper fan: 4 maximal cones fail to map into S_7, all containing C_1 = D6.D7:", bad)

# flop C_1: bistellar flip of the square {v2,v3,v6,v7} from diagonal v6v7 to v2v3
with67 = [s for s in simps if {6, 7} <= set(s)]
flopped = simps - set(with67)
apexes = set()
for s in with67:
    rest = set(s) - {6, 7}
    assert len(rest & {2, 3}) == 1
    apexes |= rest - {2, 3}
for x in apexes:
    flopped.add(tuple(sorted((2, 3, 6, x)))); flopped.add(tuple(sorted((2, 3, 7, x))))
fan2 = Fan(cones=[list(s) for s in flopped], rays=pts, check=True)
assert len(flopped) == 29 and fan2.is_complete() and fan2.is_smooth()
assert bad_cones(flopped) == []
print("flopped fan (diagonal v2v3): smooth, complete, 29 cones, every cone maps into a cone of S_7")

# images of the ray divisors
img = {i: proj(pts[i]) for i in range(len(pts))}
assert img[3] == (1, 0) and img[4] == (0, 1) and img[5] == (-1, 0) and img[6] == (0, -1)
assert img[7] == img[8] == (1, 1)
kernel = sorted(i for i in range(len(pts)) if img[i] == (0, 0))
assert kernel == [0, 1, 2, 9, 10]
fib_rays = sorted((pts[i][0], pts[i][1]) for i in kernel)
assert fib_rays == sorted([(1, 0), (0, 1), (-1, -1), (-2, -1), (-1, 0)])
# the generic ambient fiber: cones of the flopped fan contained in the kernel
fib_cones = sorted(set(tuple(sorted(set(s) & set(kernel))) for s in flopped if len(set(s) & set(kernel)) == 2))
F = Fan(cones=[[kernel.index(i) for i in c] for c in fib_cones], rays=[(pts[i][0], pts[i][1]) for i in kernel], check=True)
assert F.is_complete() and F.is_smooth() and len(fib_cones) == 5
print("generic ambient fiber: smooth complete toric surface with 5 rays (dP_7); anticanonical curves are elliptic")
print("base: D3 -> (1,0) [-1], D4 -> (0,1) [-1], D7,D8 -> (1,1) [-1], D5 -> (-1,0) [0] ~ (1,0)+(1,1), D6 -> (0,-1) [0] ~ (0,1)+(1,1)")
print("ALL CHECKS PASSED")

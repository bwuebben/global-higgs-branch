"""Exact E-string root identification and finite-circle transition of X_9.

Run from the project root: sage checks/check_x9_heterotic_transition.sage
Reconstructs both toric intersection rings; checks the flop formula, the two
E8 roots (including the elliptic-fiber subtraction), the complete Shioda and
KK divisors, the transition-face coordinates in both zero-section choices,
and the bound 7/24 <= Vol(X)/(J.F)^3 <= 83/162 on that face.
The charge computation does not assert that a root is an effective curve.
"""
import contextlib
import io

with contextlib.redirect_stdout(io.StringIO()):
    load('checks/x9_toric_base.sage')

# Recompute the intersection ring after the flop, independently of its formula.
flat_ambient = ToricVariety(fan2)
flat_ring = flat_ambient.cohomology_ring()
flat_divs = [flat_ring(flat_ambient.divisor(i)) for i in range(NR)]
flat_ac = sum(flat_divs)
flat_kap = {}
for ii, jj, kk in combinations_with_replacement(range(6), 3):
    flat_kap[(ii, jj, kk)] = flat_ambient.integrate(
        flat_divs[BASIS[ii]] * flat_divs[BASIS[jj]] * flat_divs[BASIS[kk]] * flat_ac)

def flat_tri(u, v, w):
    return sum(u[ii]*v[jj]*w[kk]*flat_kap[tuple(sorted((ii,jj,kk)))]
               for ii in range(6) for jj in range(6) for kk in range(6))

def flat_curve(u, v):
    return vector([flat_tri(u, v, CLS[bn]) for bn in BN])

for ii, jj, kk in combinations_with_replacement(range(6), 3):
    assert flat_kap[(ii,jj,kk)] == kappa(BASIS[ii],BASIS[jj],BASIS[kk]) - C1[ii]*C1[jj]*C1[kk]

Db = CLS
pull_H = Db['D5'] + Db['D4']
pull_c1 = 3*pull_H - Db['D3'] - Db['D4']
base_classes = [pull_H, Db['D3'], Db['D4']]
elliptic_fiber = flat_curve(pull_H,pull_H)
print('fiber', elliptic_fiber)
assert elliptic_fiber == vector(QQ,[1,0,0,0,0,1])
roots = {}
sections = {}
for label in ['D3','D4']:
    surface = Db[label]
    o = flat_curve(Db['D2'],surface)
    s = flat_curve(Db['E'],surface)
    ff = -flat_curve(surface,surface)
    root = s-o-ff
    roots[label] = root
    sections[label] = (o,s,ff)
    assert flat_tri(surface,surface,surface) == 0
    assert ff == elliptic_fiber
    assert (flat_tri(surface,Db['D2'],Db['D2']),
            flat_tri(surface,Db['E'],Db['E']),
            flat_tri(surface,Db['D2'],Db['E'])) == (-1,-1,0)
    assert -flat_tri(surface,surface,Db['D2']) == 1
    assert -flat_tri(surface,surface,Db['E']) == 1
    # On the rational elliptic surface F = -K = -surface|surface.
    root_div = Db['E']-Db['D2']+surface
    assert flat_tri(surface,root_div,root_div) == -2
    assert flat_tri(surface,root_div,Db['D2']) == 0
    assert flat_tri(surface,root_div,-surface) == 0
    print(label, 'K^2',flat_tri(surface,surface,surface), 'o',o,'s',s,'F',ff,'root',root)
    print('section intersections',flat_tri(surface,Db['D2'],Db['D2']),
          flat_tri(surface,Db['E'],Db['E']),flat_tri(surface,Db['D2'],Db['E']))

assert sections['D3'][0] == -C1
assert sections['D4'][0] == C2
assert roots['D3'] == roots['D4'] == vector(QQ,[0,0,0,-1,1,-2])
common_root = roots['D3']
assert sections['D4'][1]-sections['D3'][1] == C1+C2

sigma = Db['E']-Db['D2']-pull_c1+Db['D8']/2
kk_div = Db['D2']+pull_c1/2
assert sigma == vector(QQ,[-1,-2,-2,1,3/2,1])
assert [flat_tri(kk_div,base_classes[i],base_classes[j]) for i in range(3) for j in range(3)] == [1,0,0,0,-1,0,0,0,-1]
assert all(flat_tri(kk_div,kk_div,v) == 0 for v in base_classes)
assert all(flat_tri(sigma,kk_div,v) == 0 and flat_tri(sigma,Db['D8'],v) == 0 for v in base_classes)
assert all(flat_tri(sigma,u,v) == 0 for u in base_classes for v in base_classes)
assert [-flat_tri(sigma,sigma,v) for v in base_classes] == [11/2,3/2,3/2]
# The flop is MW-neutral; the shifted KK divisor has half-integral pairing.
assert sigma*C1 == 0 and kk_div*C1 == 1/2
print('sigma', sigma, 'height against H,E3,E4',[-flat_tri(sigma,sigma,v) for v in base_classes])
print('KK',kk_div)
frame = matrix([kk_div]+base_classes+[sigma,Db['D8']]).transpose()
aa,bb = var('aa bb')
face_J = aa*Db['D0']+bb*Db['D1']
coords = frame.solve_right(face_J)
assert coords == vector([3*aa+2*bb,13*aa/2+3*bb,-3*aa/2-bb,-3*aa/2-bb,2*aa+bb,-bb/2])
print('face frame (fiber,H,E3,E4,MW,SU2)',coords)
assert flat_tri(face_J,face_J,Db['E']).expand() == 0
assert flat_tri(face_J,face_J,Db['D2']).expand() == 4*aa**2
assert face_J*elliptic_fiber == 3*aa+2*bb
assert all(flat_tri(face_J,Db['E'],v).expand() == 0 for v in base_classes)
assert (flat_tri(face_J,face_J,face_J)-trip(face_J,face_J,face_J)).expand() == 0
print('face volume',flat_tri(face_J,face_J,face_J).expand()/6)
print('E volume',flat_tri(face_J,face_J,Db['E']).expand()/2)
print('D2 volume',flat_tri(face_J,face_J,Db['D2']).expand()/2)
print('fiber volume',face_J*elliptic_fiber)
print('vertical surface volumes',[flat_tri(face_J,face_J,v).expand()/2 for v in base_classes])
print('restriction to E against base curves',[flat_tri(face_J,Db['E'],v).expand() for v in base_classes])
sigma_E = Db['D2']-Db['E']-pull_c1+Db['D7']/2
kk_E = Db['E']+pull_c1/2
frame_E = matrix([kk_E]+base_classes+[sigma_E,Db['D7']]).transpose()
coords_E = frame_E.solve_right(face_J)
assert coords_E == vector([3*aa+2*bb,9*aa/2+3*bb,-3*aa/2-bb,-3*aa/2-bb,aa+bb,-3*aa/2-bb/2])
assert all(flat_tri(kk_E,kk_E,v) == 0 for v in base_classes)
assert all(flat_tri(sigma_E,kk_E,v) == 0 and flat_tri(sigma_E,Db['D7'],v) == 0 for v in base_classes)
assert [-flat_tri(sigma_E,sigma_E,v) for v in base_classes] == [11/2,3/2,3/2]
print('E-zero frame',coords_E)
xx = var('xx')
ratio = (flat_tri(face_J,face_J,face_J)/6).subs({aa:xx,bb:(1-3*xx)/2}).expand()
assert ratio == 7/24+3*xx/4-xx**2/8-5*xx**3/12
assert ratio.subs(xx=0) == 7/24 and ratio.subs(xx=1/3) == 83/162
# Its derivative decreases on [0,1/3] and is positive even at the right end.
assert diff(ratio,xx,2) == -1/4-5*xx/2
assert diff(ratio,xx).subs(xx=1/3) == 19/36
print('V/t^3, x=a/t in [0,1/3]:',ratio,'derivative',diff(ratio,xx).factor())
gauge_frame = [kk_div,sigma,Db['D8']]
assert [common_root*v for v in gauge_frame] == [0,-3/2,1]
assert [(common_root+elliptic_fiber)*v for v in gauge_frame] == [1,-3/2,1]
assert all(common_root*v == 0 for v in base_classes)
assert (face_J*(common_root+elliptic_fiber)).expand() == 0
assert (coords[0]-3*coords[4]/2+coords[5]).expand() == 0
print('r charge KK,MW,Cartan', [common_root*v for v in gauge_frame])
print('r+F charge KK,MW,Cartan', [(common_root+elliptic_fiber)*v for v in gauge_frame])
print('o3 charge KK,MW,Cartan',[flat_curve(Db['D2'],Db['D3'])*v for v in [kk_div,sigma,Db['D8']]])
print('o4 charge KK,MW,Cartan',[flat_curve(Db['D2'],Db['D4'])*v for v in [kk_div,sigma,Db['D8']]])
assert 11+12+1 == 24
assert 30*(11+12)-2*248+20+1 == 215
print('ALL CHECKS PASSED')

"""Smooth-side quantum periods and finite-order neutral-period comparison.

Run: sage -c "load('checks/check_x9_neutral_quantum_periods.sage')"
All local degrees are included at each tested spectator degree. This is not
an all-orders identification of the full integral period vector. The
seven-period all-orders proof is separate; its finite algebraic data are
checked by check_x9_neutral_period_descent.sage.
"""
import contextlib
import io
import os
os.environ.setdefault('MIRROR_DEGREE','2')
os.environ.setdefault('MAGNETIC_BULK_DEGREE','4')
assert int(os.environ['MAGNETIC_BULK_DEGREE'])>=4
with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_regular_magnetic_periods.sage')

smooth_rays=[(-1,-3,-1,0),(-1,-3,0,-1),(-1,-1,0,0),
             (-1,0,0,0),(0,1,0,0),(1,-1,0,0),(1,3,1,1)]
Qsm=matrix(ZZ,[[1,1,-2,0,0,-1,1],
              [0,0,-1,2,0,1,0],
              [0,0,1,-1,1,0,0]])
assert Qsm*matrix(ZZ,smooth_rays)==0
assert Qsm.row_module()==matrix(ZZ,smooth_rays).left_kernel()
assert sum(Qsm.columns())==vector(ZZ,[0,2,1])
smooth_kappa=[matrix(QQ,3,3,lambda j,k:nef_kappa[i+1][j+1,k+1])
              for i in range(3)]

# Triangulate the single nonsimplicial cone in either direction. The
# exceptional locus lies over an ambient fixed point missed by the CY.
fixed_cones=[(1,4,5,6),(1,3,4,6),(0,4,5,6),(0,3,4,6),
             (0,1,4,5),(0,1,3,4),(1,2,3,6),(0,2,3,6),(0,1,2,3)]
singular_cone=(0,1,2,5,6)
circuit=matrix(ZZ,[smooth_rays[i] for i in singular_cone]).left_kernel().basis()[0]
print('Ambient circuit:',circuit)
ambient_mori=[]
for direction in [1,-1]:
    simplices=fixed_cones+[tuple(i for i in singular_cone if i!=omit)
        for omit,c in zip(singular_cone,circuit) if direction*c>0]
    sfan=Fan(cones=simplices,rays=smooth_rays,check=True)
    assert sfan.is_complete() and sfan.is_simplicial()
    sx=ToricVariety(sfan)
    sr=sx.cohomology_ring()
    sd=[sr(sx.divisor(i)) for i in range(7)]
    # Qsm columns 0,4,5 form an integral basis; solve for H2,H3,H4.
    basis_indices=[0,4,5]
    basis_matrix=Qsm.matrix_from_columns(basis_indices)
    assert abs(basis_matrix.det())==1
    sh=[sum(basis_matrix.inverse()[j,i]*sd[basis_indices[j]] for j in range(3))
        for i in range(3)]
    anticanonical=sum(sd)
    for i,j,k in combinations_with_replacement(range(3),3):
        assert sx.integrate(sh[i]*sh[j]*sh[k]*anticanonical)==smooth_kappa[i][j,k]
    sc2=[sx.integrate(sx.Chern_class(2)*sh[i]*anticanonical) for i in range(3)]
    assert sc2==[sum(c2D[j]*nef[j,i+1] for j in range(6)) for i in range(3)]
    assert sx.integrate(sd[2]^3*anticanonical)==9
    assert [sx.integrate(sh[i]*sd[2]^2*anticanonical) for i in range(3)]==[-3,-3,0]
    walls={}
    for simplex in simplices:
        for wall in combinations(simplex,3):
            walls.setdefault(wall,[]).append(next(i for i in simplex if i not in wall))
    curve_rows=[]
    for wall,opposite in walls.items():
        indices=list(wall)+opposite
        relation=matrix(QQ,[smooth_rays[i] for i in indices]).T.right_kernel().basis()[0]
        if relation[3]<0:
            relation=-relation
        assert relation[3]>0 and relation[4]>0
        extended=zero_vector(QQ,7)
        for i,value in zip(indices,relation):
            extended[i]=value
        curve_rows.append(Qsm.T.solve_right(extended))
    mori=Polyhedron(rays=curve_rows)
    ambient_mori.append(mori)
    print('Ambient Mori rays:',[list(ray) for ray in mori.rays()])
    print('Smooth-side c2:',sc2)
common_mori=ambient_mori[0].intersection(ambient_mori[1])
assert ambient_mori[0]==Polyhedron(rays=[(1,0,0),(0,1,0),(0,0,1)])
assert ambient_mori[1]==Polyhedron(rays=[(-1,0,0),(1,1,0),(0,0,1)])
assert common_mori==Polyhedron(rays=[(1,1,0),(0,1,0),(0,0,1)])
print('Common Mori rays:',[list(ray) for ray in common_mori.rays()])
# All-degree pure-v vanishing: for (ell,0,0), ell>0, exactly columns
# 2 and 5 have negative arguments. Their divisor intersection vanishes
# against every K_i, so the contracted rational Hessian is zero as well.
assert all(Qsm.column(2)*smooth_kappa[i]*Qsm.column(5)==0 for i in range(3))
assert all(ray[1]+ray[2]>0 for ray in common_mori.rays())

# Independent seven-monomial Gamma coefficients versus the old ten-monomial
# local-degree sums, using precisely the same rational Hessian convention.
for n,(old_c,old_g,old_d) in endpoint_data.items():
    c,g,h=gamma_coefficients(Qsm,n)
    d=vector(QQ,[sum(smooth_kappa[i][j,k]*h[j,k] for j in range(3) for k in range(3))/2
                 for i in range(3)])
    assert c==old_c, ('fundamental',n)
    assert g==vector(old_g[1:4]), ('electric',n,g,old_g)
    assert d==vector(old_d[1:4]), ('neutral magnetic',n,d,old_d)
assert all(pair_charges[pair][i]==0 for pair in active_pairs for i in [1,2,3])
print('Exact neutral comparison:',len(endpoint_data),'spectator coefficients, all local degrees.')
print('Smooth Gamma recurrence tests:',check_contiguity(Qsm,4))
smooth_degree=int(os.environ.get('SMOOTH_GV_DEGREE','8'))
assert smooth_degree>=8
smooth_data=periods_and_gv(Qsm,smooth_kappa,smooth_degree)
nonzero={n:value for n,value in smooth_data['gv'].items() if value}
assert all(common_mori.contains(n) for n in nonzero)
expected_low={(0,0,1):10,(1,1,0):3,(0,1,1):140,(1,1,1):-20,(0,1,2):140}
assert {n:v for n,v in nonzero.items() if sum(n)<=3}==expected_low
assert sum(sum(n)<=8 for n in nonzero)==31
# Independent noncompact local-P2 one-parameter Gamma calculation.
# Its divisor dual to a line has cubic -1/3. Restricting the compact
# mirror logarithm to the extremal ray means selecting (d,d,0), not
# setting v=a and inadvertently retaining transverse curve degrees.
local_plane=periods_and_gv(matrix(ZZ,[[1,1,1,-3]]),[matrix(QQ,[[-QQ(1)/3]])],4)
assert [local_plane['gv'][(d,)] for d in range(1,5)]==[3,-6,27,-192]
for d in range(1,5):
    assert smooth_data['gv'][(d,d,0)]==local_plane['gv'][(d,)]
    assert (smooth_data['mirror_log'][0]+smooth_data['mirror_log'][1])[(d,d,0)]==local_plane['mirror_log'][0][d]
print('Smooth genus-zero GV:',len(smooth_data['gv']),'classes;',len(nonzero),'nonzero.')
print('Nonzero GV through degree three:',[(n,v) for n,v in nonzero.items() if sum(n)<=3])
print('Full nonzero GV:',sorted(nonzero.items()))
print('Plane sequence:',[smooth_data['gv'][(d,d,0)] for d in range(1,smooth_degree//2+1)])
print('ALL NEUTRAL QUANTUM CHECKS PASSED; finite-order comparison, not a full period-frame proof.')

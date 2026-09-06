"""Exact data for the all-orders neutral-period descent argument.

Run: sage -c "load('checks/check_x9_neutral_period_descent.sage')"
The proof additionally uses conifold Hodge-theoretic descent, the existing
residue-preserving mutation, and the large-volume Frobenius basis. These
geometric inputs are not replaced by finite coefficient comparisons.
"""
import contextlib
import io
with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_neutral_quantum_periods.sage')

# Charge ordering: (p0,p1,...,p6;q0,q1,...,q6), in a symplectic
# completion of the primitive electric D0/nef-dual basis. The magnetic
# slots are abstract duals, not assertions about unshifted integral Mukai
# tuples. The four charges are C1,C2,e1,e2, with zero D0 part.
I7=identity_matrix(ZZ,7)
Z7=zero_matrix(ZZ,7)
J14=block_matrix([[Z7,I7],[-I7,Z7]])
local_rows=matrix(ZZ,[[0,1,1,0],[0,0,0,0],[0,0,0,0],
                     [0,0,0,0],[1,0,1,0],[0,0,1,1]])
actual_charges=matrix(ZZ,[C1,C2,curve(4,9),curve(3,9)]).T
assert local_rows==nef.T*actual_charges
vanishing=matrix(ZZ,14,4,lambda i,j:local_rows[i-8,j] if i>=8 else 0)
assert vanishing.rank()==3
assert vanishing.T*J14*vanishing==0
assert [x for x in vanishing.elementary_divisors() if x]==[1,1,1]
assert vanishing*vector(ZZ,[-1,-1,1,-1])==0

identity14=identity_matrix(ZZ,14)
vanishing_basis=identity14.matrix_from_columns([8,12,13])
assert vanishing_basis.column_module()==vanishing.column_module()
surviving=identity14.matrix_from_columns([0,2,3,4,7,9,10,11])
annihilator=(vanishing.T*J14).right_kernel()
assert surviving.augment(vanishing_basis).column_module()==annihilator
assert surviving.T*J14*vanishing_basis==0
J8=block_matrix([[zero_matrix(ZZ,4),identity_matrix(ZZ,4)],
                 [-identity_matrix(ZZ,4),zero_matrix(ZZ,4)]])
assert surviving.T*J14*surviving==J8
assert J8.det()==1
N=sum((d.column()*(d.row()*J14) for d in vanishing.columns()),zero_matrix(ZZ,14))
assert N^2==0 and N.rank()==3
assert N.right_kernel()==annihilator
assert N.column_module().saturation()==vanishing_basis.column_module()
assert [x for x in N.elementary_divisors() if x]==[1,1,4]
print('V has rank 3 and is primitive; V-perp/V has rank 8 and unimodular symplectic form.')
print('The combined monodromy image has index 4 in V; the proof uses V, not the unsaturated image.')

# The logarithmic polynomials of a full Frobenius basis are independent.
# Only the leading cubic term of the eighth period is needed; its lower
# polynomial terms may be changed without affecting this argument.
log_ring=PolynomialRing(QQ,names=('l1','l2','l3'))
ls=vector(log_ring,log_ring.gens())
quadratics=[(ls*k*ls)/2 for k in smooth_kappa]
cubic=sum(ls[i]*quadratics[i] for i in range(3))/3
polynomials=[log_ring(1)]+list(ls)+quadratics+[cubic]
monomials=sorted(set(mon for poly in polynomials for mon in poly.monomials()),key=str)
coefficient_matrix=matrix(QQ,[[poly.monomial_coefficient(mon) for mon in monomials]
                             for poly in polynomials])
assert coefficient_matrix.rank()==8
quadratic_minor=matrix(QQ,[[k[0,0],k[0,1],k[0,2]] for k in smooth_kappa])
assert quadratic_minor==matrix(QQ,[[0,2,3],[2,6,11],[3,11,17]])
assert quadratic_minor.det()==10
assert cubic.monomial_coefficient(ls[2]^3)==QQ(83)/6
assert sum(smooth_kappa).det()==700
print('Eight logarithmic polynomials have rank 8; quadratic minor determinant 10; cubic l3^3 coefficient 83/6.')

# Exact endpoint constants, including every local degree and every zeta(2)
# tail. The six-variable rational Frobenius normalization is retained.
assert endpoint_data[(0,0,0)][0]==1
assert endpoint_data[(0,0,0)][1]==zero_vector(QQ,6)
assert vector(endpoint_data[(0,0,0)][2][1:4])==zero_vector(zeta_ring,3)
assert all(pair_charges[pair][i]==0 for pair in active_pairs for i in [1,2,3])
assert all(g[(0,0,0)]==0 for g in smooth_data['mirror_log'])
assert all(g[(0,0,0)]==0 for g in smooth_data['magnetic'])

# Restricting the full six-variable quadratic logarithmic terms leaves
# exactly the smooth-side cubic contractions, not an arbitrary frame.
old_logs=vector(log_ring,[0,ls[0],ls[1],ls[2],0,0])
for i in range(3):
    assert (old_logs*nef_kappa[i+1]*old_logs)/2==quadratics[i]
print('All old and smooth electric/neutral-magnetic constant corrections vanish; leading log polynomials match exactly.')

# Master-equation identity in abstract analytic jets: after matching
# scalar, first jets and contracted second jets, all three instanton
# gradients agree, without equating uncontracted Hessian entries.
master_ring=PolynomialRing(QQ,names=('W','g1','g2','g3','d1','d2','d3'))
W,*rest=master_ring.gens()
gs=vector(master_ring,rest[:3])
ds=vector(master_ring,rest[3:])
assert all(W*ds[i]-(gs*smooth_kappa[i]*gs)/2==
           W*ds[i]-sum(smooth_kappa[i][j,k]*gs[j]*gs[k] for j in range(3) for k in range(3))/2
           for i in range(3))
print('ALL NEUTRAL-DESCENT DATA CHECKS PASSED; the text proves all-orders equality of seven periods, not a full integral D6 marking.')

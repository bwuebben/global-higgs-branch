"""Full six-vector current exchange on three X_9 decoupling rays.

Run from the project root: sage -c "load('checks/check_x9_current_exchange.sage')"
The Hodge metric is in fixed eleven-dimensional/local-scale units. This is a
two-derivative Coulomb-chamber computation, not the unknown finite-Planck-mass
quaternionic metric of the interacting transition sector.
"""
import contextlib
import io
with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_prepotential.sage')

pol = PolynomialRing(QQ,'scl')
scl = pol.gen()
ff = pol.fraction_field()
charge = matrix(QQ,[C1,C2]).transpose()
embedding = matrix(QQ, [[1,0,1],[0,1,1]])
cooperative = vector(QQ,[-1,-1,1])
assert charge.rank() == 2 and embedding*cooperative == 0
assert charge.transpose()*J0 == vector(QQ,[1,1])

def at_infinity(value):
    if value == 0:
        return QQ(0)
    numerator, denominator = value.numerator(), value.denominator()
    assert numerator.degree() <= denominator.degree(), value
    if numerator.degree() < denominator.degree():
        return QQ(0)
    return numerator.leading_coefficient()/denominator.leading_coefficient()

def positive_on_nonnegative_axis(value):
    numerator, denominator = value.numerator(), value.denominator()
    if denominator(0) < 0:
        numerator, denominator = -numerator, -denominator
    return (numerator(0)>0 and denominator(0)>0 and
            all(c>=0 for c in numerator.list()+denominator.list()))

limits = {
    'D0': zero_matrix(QQ,2),
    'D1': matrix(QQ,[[1,-1],[-1,1]])/14,
    'interior': zero_matrix(QQ,2),
}
coefficients = {
    'D0': matrix(QQ,[[31,-9],[-9,31]])/66,
    'D1': matrix(QQ,[[96,-5],[-5,110]])/196,
    'interior': matrix(QQ,[[209,-91],[-91,209]])/590,
}
root_coefficients = {'D0': QQ(2)/3, 'D1': QQ(1), 'interior': QQ(2)/5}
print('Basis:', BASIS, '; J0 =',J0)
for tag,ww in [('D0',CLS['D0']),('D1',CLS['D1']),('interior',CLS['D0']+CLS['D1'])]:
    tt = vector(ff,J0)+scl*vector(ff,ww)
    ll = matrix(ff,6,6,lambda i,j:sum(kappa(BASIS[i],BASIS[j],BASIS[k])*tt[k] for k in range(6)))
    vv=ll*tt
    volume=tt*vv/6
    gg=-ll+vv.column()*vv.row()/(4*volume)
    gi=-ll.inverse()+tt.column()*tt.row()/(2*volume)
    assert gg*gi==identity_matrix(ff,6)
    assert gg*tt == vv/2 and tt*gg*tt == 3*volume
    assert matrix(QQ,gg.apply_map(lambda x:x(0))).is_positive_definite()
    hh=charge.transpose()*gi*charge
    assert positive_on_nonnegative_axis(hh[0,0])
    assert positive_on_nonnegative_axis(hh.det())
    hlim=hh.apply_map(at_infinity)
    correction=scl*(hh-hlim)
    leading=correction.apply_map(at_infinity)
    assert hlim == limits[tag] and leading == coefficients[tag]
    assert at_infinity(scl*sum(hh.list())) == root_coefficients[tag]
    mm = embedding.transpose()*hh*embedding
    assert mm.rank()==2 and mm*cooperative==0
    assert mm.right_kernel().dimension()==1
    mlim=embedding.transpose()*hlim*embedding
    assert mlim.rank() == (1 if tag=='D1' else 0)

    # Remove the unit graviphoton direction before contracting the charges.
    gi_vm = gi-tt.column()*tt.row()/(3*volume)
    hh_vm = charge.transpose()*gi_vm*charge
    assert hh-hh_vm == matrix(ff,2,2,[1,1,1,1])/(3*volume)
    assert (scl*(hh-hh_vm)).apply_map(at_infinity)==zero_matrix(QQ,2)

    # A non-orthogonal integral divisor basis change leaves current exchange fixed.
    change=identity_matrix(QQ,6)
    change[0,4]=2
    change[2,1]=-1
    transformed_g=change.transpose()*gg*change
    transformed_gi=change.inverse()*gi*change.transpose().inverse()
    transformed_charge=change.transpose()*charge
    assert transformed_g*transformed_gi==identity_matrix(ff,6)
    assert transformed_charge.transpose()*transformed_gi*transformed_charge==hh

    print(tag, ': V =', volume)
    print('  H(infinity) =', hlim.list())
    print('  coefficient of 1/Lambda =', leading.list())
    print('  root exchange: Lambda * (1,1) H (1,1)^T ->', root_coefficients[tag])
    print('  kernel dimensions: finite 1, limiting', mlim.right_kernel().dimension())
print('ALL CHECKS PASSED: full inverse, positivity for Lambda >= 0, asymptotics,')
print('graviphoton subtraction, basis covariance, and noncommuting kernel limits.')

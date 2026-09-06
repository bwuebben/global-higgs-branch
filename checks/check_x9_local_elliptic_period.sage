"""Continue one local del Pezzo A-period on the mass slice u=1, v=w.

This is the z3=z4=0 boundary problem, not the finite-volume compact
magnetic alignment. Numerical integration is a diagnostic, not an
interval-arithmetic error bound. Run from the project root with
sage -c "load('checks/check_x9_local_elliptic_period.sage')".
"""
import contextlib
import io
with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_mirror_continuation.sage')

PW = PolynomialRing(QQ,'w')
w = PW.gen()
KW = PW.fraction_field()
PX = PolynomialRing(KW,'x')
x = PX.gen()
quartic_local = (w*x^2+x+w)^2-4*w^2*x^2*(1+x)
aa,bb,cc,dd,ee = [quartic_local[i] for i in [4,3,2,1,0]]
II = PW(12*aa*ee-3*bb*dd+cc^2)
JJ = PW(72*aa*cc*ee+9*bb*cc*dd-27*aa*dd^2-27*bb^2*ee-2*cc^3)
assert II == 1-16*w^2+24*w^3+16*w^4
assert PW(quartic_local.discriminant()) == -256*w^7*(w-1)^2*(43*w^3-18*w^2-w+1)
assert PW(quartic_local.discriminant()) == (4*II^3-JJ^2)/27
print('Local elliptic I:',II)
print('Local elliptic J:',JJ)
singular_cubic = 43*w^3-18*w^2-w+1
assert not [r for r,multiplicity in singular_cubic.roots(AA) if 0<r<=1]
assert II(1)==25 and JJ(1)==-250
assert PW(quartic_local.discriminant()).valuation(w-1)==2
assert not [r for r,multiplicity in II.roots(AA) if 0<r<1]
print('The local fiber at w=1 is I2; no discriminant zero lies in 0<w<1.')

# The scalar A-period derivative omega=1-w*h'(w) at u=1, v=w.
# All powers of the mass coordinate u are summed at each w degree.
series_degree = 24
PS = PowerSeriesRing(QQ,'w',default_prec=series_degree+1)
ws = PS.gen()
h = PS(0).add_bigoh(series_degree+1)
for ell in range(1,series_degree//2+1):
    for k in range(ell,series_degree-ell+1):
        for m in range(k-ell,ell+1):
            multiplicities = (ell-m,ell+m-k,k-m,m,k-ell)
            assert min(multiplicities)>=0 and sum(multiplicities)==ell+k
            assert multiplicities[0]-multiplicities[2]+multiplicities[4]==0
            assert multiplicities[1]-multiplicities[3]+multiplicities[4]==0
            coefficient = (-1)^(ell+k-1)*factorial(ell+k-1)/(
                factorial(m)*factorial(ell+m-k)*factorial(ell-m)*
                factorial(k-m)*factorial(k-ell))
            assert coefficient == (-1)^(ell+k-1)*factorial(ell+k-1)/prod(
                factorial(value) for value in multiplicities)
            jet = gamma_coefficients(glsm,(m,ell,0,0,m,k))
            assert jet[1] == coefficient*glsm.column(9)
            h += coefficient*ws^(ell+k)
omega_series = 1-ws*h.derivative()
I_series,J_series = PS(II),PS(JJ)
argument = 1-J_series^2/(4*I_series^3)
hypergeom = sum(rising_factorial(QQ(1)/12,k)*rising_factorial(QQ(5)/12,k)/factorial(k)^2*
               argument^k for k in range(series_degree+1))
assert omega_series == I_series^(-QQ(1)/4)*hypergeom
print('Frobenius/hypergeometric A-period check through degree',series_degree)

# Derive the Gauss-Manin connection by exact reduction of rational forms.
g2 = KW(4*II/3)
g3 = KW(4*JJ/27)
weierstrass = 4*x^3-g2*x-g3
forms = [weierstrass,x*weierstrass]
for k in range(3):
    exact_numerator = (x^k).derivative()*weierstrass-x^k*weierstrass.derivative()/2
    forms.append(exact_numerator)
reduction = matrix(KW,[[form[i] for form in forms] for i in range(5)])
connection = []
for k in range(2):
    numerator = x^k*(g2.derivative()*x+g3.derivative())/2
    solution = reduction.solve_right(vector(KW,[numerator[i] for i in range(5)]))
    assert sum(solution[i]*forms[i] for i in range(5)) == numerator
    connection.append(list(solution[:2]))
connection = matrix(KW,connection)
weierstrass_disc = g2^3-27*g3^2
rotation = 2*g2*g3.derivative()-3*g3*g2.derivative()
assert connection == matrix(KW,[[-weierstrass_disc.derivative()/12,3*rotation/2],
                               [-g2*rotation/8,weierstrass_disc.derivative()/12]])/weierstrass_disc
print('Gauss-Manin connection:',connection)

# The exact series supplies initial conditions on the large-volume branch.
from scipy.integrate import solve_ivp
import numpy as np
w_start = QQ(1)/100
omega0 = omega_series.polynomial()(w_start)
derivative0 = omega_series.derivative().polynomial()(w_start)
eta0 = (derivative0-connection[0,0](w_start)*omega0)/connection[0,1](w_start)
h0 = float(h.polynomial()(w_start))
evaluators = [[fast_callable(SR(entry),vars=[SR(w)],domain=RDF) for entry in row]
              for row in connection.rows()]
def rhs(value,state):
    matrix_value = np.array([[entry(value) for entry in row] for row in evaluators])
    return list(matrix_value@state[:2])+[(1-state[0])/value]
def q_one(value,state):
    return np.log(value)-state[2]
q_one.terminal = False
q_one.direction = 1
results = []
for tolerance in [2e-11,2e-13]:
    solution = solve_ivp(rhs,[float(w_start),0.9999],[float(omega0),float(eta0),h0],
                         method='DOP853',rtol=tolerance,atol=tolerance/100,
                         dense_output=True,events=q_one,max_step=0.005)
    assert solution.success
    print('Tolerance',tolerance,'q=1 events:',solution.t_events[0])
    for value in [0.1,0.25,0.5,0.75,0.9,0.95,0.9999]:
        state = solution.sol(value)
        print('w, omega, log(q):',value,state[0],np.log(float(value))-state[2])
    results.append(solution)
assert max(abs(results[0].sol(value)-results[1].sol(value)).max()
           for value in [0.1,0.25,0.5,0.75,0.9,0.95,0.9999])<1e-8

# An independent all-orders evaluation uses the real hypergeometric branch.
# I>0 and Delta<0 on (0,1), so its argument stays negative: no branch
# matching is hidden in this comparison. The initial-series truncation
# and ODE integration are not used by this quadrature.
import mpmath as mp
with mp.workdps(50):
    icoeff = [mp.mpf(str(c)) for c in reversed(II.list())]
    jcoeff = [mp.mpf(str(c)) for c in reversed(JJ.list())]
    def exact_omega(value):
        ival,jval = mp.polyval(icoeff,value),mp.polyval(jcoeff,value)
        return ival^(-mp.mpf(1)/4)*mp.hyp2f1(mp.mpf(1)/12,mp.mpf(5)/12,1,
                                          1-jval*jval/(4*ival^3))
    def h_derivative(value):
        return (1-exact_omega(value))/value if value else mp.mpf(0)
    for value in ['0.25','0.5','0.95','0.9999','1']:
        endpoint = mp.mpf(value)
        hvalue = mp.quad(h_derivative,[0,endpoint/4,endpoint/2,3*endpoint/4,endpoint])
        logqvalue = mp.log(endpoint)-hvalue
        print('Independent 50-digit quadrature w, omega, log(q):',value,
              mp.nstr(exact_omega(endpoint),25),mp.nstr(logqvalue,25))
        if endpoint<1:
            ode = results[1].sol(float(endpoint))
            assert abs(float(exact_omega(endpoint))-ode[0])<1e-10
            assert abs(float(logqvalue)-(np.log(float(endpoint))-ode[2]))<1e-10
        else:
            assert abs(logqvalue)<mp.mpf('1e-45')
            assert abs(exact_omega(endpoint)-1/mp.sqrt(5))<mp.mpf('1e-45')
    print('Quadrature checks log(q(1))=0; the separate Chu--Vandermonde proof establishes exactness.')
print('LOCAL A-PERIOD CONTINUATION PASSED; compact magnetic periods not continued.')

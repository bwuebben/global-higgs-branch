"""Exact fixed-spectator resummation of magnetic periods on the four-node locus.

Run: sage -c "load('checks/check_x9_regular_magnetic_periods.sage')"
This sums every local degree, not a six-variable total-degree truncation.
The convergence and physical-stability scope are discussed in the manuscript.
"""
import contextlib
import io
import os
os.environ.setdefault('MIRROR_DEGREE','2')
with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_mirror_periods.sage')

univariate=PolynomialRing(QQ,'u')
u=univariate.gen()
rationals=univariate.fraction_field()
zeta_ring=PolynomialRing(QQ,'zeta2')
zeta2=zeta_ring.gen()

def arguments(i,L,A,B,j,k):
    return [B,A,j-i-B,L+j-k,L-i,i,k-j,i+2*A-j,k-L-2*A,A+B-L-k]

def factorial_shift(offset):
    """(u+offset)!/u!, a rational function at integer offset."""
    if offset>=0:
        return rationals(prod(u+j for j in range(1,offset+1)))
    return 1/rationals(prod(u+j for j in range(offset+1,1)))

def tail_prefactor(args,L,A,B,pair):
    """The scalar prefactor of a two-zero Gamma Hessian along its ray."""
    result=rationals(factorial(2*A+B))
    moving_numerator=moving_denominator=0
    sign_exponent=0
    for index,arg in enumerate(args):
        arg=univariate(arg)
        slope=arg[1]
        offset=ZZ(arg[0])
        if index in pair:
            assert slope==-1
            moving_numerator+=1
            result*=factorial_shift(-offset-1)
            sign_exponent+=-offset-1
        elif slope:
            assert slope==1
            moving_denominator+=1
            result/=factorial_shift(offset)
        else:
            assert offset>=0
            result/=factorial(offset)
    assert moving_numerator==moving_denominator==2
    return (-1)^sign_exponent*result

def lines_for_pair(L,A,B,pair):
    """(minimum ray parameter, symbolic arguments) for all transverse indices."""
    if pair==(2,4):
        for k in range(L+2*A,A+B-L+1):
            for j in range(k-L,k+1):
                yield max(L+1,j-B+1,j-2*A,0),arguments(u,L,A,B,j,k)
    elif pair==(6,7):
        for i in range(L+1):
            for k in range(L+2*A,A+B-L+1):
                yield max(k+1,i+2*A+1,i+B),arguments(i,L,A,B,u,k)
    elif pair==(3,9):
        for i in range(L+1):
            for j in range(i+B,i+2*A+1):
                yield max(L+j+1,A+B-L+1,L+2*A),arguments(i,L,A,B,j,u)
    elif pair==(4,9):
        for p in range(B,2*A+1):
            for q in range(L+1):
                yield max(L+1,A+B-L-p-q+1,L+2*A-p-q,0),arguments(u,L,A,B,u+p,u+p+q)
    else:
        raise ValueError(pair)

def partial_fraction_terms(function):
    polynomial,parts=rationals(function).partial_fraction_decomposition()
    terms=[]
    for part in parts:
        factors=list(part.denominator().factor())
        assert len(factors)==1
        factor,power=factors[0]
        assert factor.degree()==1 and power in [1,2]
        pole=-factor[0]/factor[1]
        assert pole in ZZ
        numerator=part.numerator()(u+pole)/factor[1]^power
        for exponent,residue in enumerate(numerator.list()):
            if residue:
                terms.append((residue,pole,power-exponent))
    return polynomial,terms

def rational_tail_sum(function,start):
    """Exact sum from start to infinity as rational + rational*zeta(2).

    Requires O(u^-2); simple-pole residues cancel before summing. Every
    pole is integral of order at most two and strictly before start.
    """
    function=rationals(function)
    if not function:
        return zeta_ring(0)
    assert function.denominator().degree()-function.numerator().degree()>=2, ('nonconvergent tail',function)
    polynomial,parts=partial_fraction_terms(function)
    assert not polynomial
    total=zeta_ring(0)
    simple_residue=QQ(0)
    for residue,pole,power in parts:
        assert pole<start
        if power==1:
            simple_residue+=residue
            total-=residue*harmonic(start-pole-1)
        else:
            total+=residue*(zeta2-harmonic(start-pole-1,2))
    assert simple_residue==0
    return total

def rational_abel_sum(function,start,ratio,mp):
    """Independent Abel sum before the local variables are specialized to 1."""
    polynomial,parts=partial_fraction_terms(function)
    to_mp=lambda value: mp.mpf(str(QQ(value).numerator()))/mp.mpf(str(QQ(value).denominator()))
    answer=mp.mpf(0)
    for power,coefficient in enumerate(polynomial.list()):
        answer+=to_mp(coefficient)*(mp.polylog(-power,ratio)-sum(
            ratio^k*mp.mpf(k)^power for k in range(1,start)))
    assert start>=1
    for residue,pole,power in parts:
        assert pole<start and power in [1,2]
        answer+=to_mp(residue)*ratio^int(pole)*(mp.polylog(power,ratio)-sum(
            ratio^k/mp.mpf(k)^power for k in range(1,start-int(pole))))
    return answer

def falling_integer(value,length):
    return prod(ZZ(value)-j for j in range(length))

def closed_tail(L,A,B,pair):
    """Finite-difference evaluation of all transverse indices, to all degrees."""
    common=QQ(factorial(2*A+B))/(factorial(A)*factorial(B))
    if pair in [(2,4),(6,7)]:
        N=B-A-2*L
        if N<0:
            return rationals(0)
        if pair==(2,4):
            coefficient=common*(-1)^(B+L)*falling_integer(B-2*A-1,L+N)/(factorial(L)*factorial(N))
            return coefficient*factorial_shift(-L-1)*factorial_shift(L-A-1)
        coefficient=common*(-1)^L*falling_integer(-L-1,N)*falling_integer(B-2*A-1,L)/(factorial(L)*factorial(N))
        return coefficient*factorial_shift(-A-B+L-1)*factorial_shift(-L-2*A-1)/(
            factorial_shift(-2*A)*factorial_shift(-B))
    N=2*A-B
    if N<0:
        return rationals(0)
    if pair==(4,9):
        coefficient=common*(-1)^(A+L+N)*falling_integer(A+2*L-B-1,L+N)/(factorial(L)*factorial(N))
        return coefficient*factorial_shift(-L-1)*factorial_shift(L-A-1)
    assert pair==(3,9)
    coefficient=common*(-1)^A*falling_integer(-L-1,L+N)/(factorial(L)*factorial(N))
    return coefficient*factorial_shift(-2*A-2*L-1)*factorial_shift(-A-B+L-1)/(
        factorial_shift(-B)*factorial_shift(-L-2*A))

active_pairs=[(2,4),(6,7),(4,9),(3,9)]  # C2,C1,e1,e2
pair_charges={pair:vector(QQ,[trip(nef.column(a),CLS[names[pair[0]]],CLS[names[pair[1]]])
                                      for a in range(6)]) for pair in combinations(range(10),2)}
assert [pair_charges[pair] for pair in active_pairs]==[
    vector(QQ,[1,0,0,0,0,0]),vector(QQ,[0,0,0,0,1,0]),
    vector(QQ,[1,0,0,0,1,1]),vector(QQ,[0,0,0,0,0,1])]
for pair in [(4,6),(7,9)]:
    assert pair_charges[pair]==0

# Bound every other two-zero sector in all six nonnegative degrees.
# A rational polyhedral cone is generated by its rays, so these inequalities
# are all-degree checks, not a finite lattice-point scan.
coordinate_rows=identity_matrix(QQ,6).rows()
spectator_row=coordinate_rows[1]+coordinate_rows[2]+coordinate_rows[3]
for pair in combinations(range(10),2):
    if pair in active_pairs or not pair_charges[pair]:
        continue
    inequalities=[[0]+list(row) for row in coordinate_rows]
    inequalities += [[0]+list(glsm.column(index)) for index in range(10) if index not in pair]
    cone=Polyhedron(ieqs=inequalities)
    assert not cone.lines()
    for ray in cone.rays():
        ray=vector(ray)
        assert all(coordinate_rows[index]*ray<=4*(spectator_row*ray) for index in [0,4,5])
print('All remaining nonzero second-derivative sectors have bounded local degrees.')

def endpoint_coefficient(L,A,B,verbose=False,abel_radius=None):
    """W0, C_a and contracted rational Hessians, including exact local tails."""
    L,A,B=map(ZZ,[L,A,B])
    cutoff=8*(L+A+B)+2
    scalar=QQ(0)
    gradient=zero_vector(QQ,6)
    contracted=zero_vector(zeta_ring,6)
    if abel_radius is not None:
        import mpmath as mp
        rho=mp.mpf(abel_radius)
        abel=[mp.mpf(0) for _ in range(6)]
        to_mp=lambda value: mp.mpf(str(QQ(value).numerator()))/mp.mpf(str(QQ(value).denominator()))
    # All nonactive sectors are bounded. Check a generous finite box;
    # active sectors are covered separately by their complete ray lines.
    box=4*(L+A+B)+1
    for i in range(box+1):
      for j in range(box+1):
        for k in range(box+1):
            args=arguments(i,L,A,B,j,k)
            negative=tuple(index for index,arg in enumerate(args) if arg<0)
            if len(negative)>2 or negative in active_pairs:
                continue
            if len(negative)==2 and not pair_charges[negative]:
                continue
            c,g,h=gamma_coefficients(glsm,(i,L,A,B,j,k))
            scalar+=c
            gradient+=g
            coefficient=vector(QQ,[sum(nef_kappa[a][r,s]*h[r,s] for r in range(6) for s in range(6))/2
                                        for a in range(6)])
            contracted+=coefficient
            if abel_radius is not None:
                abel=[abel[a]+to_mp(coefficient[a])*rho^(i+j+k) for a in range(6)]
    for pair in active_pairs:
        finite=QQ(0)
        tail=rationals(0)
        for minimum,args in lines_for_pair(L,A,B,pair):
            assert 0<=minimum<cutoff
            function=tail_prefactor(args,L,A,B,pair)
            for h in range(minimum,cutoff):
                finite+=function(h)
            tail+=function
            if abel_radius is not None:
                # args[5]=i, args[2]=j-i-B, args[6]=k-j.
                local_degree=3*univariate(args[5])+2*univariate(args[2])+2*B+univariate(args[6])
                abel_scalar=rho^int(local_degree[0])*rational_abel_sum(
                    function,minimum,rho^int(local_degree[1]),mp)
                abel=[abel[a]+to_mp(pair_charges[pair][a])*abel_scalar for a in range(6)]
        assert tail==closed_tail(L,A,B,pair), ('finite-difference formula',L,A,B,pair)
        summed=zeta_ring(finite)+rational_tail_sum(tail,cutoff)
        contracted+=summed*pair_charges[pair]
        if verbose:
            print('  pair',pair,'tail',tail,'sum',summed)
    assert all(gradient[a]==0 for a in [0,4,5])
    if abel_radius is not None:
        exact=[to_mp(entry[0])+to_mp(entry[1])*mp.zeta(2) for entry in contracted]
        error=max(abs(left-right) for left,right in zip(abel,exact))
        assert error<mp.mpf('1e-6'), ('Abel endpoint discrepancy',L,A,B,error)
    return scalar,gradient,contracted

endpoint_data={}
for d in range(int(os.environ.get('MAGNETIC_BULK_DEGREE','2'))+1):
  for spectator_L,spectator_A,spectator_B in IntegerVectors(d,3):
    result=endpoint_coefficient(spectator_L,spectator_A,spectator_B,verbose=(d==0))
    endpoint_data[(spectator_L,spectator_A,spectator_B)]=result
    print('Spectator degree',(spectator_L,spectator_A,spectator_B),'W0',result[0],'C',result[1],'D',result[2])
assert endpoint_data[(0,0,0)][2]==vector(zeta_ring,[2*zeta2,0,0,0,2*zeta2,2*zeta2])
assert endpoint_data[(1,0,0)][2]==endpoint_data[(0,1,0)][2]==vector(zeta_ring,[1,0,0,0,1,2])
assert endpoint_data[(0,0,1)][2]==vector(zeta_ring,[0,0,0,10,10,10])

# Independent radial Abel continuation retains all transverse local weights
# until the final limit; it checks that early specialization has lost no term.
import mpmath as mp
with mp.workdps(100):
    for spectator in [(0,0,0),(1,0,0),(0,1,1),(1,0,2)]:
        endpoint_coefficient(*spectator,abel_radius=mp.mpf(1)-mp.mpf('1e-12'))
print('Independent radial Abel sums agree with the exact endpoint values in four test sectors.')

# Flat-coordinate magnetic corrections and the exact real-branch phase test.
spec_ring=PolynomialRing(QQ,names=('y2','y3','y4'))
y2,y3,y4=spec_ring.gens()
Jspec=y2*nef.column(1)+y3*nef.column(2)+y4*nef.column(3)
PA=CLS['D4']
PV=CLS['D7']+CLS['D8']
PT=CLS['D6']
qG=L(PA)*(CLS['D2']+PA)-L(PA)*PA/2
qH=-L(PV)*PV/2
assert Jspec*qG==-y3-3*y4/2 and Jspec*qH==y3+3*y4/2
assert trip(Jspec,Jspec,PA)/2==y3^2+3*y3*y4+3*y4^2/2
assert trip(Jspec,Jspec,PV)/2==2*y2*y3+y3^2+3*y2*y4+5*y3*y4+11*y4^2/2
volumeT=2*y2*y3+2*y3^2+3*y2*y4+8*y3*y4+7*y4^2
assert trip(Jspec,Jspec,PT)/2==volumeT
assert qG+qH==C2-C1/2
assert L(PV)*PT==C1
for p,constant,linear in [(PA,-2,-1),(PV,2,2),(PT,0,1)]:
    weights=nef.inverse()*p
    assert weights*endpoint_data[(0,0,0)][2]==constant*zeta2
    assert weights*endpoint_data[(1,0,0)][2]==linear
    assert weights*endpoint_data[(0,1,0)][2]==linear
    assert weights*endpoint_data[(0,0,1)][2]==0
phase_ring=PolynomialRing(QQ,names=('realG','realH','imagG'))
realG,realH,imagG=phase_ring.gens()
assert imagG*realH-realG*(-imagG)==imagG*(realG+realH)
print('Psi_A=-1/12-(v+a)/(4*pi^2)+O(bulk^2); Psi_V=1/12+2*(v+a)/(4*pi^2)+O(bulk^2).')
print('Exact endpoint alignment: Im(ZG*conj(ZH))=-(y3+3*y4/2)*(Vol_T+Psi_T+1/12).')
print('ALL REGULAR MAGNETIC CHECKS PASSED; constituent survival and other decay channels remain open.')

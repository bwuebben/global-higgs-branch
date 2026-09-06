"""Exact weighted finite-difference checks for Appendix C of paper 2.

Run: sage -c "load('checks/check_x9_weighted_abel.sage')"
Independent of the period implementation. The finite tests check the
weighted identity used in the uniform Abel estimate, not convergence
by finite truncation. The all-degrees bound is proved in the appendix.
"""
R = PolynomialRing(QQ, 'rho')
rho = R.gen()
falling = lambda x,k: prod(ZZ(x)-j for j in range(k))
G = lambda x: factorial(ZZ(x)-1)

def raw(L,m,n,h,kind):
    N = n-m-2*L if kind in ['C1','C2'] else 2*m-n
    if N < 0:
        return R(0)
    ans=R(0)
    for s in range(L+1):
      for t in range(N+1):
        if kind=='C2': i,j,k=h,2*m+t+s,L+2*m+t
        if kind=='C1': i,j,k=s,h,L+2*m+t
        if kind=='e1': i,j,k=h,h+n+t,h+n+t+s
        if kind=='e2': i,j,k=s,s+n+t,h
        args=[n,m,j-i-n,L+j-k,L-i,i,k-j,i+2*m-j,k-L-2*m,m+n-L-k]
        negative=[x for x in args if x<0]
        assert len(negative)==2
        pref=QQ(factorial(2*m+n))*prod((-1)^(-x-1)*factorial(-x-1) for x in negative)
        pref/=prod(factorial(x) for x in args if x>=0)
        ans+=pref*rho^(i+j+k)
    return ans

def expanded(L,m,n,h,kind):
    N = n-m-2*L if kind in ['C1','C2'] else 2*m-n
    if N < 0:
        return R(0)
    K=QQ(factorial(2*m+n))/(factorial(m)*factorial(n)*factorial(L)*factorial(N))
    ans=R(0)
    for x in range(L+1):
      for y in range(N+1):
        lost=x+y
        count=L+N-lost
        choose=binomial(L,x)*binomial(N,y)
        if kind=='C2':
            weight=rho^(h+L+4*m+L-x+2*(N-y))*(1-rho)^x*(1-rho^2)^y
            amplitude=(-1)^(n+L)*falling(n-2*m-1,count)
            ratio=QQ(G(h-L)*G(h+L-m+lost))/G(h+1)^2
        if kind=='C1':
            weight=rho^(h+L+2*m+L-x+N-y)*(1-rho)^(x+y)
            amplitude=(-1)^L*falling(-L-1,N-y)*falling(n-2*m-1,L-x)
            ratio=QQ(G(h-m-n+L+y)*G(h-L-2*m+x))/(G(h-2*m+1)*G(h-n+1))
        if kind=='e1':
            weight=rho^(3*h+2*n+L-x+2*(N-y))*(1-rho)^x*(1-rho^2)^y
            amplitude=(-1)^(m+count)*falling(m+2*L-n-1,count)
            ratio=QQ(G(h-L)*G(h+L-m))/(G(h+1)*G(h-lost+1))
        if kind=='e2':
            weight=rho^(h+n+2*(L-x)+N-y)*(1-rho^2)^x*(1-rho)^y
            amplitude=(-1)^m*falling(-L-1,count)
            ratio=QQ(G(h-2*m-2*L+lost)*G(h-m-n+L))/(G(h-n+1)*G(h-L-2*m+1))
        ans+=K*choose*weight*amplitude*ratio
    return ans

checks=0
for d in range(7):
  for L,m,n in IntegerVectors(d,3):
    for kind in ['C1','C2','e1','e2']:
      for h in [8*d+2,8*d+9]:
        assert raw(L,m,n,h,kind)==expanded(L,m,n,h,kind),(L,m,n,h,kind)
        checks+=1
print('Weighted identities checked exactly as polynomials in rho:',checks)
print('Expansion uses (1-uE)^r = sum binom(r,s)*(1-u)^s*u^(r-s)*(1-E)^(r-s).')
print('At lost order s the h exponent is -m-2+s, as printed.')
print('For s>m+2, p=s-m-2 <= 3d and p!/h0^p <= 1 at h0=8d+2.')
print('This controls the otherwise growing Abel moments uniformly in degree.')

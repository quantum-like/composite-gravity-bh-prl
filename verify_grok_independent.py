# Grok independent check (does not reuse Claude's script). r0=1, c=1, Cartesian background.
# Route: A := Levi-Civita spin connection of tetrad (45)+(51) (symmetric gauge), matrix form.
# F_mn = d_m w_n - d_n w_m + [w_m,w_n]; vacuum YM: sum_m d_m F^{mn} + [w_m, F^{mn}] = 0, spacetime idx raised with eta.
# High-precision nested numerical differentiation (mpmath). Part 2: metric invariants with sympy (spherical).
import mpmath as mp, sympy as sp
mp.mp.dps=50
eta=mp.diag([-1,1,1,1])
def E(p):
    x=p[1:]; r=mp.sqrt(sum(v**2 for v in x)); a=1+0.5/r; b=1-1/r; q=1+1/r
    M=mp.zeros(4,4); M[0,0]=b
    for i in range(3):
        for j in range(3): M[i+1,j+1]=a*(i==j)+(q-a)*x[i]*x[j]/r**2
    return M
def d(f,p,m):  # partial derivative along coord m (static: t-derivative 0)
    if m==0: return f(p)*0
    def g(s):
        pp=list(p); pp[m]=s; return f(pp)
    return mp.diff(g,p[m])
def dmat(f,p,m):
    if m==0: return f(p)*0
    h=mp.mpf('1e-12')
    # 6th-order central difference, elementwise
    def at(s):
        pp=list(p); pp[m]=p[m]+s*h; return f(pp)
    return (-at(-3)+9*at(-2)-45*at(-1)+45*at(1)-9*at(2)+at(3))/(60*h)
def omega(p):
    e=E(p); ei=e**-1; dE=[dmat(E,p,m) for m in range(4)]
    # w_m = e * d_m(e^-1) + e*Gamma_m*e^-1 ; equivalently w_m = -(d_m e) e^-1 + e Gam_m e^-1
    g=e.T*eta*e; gi=g**-1
    dg=[dE[m].T*eta*e+e.T*eta*dE[m] for m in range(4)]
    W=[]
    for m in range(4):
        Gm=mp.zeros(4,4)  # Gm[n,s]=Gamma^n_{m s}
        for n in range(4):
            for s in range(4):
                Gm[n,s]=sum(gi[n,l]*(dg[m][l,s]+dg[s][l,m]-dg[l][m,s]) for l in range(4))/2
        W.append(-dE[m]*ei+e*Gm*ei)
    return W
def Fmn(p):
    W=omega(p)
    dW=[dmat(lambda pp,k=k:omega(pp)[k],p,m) for m in range(4) for k in range(4)]
    dW=[[dW[4*m+k] for k in range(4)] for m in range(4)]
    return [[dW[m][n]-dW[n][m]+W[m]*W[n]-W[n]*W[m] for n in range(4)] for m in range(4)], W
mp.mp.dps=80
def raiseF(F):
    return [[eta[m,m]*eta[n,n]*F[m][n] for n in range(4)] for m in range(4)]
pts=[(1.5,0,0),(2,0,0),(3,0,0),(2,2,1),(0.6,0.3,0.4),(2,3,6)]
for P in pts:
    p=[mp.mpf(0)]+[mp.mpf(v) for v in P]; r=mp.sqrt(sum(v**2 for v in p[1:]))
    F,W=Fmn(p); Fu=raiseF(F)
    pref=1/(r**2*(r+1))
    # compare with (49): boost A_(0l)0 = pref*x_l ; rotation A_(23)? compare W_y[1,2] style
    print(f"\n== point {P}, r={mp.nstr(r,8)}  (49)-prefactor*x={mp.nstr(pref*p[1],12)}")
    print("  w_t[0,1] (boost x, time comp) =",mp.nstr(W[0][0,1],12),"  w_t[1,0] =",mp.nstr(W[0][1,0],12))
    print("  w_y[1,2] =",mp.nstr(W[2][1,2],12)," w_x[1,2] =",mp.nstr(W[1][1,2],12), " pref*x_y-ish:",mp.nstr(pref*p[2],12))
    FF=-sum((F[m][n]*Fu[m][n]).apply(lambda v:v)[i,i] for m in range(4) for n in range(4) for i in range(4))/2
    # tr convention: F.F(Killing-type, labels eta-raised) = sum_{a} F^a F_a = -(1/2) tr(F F) for so(1,3) matrices
    pred=4*(6*r**2+4*r+1)/(r**4*(r+1)**4)
    print("  F.F = -1/2 sum tr(F_mn F^mn) =",mp.nstr(FF,15)," draft formula =",mp.nstr(pred,15))
    res=[]
    for n in range(4):
        def Fcol(pp,n=n):
            Fp,_=Fmn(pp); Fpu=raiseF(Fp); return Fpu
        M=mp.zeros(4,4)
        for m in range(1,4):
            M+=dmat(lambda pp,m=m,n=n: raiseF(Fmn(pp)[0])[m][n],p,m)
        for m in range(4):
            M+=W[m]*Fu[m][n]-Fu[m][n]*W[m]
        res.append(mp.nstr(max(abs(v) for v in M),5))
        # opposite commutator sign, to show discrimination
    Mo=mp.zeros(4,4)
    for m in range(1,4): Mo+=dmat(lambda pp,m=m: raiseF(Fmn(pp)[0])[m][1],p,m)
    for m in range(4): Mo-=W[m]*Fu[m][1]-Fu[m][1]*W[m]
    print("  max|D_m F^{m nu}| nu=0..3:",res,"  [opposite-sign D, nu=1]:",mp.nstr(max(abs(v) for v in Mo),5))

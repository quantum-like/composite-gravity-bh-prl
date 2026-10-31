"""Independent checks for the PRL draft 'Gauge-regular but metric-singular' (Öttinger 2026 static solution).
Units: r0 = 1 (geometrized; r0 = G m / c^2). Exact rational arithmetic at Pythagorean lattice points.
Conventions follow Öttinger, J. Phys. Commun. 10, 105003 (2026): tetrad (45), metric (46), solution (51),
composition rule (4), field strength (18) with structure constants of Sec. III.A / (A2), field equations (19).
Lie-algebra labels are raised/lowered with eta on the index pairs of Table I (as stated below Eq. (18))."""
import sympy as sp, itertools
t,x,y,z,r = sp.symbols('t x y z r', real=True); X=[t,x,y,z]
def D(f,mu):                       # partial derivative treating r=|x| as a function of x,y,z
    return 0 if mu==0 else sp.diff(f,X[mu]) + (X[mu]/r)*sp.diff(f,r)
a=1+1/(2*r); b=1-1/r; q=1+1/r; eta=sp.diag(-1,1,1,1)
bt=sp.zeros(4,4); bt[0,0]=b
bbar=sp.zeros(4,4); bbar[0,0]=1/b
for k in range(1,4):
    for m in range(1,4):
        bt[k,m]=a*(sp.KroneckerDelta(k,m)-X[k]*X[m]/r**2)+q*X[k]*X[m]/r**2
        bbar[m,k]=(1/a)*(sp.KroneckerDelta(k,m)-X[k]*X[m]/r**2)+(1/q)*X[k]*X[m]/r**2
g=bt.T*eta*bt; blow=eta*bt
def A_comp(k,l,rho):               # composition rule (4)
    e=sp.Rational(1,2)*sum(bbar[mu,k]*(D(g[nu,rho],mu)-D(g[mu,rho],nu))*bbar[nu,l] for mu in range(4) for nu in range(4))
    return e+sp.Rational(1,2)*sum(D(blow[k,mu],rho)*bbar[mu,l]-bbar[mu,k]*D(blow[l,mu],rho) for mu in range(4))
P=1/(r**2*(r+1))
def A49(k,l,rho):                  # Eq. (49)
    if k==0 and l>0: return P*X[l]*sp.KroneckerDelta(rho,0)
    if l==0 and k>0: return -P*X[k]*sp.KroneckerDelta(rho,0)
    if k>0 and l>0: return P*(X[l]*sp.KroneckerDelta(k,rho)-X[k]*sp.KroneckerDelta(l,rho))
    return sp.Integer(0)
def ev(e,p,s):
    xv,yv,zv=[sp.Rational(c)*s for c in p]; return sp.nsimplify(e.subs({x:xv,y:yv,z:zv,r:sp.sqrt(xv**2+yv**2+zv**2)}))
pts=[((1,2,2),sp.Rational(1,2)),((1,2,2),sp.Rational(2,3)),((1,2,2),1),((2,3,6),sp.Rational(1,7)),((2,3,6),1),((3,4,12),sp.Rational(1,13))]
print("(32) coordinate condition residuals:",[sp.simplify(sum(D(bt[k,m],m) for m in range(4)).subs(z,sp.sqrt(r**2-x**2-y**2))) for k in range(4)])
print("max|A(4)-A(49)| over points (r=1.5,2,3,1,7,1):",max(abs(ev(A_comp(k,l,rho)-A49(k,l,rho),p,s)) for (p,s) in pts for k in range(4) for l in range(4) for rho in range(4) if ev(r,p,s)!=1))
print("A_(01)0 from (4) on x-axis (limit r->r0 exists):",sp.factor(sp.simplify(A_comp(0,1,0).subs({y:0,z:0,x:r}))))
# matrix form (A_rho)^k_l = eta^{kk} A_(kl)rho ; from (6)-(7): A -> L A L^-1 + dL L^-1  =>  F = dA - [A,A], D = d - [A,.]
Am=[sp.Matrix(4,4,lambda k,l: eta[k,k]*A49(k,l,rho)) for rho in range(4)]
def comm(M,N): return M*N-N*M
def F(s): return {(mu,nu): Am[nu].applyfunc(lambda e:D(e,mu))-Am[mu].applyfunc(lambda e:D(e,nu))+s*comm(Am[mu],Am[nu]) for mu in range(4) for nu in range(4)}
pairs=[(0,1),(0,2),(0,3),(2,3),(3,1),(1,2)]
def FF(Fd):   # sum_a F_a.F^a with Lie labels raised by eta on pairs, spacetime indices by eta
    return sum(eta[k,k]*eta[l,l]*eta[mu,mu]*eta[nu,nu]*(eta[k,k]*Fd[(mu,nu)][k,l])**2 for (k,l) in pairs for mu in range(4) for nu in range(4))
def YM(Fd,s):  # eta^{mu mu}( d_mu F_{mu nu} + s [A_mu, F_{mu nu}] )
    return [sum((eta[mu,mu]*(Fd[(mu,nu)].applyfunc(lambda e:D(e,mu))+s*comm(Am[mu],Fd[(mu,nu)])) for mu in range(4)),sp.zeros(4,4)) for nu in range(4)]
Fm,Fp=F(-1),F(+1)
paper=lambda rv: 4*(6*rv**2+4*rv+1)/(rv**4*(rv+1)**4)
for (p,s) in pts[:3]:
    rv=ev(r,p,s); print(f"r={rv}: F.F(consistent sign)={ev(FF(Fm),p,s)}  draft formula={paper(rv)}  F.F(opposite sign)={ev(FF(Fp),p,s)}")
Em,Ep=YM(Fm,-1),YM(Fp,+1)
for (p,s) in pts:
    print(f"r={ev(r,p,s)}: max|(19) residual| consistent={max(abs(ev(Em[nu][i,j],p,s)) for nu in range(4) for i in range(4) for j in range(4))}  opposite={max(abs(ev(Ep[nu][i,j],p,s)) for nu in range(4) for i in range(4) for j in range(4))}")
zs=sp.sqrt(r**2-x**2-y**2)
print("symbolic (19) residual identically zero:",all(sp.simplify(sp.together(Em[nu][i,j]).subs(z,zs))==0 for nu in range(4) for i in range(4) for j in range(4)))
# metric side (spherical), areal radius, PPN, photon sphere, surface gravity
R_=sp.symbols('R',positive=True); u=sp.symbols('u',positive=True)
B=((2*R_-3)/(2*R_-1))**2; gRR=((2*R_+1)/(2*R_-1))**2
print("B(u)=",sp.expand(sp.series(B.subs(R_,1/u),u,0,4).removeO())," g_RR(u)=",sp.expand(sp.series(gRR.subs(R_,1/u),u,0,3).removeO()))
print("photon sphere R/r0 =",[sp.nsimplify(s_) for s_ in sp.solve(sp.diff(R_**2/B,R_),R_)])
print("Killing kappa(R) = d sqrt(B)/dR / sqrt(g_RR) =",sp.factor(sp.diff((2*R_-3)/(2*R_-1),R_)/((2*R_+1)/(2*R_-1)))," at R=3/2:",sp.factor(sp.diff((2*R_-3)/(2*R_-1),R_)/((2*R_+1)/(2*R_-1))).subs(R_,sp.Rational(3,2)))

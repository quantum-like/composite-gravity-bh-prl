"""Pre-submit check: does Öttinger static solution (49) satisfy YM eq. (19) for r>0?
Source: arXiv:2604.07008 / DOI 10.1088/2399-6528/aeb16c. Uses only (18),(19),(A2),(49), coupling=1.
Lie labels as lowered antisymmetric pair matrices X_(kl); bracket from (A2):
  f_(kl)^{bc} B_b C_c = eta^{k'l'}[B_(k'l)C_(kl') - C_(k'l)B_(kl')] = (C eta B - B eta C)_(kl).
sign s=+1 is paper's f; s=-1 diagnostic. Also checks vacuum source: point mass => T=0 for r>0 => J (23) = 0.
"""
import os
import sympy as sp, mpmath
t,x,y,z=sp.symbols('t x y z',real=True); r0=sp.symbols('r0',positive=True)
X=[t,x,y,z]; r=sp.sqrt(x**2+y**2+z**2); eta=sp.diag(-1,1,1,1)
P=r0/(r**2*(r+r0)); xs=[None,x,y,z]
def A(rho):  # (48)/(49): A_(0l)rho = P x_l d_rho0 ; A_(kl)rho = P(x_l d_k rho - x_k d_l rho)
    M=sp.zeros(4,4)
    for l in (1,2,3):
        if rho==0: M[0,l]=P*xs[l]; M[l,0]=-M[0,l]
    for k in (1,2,3):
        for l in (1,2,3):
            if k!=l: M[k,l]=P*(xs[l]*(1 if k==rho else 0)-xs[k]*(1 if l==rho else 0))
    return M
def br(B,C,s): return s*(C*eta*B-B*eta*C)
out=[]
for s in (1,-1):
    As=[A(m) for m in range(4)]
    F=[[sp.diff(As[n],X[m])-sp.diff(As[m],X[n])+br(As[m],As[n],s) for n in range(4)] for m in range(4)]
    Fu=[[F[m][n]*eta[m,m]*eta[n,n] for n in range(4)] for m in range(4)]
    worst=0
    for nu in range(4):
        L=sp.zeros(4,4)
        for m in range(4):
            L+=sp.diff(Fu[m][nu],X[m])+br(As[m],Fu[m][nu],s)  # A_mu (lower) F^{mu nu}: one eta per index, not two
        f_=sp.lambdify((x,y,z,r0),L,'mpmath')
        # numeric check at points incl r=r0
        for pt in [(0.3,0.4,0.5),(1,0,0),(0.6,-0.8,0.0),(2.0,1.0,-1.0)]:
            v=f_(*pt,1)
            worst=max(worst,max(abs(e) for e in v))
            if nu<2 and pt==(0.3,0.4,0.5): out.append(f's={s:+d} nu={nu} pt={pt}: '+str([[float(e) for e in row] for row in v.tolist()]))
    out.append(f"s={s:+d}: max |LHS(19)| at sample points (r0=1, incl r=r0) = {worst}")
    # Lorenz gauge (38)
    lor=sp.simplify(sum((sp.diff(As[m],X[m])*eta[m,m] for m in range(4)),sp.zeros(4,4)).subs({x:0.3,y:0.4,z:0.5,r0:1}))
    out.append(f"Lorenz (38) d^mu A_mu = {0 if lor==sp.zeros(4,4) else lor}")
txt="\n".join(out); print(txt)
os.makedirs('results', exist_ok=True)
with open('results/verify_49_vs_19_out.txt','w') as f:
    f.write(txt+"\n")

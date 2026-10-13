import sys, json; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))
from sympy.polys.polytools import rem as captured_rem
import qguard, sympy as sp, sympy.polys as spp
d=qguard.d; y=sp.Symbol('y'); Dm=sp.Dummy('D'); X=sp.Symbol('X')
E={"lin":d*X+1/(d-1),"pow":1/(d**2-1),"alone":d,"mixed":(d**2-2*d-2)/(32*(d**2-1))}
num=lambda e: sp.numer(sp.together(e))
P={"polys_ns_rem":lambda e: spp.rem(num(e),d-4,d),"captured_rem":lambda e: captured_rem(num(e),d-4,d),
   "reduced":lambda e: sp.reduced(num(e),[d-4],d),"groebner":lambda e: sp.groebner([num(e)-X, d-4],d,X),
   "Matrix_Dummy":lambda e: sp.Matrix([e]).subs(d,Dm).subs(Dm,4),"rep_eval":lambda e: sp.Poly(num(e),d).rep.eval(4),
   "rem":lambda e: sp.rem(num(e),d-4,d),"resultant":lambda e: sp.resultant(num(e),d-4,d),"compose":lambda e: sp.compose(num(e),4,d),
   "Poly_shift":lambda e: sp.Poly(num(e),d).shift(-4),"Matrix_y":lambda e: sp.Matrix([e]).subs(d,y).subs(y,4),
   "subs":lambda e: e.subs(d,4),"xreplace":lambda e: e.xreplace({d:4}),"evalf":lambda e: e.evalf(subs={d:4}),
   "lambdify":lambda e: sp.lambdify(d,e)(4),"limit":lambda e: sp.limit(e,d,4),"apart_residue":lambda e: sp.apart(num(e)/(d-4),d),"apart_list":lambda e: sp.apart_list(num(e)/(d-4),d)}
miss=[]
for k,f in P.items():
    for n,e in E.items():
        try: f(e); miss.append((k,n,"NOT_RAISED"))
        except qguard.PrematureD4Substitution: pass
        except Exception as x: miss.append((k,n,type(x).__name__))
e = E["mixed"]
alg = {"cancel_factor": str(sp.factor(sp.cancel(e + 1 / (d + 1)))), "domain_only_poly": str(sp.Poly(X**2*e + X, X).as_expr().equals(X**2*e + X)),
       "authorised_value": str(qguard.dred.substitute_d(e, 4, "Q01 probe")), "d_free_rem": str(sp.rem(X**3 + 1, X - 1, X))}
print(json.dumps({"forms": sorted(P), "n_forms": len(P), "n_exprs": len(E), "missed": miss, "algebra": alg,
                  "entry_points": qguard.ENTRY_POINTS}))

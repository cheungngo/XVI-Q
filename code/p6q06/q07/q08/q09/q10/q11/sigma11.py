#!/usr/bin/env python3
"""XVI-Q Q11 goal 1: BRST-trivial directions as genuine sigma-images in the gauge-fixed BV variables.
Gauge fixing by the canonical shift with Psi = cbar.F(h) (B integrated out) means the antifield of h enters only
through Ktilde = K + cbar.F'(h), and the ghost equation (Q08, rho = -2) is exactly this statement on the amplitudes:
Gamma_{cbar..} = rho * Gamma_{K..}[K -> T(cbar)].  A local gh -1 functional must therefore be built from Ktilde:
  Y_F = Ktilde.F(h,h) :  (Y_K, Y_cbar) = (K.F, rho*T(K.F))   -- Q08 kept only K.F (no cbar part)
  Y_G = c*.G(c,h)     :  sigma components X_Khc (c* -> D^T Ktilde), X_cbar ch = rho*T(X_Khc), X_c*cc (h -> s0 h)
                        -- Q08 had no X_c*cc piece and a hand-made cbar piece
sigma = (S0, .) is then applied with the SAME channel operators that build the rows (D-op: K -> dS2/dh; A-op:
h -> s0 h; E-op: cbar -> sigma_F F(h)), with the antibracket's relative sign -1 between antifield-type (D, E) and
field-type (A) pieces (Q09 derivation).  No per-generator fitting: one test per generator of
   rows(X) = 0, Bose(X) = 0, ghost-equation(X) = 0.
The only convention not fixed a priori is the sign eps_S of the c*cc piece relative to the c-ordering of the
rows (two identical Grassmann c's); both signs are reported."""
import json, pickle, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
HERE = Path(__file__).resolve().parent
Q = HERE.parents[6]
sys.path[:0] = [str(HERE.parents[5] / "p6q05"), str(HERE.parents[4]), str(HERE.parents[3]), str(HERE.parents[2])]
import rows as R  # noqa: E402
import gravity as G  # noqa: E402
import uvpole as U  # noqa: E402
import hhc  # noqa: E402
from evfit import monomials  # noqa: E402
from geq import Tmap  # noqa: E402
import triv  # noqa: E402
LOG = Q / "build" / "logs"
I = sp.I
RHO = -2
MS = {"K": monomials({"aK": 2, "gc1": 1, "a3": 2}, 3), "C": monomials({"gb2": 1, "gc1": 1, "a3": 2}, 4),
      "S": monomials({"z": 1, "gc1": 1, "gc2": 1}, 3), "H": monomials({"a1": 2, "a2": 2, "a3": 2}, 4)}
NK, NC, NS, NH = (len(MS[k]) for k in "KCSH")
OFF = {"K": 0, "C": NK, "S": NK + NC, "H": NK + NC + NS}
NDIM = OFF["H"] + NH
SW23 = {"p1": R.P1, "p2": R.P3, "a2": "a3", "a3": "a2"}
REL = {"1": {}, "2": {"p1": R.P2, "p2": R.P1, "a2": "a1"}, "3": {"p1": R.P3, "p2": R.P1, "a2": "a1", "a3": "a2"}}
MOM = {"1": hhc.p1, "2": hhc.p2, "3": hhc.p3}
SIG = None


def vec(sector, e):
    v, out = triv.coords(sp.expand(e), MS[sector], sector)
    w = [0] * NDIM
    for i, c in enumerate(v):
        w[OFF[sector] + i] = c
    return w, out != 0


def add(*vs):
    return [sum(x) for x in zip(*vs)]


def scale(c, v):
    return [c * x for x in v]


def D_op(f):
    """X_hhh from Y = K.f (K at p1 with aK; h a2 at p2, a3 at p3; f symmetric in the h legs)."""
    xh = 0
    for i in "123":
        fi = R.remap(f, REL[i]) if REL[i] else f
        fi = R.polarise(fi, "aK", "KU", "KW")
        eh, gf = G.S2_EH_gf(hhc.ext(f"a{i}"), G.probe_mat("E"), MOM[i])
        xh += R.strip_t(R.contract_id(fi, I * (eh + gf), "K", "E", "h"))
    return xh


def A_on_K(f):
    """X_Khc from Y = K.f : h -> s0 h on the h legs, layout K(aK) p2, c(gc1) p1, h(a3) p3."""
    fk = R.remap(f, {"p1": R.P2, "p2": R.P1})
    return R.strip_t(-I * R.polarise(fk, "a2", "p1", "gc1"))


def E_op(phi):
    """X_hhh from Y = cbar.phi (cbar gv at p1; h a2 p2, a3 p3): cbar -> sigma_F F(h)."""
    xh = 0
    for i in "123":
        fi = R.remap(phi, REL[i]) if REL[i] else phi
        D = G.Dvec(hhc.ext(f"a{i}"), U.V({"p1": 1}) if i == "1" else (hhc.F2 if i == "2" else hhc.F3))
        xh += SIG * hhc.subst_vec(fi, "gv", D)
    return R.strip_t(xh)


def A_on_C(phi):
    fc = R.remap(phi, {"p1": R.P2, "p2": R.P1})
    fc = hhc.subst_vec(fc, "gv", {"gb2": 1})
    return R.strip_t(-I * R.polarise(sp.expand(fc), "a2", "p1", "gc1"))


def cbar_partner(f):
    """phi with Y_cbar = cbar.phi = rho*T(K.f): T acts at the K momentum (moved to p2 for geq.Tmap), back to
    the cbar(gv)@p1 layout."""
    f2 = R.remap(f, {"p1": R.P2, "p2": R.P1, "a2": "a1"})            # K@p2, h a1@p1, a3@p3
    t = sp.expand(RHO * Tmap(f2))                                       # cbar(gb2)@p2
    back = R.remap(t, {"p1": R.P2, "p2": R.P1, "a1": "a2", "gb2": "gv"})  # cbar(gv)@p1, a2@p2, a3@p3
    return sp.expand(back)


def family_F():
    out = []
    for m in monomials({"aK": 2, "a2": 2, "a3": 2}, 2):
        f = sp.expand(m + R.remap(m, SW23))
        if f == 0:
            continue
        phi = cbar_partner(f)
        vh1, l1 = vec("H", D_op(f))
        vk, l2 = vec("K", A_on_K(f))
        vh2, l3 = vec("H", E_op(phi)) if phi != 0 else ([0] * NDIM, False)
        vc, l4 = vec("C", A_on_C(phi)) if phi != 0 else ([0] * NDIM, False)
        out.append({"anti": add(vh1, vh2), "field": add(vk, vc), "Q08_anti": vh1, "Q08_field": vk,
                    "leak": l1 or l2 or l3 or l4})
    return out


def family_G():
    out = []
    for m in monomials({"gv": 1, "gc1": 1, "a3": 2}, 2):            # c*(gv)@p2, c(gc1)@p1, h(a3)@p3
        xk = sp.expand(-I * U.H("aK", "p2") * hhc.subst_vec(m, "gv", {"aK": 1}))   # c* -> D^T Ktilde (K part)
        xc = sp.expand(RHO * Tmap(xk))                                             # its cbar partner
        s = R.remap(m, {"p1": R.P1, "p2": R.P3, "gv": "z"})                       # c*(z)@p3, c@p1, h(a3)@p2
        s = R.strip_t(-I * R.polarise(sp.expand(s), "a3", "p2", "gc2"))            # h -> s0 h = D c(gc2)@p2
        s_sw = R.remap(s, {"p1": R.P2, "p2": R.P1, "gc1": "gc2", "gc2": "gc1"})
        xs_anti = sp.expand(s - s_sw)
        vk, l1 = vec("K", xk)
        vc, l2 = vec("C", xc)
        vs, l3 = vec("S", xs_anti)
        out.append({"anti": add(vk, vc), "field": vs, "leak": l1 or l2 or l3})
    return out


def main():
    global SIG
    t0 = time.time()
    SIG = sp.sympify(json.loads((LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    d = pickle.loads((LOG / "q08_triv.pkl").read_bytes())
    A, Ag = d["A"], d["Ag"]
    B = pickle.loads((LOG / "q11_blocks.pkl").read_bytes()) if (LOG / "q11_blocks.pkl").is_file() else None

    def zero(M, v):
        return all(sp.expand(x) == 0 for x in M * sp.Matrix(v))

    res = {}
    keep = {}
    for name, fam in (("F", family_F()), ("G", family_G())):
        st = {"n": len(fam), "leak": sum(g["leak"] for g in fam), "zero_vector": 0}
        for lab, lam in (("derived_minus1", -1), ("plus1", 1)):
            ok, okg, vecs, fails = 0, 0, [], {}
            for g in fam:
                v = add(g["anti"], scale(lam, g["field"]))
                if not any(sp.expand(x) != 0 for x in v):
                    continue
                if zero(A, v):
                    ok += 1
                    vecs.append(v)
                    okg += zero(Ag, v)
                elif B is not None:
                    k = " ".join(b for b, M in B.items() if not zero(M, v))
                    fails[k] = fails.get(k, 0) + 1
            st[lab] = {"in_N": ok, "in_N_and_ghost_eq": okg, "fail_blocks": fails}
            keep[(name, lab)] = vecs
        st["zero_vector"] = sum(1 for g in fam if not any(sp.expand(x) != 0 for x in add(g["anti"], g["field"]))
                                and not any(sp.expand(x) != 0 for x in g["anti"]))
        res[name] = st
    pickle.dump(keep, open(LOG / "q11_sigma_vecs.pkl", "wb"))
    res["seconds"] = round(time.time() - t0, 1)
    (LOG / "q11_sigma.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()

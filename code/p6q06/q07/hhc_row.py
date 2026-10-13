#!/usr/bin/env python3
"""XVI-Q Q07: assemble row h.h.c from the hhc.py channel files; residual, 4D part, matter-sector check, unconstrained fit,
mutations."""
from __future__ import annotations
import json, sys  # noqa: E401
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
Q = HERE.parents[2]
sys.path.insert(0, str(HERE.parents[1] / "p6q05"))
sys.path.insert(0, str(HERE))
import rows as R  # noqa: E402
import hhc  # noqa: E402
from evct import four_d  # noqa: E402
LOG = Q / "build" / "logs"


def rd(k):
    return json.loads((LOG / f"q07_hhc_{k}.json").read_text())


def parts():
    sig = sp.sympify(rd("tree")["sigma_F"])
    P = {"A": sp.sympify(rd("A")["value"]) if (LOG / "q07_hhc_A.json").is_file() else hhc.chan_A()}
    for i in ("2", "3"):
        b = rd("B" + i)
        for s in ("grav", "ghost", "scalar"):
            P[f"B{i}_{s}"] = sp.sympify(b[s])
        P["D" + i] = sp.sympify(rd("D" + i)["value"])
        P["E" + i] = sig * sp.sympify(rd("E" + i)["value"])
    P["C"] = sp.sympify(rd("C")["value"])
    return P, sig


def row():
    P, sig = parts()
    Rr = sp.expand(sum(P.values()))
    R4 = four_d(Rr)
    # matter (scalar-loop) sector alone: A_scalar + B_scalar
    As = R.strip_t(-sp.I * R.polarise(sp.expand(hhc.load_any("G14") + hhc.load_any("G15")), "a1", "p1", "gc1"))
    matter = sp.expand(As + P["B2_scalar"] + P["B3_scalar"])
    groups = {"A": ["A"], "B": [k for k in P if k.startswith("B")], "C": ["C"], "D": ["D2", "D3"], "E": ["E2", "E3"]}
    muts = {f"drop_{g}": four_d(sp.expand(Rr - sum(P[k] for k in ks))) != 0 for g, ks in groups.items()}
    muts["flip_sigma_F"] = four_d(sp.expand(Rr - 2 * (P["E2"] + P["E3"]))) != 0
    cs = sp.symbols("y0:5")
    G4 = {g: four_d(sp.expand(sum(P[k] for k in ks))) for g, ks in groups.items()}
    tot = sp.expand(sum(c * G4[g] for c, g in zip(cs, groups)))
    gens = sorted([s for s in tot.free_symbols if s not in cs], key=str)
    fit = sp.linsolve(sp.Poly(tot, *gens).coeffs(), cs)
    return {"sigma_F": str(sig), "channel_terms": {k: len(sp.Add.make_args(v)) for k, v in P.items()},
            "residual": str(Rr), "residual_zero": Rr == 0, "residual_4D": str(R4), "residual_4D_zero": R4 == 0,
            "residual_n_terms": len(sp.Add.make_args(Rr)) if Rr != 0 else 0,
            "residual_4D_n_terms": len(sp.Add.make_args(R4)) if R4 != 0 else 0,
            "matter_sector_closes": matter == 0, "matter_sector_4D_closes": four_d(matter) == 0,
            "unconstrained_fit_4D": {"groups": list(groups), "solution": [str(x) for x in list(fit)[0]] if fit else None},
            "mutations_detected": muts}


if __name__ == "__main__":
    r = row()
    (LOG / "q07_hhc_row.json").write_text(json.dumps(r, indent=1) + "\n")
    print({k: (v if not isinstance(v, str) or len(v) < 150 else v[:150]) for k, v in r.items()})

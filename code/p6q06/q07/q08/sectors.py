#!/usr/bin/env python3
"""XVI-Q Q08: sector decomposition of the homogeneous space N (from triv.py's pickled matrices): kernel dims of the
rows restricted to single sectors, and the trivial span (best fitted lambda) intersected with the ghost-equation
subspace."""
import json, pickle, sys  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix
Q = Path(__file__).resolve().parents[4]
LOG = Q / "build" / "logs"
d = pickle.loads((LOG / "q08_triv.pkl").read_bytes())
A, Ag, fams = d["A"], d["Ag"], d["fams"]
tr = json.loads((LOG / "q08_triv.json").read_text())
nK, nC, nS = 134, 81, 18
cuts = {"Kch": (0, nK), "cbar_ch": (nK, nK + nC), "cstar_cc": (nK + nC, nK + nC + nS), "hhh": (nK + nC + nS, A.shape[1])}
out = {}
for k, (a, b) in cuts.items():
    sub = A[:, a:b]
    out[f"dim_N_{k}_only"] = (b - a) - DomainMatrix.from_Matrix(sub).rank()
lams = {n: sp.sympify(v["best_uniform_lambda"]) for n, v in tr["families"].items()}
B = sp.Matrix.hstack(*[sp.Matrix(p) + lams[n] * sp.Matrix(q) for n, pr in fams.items() for p, q in pr])
AA = A.col_join(Ag)
rB = DomainMatrix.from_Matrix(B).rank()
out["dim_trivial_best_lambda_cap_N_and_ghost_eq"] = rB - DomainMatrix.from_Matrix(AA * B).rank()
(LOG / "q08_sectors.json").write_text(json.dumps(out, indent=1) + "\n")
print(out)

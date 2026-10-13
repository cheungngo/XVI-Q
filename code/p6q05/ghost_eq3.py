#!/usr/bin/env python3
"""XVI-Q Q05 diagnostic: three-point ghost equation (G04+G09) vs (G03+G07)[K -> T(cbar)], two-point D06 vs D13[K->T]."""
import json, sys  # noqa: E401
from pathlib import Path
import sympy as sp
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rows as R, uvpole as U, tensor as TT, graphs as GR  # noqa: E401,E402
k = U.V({"p2": 1})   # untracked by the Taylor parameter t (G03/G07 extract t^3)
gb = U.vec("gb2")
Tm = TT.sym(k, gb)
Tm.eta = -sp.Rational(1, 2) * U.dot(k, gb)
orig = GR.ext
GR.ext = lambda a: Tm if a == "aK" else orig(a)
g3, g7 = R.strip_t(GR.G03()), R.strip_t(GR.G07())
GR.ext = orig
g4, g9 = R.strip_t(R.load("G04")), R.strip_t(R.load("G09"))
out = {"ratio_total": str(sp.cancel((g4 + g9) / (g3 + g7))) if (g3 + g7) != 0 else "den0",
       "ratio_G04_G03T": str(sp.cancel(g4 / g3)) if g3 != 0 else "den0",
       "ratio_G09_G07T": str(sp.cancel(g9 / g7)) if g7 != 0 else "den0"}
(HERE.parents[1] / "build/logs/q05_ghost_eq3.json").write_text(json.dumps(out, indent=1) + "\n")
print(out)

#!/usr/bin/env python3
"""XVI-Q Q06 independent second route for P6-G10 (c* c c): loop momentum relabelled so that the internal graviton
carries l - p1 (ghost lines l and p1 + p2 - l), and the V_{c* c c} vertex is entered with its two c slots in the
opposite order (compensated by the Grassmann sign -1).  Must equal the Q05 G10 pole exactly."""
from __future__ import annotations
import sys
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "p6q05"))
import graphs as GR  # noqa: E402
import gravity as G  # noqa: E402
import uvpole as U  # noqa: E402

p1, p2, l = GR.p1, GR.p2, GR.l


def G10_route2():
    vA = G.V_cb_c_h(G.probe_vec("xAb"), U.vec("gc1"), p1, G.probe_mat("hA"), l - p1)
    vB = G.V_cb_c_h(G.probe_vec("xBb"), U.vec("gc2"), p2, G.probe_mat("hB"), p1 - l)
    vC = G.V_cs_c_c(U.vec("z"), G.probe_vec("xCb"), p1 + p2 - l, G.probe_vec("xCa"), l)
    n = G.contract_vertices(vA, vB, [("h", "hA", "hB")])
    n = G.contract_vertices(n, vC, [("c", "xAb", "xCa", GR.GH), ("c", "xBb", "xCb", GR.GH)])
    return GR.r4(U.loop_pole(-n, [U.V(), -p1, -(p1 + p2)], 3, 1), "G10 route2")


if __name__ == "__main__":
    r = G10_route2()
    print(sp.expand(r - GR.G10()) == 0)

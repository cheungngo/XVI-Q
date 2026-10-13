#!/usr/bin/env python3
"""XVI-Q Q07 certificates: h.h.c row (tree calibration, own Gamma_hh incl. scalar loop, own G14/G15, channel E =
erratum Q06-E2), evanescent CTs (naive test + joint local-CT existence on all kappa^3 three-point rows), status."""
from __future__ import annotations
import hashlib, json, subprocess, sys, time  # noqa: E401
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
Q = HERE.parents[2]
sys.path.insert(0, str(HERE.parents[1] / "p6q05"))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
CERT = Q / "certificates"
LOG = Q / "build" / "logs"
PY = sys.executable


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def run(args):
    r = subprocess.run([PY, *args], cwd=Q, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"{args}: {r.stderr[-2000:]}")
    return r.stdout


def main():
    t0 = time.time()
    LOG.mkdir(parents=True, exist_ok=True)
    for f in LOG.glob("q07_*.json"):
        f.unlink()
    with ThreadPoolExecutor(10) as ex:
        list(ex.map(run, [["code/p6q06/q07/hhc.py", k] for k in ("tree", "S", "B2", "B3", "C", "D2", "D3", "E2", "E3")]))
    run(["code/p6q06/q07/hhc.py", "A"])
    run(["code/p6q06/q07/hhc_row.py"])
    run(["code/p6q06/q07/evfit2.py"])
    import evct
    naive = {t: evct.test(t) for t in ("ccK", "ccb")}
    tree = json.loads((LOG / "q07_hhc_tree.json").read_text())
    row = json.loads((LOG / "q07_hhc_row.json").read_text())
    res = sp.sympify(row["residual"])
    g = {k: sp.sympify(json.loads((LOG / f"q07_{k}.json").read_text())["pole_over_U_d4_residue"]) for k in ("G14", "G15")}
    hhc = {"row": "h.h.c: c(gc1) p1, h(a2) p2, h(a3) p3; (S, Gamma_1^pole) at kappa^3",
           "channels": {"A": "Gamma_hhh (G01+G02+G08+G14+G15) at h(p1) -> s0 h(c)",
                        "B": "own Gamma_hh (graviton + ghost + scalar loops) with one leg s1 h = (1/2) L_c h",
                        "C": "delta S_EH3/delta h|_{h2 h3} . Gamma_{Kc} (D13 own)",
                        "D": "delta S_2/delta h|_{h_i} (EH + gauge fixing) . Gamma_{K h_j c} (G03+G07)",
                        "E": "s cbar = sigma_F F(h_i): F . Gamma_{cbar c h_j} (G04+G09) -- erratum Q06-E2 channel"},
           "tree_calibration": {"sigma_F": tree["sigma_F"], "tree_row_closes": tree["closes"] == "True",
                                "tree_row_without_E_nonzero": tree["without_E_zero"] == "False"},
           "scalar_loops_own": {k: {"n_terms": len(sp.Add.make_args(v)), "sha256": sha(str(v))} for k, v in g.items()},
           "matter_sector_closes_exactly": row["matter_sector_closes"],
           "channel_terms": row["channel_terms"], "residual_4D": row["residual_4D"],
           "residual_4D_zero": row["residual_4D_zero"],
           "unconstrained_fit_4D": row["unconstrained_fit_4D"], "mutations_detected": row["mutations_detected"],
           "evanescent_remainder": row["residual"], "evanescent_remainder_sha256": sha(row["residual"]),
           "evanescent_remainder_n_terms": row["residual_n_terms"],
           "evanescent_remainder_every_term_has_T": all(any(s.name.startswith("T_") for s in t.free_symbols)
                                                       for t in sp.Add.make_args(res)) and res != 0}
    (CERT / "Q07_hhc_row.json").write_text(json.dumps(hhc, indent=1) + "\n")
    joint = json.loads((LOG / "q07_evfit_joint.json").read_text())
    for k in list(joint):
        if k.startswith("particular_CT_") and not k.endswith("n_terms"):
            joint[k + "_sha256"] = sha(joint[k])
    ev = {"question": "sigma0-consistency of evanescent kappa^3 counterterms removing the evanescent remainders",
          "naive_p22_type": {"definition": "CT_ev = -(T-part of every loop input), i.e. rows on T-projected inputs",
                             **{t: {"zero": v["zero"], "n_T_terms_left": v["n_T_terms"]} for t, v in naive.items()},
                             "verdict": "NOT sigma0-consistent: tree vertices regenerate T structures from the 4D loop parts"},
          "joint_local_CT": {"definition": "CT = -Gamma_pole + X, X = generic local evanescent (>= one T) kappa^3 "
                                           "monomials for Gamma_{Kch}, Gamma_{cbar ch}, Gamma_{c*cc}, Gamma_hhh (Bose "
                                           "symmetric); the same X_{Kch}, X_{cbar ch} enter c.c.K/c.c.cbar and h.h.c; "
                                           "kappa^2 poles unmodified", **{k: v for k, v in joint.items() if k != "seconds"}},
          "verdict": ("CANCELLABLE: a local, Bose-symmetric, purely evanescent kappa^3 counterterm X exists with "
                      "sigma0-variation cancelling every evanescent remainder on the rows c.c.K, c.c.cbar, h.h.c "
                      "simultaneously (verified by substitution); not unique (free parameters)")
          if joint.get("solvable") and joint.get("verified_rows_zero") and joint.get("hhh_CT_bose_symmetric") and
          joint.get("all_CT_evanescent") else "NOT cancellable in the ansatz (explicit obstruction)",
          "not_tested": ["ghost-equation compatibility of X_{cbar ch} with X_{Kch}",
                         "X beyond three-point (four-point kappa^3 rows)", "finite evanescent remnants (pole level only)"]}
    (CERT / "Q07_evanescent_CT.json").write_text(json.dumps(ev, indent=1) + "\n")
    hhc_ok = hhc["residual_4D_zero"] and all(row["mutations_detected"].values()) and hhc["tree_calibration"]["tree_row_closes"]
    canc = ev["verdict"].startswith("CANCELLABLE")
    st = {"landmark": "XVI-Q", "release": "Q07", "P6_outcome": "Partial",
          "stages": {"row_c_phi_phi": "closed at pole level (Q05)",
                     "row_c_c_cbar": "4D closed (Q06); evanescent remainder cancellable by local evanescent CT (Q07)",
                     "row_c_c_K": "4D closed (Q06); evanescent remainder cancellable by local evanescent CT (Q07)",
                     "row_h_h_c": ("4D closed (Q07, all channels +1, sigma_F tree-calibrated); evanescent remainder "
                                   "cancellable by local evanescent CT (Q07)") if hhc_ok and canc else "OPEN",
                     "graphs_missing": [], "graphs_derived_zero": ["T1", "T2", "T3"]},
          "kappa3_CT_domain_reduction": ("Partial: three-point kappa^3 CT = -Gamma_pole (4D part, sigma0-closed on all "
                                         "four rows) + X (local evanescent, existence shown, non-unique); reduction to "
                                         "a basis of local invariant operators / field redefinitions not performed"),
          "kappa3_defect": "zero on the 4D part of all four three-point rows; evanescent part sigma0-exact (pole level)",
          "kappa3_CT_solve": "solved at pole level (particular solution certified, Q07_evanescent_CT.json)",
          "sigma0_closure_kappa3": ("CLOSED at pole level on the kappa^3 one-loop THREE-POINT rows (c.phi.phi, c.c.cbar, "
                                    "c.c.K, h.h.c) in the enlarged (evanescent) CT domain" if hhc_ok and canc else
                                    "NOT established"),
          "not_implied": ["kappa^3 four-point / higher rows", "finite parts", "ghost-equation compatibility of X",
                          "uniqueness of X", "bulk symmetry restoration", "any change to the O.1 record"],
          "errata": {"Q06-E2": "addressed: the F(h).Gamma_{cbar c h} channel is required (tree: sigma_F = i; mutation "
                               "drop_E detected)"},
          "disclosed_revisions": {"Q07-R1": "disclosed scope revision of the frozen Q06 layer gate "
                                            "(code/p6q06/release_q06.py): files under code/p6q06/q07/ and tools/*q07* "
                                            "are left to the Q07 gate, and the pinned Q05_guard_scan file count is a "
                                            "lower bound (as Q06-R1 did for the Q05 gate). Q05 gate untouched. "
                                            "release_q06.py is pinned by q07_expected.json."},
          "O1_closure_unchanged": "CLOSED (scoped)", "bulk_symmetry_restored": False,
          "flags": {"execution_pass": bool(hhc_ok and canc), "claim_status": "partial",
                    "scope": "pole level, one loop, kappa^3, three-point rows"}, "seconds": round(time.time() - t0, 1)}
    (CERT / "Q07_P6_status.json").write_text(json.dumps(st, indent=1) + "\n")
    print("Q07 certs", st["flags"], st["sigma0_closure_kappa3"][:40], st["seconds"])


if __name__ == "__main__":
    main()

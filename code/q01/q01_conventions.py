"""XVI-Q Q01 conventions: pins and EVERY certified number the gate value-checks."""
LANDMARK, RELEASE, VERSION_LABEL, DATE = "XVI-Q", "Q01", "XVI-Q Q01", "2026-10-10"
XVI_P_P23_ZIP_SHA256 = "3f4a88db462923553272e8a41407c0c06aacbc431736efdfeda10c27aae8a177"
P23_INDEPENDENT_REPORT_SHA256 = "2e3e5bf325c9a0e7"            # prefix, full sha in intake/imports.json
O1_LEDGER_SHA256 = "3cd3d1ba491381ba0a73d1677bfed10d3743e0575f61c71c325e626e4a78a645"
EXPORT_SHA256 = {   # FILE sha256 of the exported .srepr (byte-identical to the published p23 exports; P13 records polynomial content shas)
    2: "8ef4395c51434fe073109065e31773460c3708937a16ea29a036abf6c72c5030",
    3: "183798063fa89839e9a80ae4a81e6c2f9c491a7ad386540d1a4e33382b0686de",
    4: "d71d9a839b5e6388abe32ccb6bb334b87cca8e104f98a99ec198b923dde099fc"}
EXPECTED = {
    "ell4_n_terms": 3360, "ell4_normalisation_from_XVI_O_probe": "-2", "ell4_full_polynomial_equal": True,
    "ell4_probe_equal": True,
    "W_terms": {"kappa0": 132, "kappa1": 2148, "kappa2": 15428},
    "E_c_W_all_zero": {"kappa0": True, "kappa1": True, "kappa2": True},
    "master_mutation_residuals": {"flip_sigma_c": 6, "drop_cstar": 6, "alter_gf_shift": 2,
                                  "flip_odd_derivative_sign": 6, "corrupt_a2": 1, "corrupt_ell4": 1,
                                  "drop_s1_transport": 2},
    "classical_certs_content_identical_to_p23": ["P13_quartic_verified", "P13_brst_derivation", "P13_master_kappa2",
                                                 "P13_boundary_identity", "P13_gauge_fixing_map"],
    "guard_forms": 18, "guard_exprs": 4,   # Q02 revision: guard_forms 16 -> 18 (+ apart_residue, apart_list)
    "O2_items": {"a_ell4": "Constructed", "b_antibracket_expansion": "Constructed", "c_ghost_antifield": "Constructed"},
}
HARD_STOPS = {"O1_closure": "NOT_CLOSED", "bulk_symmetry_restored": False, "insertion": "frozen",
              "XV_B_opened": False, "XIV_J_opened": False, "Route_B_discharged": False, "Missing_never_zero": True}

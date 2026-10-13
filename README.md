# XVI-Q — release Q01 (O.2 classical κ² BV master on S_inv Constructed; O.1 NOT_CLOSED; bulk symmetry not claimed) — 2026-10-10

XVI-Q is a new landmark, registered locally only (registration/local_registration.json; not published, no GitHub
edit). XVI-O, XVI-N and XVI-P are read-only and imported by sha256 (intake/imports.json); XVI-P is frozen at p23,
zip 3f4a88db…. Errata from the p23 review: intake/errata_p23.json.

O.2 (2(S0,S2)+(S1,S1)=0):
- (a) full ℓ₄: Constructed. The export is independently reconstructed (3360/3360 terms), with the normalisation
  matching XVI-O's probe and convention.
- (b) explicit antibracket expansion: Constructed. E_c(W₂) = 0 exactly; W₂ has 15428 terms; independent engine.
- (c) ghost/antifield sector: Constructed. Evidence is the re-executed exact P13 master with s₂ derived zero; there
  is no second independent engine for this sector.

Guard: code/q01/qguard.py (p23 bypasses closed). Gate: code/q01/release_q01.py value-checks every certified number.
Build: `python3 tools/build.py`. Rebuild: `python3 tools/rebuild_from_source.py` (REBUILD_INSTRUCTIONS.md).
Hard stops: O.1 NOT_CLOSED; bulk_symmetry_restored false; insertion frozen; XV-B/XIV-J unopened; Route B not
discharged; XVI-M pin unchanged; Missing never zero.


## Q04 (current; P.6 entry, honest partial)

One-loop kappa^3 three-point sector: inventory 15 (== XVI-P P23) + 3 scaleless kappa^3 tadpoles; kappa^3 identity ledger (14 gh1 rows, 4 need P.6 inputs); own UV-pole engine validated vs XVI-P D11 (all 14 structures); scalar-loop hhh (2 graphs) at pole level with the matter Ward identity verified; 13 graphs Missing; kappa^3 CT/closure Missing. O.1 CLOSED (scoped) unchanged.

## Q03 (correction release)

O.1 CLOSED (scoped), independently confirmed by the Q02 review (report sha e4b54702...). Q03 records the enlarged counterterm domain (p14r 11 4D directions + p22 evanescent operators), states that the XVI-P prescription (not XVI-O o70-sync kernel; reason 8 not answered) is used, lists the scope, what it does not imply, and the caveats (certificates/Q03_O1_record.json). Gate re-hashes provenance/code/report and checks gate reads are in the source set; d guard closes the PolyRing/FracField bypass. No new science.

## Q02
Bulk O.5 graphs D04 (=0), D05a (=0), D06 (= XVI-P Sigma), D12 reconciled; ST4a/ST4b on 23 bulk two-point rows, outcome cancellable, full sigma_0 closure; O.1 walk: CLOSED (scoped: one-loop kappa^2 bulk two-point, m=0, declared prescription, external momenta in QdS). bulk_symmetry_restored false. Build: `python3 tools/build.py` (no sync by default; `--sync <dir>` to opt in). Q01 archived in archival/q01.

## Q06 (P.6 continued; honest partial)
- Ghost rows c.c.cbar / c.c.K: Q05's OPEN residual traced to a Grassmann sign in the K-channel assembly; with the fix the
  4D part of both rows vanishes (unique fit, all mutations detected); a purely evanescent remainder is recorded (OPEN).
- G10 recomputed by a second route: identical.  V4_grav derived and validated by the exact EH gauge identity.
- G01 (V3^3) and G02 (V3-V4) evaluated at pole level with a ring fast path (equal to the Expr path on G05/G06); full
  Bose symmetry.  All 15 P.6 graphs now have poles.
- Missing: h.h.c row, kappa^3 counterterm domain / defect / CT solve, sigma0 closure.  O.1 CLOSED (scoped) unchanged.
- Build: `python3 tools/build_q06.py` (no sync by default); rebuild: `python3 tools/rebuild_q06.py`.

## Q05 (P.6 continued; honest partial)
Own gravity and ghost rules, validated against XVI-P D01/D03/D06/D13 (constant ratios I, I, I, -I; the D13 orientation is
disclosed). G03-G13 evaluated at pole level. Row c.phi.phi closes exactly at the pole. Rows c.c.cbar and c.c.K are OPEN
(unresolved assembly; not an obstruction). G01 and G02 are Missing. sigma0 closure at kappa^3 is not established.
O.1 CLOSED (scoped) is unchanged. Build: `python3 tools/build_q05.py [--sync DIR]`; rebuild: `python3 tools/rebuild_q05.py`.

## Q07 (P.6 three-point rows)
h.h.c row assembled (own Gamma_hh incl. scalar loop, own G14/G15, F(h) channel with sigma_F = i tree-calibrated): 4D part
closed. Evanescent remainders of c.c.K, c.c.cbar, h.h.c cancelled jointly by a local evanescent kappa^3 counterterm
(existence, non-unique). Code: code/p6q06/q07/; gate code/p6q06/q07/release_q07.py; build tools/build_q06_q07.py.
Disclosed revision Q07-R1 (Q06 gate scope). Status: honest partial.

## Q08 (evanescent counterterm structure)
X ghost-equation compatible (rho = -2); free directions 152 -> 104; >= 97 (fitted) BRST-trivial; basis, physicality,
p22 consistency open. Code code/p6q06/q07/q08/; gate code/p6q06/q07/q08/release_q08.py; build tools/build_q06_q07_q08.py.
Disclosed revision Q08-R1 (Q07 gate scope).

## Q09 (review fixes, derived trivial directions)
sigma0 status relabelled; Q04 report bound by hash; guard_scan3 extension. Code code/p6q06/q07/q08/q09/; build tools/build_q06_q07_q08_q09.py. Disclosed revision Q09-R1 (Q08 gate scope).

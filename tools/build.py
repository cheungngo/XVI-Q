#!/usr/bin/env python3
"""XVI-Q Q04 full build (Q04 = P.6 entry layer: q04 certs, Q03 gate as layer gate, Q04 gate owns RELEASE_STATUS).  Q03: (Q03 = correction layer on Q02: Q03 O.1 record + ring/field guard certs, Q02 gate run as layer gate,
Q03 gate owns RELEASE_STATUS).  Q02 description: (from scratch): clean -> intake/extract -> classical re-execution of XVI-P p23 stages ->
7 master mutations -> Q01 certificates (o2_q01) -> Q02 certificates (bulk graphs, bulk identity, O.1 walk) ->
Q01 layer gate -> Q02 gate (value-checks every certified number) -> tests -> PDFs (science/appendix/combined) ->
zips (source-only/review-focus/full) -> external manifest LAST -> OPTIONAL sync (Q02: default NONE; release opt-in: --sync <dir> or Q_SYNC_DIR)."""
from __future__ import annotations
import hashlib, json, os, shutil, subprocess, sys, time, zipfile  # noqa: E401
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

Q = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Q / "tools"))
import source_set as SS  # noqa: E402

PY = sys.executable
SYNC = (sys.argv[sys.argv.index("--sync") + 1] if "--sync" in sys.argv else os.environ.get("Q_SYNC_DIR", "none"))  # Q02: NO sync by default
RUN = Q / "build" / "xvip_run"
STAGES = ["code/classical/inherited.py", "code/classical/boundary_identity.py", "code/classical/quartic.py",
          "code/classical/brst_from_geometry.py", "code/classical/gauge_fixing_map.py", "code/classical/master_kappa2.py"]
MUTS = ["flip_sigma_c", "drop_cstar", "alter_gf_shift", "flip_odd_derivative_sign", "corrupt_a2", "corrupt_ell4",
        "drop_s1_transport"]
PDFS = ["XVI-Q-science", "XVI-Q-appendix", "XVI-Q"]


def sh(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def step(name, cmd, cwd=Q, env=None):
    t0 = time.time()
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, env=env)
    print(f"[{name}] rc={r.returncode} {time.time() - t0:.1f}s {r.stdout.strip()[-300:]}", flush=True)
    if r.returncode != 0:
        print(r.stdout[-3000:], r.stderr[-3000:])
        raise SystemExit(f"stage {name} failed")
    return r


def clean():
    for d in ("certificates", "ledgers", "dist", "imports"):
        shutil.rmtree(Q / d, ignore_errors=True)
    if (Q / "build").is_dir():                       # keep build/logs (persistent logs), wipe everything else
        for x in (Q / "build").iterdir():
            if x.name != "logs":
                shutil.rmtree(x) if x.is_dir() else x.unlink()
    for f in ["RELEASE_STATUS.json", "RELEASE_MANIFEST_EXTERNAL.txt"] + [f"manuscript/{p}.pdf" for p in PDFS]:
        (Q / f).unlink(missing_ok=True)
    for d in ("certificates", "ledgers", "dist", "build", "imports"):
        (Q / d).mkdir(exist_ok=True)


def extract():
    z = zipfile.ZipFile(Q / "provenance/XVI-P/XVI-P-p23-source-only.zip")
    assert sh(Q / "provenance/XVI-P/XVI-P-p23-source-only.zip") == \
        "3f4a88db462923553272e8a41407c0c06aacbc431736efdfeda10c27aae8a177"
    z.extractall(Q / "imports")
    (Q / "imports" / "XVI-P").rename(Q / "imports" / "XVI-P-p23")
    shutil.copytree(Q / "imports" / "XVI-P-p23", RUN)
    for d in ("certificates", "ledgers", "exports", "build/p14"):
        (RUN / d).mkdir(parents=True, exist_ok=True)


def classical():
    for s in STAGES:
        step(Path(s).stem, [PY, s], cwd=RUN)


def mutations():
    def one(m):
        d = Q / "build" / f"mut_{m}"
        shutil.copytree(RUN, d)
        return m, subprocess.run([PY, "code/classical/master_kappa2.py"], cwd=d, capture_output=True, text=True,
                                 env=dict(os.environ, P13_MUTATIONS=m)).returncode
    with ThreadPoolExecutor(7) as ex:
        res = list(ex.map(one, MUTS))
    print("[master_mutations]", res, flush=True)


def tex_escape_path(s):
    return s.replace("_", r"\_")


def manuscript():
    ms = Q / "manuscript"
    status = json.loads((Q / "RELEASE_STATUS.json").read_text())
    pre = (r"\documentclass[10pt]{article}\usepackage[margin=2cm]{geometry}\usepackage{amsmath,amssymb}"
           r"\usepackage{fancyhdr}\pagestyle{fancy}\fancyfoot[C]{XVI-Q Q04 \textbar{} page \thepage}"
           r"\usepackage{hyperref}\begin{document}")
    title = (r"\title{Landmark XVI-Q --- Q04 (P.6 entry, honest partial)\\\large One-loop $\kappa^3$ three-point sector: "
             r"inventory (15 + 3 scaleless tadpoles), $\kappa^3$ identity ledger, validated UV-pole engine, scalar-loop hhh "
             r"evaluated with the matter Ward identity verified at the pole; 13 graphs Missing; $\kappa^3$ CT/closure Missing. "
             r"O.1 CLOSED (scoped, $\kappa^2$ two-point) unchanged}"
             r"\date{2026-10-11}\maketitle\thispagestyle{fancy}\tableofcontents\newpage")
    sci = r"\input{notes/Q04_science.tex}\newpage\section*{Q03 layer}\addcontentsline{toc}{section}{Q03 layer}\input{notes/Q03_science.tex}\newpage\section*{Q02 layer (science; O.1 record superseded by Q03)}\addcontentsline{toc}{section}{Q02 layer}\input{notes/Q02_science.tex}\newpage\section*{Q01 layer (unchanged science)}\addcontentsline{toc}{section}{Q01 layer}\input{notes/Q01_science.tex}"
    app = [r"\section{Appendix: certificates, status, intake}"]
    for f in ["RELEASE_STATUS.json", "certificates/Q04_P6_status.json", "certificates/Q04_P6_inventory.json",
              "certificates/Q04_P6_identity.json", "certificates/Q04_P6_eval.json", "certificates/Q03_O1_record.json", "certificates/Q03_guard_ring.json",
              "ledgers/Q03_milestones.json", "certificates/Q02_O1.json", "certificates/Q02_bulk_identity.json",
              "certificates/Q02_bulk_graphs.json", "certificates/Q02_D12_reconciliation.json", "certificates/Q02_D04_scrutiny.json", "ledgers/Q02_milestones.json",
              "certificates/Q01_O2.json", "certificates/Q01_guard.json",
              "certificates/Q01_intake.json", "ledgers/Q01_milestones.json", "intake/errata_p23.json",
              "registration/local_registration.json"]:
        app += [r"\subsection*{" + tex_escape_path(f) + "}", r"{\scriptsize\begin{verbatim}",
                SS.wrap((Q / f).read_text()), r"\end{verbatim}}"]
    app += [r"\section{Appendix: source listings}"]
    for f in SS.source_files(Q):
        if f.endswith(".py"):
            app += [r"\subsection*{" + tex_escape_path(f) + "}", r"{\scriptsize\begin{verbatim}",
                    SS.wrap((Q / f).read_text()), r"\end{verbatim}}"]
    (ms / "body_appendix.tex").write_text("\n".join(app))
    docs = {"XVI-Q-science": sci, "XVI-Q-appendix": r"\input{body_appendix.tex}",
            "XVI-Q": sci + r"\newpage" + r"\input{body_appendix.tex}"}
    for name, body in docs.items():
        (ms / f"{name}.tex").write_text(pre + title + body + r"\end{document}")
        for _ in range(2):
            r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", f"{name}.tex"], cwd=ms,
                               capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout[-3000:])
            raise SystemExit(f"pdflatex {name} failed")
    for p in ms.glob("*.aux"):
        p.unlink()
    print("[pdfs]", {p: subprocess.run(["pdfinfo", f"manuscript/{p}.pdf"], cwd=Q, capture_output=True,
                                       text=True).stdout.split("Pages:")[1].split()[0] for p in PDFS}, flush=True)
    return status


def zips():
    dist = Q / "dist"
    SS.make_zip(dist / "XVI-Q-Q04-source-only.zip", Q)
    gen = ["RELEASE_STATUS.json"] + [p.relative_to(Q).as_posix() for p in sorted((Q / "certificates").glob("*.json"))] + \
          [p.relative_to(Q).as_posix() for p in sorted((Q / "ledgers").glob("*.json"))] + [f"manuscript/{p}.pdf" for p in PDFS]
    with zipfile.ZipFile(dist / "XVI-Q-review-focus.zip", "w", zipfile.ZIP_DEFLATED) as z:
        for f in gen + ["README.md", "intake/errata_p23.json", "registration/local_registration.json"]:
            z.write(Q / f, f"XVI-Q/{f}")
    with zipfile.ZipFile(dist / "XVI-Q-full.zip", "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(set(SS.source_files(Q)) | set(gen)):
            z.write(Q / f, f"XVI-Q/{f}")


def manifest():
    files = [f"manuscript/{p}.pdf" for p in PDFS] + ["dist/XVI-Q-Q04-source-only.zip", "dist/XVI-Q-review-focus.zip",
                                                     "dist/XVI-Q-full.zip", "RELEASE_STATUS.json"] + \
            [p.relative_to(Q).as_posix() for p in sorted((Q / "certificates").glob("Q0*.json"))]
    lines = ["# XVI-Q Q04 external manifest (written LAST; PDFs never contain these hashes)"]
    for f in files:
        pages = ""
        if f.endswith(".pdf"):
            pages = " pages=" + subprocess.run(["pdfinfo", f], cwd=Q, capture_output=True, text=True).stdout.split("Pages:")[1].split()[0]
        lines.append(f"{sh(Q / f)}  {f}{pages}")
    lines.append(f"source_tree_sha256 {SS.tree_sha(Q)}")
    (Q / "RELEASE_MANIFEST_EXTERNAL.txt").write_text("\n".join(lines) + "\n")
    for f in [f"manuscript/{p}.pdf" for p in PDFS]:            # no manifest hash inside any PDF
        txt = subprocess.run(["pdftotext", f, "-"], cwd=Q, capture_output=True, text=True).stdout
        assert not any(l.split()[0] in txt for l in lines[1:-1]), f


def sync():
    if SYNC == "none":
        return
    top = Path(SYNC)
    for p in PDFS:
        shutil.copy(Q / "manuscript" / f"{p}.pdf", top / f"{p}.pdf")
    for z in ("XVI-Q-Q04-source-only.zip", "XVI-Q-review-focus.zip", "XVI-Q-full.zip"):
        shutil.copy(Q / "dist" / z, top / z)
    print(f"[sync] {top}", flush=True)


def main():
    clean()
    extract()
    classical()
    mutations()
    step("certificates", [PY, "code/q01/o2_q01.py"])
    step("certificates_q02", [PY, "code/q02/q02_certs.py"])
    step("gate_q01_layer", [PY, "code/q01/release_q01.py"])
    step("certificates_q03", [PY, "code/q03/q03_certs.py"])
    step("gate_q02_layer", [PY, "code/q02/release_q02.py", "--layer"])
    step("certificates_q04", [PY, "code/q04/q04_certs.py"])
    step("gate_q03_layer", [PY, "code/q03/release_q03.py", "--layer"])
    step("gate", [PY, "code/q04/release_q04.py"])
    step("tests", [PY, "-m", "pytest", "-q", "-p", "no:cacheprovider", "tests"])
    manuscript()
    zips()
    manifest()
    sync()
    print("BUILD_OK", flush=True)


if __name__ == "__main__":
    main()

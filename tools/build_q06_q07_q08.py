#!/usr/bin/env python3
"""XVI-Q Q08 full build (from scratch).  Reuses the helpers of tools/build.py (pinned by the Q04 layer gate, unchanged):
clean -> extract -> classical -> mutations -> Q01..Q04 certificates -> Q01..Q03 layer gates -> Q04 gate (layer record)
-> Q05 certificates -> Q05 layer gate -> Q06 certificates -> Q06 layer gate -> Q07 certificates -> Q07 layer gate -> Q08 certificates -> Q08 gate (owns RELEASE_STATUS) -> tests -> PDFs -> zips -> manifest LAST -> optional sync
(default none; release opt-in: --sync <dir> or Q_SYNC_DIR)."""
from __future__ import annotations
import hashlib, json, os, shutil, subprocess, sys, zipfile  # noqa: E401
from pathlib import Path

Q = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Q / "tools"))
import build as B  # noqa: E402
import source_set_q06_q07_q08 as SS  # noqa: E402

PY = sys.executable
SYNC = (sys.argv[sys.argv.index("--sync") + 1] if "--sync" in sys.argv else os.environ.get("Q_SYNC_DIR", "none"))
PDFS = ["XVI-Q-science", "XVI-Q-appendix", "XVI-Q"]
ZIP = "XVI-Q-Q08-source-only.zip"
NOTES = ["Q08", "Q07", "Q06", "Q05", "Q04", "Q03", "Q02", "Q01"]
APPX = ["RELEASE_STATUS.json", "certificates/Q08_P6_status.json", "certificates/Q08_ghost_equation.json",
        "certificates/Q08_free_directions.json", "ledgers/Q08_milestones.json", "certificates/Q07_P6_status.json", "certificates/Q07_hhc_row.json",
        "certificates/Q07_evanescent_CT.json", "ledgers/Q07_milestones.json", "certificates/Q06_P6_status.json", "certificates/Q06_ghost_rows.json",
        "certificates/Q06_G10_crosscheck.json", "certificates/Q06_V4_validation.json", "certificates/Q06_fastpath.json",
        "certificates/Q06_P6_hhh.json", "ledgers/Q06_milestones.json", "certificates/Q05_P6_status.json", "certificates/Q05_P6_rows.json",
        "certificates/Q05_rules_validation.json", "certificates/Q05_P6_graphs.json", "certificates/Q05_guard_scan.json",
        "ledgers/Q05_milestones.json", "certificates/Q04_release.json", "certificates/Q04_P6_status.json",
        "certificates/Q04_P6_inventory.json", "certificates/Q04_P6_identity.json", "certificates/Q04_P6_eval.json",
        "certificates/Q03_O1_record.json", "certificates/Q02_bulk_identity.json", "certificates/Q01_O2.json",
        "intake/errata_p23.json", "registration/local_registration.json"]


def sh(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def manuscript():
    ms = Q / "manuscript"
    pre = (r"\documentclass[10pt]{article}\usepackage[margin=2cm]{geometry}\usepackage{amsmath,amssymb}"
           r"\usepackage{fancyhdr}\pagestyle{fancy}\fancyfoot[C]{XVI-Q Q08 \textbar{} page \thepage}"
           r"\usepackage{hyperref}\begin{document}")
    title = (r"\title{Landmark XVI-Q --- Q08 (P.6 evanescent counterterm structure, honest partial)\\\large "
             r"One-loop $\kappa^3$: the local evanescent counterterm X can be chosen ghost-equation compatible; "
             r"free directions 152 to 104, partly BRST-trivial (fitted families); operator basis, physical ambiguity "
             r"and p22 consistency open. O.1 CLOSED (scoped, $\kappa^2$ two-point) unchanged}"
             r"\date{2026-10-11}\maketitle\thispagestyle{fancy}\tableofcontents\newpage")
    parts = [r"\input{notes/Q08_science.tex}"]
    for n in NOTES[1:]:
        parts.append(r"\newpage\section*{" + n + r" layer}\addcontentsline{toc}{section}{" + n + r" layer}\input{notes/"
                     + n + r"_science.tex}")
    sci = "".join(parts)
    app = [r"\section{Appendix: certificates, status, intake}"]
    for f in APPX:
        app += [r"\subsection*{" + B.tex_escape_path(f) + "}", r"{\scriptsize\begin{verbatim}",
                SS.wrap((Q / f).read_text()), r"\end{verbatim}}"]
    app += [r"\section{Appendix: source listings}"]
    for f in SS.source_files(Q):
        if f.endswith(".py"):
            app += [r"\subsection*{" + B.tex_escape_path(f) + "}", r"{\scriptsize\begin{verbatim}",
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


def zips():
    dist = Q / "dist"
    for old in dist.glob("*.zip"):
        old.unlink()
    SS.make_zip(dist / ZIP, Q)
    gen = ["RELEASE_STATUS.json"] + [p.relative_to(Q).as_posix() for p in sorted((Q / "certificates").glob("*.json"))] + \
          [p.relative_to(Q).as_posix() for p in sorted((Q / "ledgers").glob("*.json"))] + [f"manuscript/{p}.pdf" for p in PDFS]
    with zipfile.ZipFile(dist / "XVI-Q-review-focus.zip", "w", zipfile.ZIP_DEFLATED) as z:
        for f in gen + ["README.md", "intake/errata_p23.json", "registration/local_registration.json"]:
            z.write(Q / f, f"XVI-Q/{f}")
    with zipfile.ZipFile(dist / "XVI-Q-full.zip", "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(set(SS.source_files(Q)) | set(gen)):
            z.write(Q / f, f"XVI-Q/{f}")


def manifest():
    files = [f"manuscript/{p}.pdf" for p in PDFS] + [f"dist/{ZIP}", "dist/XVI-Q-review-focus.zip", "dist/XVI-Q-full.zip",
                                                     "RELEASE_STATUS.json"] + \
            [p.relative_to(Q).as_posix() for p in sorted((Q / "certificates").glob("Q0*.json"))]
    lines = ["# XVI-Q Q08 external manifest (written LAST; PDFs never contain these hashes)"]
    for f in files:
        pages = ""
        if f.endswith(".pdf"):
            pages = " pages=" + subprocess.run(["pdfinfo", f], cwd=Q, capture_output=True, text=True).stdout.split("Pages:")[1].split()[0]
        lines.append(f"{sh(Q / f)}  {f}{pages}")
    lines.append(f"source_tree_sha256 {SS.tree_sha(Q)}")
    (Q / "RELEASE_MANIFEST_EXTERNAL.txt").write_text("\n".join(lines) + "\n")
    for f in [f"manuscript/{p}.pdf" for p in PDFS]:
        txt = subprocess.run(["pdftotext", f, "-"], cwd=Q, capture_output=True, text=True).stdout
        assert not any(ln.split()[0] in txt for ln in lines[1:-1]), f


def sync():
    if SYNC == "none":
        return
    top = Path(SYNC)
    for p in PDFS:
        shutil.copy(Q / "manuscript" / f"{p}.pdf", top / f"{p}.pdf")
    for z in (ZIP, "XVI-Q-review-focus.zip", "XVI-Q-full.zip"):
        shutil.copy(Q / "dist" / z, top / z)
    print(f"[sync] {top}", flush=True)


def main():
    B.clean()
    B.extract()
    B.classical()
    B.mutations()
    st = B.step
    st("certificates", [PY, "code/q01/o2_q01.py"])
    st("certificates_q02", [PY, "code/q02/q02_certs.py"])
    st("gate_q01_layer", [PY, "code/q01/release_q01.py"])
    st("certificates_q03", [PY, "code/q03/q03_certs.py"])
    st("gate_q02_layer", [PY, "code/q02/release_q02.py", "--layer"])
    st("certificates_q04", [PY, "code/q04/q04_certs.py"])
    st("gate_q03_layer", [PY, "code/q03/release_q03.py", "--layer"])
    st("gate_q04_layer", [PY, "code/q04/release_q04.py"])
    st("certificates_q05", [PY, "code/p6q05/q05_certs.py"])
    st("gate_q05_layer", [PY, "code/p6q05/release_q05.py"])
    st("certificates_q06", [PY, "code/p6q06/q06_certs.py"])
    st("gate_q06_layer", [PY, "code/p6q06/release_q06.py"])
    st("certificates_q07", [PY, "code/p6q06/q07/q07_certs.py"])
    st("gate_q07_layer", [PY, "code/p6q06/q07/release_q07.py"])
    st("certificates_q08", [PY, "code/p6q06/q07/q08/q08_certs.py"])
    st("gate", [PY, "code/p6q06/q07/q08/release_q08.py"])
    st("tests", [PY, "-m", "pytest", "-q", "-p", "no:cacheprovider", "tests", "code/p6q06/q07/test_q07.py", "code/p6q06/q07/q08/test_q08.py"])
    manuscript()
    zips()
    manifest()
    sync()
    print("BUILD_OK", flush=True)


if __name__ == "__main__":
    main()

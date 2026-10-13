# XVI-Q Q01 rebuild instructions
1. Unpack XVI-Q-Q04-source-only.zip (source only; the XVI-P p23 zip is included under provenance/XVI-P and is
   re-extracted and sha-checked by the build).
2. Run `python3 tools/rebuild_from_source.py --report <file>`. It runs `tools/build.py` from scratch in a temporary
   copy and compares every certificate (content sha with `seconds` removed) and the ℓ₂..ℓ₄ export shas with
   reproducibility/reference_values_q04.json.
3. Environment: python3, sympy, pytest, pdflatex, poppler-utils.

Q04: tools/build.py does NOT sync into /workspace unless `--sync <dir>` (or Q_SYNC_DIR) is given.

## Q06
From the Q06 source-only zip: `python3 tools/rebuild_q06.py` builds a copy of the source tree in a temporary directory
(tools/build_q06.py, no sync) and compares it with reproducibility/reference_values_q06.json.

## Q05
From the Q05 source-only zip: `python3 tools/rebuild_q05.py` builds a copy of the source tree in a temporary directory
(tools/build_q05.py, no sync) and compares it with reproducibility/reference_values_q05.json.

## Q07
python3 tools/build_q06_q07.py            (no sync by default; --sync <dir> opt-in)
python3 tools/rebuild_q06_q07.py          (rebuild from the source-only zip and compare with reproducibility/reference_values_q07.json)

## Q08
python3 tools/build_q06_q07_q08.py   (no sync by default)
python3 tools/rebuild_q06_q07_q08.py

## Q09
python3 tools/build_q06_q07_q08_q09.py   (no sync by default)
python3 tools/rebuild_q06_q07_q08_q09.py

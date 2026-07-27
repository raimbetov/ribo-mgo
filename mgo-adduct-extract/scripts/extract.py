#!/usr/bin/env python3
"""Spec-driven extractor for proteomically-detected methylglyoxal adduct tables.

Every dataset is described by one YAML file in ../datasets/. This script contains
no per-paper knowledge: the page ranges, column geometry, expected counts and
acceptance checks all live in the spec, so the spec is the reviewable artifact and
this file is just the mechanism.

    python3 extract.py ../datasets/irshad-2019-hmec1.yaml --repo-root /path/to/repo
    python3 extract.py --all --repo-root /path/to/repo

A dataset is a list of `sources`, each parsed by one `model`, concatenated into a
single output TSV. Nothing is deduplicated: a protein listed in two source tables
appears twice, distinguishable by its provenance columns.

Row models
----------
xlsx            one worksheet; header row located by a marker cell
pdf_anchor      PDF table whose rows are anchored on a number in the left margin,
                with columns cut by x-position. `assign: nearest` attaches a word
                to the closest anchor (needed when a row's number is printed
                *below* its name); `assign: below` attaches it to the nearest
                anchor at or above (the usual case, and the safe one when rows
                span several baselines).
pdf_text_rows   regex over `pdftotext -layout` output, for pages whose embedded
                text has no usable word spacing.
literal         rows transcribed from prose or from a figure, carried in the spec.
label_union     union of several named label lists (figure panels), emitting
                membership columns.

DESIGN RULE: nothing may be dropped silently. Every code path that discards input
-- a blank spreadsheet row, a repeated anchor number, a word that falls in a gap
between declared columns -- counts what it discarded and fails the build unless
the spec explicitly acknowledges it. A check naming a column that does not exist
is an error, not a pass.

Exit status is non-zero if any check fails; no file is written in that case.
"""
import argparse
import csv
import os
import re
import subprocess
import sys

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required: pip install -r ../requirements.txt")

PROVENANCE = ["source_file", "source_table", "source_row"]


# ---------------------------------------------------------------- helpers

def norm_ws(s):
    return re.sub(r"\s+", " ", s or "").strip()


def in_x(word, rng):
    return rng[0] <= word["x0"] < rng[1]


def as_list(value, field):
    """YAML lets a single string stand where a list is meant; that would silently
    become a per-character match. Normalise, and reject anything else."""
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, (list, tuple)):
        return list(value)
    raise SpecError(f"{field} must be a string or list, got {type(value).__name__}")


class SpecError(Exception):
    pass


def load_pdfplumber():
    try:
        import pdfplumber
    except ImportError:
        sys.exit("pdfplumber is required: pip install -r ../requirements.txt")
    return pdfplumber


# ---------------------------------------------------------------- row models
# Each model returns (header_or_None, records, diagnostics).
# diagnostics is a list of (severity, message) where severity is 'note' or 'problem'.

def model_xlsx(src, root):
    """One worksheet. Header row is the first whose first cell matches `header_marker`."""
    try:
        import openpyxl
    except ImportError:
        sys.exit("openpyxl is required: pip install -r ../requirements.txt")

    path = os.path.join(root, src["file"])
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[src["sheet"]]
    rows = list(ws.iter_rows(values_only=True))
    marker = src.get("header_marker", "Accession")
    try:
        hidx = next(i for i, r in enumerate(rows)
                    if r and str(r[0]).strip() == marker)
    except StopIteration:
        raise SpecError(f"no header row whose first cell is {marker!r} "
                        f"in sheet {src['sheet']!r}")
    header = [str(c).strip() if c is not None else "" for c in rows[hidx]]

    diag = []
    blank_hdr = [i for i, h in enumerate(header) if not h]
    dupe_hdr = sorted({h for h in header if h and header.count(h) > 1})
    if blank_hdr:
        diag.append(("problem", f"blank header cell(s) at index {blank_hdr}"))
    if dupe_hdr:
        diag.append(("problem", f"duplicate header name(s) {dupe_hdr} -- "
                                f"columns would silently collapse"))

    out, skipped = [], []
    for i, r in enumerate(rows[hidx + 1:], start=hidx + 2):
        if not r or all(c in (None, "") for c in r[:len(header)]):
            continue                       # wholly empty row: genuine padding
        if r[0] in (None, ""):
            skipped.append(i)              # has content but no key: never silent
            continue
        cells = ["" if c is None else str(c).strip() for c in r[:len(header)]]
        cells += [""] * (len(header) - len(cells))
        rec = dict(zip(header, cells))
        rec["__row__"] = str(i)
        out.append(rec)

    wb.close()
    allowed = src.get("allow_rows_without_key", 0)
    if skipped:
        sev = "note" if len(skipped) == allowed else "problem"
        diag.append((sev, f"{len(skipped)} row(s) have content but an empty "
                          f"first column and were not emitted: {skipped[:10]}"))
    return header, out, diag


def model_pdf_anchor(src, root):
    """PDF table anchored on a row number in the left margin."""
    pdfplumber = load_pdfplumber()
    path = os.path.join(root, src["file"])
    cols = [(c["name"], tuple(c["x"])) for c in src["columns"]]
    anchor_x = tuple(src["anchor_x"])
    lo, hi = src["pages"]
    band = float(src.get("band", 3.0))
    y_tol = float(src["y_tol"])
    same_line = float(src.get("same_line_tol", band))
    assign = src.get("assign", "below")
    max_no = int(src["anchor_max"])
    drop_words = as_list(src.get("drop_baselines_containing"),
                         "drop_baselines_containing")
    drop_numeric = src.get("drop_numeric_baselines", None)
    once = src.get("anchor_once_per_table", False)

    out, seen, diag = [], set(), []
    suppressed, unassigned = [], []

    with pdfplumber.open(path) as pdf:
        if hi >= len(pdf.pages):
            diag.append(("problem", f"pages [{lo},{hi}] exceeds document "
                                    f"({len(pdf.pages)} pages)"))
        for pno in range(lo, min(hi, len(pdf.pages) - 1) + 1):
            page = pdf.pages[pno]
            words = page.extract_words(use_text_flow=False, keep_blank_chars=False)

            bands = {}
            for w in words:
                bands.setdefault(round(w["top"] / band), []).append(w)
            skip = set()
            if drop_words:
                skip |= {k for k, ws in bands.items()
                         if any(any(d in w["text"] for d in drop_words) for w in ws)}
            if drop_numeric is not None:
                skip |= {k for k, ws in bands.items()
                         if all(re.fullmatch(r"[\d,()%]+", w["text"])
                                and w["x0"] >= drop_numeric for w in ws)}

            anchors, taken = [], set()
            for w in sorted(words, key=lambda w: w["top"]):
                if not in_x(w, anchor_x) or round(w["top"] / band) in skip:
                    continue
                if not re.fullmatch(r"\d{1,4}", w["text"]):
                    continue
                n = int(w["text"])
                if not 1 <= n <= max_no:
                    continue
                if n in taken or (once and n in seen):
                    suppressed.append((pno + 1, n))     # never silent
                    continue
                taken.add(n)
                anchors.append(w)
            seen |= taken
            if not anchors:
                continue
            anchors.sort(key=lambda w: w["top"])

            buckets = {id(a): [] for a in anchors}
            xmin, xmax = cols[0][1][0], cols[-1][1][1]
            for w in words:
                if in_x(w, anchor_x) or round(w["top"] / band) in skip:
                    continue
                if w["x0"] < xmin or w["x0"] >= xmax:
                    continue
                if assign == "nearest":
                    a = min(anchors, key=lambda a: abs(a["top"] - w["top"]))
                    if abs(a["top"] - w["top"]) > y_tol:
                        continue
                else:
                    above = [a for a in anchors if a["top"] <= w["top"] + same_line]
                    if not above:
                        continue
                    a = above[-1]
                    if w["top"] - a["top"] > y_tol:
                        continue
                # A word inside the table's x-span that matches no declared column
                # falls in a gap between column ranges and would vanish.
                if not any(in_x(w, rng) for _, rng in cols):
                    unassigned.append((pno + 1, round(w["x0"]), w["text"][:20]))
                    continue
                buckets[id(a)].append(w)

            for a in anchors:
                ws = sorted(buckets[id(a)],
                            key=lambda w: (round(w["top"] / band), w["x0"]))
                rec = {src.get("anchor_name", "No"): a["text"]}
                for name, rng in cols:
                    rec[name] = norm_ws(" ".join(w["text"] for w in ws
                                                 if in_x(w, rng))).rstrip(",")
                rec["__row__"] = f"p{pno + 1} row {a['text']}"
                rec["__no__"] = int(a["text"])
                out.append(rec)

    if suppressed:
        diag.append(("problem",
                     f"{len(suppressed)} repeated anchor number(s) suppressed "
                     f"(page, number): {suppressed[:8]} -- check the page range"))
    if unassigned:
        allowed = src.get("allow_unassigned_words", 0)
        sev = "note" if len(unassigned) <= allowed else "problem"
        diag.append((sev,
                     f"{len(unassigned)} word(s) fell in a gap between declared "
                     f"columns and were discarded (page, x, text): "
                     f"{unassigned[:8]}"))
    return None, out, diag


def model_pdf_text_rows(src, root):
    """Regex over `pdftotext -layout` output."""
    path = os.path.join(root, src["file"])
    lo, hi = src["pages"]                      # 1-based, as pdftotext expects
    diag = []
    out, seen = [], set()
    for pno in range(lo, hi + 1):
        txt = subprocess.run(["pdftotext", "-layout", "-f", str(pno), "-l",
                              str(pno), path, "-"],
                             capture_output=True, text=True, check=True).stdout
        rx = re.compile(src["regex"])
        groups = src["groups"]
        max_no = int(src.get("max_no", 10 ** 6))
        for line in txt.splitlines():
            m = rx.match(line)
            if not m:
                continue
            rec = {g: m.group(i + 1).strip() for i, g in enumerate(groups)}
            n = int(rec[groups[0]])
            if not 1 <= n <= max_no:
                continue
            if n in seen:
                diag.append(("problem", f"row number {n} matched more than once"))
                continue
            seen.add(n)
            rec["__row__"] = src["row_ref"].format(no=n, page=pno)
            rec["__no__"] = n
            out.append(rec)
    return None, out, diag


def model_literal(src, root):
    """Rows transcribed from prose or a figure, carried verbatim in the spec."""
    names = src["columns"]
    out = []
    for i, row in enumerate(src["rows"], start=1):
        if len(row) != len(names):
            raise SpecError(f"literal row {i} has {len(row)} values "
                            f"but {len(names)} columns are declared")
        rec = dict(zip(names, [str(v) for v in row]))
        rec["__row__"] = src.get("row_ref", "prose")
        rec["__no__"] = i
        out.append(rec)
    return None, out, []


def model_label_union(src, root):
    """Union of named label lists (figure panels), with membership columns."""
    panels = src["panels"]
    diag = []
    for p in panels:
        labels = p["labels"]
        exp = p.get("expected_count")
        dupes = sorted({x for x in labels if labels.count(x) > 1})
        bad = (exp is not None and len(labels) != exp) or dupes
        if exp is not None and len(labels) != exp:
            diag.append(("problem", f"panel {p['key']}: {len(labels)} labels, "
                                    f"expected {exp}"))
        if dupes:
            diag.append(("problem", f"panel {p['key']}: duplicate labels {dupes}"))
        if exp is not None and not bad:
            diag.append(("note", f"panel {p['key']}: {len(labels)} labels "
                                 f"(expected {exp}), no duplicates"))

    union = sorted({g for p in panels for g in p["labels"]})
    out = []
    for i, g in enumerate(union, start=1):
        member = [p for p in panels if g in p["labels"]]
        rec = {
            src["label_column"]: g,
            src["types_column"]: "; ".join(p["type"] for p in member),
            src["panels_column"]: "".join(p["key"] for p in member),
        }
        rec["__row__"] = src.get("row_ref", "")
        rec["__no__"] = i
        out.append(rec)
    return None, out, diag


MODELS = {
    "xlsx": model_xlsx,
    "pdf_anchor": model_pdf_anchor,
    "pdf_text_rows": model_pdf_text_rows,
    "literal": model_literal,
    "label_union": model_label_union,
}


# ---------------------------------------------------------------- checks

def tokens(value, pattern):
    return re.findall(pattern, value or "")


def need_column(col, cols, kind):
    if col not in cols:
        raise SpecError(f"check '{kind}' names column {col!r}, which this "
                        f"dataset does not produce. Columns are: {sorted(cols)}")


def run_checks(spec, header, cols, rows, report, registry):
    """Returns (problems, columns_to_drop)."""
    problems, drop = [], []
    for chk in spec.get("checks", []):
        kind = chk["type"]

        if kind == "contiguous_numbering":
            nums = [r["__no__"] for r in rows if r.get("__no__") is not None]
            if len(nums) != len(rows):
                raise SpecError("check 'contiguous_numbering' requires every "
                                "source to number its rows; this dataset has "
                                f"{len(rows) - len(nums)} unnumbered row(s)")
            exp = list(range(1, len(nums) + 1))
            got = sorted(nums)
            if got != exp:
                missing = sorted(set(exp) - set(got))
                extra = sorted(n for n in set(got) if got.count(n) > 1)
                problems.append(f"numbering not contiguous 1..{len(nums)}; "
                                f"missing {missing[:10]} repeated {extra[:10]}")
            report.append(f"  numbering contiguous 1..{len(nums)}: "
                          f"{'yes' if got == exp else 'NO'}")

        elif kind == "column_sum":
            need_column(chk["column"], cols, kind)
            col = chk["column"]
            bad = [r["__no__"] for r in rows if not str(r.get(col, "")).isdigit()]
            if bad:
                problems.append(f"non-numeric '{col}' in rows {bad[:10]}")
            total = sum(int(r[col]) for r in rows if str(r.get(col, "")).isdigit())
            if total != chk["equals"]:
                problems.append(f"sum of '{col}' = {total}, expected {chk['equals']}")
            report.append(f"  sum of '{col}' = {total} "
                          f"(expected {chk['equals']}) -- {chk.get('evidence','')}")

        elif kind == "token_sum":
            need_column(chk["column"], cols, kind)
            col, pat = chk["column"], chk["pattern"]
            total = sum(len(tokens(r.get(col, ""), pat)) for r in rows)
            if total != chk["equals"]:
                problems.append(f"token sum of '{col}' = {total}, "
                                f"expected {chk['equals']}")
            report.append(f"  token sum of '{col}' matching /{pat}/ = {total} "
                          f"(expected {chk['equals']}) -- {chk.get('evidence','')}")

        elif kind == "declared_vs_listed":
            need_column(chk["count_column"], cols, kind)
            need_column(chk["token_column"], cols, kind)
            col, tcol, pat = chk["count_column"], chk["token_column"], chk["pattern"]
            off = [r["__no__"] for r in rows
                   if str(r.get(col, "")).isdigit()
                   and len(tokens(r.get(tcol, ""), pat)) != int(r[col])]
            # Known publication typos are pinned by row number, so a NEW mismatch
            # -- or the disappearance of a known one -- fails the build.
            expected = chk.get("expect_offenders")
            if expected is None:
                if off:
                    problems.append(f"declared/listed mismatch in rows {off}")
            elif sorted(off) != sorted(expected):
                problems.append(f"declared/listed offenders {off} != "
                                f"pinned {sorted(expected)}")
            report.append(f"  rows whose declared count != printed list: {off}"
                          + (f" (pinned, {chk.get('reason','')})" if expected else ""))

        elif kind == "subset":
            need_column(chk["subset_column"], cols, kind)
            need_column(chk["superset_column"], cols, kind)
            sub, sup, pat = chk["subset_column"], chk["superset_column"], chk["pattern"]
            bad = [r["__no__"] for r in rows
                   if not set(tokens(r.get(sub, ""), pat))
                   <= set(tokens(r.get(sup, ""), pat))]
            report.append(f"  '{sub}' subset of '{sup}': "
                          f"{len(bad)} violation(s) {bad[:10]}")
            if bad:
                if chk.get("on_fail") == "drop_column":
                    drop.append(sub)
                    report.append(f"  -> dropping '{sub}' ({chk.get('reason','')})")
                else:
                    problems.append(f"'{sub}' not a subset of '{sup}' in rows {bad}")

        elif kind == "regex_column":
            need_column(chk["column"], cols, kind)
            col, pat = chk["column"], chk["pattern"]
            bad = [r["__no__"] for r in rows
                   if not re.fullmatch(pat, r.get(col, ""))]
            if bad:
                problems.append(f"'{col}' fails /{pat}/ in rows {bad[:10]}")
            report.append(f"  '{col}' matches /{pat}/: {len(rows) - len(bad)}/{len(rows)}")

        elif kind == "unique_column":
            need_column(chk["column"], cols, kind)
            col = chk["column"]
            vals = [r.get(col, "") for r in rows]
            dupes = sorted({v for v in vals if vals.count(v) > 1})
            expected = chk.get("expect_duplicates")
            if expected is None:
                if dupes:
                    problems.append(f"duplicate '{col}': {dupes[:10]}")
            elif sorted(dupes) != sorted(expected):
                problems.append(f"duplicates in '{col}' {dupes[:10]} != "
                                f"pinned {sorted(expected)[:10]}")
            report.append(f"  '{col}' unique: {len(set(vals))}/{len(vals)}"
                          + (f" duplicates {dupes[:6]}" if dupes else ""))

        elif kind == "no_blank":
            allbad = {}
            for col in chk["columns"]:
                need_column(col, cols, kind)
                bad = [r["__no__"] for r in rows if not r.get(col, "").strip()]
                if bad:
                    allbad[col] = bad
                    problems.append(f"blank '{col}' in rows {bad[:10]}")
            report.append(f"  no blanks in {chk['columns']}: "
                          + ("ok" if not allbad else f"FAILED {allbad}"))

        elif kind == "blank_count":
            # Some columns are legitimately blank for a known subset of rows
            # (a source table the paper gives no gene symbol for). Pin the
            # count so a parse that starts losing values is caught.
            need_column(chk["column"], cols, kind)
            col = chk["column"]
            blanks = [r["__no__"] for r in rows if not r.get(col, "").strip()]
            if len(blanks) != chk["equals"]:
                problems.append(f"{len(blanks)} blank '{col}' cells, "
                                f"expected {chk['equals']}")
            report.append(f"  blank '{col}' cells: {len(blanks)} "
                          f"(expected {chk['equals']}) -- {chk.get('reason','')}")

        elif kind == "any_nonempty_matching":
            pat = re.compile(chk["column_pattern"])
            matching = [c for c in cols if pat.search(c)]
            if not matching:
                raise SpecError(f"check 'any_nonempty_matching': no column "
                                f"matches /{chk['column_pattern']}/")
            bad = [i for i, r in enumerate(rows, 1)
                   if not any(r.get(c, "").strip() for c in matching)]
            if bad:
                problems.append(f"{len(bad)} row(s) with no adduct in {matching}")
            report.append(f"  every row has an adduct in {matching}: "
                          f"{len(rows) - len(bad)}/{len(rows)}")

        elif kind == "cross_method_identifiers":
            problems += cross_method_check(spec, chk, rows, cols, report)

        elif kind == "covers_all_of":
            # Every identifier in a dependency dataset must appear here. This is
            # the external reconciliation for figure-derived data.
            need_column(chk["column"], cols, kind)
            dep = registry.get(chk["from_file"])
            if dep is None:
                raise SpecError(f"check 'covers_all_of' needs {chk['from_file']!r} "
                                f"to be built in the same run")
            want = {r[chk["from_column"]] for r in dep}
            have = {r.get(chk["column"], "") for r in rows}
            missing = sorted(want - have)
            if missing:
                problems.append(f"{len(missing)} value(s) from "
                                f"{chk['from_file']} missing here: {missing[:10]}")
            report.append(f"  covers all {len(want)} '{chk['from_column']}' "
                          f"from {chk['from_file']}: "
                          f"{len(want) - len(missing)}/{len(want)}")

        else:
            raise SpecError(f"unknown check type {kind!r}")
    return problems, drop


def cross_method_check(spec, chk, rows, cols, report):
    """Re-parse the source with pdftotext and compare the identifier sets.

    Catches identifiers invented by coordinate assignment and whole rows dropped
    by it, without depending on row alignment between the two methods. Compared
    only against rows from the source being re-parsed.
    """
    need_column(chk["column"], cols, "cross_method_identifiers")
    root = spec["__root__"]
    idx = chk.get("source_index", 0)
    src = spec["sources"][idx]
    table = src.get("table", "")
    path = os.path.join(root, src["file"])
    lo, hi = src["pages"]
    if src["model"] == "pdf_anchor":
        lo, hi = lo + 1, hi + 1                 # 0-based spec -> 1-based tool
    txt = subprocess.run(["pdftotext", "-layout", "-f", str(lo), "-l", str(hi),
                          path, "-"],
                         capture_output=True, text=True, check=True).stdout
    pat = re.compile(chk["pattern"])
    found = {t for t in txt.split() if pat.fullmatch(t)}

    # Only rows that came from this source may be compared.
    scoped = [r for r in rows if r.get("source_table") == table
              and r.get("source_file") == src["file"]]
    if not scoped:
        raise SpecError("check 'cross_method_identifiers': no rows came from "
                        f"source_index {idx} (table {table!r})")
    mine = {r.get(chk["column"], "") for r in scoped}

    excuse = set()
    for col in chk.get("also_accounted_by", []):
        need_column(col, cols, "cross_method_identifiers.also_accounted_by")
        excuse |= {r.get(col, "") for r in scoped}

    invented = sorted(mine - found)
    unexplained = sorted(found - mine - excuse)
    report.append(f"  cross-method on {table}: {len(scoped)} rows, "
                  f"{len(invented)} invented, {len(unexplained)} unexplained")
    problems = []
    if invented:
        problems.append(f"identifiers absent from check parse: {invented[:8]}")
    if unexplained:
        problems.append(f"unaccounted identifiers in check parse: {unexplained[:8]}")
    return problems


# ---------------------------------------------------------------- assembly

def build(spec, root, registry):
    spec["__root__"] = root
    out_cols, rows, report, problems = [], [], [], []

    for si, src in enumerate(spec["sources"]):
        if src["model"] not in MODELS:
            raise SpecError(f"unknown model {src['model']!r}")
        header, recs, diag = MODELS[src["model"]](src, root)

        for sev, msg in diag:
            report.append(f"  [{src.get('table', src['model']).split(' (')[0]}] {msg}")
            if sev == "problem":
                problems.append(msg)

        if src.get("emit") == "passthrough":
            emit = {c: {"column": c} for c in header}
        else:
            emit = src["emit"]
            bad_keys = [k for k in emit if not isinstance(k, str)]
            if bad_keys:
                raise SpecError(f"source {si}: column names must be strings; "
                                f"{bad_keys!r} were parsed as another YAML type "
                                f"-- quote them")
        # The output schema is the union over all sources, in first-seen order,
        # so a column only one source emits is not lost.
        for c in emit:
            if c not in out_cols:
                out_cols.append(c)

        table = src.get("table", "")
        sfile = src.get("provenance_file", src["file"])
        for rec in recs:
            row = {c: "" for c in out_cols}
            for col, rule in emit.items():
                if "column" in rule:
                    row[col] = rec.get(rule["column"], "")
                elif "const" in rule:
                    row[col] = str(rule["const"])
            row["source_file"] = sfile
            row["source_table"] = table.format(**rec) if "{" in table else table
            row["source_row"] = rec.get("__row__", "")
            row["__no__"] = rec.get("__no__")
            rows.append(row)

    # Backfill any column introduced by a later source onto earlier rows.
    for r in rows:
        for c in out_cols:
            r.setdefault(c, "")

    for lk in spec.get("lookups", []):
        for key in ("into", "into_flag"):
            if lk.get(key) and lk[key] not in out_cols:
                raise SpecError(f"lookup target {lk[key]!r} is not in any "
                                f"source's emit, so it would be discarded")
    return out_cols, rows, problems, report


def resolve_lookups(spec, rows, registry):
    """Fill a column from another dataset built in THIS run.

    Reading a previously written file would let a stale or failed dependency
    silently supply values, so the dependency must be present in the run registry.
    """
    for lk in spec.get("lookups", []):
        dep = registry.get(lk["from_file"])
        if dep is None:
            raise SpecError(f"lookup source {lk['from_file']!r} was not built in "
                            f"this run -- build it first, or use --all")
        table = {r[lk["match_on"]]: r[lk["take"]] for r in dep}
        for r in rows:
            key = r.get(lk["key_column"], "")
            if lk.get("into_flag"):
                r[lk["into_flag"]] = "yes" if key in table else "no"
            if lk.get("into"):
                r[lk["into"]] = table.get(key, "")
    return rows


def write_tsv(path, cols, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "") for c in cols])


def process(path, root, outdir, registry):
    with open(path, encoding="utf-8") as fh:
        spec = yaml.safe_load(fh)

    name = spec["output"]
    print(f"\n{spec['id']}  ->  {name}")
    try:
        out_cols, rows, problems, report = build(spec, root, registry)
        rows = resolve_lookups(spec, rows, registry)
        cols = out_cols + PROVENANCE
        p2, drop = run_checks(spec, None, set(cols), rows, report, registry)
        problems += p2
        cols = [c for c in out_cols if c not in drop] + PROVENANCE
    except SpecError as e:
        print(f"  FAIL  spec error: {e}")
        return 1

    exp = spec.get("expected_rows")
    if exp is not None and len(rows) != exp:
        problems.append(f"produced {len(rows)} rows, expected {exp}")
    print(f"  rows = {len(rows)}" + (f" (expected {exp})" if exp else ""))
    for line in report:
        print(line)

    if problems:
        print("  FAIL  " + "; ".join(problems))
        return 1

    os.makedirs(outdir, exist_ok=True)
    write_tsv(os.path.join(outdir, name), cols, rows)
    registry[name] = rows
    print(f"  PASS  wrote {name} ({len(rows)} rows)")
    for note in spec.get("caveats", []):
        print(f"  NOTE  {note}")
    return 0


def order_specs(paths):
    """Topologically order so a dataset's lookup dependencies build first."""
    meta = {}
    for p in paths:
        with open(p, encoding="utf-8") as fh:
            spec = yaml.safe_load(fh)
        meta[p] = (spec["output"],
                   {lk["from_file"] for lk in spec.get("lookups", [])})
    produced = {out: p for p, (out, _) in meta.items()}

    ordered, done, guard = [], set(), 0
    while len(ordered) < len(paths):
        progressed = False
        for p in paths:
            if p in done:
                continue
            _, deps = meta[p]
            if all(d not in produced or produced[d] in done for d in deps):
                ordered.append(p)
                done.add(p)
                progressed = True
        guard += 1
        if not progressed or guard > len(paths) + 1:
            remaining = [p for p in paths if p not in done]
            sys.exit(f"circular or unresolvable lookup dependency among "
                     f"{[os.path.basename(p) for p in remaining]}")
    return ordered


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("spec", nargs="*", help="dataset spec YAML file(s)")
    ap.add_argument("--all", action="store_true", help="build every spec")
    ap.add_argument("--repo-root",
                    default=os.path.dirname(os.path.dirname(here)),
                    help="repository root that source paths are relative to "
                         "(default: the directory containing this skill)")
    ap.add_argument("--out", default=None, help="output directory")
    args = ap.parse_args()

    root = os.path.abspath(args.repo_root)
    outdir = os.path.abspath(args.out or os.path.join(root, "01-data-extracted"))
    dsdir = os.path.join(here, os.pardir, "datasets")

    specs = sorted(os.path.join(dsdir, f) for f in os.listdir(dsdir)
                   if f.endswith(".yaml")) if args.all else args.spec
    if not specs:
        ap.error("give one or more spec files, or --all")

    registry, rc = {}, 0
    for s in order_specs(specs):
        rc |= process(s, root, outdir, registry)
    print()
    print("all specs built" if rc == 0 else "BUILD FAILED")
    return rc


if __name__ == "__main__":
    sys.exit(main())

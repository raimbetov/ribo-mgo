#!/usr/bin/env python3
"""Turn 01-data-extracted/ TSVs into uniprot-name-match input files (stage 01 -> 02).

Each paper identifies proteins differently, so each dataset needs its own rule.
This script applies them and writes one plain-text file per dataset, in the
format uniprot-name-match expects (one clue-set per line, TAB-separated fields).

    python3 make_uniprot_input.py                 # write to 02-data-processed/uniprot-input/
    python3 make_uniprot_input.py --out DIR
    python3 make_uniprot_input.py --print-plan    # show the rules without writing

Why each rule is what it is (all measured against live UniProt):

* ALWAYS PAIR AN IDENTIFIER WITH A NAME where the source has one. An accession
  alone has nothing to fall back on if it is ever retired; `accession<TAB>name`
  degrades gracefully to the name. Donnellan's `Description` column carries the
  full protein name, so it costs nothing to include.

* STRIP ISOFORM SUFFIXES. 53 Donnellan accessions look like `P04075-2`. UniProt
  does not resolve the isoform suffix, so the row scores `low` even though the
  canonical entry is found. The suffix is kept in the `From` column via the
  paired name, and recorded in the isoform report.

* SPLIT IRSHAD'S PARENTHETICAL. Its `Protein` column reads `Pyruvate kinase-M
  (PKM)`. Fed raw, the hint is invisible to the parser and the sample mapped
  2/20; split into `Pyruvate kinase-M<TAB>PKM` it mapped 19/20. Note the hints
  are a mix of gene symbols and UniProt mnemonics (`PROF1` for PFN1, `AATM`,
  `G3P`), which still resolve, sometimes at `medium`.

* BARE GENE SYMBOLS ARE FINE NOW. Zheng's Fig. 4 rows without an accession are
  emitted as the lone symbol; uniprot-name-match treats a single-field line that
  looks like a gene symbol as a gene symbol.
"""
import argparse
import csv
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, "01-data-extracted")

ISOFORM_RE = re.compile(r"^([A-Z0-9]+)-\d+$")
# "Pyruvate kinase-M (PKM)" -> ("Pyruvate kinase-M", "PKM")
HINT_RE = re.compile(r"^(.*?)\s*\(([A-Za-z0-9_]+)\)\s*$")
# Donnellan Description: "... OS=Homo sapiens OX=9606 GN=HSPE1 PE=1 SV=2"
DESC_TAIL_RE = re.compile(r"\s+OS=.*$")


def clean_desc(desc):
    """Drop the OS=/OX=/GN=/PE=/SV= tail from a Proteome Discoverer description."""
    return DESC_TAIL_RE.sub("", desc or "").strip()


def strip_isoform(acc):
    m = ISOFORM_RE.match(acc or "")
    return (m.group(1), True) if m else (acc, False)


def read(name):
    with open(os.path.join(SRC, name), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def line(*fields):
    """One clue-set. Empty fields are dropped so a missing name never leaves a
    dangling tab, which the parser would read as an empty second field."""
    return "\t".join(f for f in (x.strip() for x in fields) if f)


# --------------------------------------------------------------- per-dataset rules

def donnellan(rows):
    out, iso = [], []
    for r in rows:
        acc, was_iso = strip_isoform(r["Accession"])
        if was_iso:
            iso.append(r["Accession"])
        out.append(line(acc, clean_desc(r["Description"])))
    return out, iso


def ashour_cellfree(rows):
    # Entry name plus the protein name, so a retired entry name still resolves.
    return [line(r["Uniprot ID"], r["Protein"]) for r in rows], []


def name_gene(rows, name_col, gene_col):
    return [line(r[name_col], r[gene_col]) for r in rows], []


def alhujaily(rows):
    # Gene first where present; 9 TABLE 2 footnote rows have no gene symbol at
    # all and fall back to the protein name alone.
    return [line(r["Gene"], r["Name of protein"]) for r in rows], []


def irshad_hmec1(rows):
    out = []
    for r in rows:
        m = HINT_RE.match(r["Protein"])
        out.append(line(m.group(1), m.group(2)) if m else line(r["Protein"]))
    return out, []


def zheng_s1(rows):
    return [line(r["Protein ID"], r["Protein name"]) for r in rows], []


def zheng_fig4(rows):
    # Accession where the Table S1 lookup supplied one, else the bare symbol.
    return [line(r["Protein ID"] or r["Gene name"]) for r in rows], []


PLAN = [
    ("donnellan-2022-caov3",   "Donnellan et al. 2022 - Caov3.tsv",              donnellan,        "Accession (isoform stripped) + Description"),
    ("donnellan-2022-ov90",    "Donnellan et al. 2022 - OV90.tsv",               donnellan,        "Accession (isoform stripped) + Description"),
    ("donnellan-2022-pbl",     "Donnellan et al. 2022 - PBL.tsv",                donnellan,        "Accession (isoform stripped) + Description"),
    ("donnellan-2022-wil2ns",  "Donnellan et al. 2022 - WIL2-NS.tsv",            donnellan,        "Accession (isoform stripped) + Description"),
    ("ashour-2020-cellfree",   "Ashour et al. 2020 - PDLF cell-free.tsv",        ashour_cellfree,  "Uniprot ID (entry name) + Protein"),
    ("ashour-2020-incell",     "Ashour et al. 2020 - PDLF in-cell.tsv",
     lambda rs: name_gene(rs, "Protein", "Gene"),                                                  "Protein + Gene"),
    ("irshad-2019-haec",       "Irshad et al. 2019 - HAEC.tsv",
     lambda rs: name_gene(rs, "Protein", "Gene"),                                                  "Protein + Gene"),
    ("irshad-2019-hmec1",      "Irshad et al. 2019 - HMEC-1.tsv",                irshad_hmec1,     "Protein split on '(HINT)' -> name + hint"),
    ("alhujaily-2021-hek293",  "Alhujaily et al. 2021 - HEK293.tsv",             alhujaily,        "Gene + Name of protein"),
    ("zheng-2024-shsy5y",      "Zheng et al. 2024 - SH-SY5Y.tsv",                zheng_s1,         "Protein ID + Protein name"),
    ("zheng-2024-shsy5y-fig4", "Zheng et al. 2024 - SH-SY5Y (Fig4 recovered).tsv", zheng_fig4,     "Protein ID where known, else bare Gene name"),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(ROOT, "02-data-processed",
                                                  "uniprot-input"))
    ap.add_argument("--print-plan", action="store_true",
                    help="show the per-dataset rule and exit")
    args = ap.parse_args()

    if args.print_plan:
        for key, src, _, rule in PLAN:
            print(f"{key:26s} {rule}\n{'':26s} <- {src}")
        return 0

    os.makedirs(args.out, exist_ok=True)
    total, total_iso = 0, 0
    print(f"{'dataset':26s} {'lines':>6}  {'isoforms':>8}  example")
    for key, src, fn, _rule in PLAN:
        if not os.path.exists(os.path.join(SRC, src)):
            print(f"{key:26s}  SOURCE MISSING: {src}")
            return 1
        rows = read(src)
        lines, iso = fn(rows)
        if len(lines) != len(rows):
            print(f"{key:26s}  FAIL: {len(lines)} lines from {len(rows)} rows")
            return 1
        blank = [i for i, l in enumerate(lines, 1) if not l.strip()]
        if blank:
            print(f"{key:26s}  FAIL: blank clue-set on row(s) {blank[:5]}")
            return 1
        path = os.path.join(args.out, f"{key}.txt")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
        total += len(lines)
        total_iso += len(iso)
        print(f"{key:26s} {len(lines):6d}  {len(iso):8d}  {lines[0][:58]!r}")

    print(f"\n{total} clue-sets written to {args.out}")
    print(f"{total_iso} isoform suffixes stripped (canonical entry is the target;"
          f" the original accession is still in 01-data-extracted/)")
    print("\nNext: run uniprot-name-match over each file, e.g.")
    print(f"  for f in {args.out}/*.txt; do")
    print('    python3 uniprot-name-match/scripts/run_uniprot_name_match.py "$f" \\')
    print('      -o "${f%.txt}.uniprot.tsv"')
    print("  done")
    print("\nParse the results with a real CSV reader: when a clue-set has two "
          "fields the `From` column contains a TAB and is quoted, so `cut -f` "
          "mis-splits it. The skill also writes CRLF line endings.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

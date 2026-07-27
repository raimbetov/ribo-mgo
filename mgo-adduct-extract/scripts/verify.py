#!/usr/bin/env python3
"""Cross-cutting verification of extracted TSVs against their dataset specs.

Independent of how each file was produced. Per-table checks (site-sum
reconciliations, cross-method identifier agreement) run inside extract.py at build
time; this script re-checks the written files:

  * well-formed TSV -- no ragged rows, no stray tabs or newlines
  * row count matches the spec's expected_rows, with the paper's own wording as
    the stated evidence
  * every row carries traceable provenance
  * cytoplasmic ribosomal proteins present, as a smoke test that the intended
    table was extracted rather than a neighbouring one

    python3 verify.py --repo-root /path/to/repo
"""
import argparse
import csv
import os
import re
import sys

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required: pip install -r ../requirements.txt")

# Papers identify proteins inconsistently -- gene symbols (Donnellan, Irshad,
# Zheng) or UniProt entry names (Ashour) -- so match both forms.
RP_SYMBOL = re.compile(
    r"\b(?:RP[LS]\d+[A-Z]?\d*|RPLP[0-2]|RPSA|FAU|UBA52|RACK1)\b")
RP_ENTRY = re.compile(
    r"\b(?:RS[0-9]+[A-Z]?|RL[0-9]+[A-Z]?|RLA[0-9]|RSSA|GBLP)_HUMAN\b")


def rp_reference(root):
    """Cytoplasmic ribosomal protein symbols from the repo's reference set."""
    rpdir = os.path.join(root, "rp-script")
    syms = set()
    if not os.path.isdir(rpdir):
        return syms
    for fn in os.listdir(rpdir):
        if fn.endswith((".tsv", ".csv", ".md")):
            with open(os.path.join(rpdir, fn), encoding="utf-8",
                      errors="replace") as fh:
                syms |= set(RP_SYMBOL.findall(fh.read()))
    return syms


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo-root",
                    default=os.path.dirname(os.path.dirname(here)))
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    root = os.path.abspath(args.repo_root)
    outdir = os.path.abspath(args.out or os.path.join(root, "01-data-extracted"))
    dsdir = os.path.join(here, os.pardir, "datasets")

    rps = rp_reference(root)
    print(f"reference cytoplasmic RP symbols loaded: {len(rps)}\n")

    failures, total = [], 0
    for fn in sorted(f for f in os.listdir(dsdir) if f.endswith(".yaml")):
        with open(os.path.join(dsdir, fn), encoding="utf-8") as fh:
            spec = yaml.safe_load(fh)
        path = os.path.join(outdir, spec["output"])
        if not os.path.exists(path):
            print(f"FAIL  {spec['output']}\n        not built")
            failures.append(spec["output"])
            continue

        with open(path, newline="", encoding="utf-8") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        header, data = rows[0], rows[1:]
        total += 1

        ragged = [i for i, r in enumerate(data, 2) if len(r) != len(header)]
        prov = [header.index(c) for c in ("source_file", "source_table")
                if c in header]
        noprov = [i for i, r in enumerate(data, 2)
                  if any(not r[j].strip() for j in prov)]
        exp = spec.get("expected_rows")

        hits = set()
        for r in data:
            txt = " ".join(r)
            hits |= {t for t in RP_SYMBOL.findall(txt) if t in rps}
            hits |= set(RP_ENTRY.findall(txt))

        ok = not ragged and not noprov and (exp is None or len(data) == exp)
        if not ok:
            failures.append(spec["output"])
        print(f"{'PASS' if ok else 'FAIL'}  {spec['output']}")
        print(f"        {spec['cell_line']} | {spec['system']} | {spec['condition']}")
        print(f"        rows={len(data)}" + (f" expected={exp}" if exp else "")
              + f"  cols={len(header)}  ragged={len(ragged)}"
              + f"  rows_missing_provenance={len(noprov)}")
        print(f"        count evidence: {' '.join(spec['count_evidence'].split())}")
        print(f"        cytoplasmic RPs present: {len(hits)}"
              + (f"  e.g. {sorted(hits)[:8]}" if hits else ""))

    print()
    if failures:
        print("FAILURES: " + ", ".join(failures))
        return 1
    print(f"All {total} files passed structural and count checks.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

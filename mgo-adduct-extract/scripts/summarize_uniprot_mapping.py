#!/usr/bin/env python3
"""QA the stage 01 -> 02 UniProt mapping and report the ribosomal intersection.

    python3 mgo-adduct-extract/scripts/summarize_uniprot_mapping.py

Checks the contract the mapper promises, then answers the project's actual
question: which cytoplasmic ribosomal proteins carry an MGO adduct, per paper
and cell line.

Contract checks (a failure here means the mapping is not usable as-is):
  * one output row per input line, in order -- inputs are never dropped
  * `From` is the raw input line, verbatim, so every row traces to a supplement
  * accepted rows carry a well-formed accession
"""
import csv
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IN_DIR = os.path.join(ROOT, "02-data-processed", "uniprot-input")
MAP_DIR = os.path.join(ROOT, "02-data-processed", "uniprot-mapped")
RP_CSV = os.path.join(ROOT, "rp-script", "human_ribosomal_proteins.csv")

ACC_RE = re.compile(
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}$"
)


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_rp():
    if not os.path.exists(RP_CSV):
        return None
    rp = {}
    with open(RP_CSV, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rp[r["UniProt_ID"].strip()] = (r["Gene_Symbol"].strip(),
                                           r["Unified_Nomenclature"].strip())
    return rp


def main():
    rp = load_rp()
    problems = []
    grand = {"rows": 0, "accepted": 0, "review": 0}
    rp_union = {}

    print("%-26s %6s %9s %8s   %s" % ("dataset", "rows", "accepted", "review", "ribosomal"))
    print("-" * 74)

    for src in sorted(glob.glob(os.path.join(IN_DIR, "*.txt"))):
        key = os.path.basename(src)[:-4]
        final = os.path.join(MAP_DIR, f"{key}.uniprot.tsv")
        review = os.path.join(MAP_DIR, f"{key}.uniprot.review.tsv")
        if not os.path.exists(final):
            print("%-26s  NOT MAPPED" % key)
            problems.append(f"{key}: no output")
            continue

        inputs = [l.rstrip("\n") for l in open(src, encoding="utf-8") if l.strip()]
        rows = read_tsv(final)
        rev = read_tsv(review) if os.path.exists(review) else []

        # contract: one row per input, same order, From verbatim
        if len(rows) != len(inputs):
            problems.append(f"{key}: {len(rows)} output rows for {len(inputs)} inputs")
        else:
            for i, (want, got) in enumerate(zip(inputs, rows), 1):
                if got["From"].strip() != want.strip():
                    problems.append(f"{key} row {i}: From={got['From']!r} != input {want!r}")
                    break

        accepted = [r for r in rows if r["Entry"]]
        for r in accepted:
            if not ACC_RE.fullmatch(r["Entry"]):
                problems.append(f"{key}: malformed accession {r['Entry']!r}")
                break

        hits = {}
        if rp:
            for r in accepted:
                if r["Entry"] in rp:
                    hits[r["Entry"]] = rp[r["Entry"]]
            rp_union.update(hits)

        grand["rows"] += len(rows)
        grand["accepted"] += len(accepted)
        grand["review"] += len(rev)
        pct = 100.0 * len(accepted) / len(rows) if rows else 0
        print("%-26s %6d %5d %3.0f%% %8d   %d" % (key, len(rows), len(accepted), pct, len(rev), len(hits)))

    print("-" * 74)
    pct = 100.0 * grand["accepted"] / grand["rows"] if grand["rows"] else 0
    print("%-26s %6d %5d %3.0f%% %8d   %d unique" % (
        "TOTAL", grand["rows"], grand["accepted"], pct, grand["review"], len(rp_union)))

    # Why did the review cases land there? Confidence alone is not actionable;
    # "no candidates at all" needs a different fix from "two plausible entries".
    buckets = {}
    zero_cand = []
    for src in sorted(glob.glob(os.path.join(IN_DIR, "*.txt"))):
        key = os.path.basename(src)[:-4]
        review = os.path.join(MAP_DIR, f"{key}.uniprot.review.tsv")
        if not os.path.exists(review):
            continue
        for r in read_tsv(review):
            buckets[(r["Confidence"], r["Reason"])] = buckets.get(
                (r["Confidence"], r["Reason"]), 0) + 1
            if not r["Selected Entry"]:
                zero_cand.append((key, r["From"]))
    if buckets:
        print("\nReview cases by reason:")
        for (conf, reason), n in sorted(buckets.items(), key=lambda kv: -kv[1]):
            print(f"  {n:5d}  [{conf}] {reason}")
    if zero_cand:
        print(f"\nNo candidate retrieved at all ({len(zero_cand)}) -- these are the "
              f"only inputs UniProt could not answer:")
        for key, frm in zero_cand[:20]:
            print(f"  {key}: {frm!r}")
        if len(zero_cand) > 20:
            print(f"  ... and {len(zero_cand) - 20} more")

    if rp:
        print(f"\nCytoplasmic ribosomal proteins with an MGO adduct "
              f"({len(rp_union)}/{len(rp)} of the reference set):")
        for acc, (gene, uni) in sorted(rp_union.items(), key=lambda kv: kv[1][0]):
            print(f"  {acc}  {gene:<10} {uni}")
    else:
        print("\nrp-script/human_ribosomal_proteins.csv not found; skipped intersection.")

    if problems:
        print(f"\nCONTRACT PROBLEMS ({len(problems)}):")
        for p in problems:
            print("  -", p)
        return 1
    print("\nContract checks passed: row counts, order, verbatim From, accession format.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Merge hand-resolved review decisions into the mapped datasets (stage 02 final).

    python3 mgo-adduct-extract/scripts/apply_user_reviews.py

Reads 02-data-processed/uniprot-mapped/user-reviewed/*.user-reviewed.tsv and
folds those decisions into the machine mapping, writing the completed datasets
to 03-data-final/.

The user-reviewed files carry, per row:
    From | Entry | Entry Name | Status | Notes
with Status either `reviewed` (use this entry) or `drop` (no UniProtKB entry
exists -- the Notes give the UniParc ID).

Three things this does that a plain join would not:

* VERIFIES every assigned accession against live UniProt instead of trusting
  the typed entry name, and fills Protein names / Gene Names / Organism from
  the record itself, so a hand-entered row is as complete as a machine one.
* HONOURS `drop`. Those rows leave the dataset, so the output is no longer one
  row per input line -- 03-data-final/ is the only stage where that holds, and
  the reason the row count is reported explicitly.
* ALLOWS ONE INPUT TO YIELD SEVERAL ROWS. Irshad prints "Nucleoside diphosphate
  kinase (NME1-NME2)", one label covering two gene products; the review assigns
  both NDKA and NDKB.

Provenance is kept: the `Source` column marks each row `auto` or `user-reviewed`,
and `From` still holds the input line verbatim.
"""
import argparse
import csv
import glob
import json
import os
import sys
import time
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MAP_DIR = os.path.join(ROOT, "02-data-processed", "uniprot-mapped")
UR_DIR = os.path.join(MAP_DIR, "user-reviewed")
OUT_DIR = os.path.join(ROOT, "03-data-final")
RP_CSV = os.path.join(ROOT, "rp-script", "human_ribosomal_proteins.csv")

COLUMNS = ["From", "Entry", "Entry Name", "Protein names", "Gene Names",
           "Organism", "Source", "Notes"]
API = "https://rest.uniprot.org/uniprotkb/search"


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_user_reviewed(path):
    """Positional, not DictReader: both entry columns are headed 'Assigned by
    user', and DictReader would silently collapse them."""
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    out = []
    for row in rows[1:]:
        if not row or not row[0].strip():
            continue
        row = row + [""] * (5 - len(row))
        out.append({"From": row[0], "entry": row[1].strip(),
                    "entry_name": row[2].strip(),
                    "status": row[3].strip().lower(), "notes": row[4].strip()})
    return out


def fetch(accessions, delay=0.2):
    """Look up records in batches; returns {accession: record}."""
    found = {}
    accessions = sorted(set(a for a in accessions if a))
    for i in range(0, len(accessions), 25):
        batch = accessions[i:i + 25]
        query = " OR ".join(f"accession:{a}" for a in batch)
        url = (f"{API}?query={urllib.parse.quote(query)}&size=500"
               "&fields=accession,id,protein_name,gene_names,organism_name")
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                data = json.load(fh)
        except Exception as exc:                       # noqa: BLE001
            print(f"  WARNING: lookup failed for {batch[0]}..: {exc}")
            continue
        for rec in data.get("results", []):
            found[rec.get("primaryAccession", "")] = rec
        time.sleep(delay)
    return found


def describe(rec):
    desc = rec.get("proteinDescription", {})
    names = []
    rn = desc.get("recommendedName", {}).get("fullName", {}).get("value")
    if rn:
        names.append(rn)
    for alt in desc.get("alternativeNames", []):
        v = alt.get("fullName", {}).get("value")
        if v:
            names.append(f"({v})")
    for sub in desc.get("submissionNames", []):
        v = sub.get("fullName", {}).get("value")
        if v:
            names.append(f"({v})" if names else v)
    genes = []
    for g in rec.get("genes", []):
        v = g.get("geneName", {}).get("value")
        if v and v not in genes:
            genes.append(v)
        for syn in g.get("synonyms", []):
            v = syn.get("value")
            if v and v not in genes:
                genes.append(v)
    org = rec.get("organism", {})
    label = org.get("scientificName", "")
    if label and org.get("commonName"):
        label = f"{label} ({org['commonName']})"
    return " ".join(names), " ".join(genes), label


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=OUT_DIR)
    args = ap.parse_args()

    ur_files = sorted(glob.glob(os.path.join(UR_DIR, "*.user-reviewed.tsv")))
    decisions = {}
    for path in ur_files:
        key = os.path.basename(path).replace(".uniprot.user-reviewed.tsv", "")
        decisions[key] = read_user_reviewed(path)

    wanted = [d["entry"] for rows in decisions.values() for d in rows
              if d["status"] != "drop"]
    print(f"Verifying {len(set(wanted))} user-assigned accessions against UniProt...")
    records = fetch(wanted)

    problems = []
    for key, rows in decisions.items():
        for d in rows:
            if d["status"] == "drop":
                continue
            rec = records.get(d["entry"])
            if rec is None:
                problems.append(f"{key}: {d['entry']} did not resolve in UniProt")
            elif d["entry_name"] and rec.get("uniProtkbId", "") != d["entry_name"]:
                problems.append(
                    f"{key}: {d['entry']} is {rec.get('uniProtkbId')}, "
                    f"review says {d['entry_name']}")
    if problems:
        print("\nASSIGNMENT PROBLEMS:")
        for p in problems:
            print("  -", p)
        return 1
    print("All user-assigned accessions resolve and their entry names agree.\n")

    os.makedirs(args.out, exist_ok=True)
    rp = {}
    if os.path.exists(RP_CSV):
        with open(RP_CSV, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                rp[r["UniProt_ID"].strip()] = (r["Gene_Symbol"].strip(),
                                               r["Unified_Nomenclature"].strip())

    print("%-26s %6s %8s %7s %7s %6s   %s" %
          ("dataset", "input", "resolved", "byhand", "added", "drop", "ribosomal"))
    print("-" * 82)
    totals = dict(inp=0, res=0, hand=0, added=0, drop=0)
    rp_union = {}
    dropped_rows = []

    for src in sorted(glob.glob(os.path.join(MAP_DIR, "*.uniprot.tsv"))):
        key = os.path.basename(src).replace(".uniprot.tsv", "")
        rows = read_tsv(src)
        by_from = {}
        for d in decisions.get(key, []):
            by_from.setdefault(d["From"], []).append(d)

        out_rows, n_hand, n_added, n_drop = [], 0, 0, 0
        for r in rows:
            ds = by_from.get(r["From"])
            if not ds:
                if r["Entry"]:
                    out_rows.append({**{c: r.get(c, "") for c in COLUMNS[:6]},
                                     "Source": "auto", "Notes": ""})
                continue
            if ds[0]["status"] == "drop":
                n_drop += 1
                dropped_rows.append((key, r["From"], ds[0]["notes"]))
                continue
            for i, d in enumerate(ds):
                rec = records[d["entry"]]
                pname, genes, org = describe(rec)
                out_rows.append({
                    "From": r["From"], "Entry": d["entry"],
                    "Entry Name": rec.get("uniProtkbId", ""),
                    "Protein names": pname, "Gene Names": genes, "Organism": org,
                    "Source": "user-reviewed", "Notes": d["notes"]})
                if i == 0:
                    n_hand += 1
                else:
                    n_added += 1

        hits = {r["Entry"]: rp[r["Entry"]] for r in out_rows if r["Entry"] in rp}
        rp_union.update(hits)

        dest = os.path.join(args.out, f"{key}.uniprot.final.tsv")
        with open(dest, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
            w.writeheader()
            w.writerows(out_rows)

        totals["inp"] += len(rows)
        totals["res"] += len(out_rows)
        totals["hand"] += n_hand
        totals["added"] += n_added
        totals["drop"] += n_drop
        print("%-26s %6d %8d %7d %7d %6d   %d" %
              (key, len(rows), len(out_rows), n_hand, n_added, n_drop, len(hits)))

    print("-" * 82)
    print("%-26s %6d %8d %7d %7d %6d   %d unique" %
          ("TOTAL", totals["inp"], totals["res"], totals["hand"],
           totals["added"], totals["drop"], len(rp_union)))

    if dropped_rows:
        print(f"\nDropped ({len(dropped_rows)}) -- no UniProtKB entry exists:")
        for key, frm, note in dropped_rows:
            print(f"  {key}: {frm.replace(chr(9), ' | ')}  [{note or 'no note'}]")

    if rp:
        print(f"\nCytoplasmic ribosomal proteins with an MGO adduct "
              f"({len(rp_union)}/{len(rp)} of the reference set):")
        for acc, (gene, uni) in sorted(rp_union.items(), key=lambda kv: kv[1][0]):
            print(f"  {acc}  {gene:<10} {uni}")

    print(f"\nWrote {len(glob.glob(os.path.join(args.out, '*.tsv')))} datasets to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

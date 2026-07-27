#!/usr/bin/env python3
"""Merge 03-data-final/ into one canonical list of MGO-glycated proteins.

    python3 mgo-adduct-extract/scripts/merge_final_datasets.py

One row per UniProt accession, with the provenance of every dataset that
reported it: study, cell line, experimental system, and the label the paper
printed. Writes 04-data-merged/mgo-glycated-proteins.tsv.

`Evidence` is the column to read before drawing any conclusion. Irshad's HMEC-1
and Ashour's PDLF cell-free datasets are *lysates spiked with exogenous MG*, not
adducts formed in living cells, and they are an order of magnitude deeper than
the in-cell lists precisely because of that (220 and 172 rows vs 2 and 5).
Pooling the two silently would let cell-free depth masquerade as in-cell
coverage, so each protein is marked `in-cell`, `cell-free`, or `both`.
"""
import argparse
import csv
import glob
import os
import sys
from collections import OrderedDict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, "03-data-final")
RP_CSV = os.path.join(ROOT, "rp-script", "human_ribosomal_proteins.csv")
OUT = os.path.join(ROOT, "04-data-merged", "mgo-glycated-proteins.tsv")

# study, cell line, system (verbatim from 01-data-extracted/MANIFEST.md),
# evidence class, caveat
PROVENANCE = {
    "donnellan-2022-wil2ns":  ("Donnellan 2022", "WIL2-NS",        "in-cell, 500 uM MGO 24 h",      "in-cell",   ""),
    "donnellan-2022-ov90":    ("Donnellan 2022", "OV90 parental",  "in-cell, basal",                "in-cell",   ""),
    "donnellan-2022-caov3":   ("Donnellan 2022", "Caov3 parental", "in-cell, basal",                "in-cell",   ""),
    "donnellan-2022-pbl":     ("Donnellan 2022", "PBL, 3 donors",  "ex vivo, basal",                "in-cell",   "primary cells, ex vivo"),
    "irshad-2019-hmec1":      ("Irshad 2019",    "HMEC-1",         "cell-free lysate + 500 uM MG",  "cell-free", ""),
    "irshad-2019-haec":       ("Irshad 2019",    "HAEC",           "in-cell, 20 mM glucose",        "in-cell",   "from main-text prose"),
    "ashour-2020-cellfree":   ("Ashour 2020",    "PDLF primary",   "cell-free lysate + MG",         "cell-free", ""),
    "ashour-2020-incell":     ("Ashour 2020",    "PDLF primary",   "in-cell, endogenous",           "in-cell",   "from main-text prose"),
    "alhujaily-2021-hek293":  ("Alhujaily 2021", "HEK293",         "in-cell, 131 uM MG 6 h",        "in-cell",   "supplement covers ~21% of the 681 reported"),
    "zheng-2024-shsy5y":      ("Zheng 2024",     "SH-SY5Y",        "in-cell, 0-1000 uM MGO",        "in-cell",   ""),
    "zheng-2024-shsy5y-fig4": ("Zheng 2024",     "SH-SY5Y",        "in-cell, 0-1000 uM MGO",        "in-cell",   "recovered from Fig. 4 heatmap labels"),
}

COLUMNS = [
    "Entry", "Entry Name", "Gene Names", "Protein names",
    "Ribosomal", "RP Gene", "RP Nomenclature", "RP Subunit",
    "Evidence", "N studies", "N datasets",
    "Studies", "Cell lines", "Systems", "Datasets",
    "Source labels", "Caveats", "Notes",
]


def join(values):
    """De-duplicated, order-preserving, semicolon-joined."""
    return "; ".join(OrderedDict.fromkeys(v for v in values if v))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()

    rp = {}
    if os.path.exists(RP_CSV):
        with open(RP_CSV, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                rp[r["UniProt_ID"].strip()] = (r["Gene_Symbol"].strip(),
                                               r["Unified_Nomenclature"].strip(),
                                               r["Subunit"].strip())

    files = sorted(glob.glob(os.path.join(SRC, "*.uniprot.final.tsv")))
    if not files:
        print(f"No datasets in {SRC}")
        return 1

    proteins = OrderedDict()
    unknown = []
    for path in files:
        key = os.path.basename(path).replace(".uniprot.final.tsv", "")
        if key not in PROVENANCE:
            unknown.append(key)
            continue
        study, line, system, evidence, caveat = PROVENANCE[key]
        with open(path, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        seen_here = set()
        for r in rows:
            acc = r["Entry"]
            if not acc:
                continue
            p = proteins.setdefault(acc, {
                "Entry": acc, "Entry Name": r["Entry Name"],
                "Gene Names": r["Gene Names"], "Protein names": r["Protein names"],
                "studies": [], "lines": [], "systems": [], "datasets": [],
                "labels": [], "caveats": [], "notes": [], "evidence": set(),
            })
            # A dataset counts once per protein even when the paper printed it
            # on several rows (Donnellan's histone accessions recur with
            # different site sets).
            if (acc, key) not in seen_here:
                seen_here.add((acc, key))
                p["studies"].append(study)
                p["lines"].append(line)
                p["systems"].append(system)
                p["datasets"].append(key)
                p["evidence"].add(evidence)
                if caveat:
                    p["caveats"].append(f"{key}: {caveat}")
            label = r["From"].replace("\t", " | ")
            if label not in p["labels"]:
                p["labels"].append(label)
            if r.get("Notes"):
                note = f"{key}: {r['Notes']}"
                if note not in p["notes"]:
                    p["notes"].append(note)

    if unknown:
        print("Datasets with no provenance entry (add them to PROVENANCE):")
        for k in unknown:
            print("  -", k)
        return 1

    rows_out = []
    for acc, p in proteins.items():
        ev = p["evidence"]
        evidence = "both" if len(ev) > 1 else next(iter(ev))
        gene, nom, sub = rp.get(acc, ("", "", ""))
        rows_out.append({
            "Entry": acc, "Entry Name": p["Entry Name"],
            "Gene Names": p["Gene Names"], "Protein names": p["Protein names"],
            "Ribosomal": "yes" if acc in rp else "",
            "RP Gene": gene, "RP Nomenclature": nom, "RP Subunit": sub,
            "Evidence": evidence,
            "N studies": len(set(p["studies"])), "N datasets": len(p["datasets"]),
            "Studies": join(p["studies"]), "Cell lines": join(p["lines"]),
            "Systems": join(p["systems"]), "Datasets": join(p["datasets"]),
            "Source labels": join(p["labels"]),
            "Caveats": join(p["caveats"]), "Notes": join(p["notes"]),
        })

    # Most-corroborated first, ribosomal ahead of the rest at equal support, then
    # accession so the file is stable across rebuilds.
    rows_out.sort(key=lambda r: (-r["N studies"], -r["N datasets"],
                                 r["Ribosomal"] != "yes", r["Entry"]))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)

    # ---------------------------------------------------------------- report
    total = len(rows_out)
    by_ev = {}
    for r in rows_out:
        by_ev[r["Evidence"]] = by_ev.get(r["Evidence"], 0) + 1
    print(f"{sum(1 for _ in files)} datasets -> {total} unique proteins\n")
    print("Evidence:")
    for k in ("in-cell", "cell-free", "both"):
        if k in by_ev:
            print(f"  {k:<10} {by_ev[k]:5d}")
    print("\nCorroboration:")
    for n in sorted({r["N studies"] for r in rows_out}, reverse=True):
        c = sum(1 for r in rows_out if r["N studies"] == n)
        print(f"  reported by {n} stud{'y' if n == 1 else 'ies'}: {c}")

    ribo = [r for r in rows_out if r["Ribosomal"] == "yes"]
    print(f"\nCytoplasmic ribosomal proteins: {len(ribo)}/{len(rp)}")
    rb = {}
    for r in ribo:
        rb[r["Evidence"]] = rb.get(r["Evidence"], 0) + 1
    for k in ("in-cell", "cell-free", "both"):
        if k in rb:
            print(f"  {k:<10} {rb[k]:5d}")
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

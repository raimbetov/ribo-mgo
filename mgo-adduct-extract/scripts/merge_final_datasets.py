#!/usr/bin/env python3
"""Merge 03-data-final/ into one canonical list of MGO-glycated proteins.

    python3 mgo-adduct-extract/scripts/merge_final_datasets.py
    python3 mgo-adduct-extract/scripts/merge_final_datasets.py --no-canonicalize

Writes 04-data-merged/mgo-glycated-proteins.tsv -- one row per protein, with the
provenance of every dataset that reported it: study, cell line, experimental
system, and the label the paper printed.

CANONICALIZATION. The papers key their tables differently: Donnellan's
supplement is Proteome Discoverer output carrying TrEMBL and isoform
accessions, while the others give gene symbols or Swiss-Prot IDs. Left alone,
the same protein appears under several accessions -- and worse, a protein
reported *only* under a TrEMBL accession never matches a reference set keyed by
Swiss-Prot, so it silently vanishes from any intersection. That is how RPS14
(A0A2R8Y811) and RPL10 (F8W7C6) went missing from the ribosomal count.

So each unreviewed accession is resolved to its reviewed Swiss-Prot counterpart
via an exact gene-symbol lookup, and rows that collapse onto one entry are
merged. Every original accession is kept in the `Accessions` column and every
substitution is logged to 04-data-merged/accession-canonicalization.tsv, so
nothing is lost and each step can be checked. An entry that cannot be resolved
-- no gene name, no reviewed entry, or an ambiguous gene -- is kept exactly as
it was rather than guessed at.

EVIDENCE. Read the `Evidence` column before drawing conclusions. Irshad HMEC-1
and Ashour PDLF cell-free are lysates spiked with exogenous MG, not adducts
formed in living cells, and are an order of magnitude deeper for that reason
(220 and 172 rows vs 2 and 5). Pooling silently would let cell-free depth
masquerade as in-cell coverage, so each protein is marked `in-cell`,
`cell-free`, or `both`.
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
from collections import OrderedDict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, "03-data-final")
RP_CSV = os.path.join(ROOT, "rp-script", "human_ribosomal_proteins.csv")
OUT_DIR = os.path.join(ROOT, "04-data-merged")
OUT = os.path.join(OUT_DIR, "mgo-glycated-proteins.tsv")
AUDIT = os.path.join(OUT_DIR, "accession-canonicalization.tsv")
API = "https://rest.uniprot.org/uniprotkb/search"
REVIEWED = "UniProtKB reviewed (Swiss-Prot)"

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
    "Entry", "Entry Name", "Gene", "Gene Names", "Protein names", "Reviewed",
    "Ribosomal", "RP Gene", "RP Nomenclature", "RP Subunit",
    "Evidence", "N studies", "N datasets",
    "Studies", "Cell lines", "Systems", "Datasets",
    "Accessions", "Source labels", "Caveats", "Notes",
]
AUDIT_COLUMNS = ["From accession", "From entry name", "Gene used",
                 "To accession", "To entry name", "Outcome"]


def dedup(values):
    return list(OrderedDict.fromkeys(v for v in values if v))


def join(values):
    return "; ".join(dedup(values))


def api_get(query, fields, size=500, timeout=60):
    url = f"{API}?query={urllib.parse.quote(query)}&fields={fields}&size={size}"
    with urllib.request.urlopen(url, timeout=timeout) as fh:
        return json.load(fh).get("results", [])


def fetch_records(accessions, delay=0.15):
    """{accession: record} for a set of accessions, in batches."""
    out = {}
    accessions = sorted(set(a for a in accessions if a))
    fields = "accession,id,protein_name,gene_names,organism_name"
    for i in range(0, len(accessions), 25):
        batch = accessions[i:i + 25]
        query = " OR ".join(f"accession:{a}" for a in batch)
        try:
            for rec in api_get(query, fields):
                out[rec.get("primaryAccession", "")] = rec
        except Exception as exc:                                # noqa: BLE001
            print(f"  WARNING: metadata lookup failed near {batch[0]}: {exc}")
        time.sleep(delay)
    return out


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
    return " ".join(names), genes


def canonicalise(records, delay=0.15):
    """{accession: canonical accession} plus an audit row per unreviewed entry.

    Only unreviewed entries move, and only on an exact gene-symbol match to a
    single reviewed human entry. Anything ambiguous stays put -- a wrong merge
    is far more damaging than an unmerged row, because it silently attributes
    one paper's evidence to a different protein.
    """
    mapping, audit, wanted = {}, [], set()
    for acc, rec in sorted(records.items()):
        if rec.get("entryType") == REVIEWED:
            continue
        _, genes = describe(rec)
        entry_name = rec.get("uniProtkbId", "")
        if not genes:
            audit.append([acc, entry_name, "", "", "", "kept: no gene name"])
            continue
        resolved = False
        for gene in genes:
            try:
                hits = api_get(
                    f"gene_exact:{gene} AND reviewed:true AND organism_id:9606",
                    "accession,id", size=10)
            except Exception:                                   # noqa: BLE001
                hits = []
            time.sleep(delay)
            if len(hits) == 1:
                target = hits[0]["primaryAccession"]
                mapping[acc] = target
                wanted.add(target)
                audit.append([acc, entry_name, gene, target,
                              hits[0].get("uniProtkbId", ""), "canonicalised"])
                resolved = True
                break
            if len(hits) > 1:
                audit.append([acc, entry_name, gene, "", "",
                              f"kept: {len(hits)} reviewed entries for this gene"])
                resolved = True
                break
        if not resolved:
            audit.append([acc, entry_name, genes[0], "", "",
                          "kept: no reviewed entry for this gene"])
    return mapping, audit, wanted


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--audit", default=AUDIT)
    ap.add_argument("--no-canonicalize", action="store_true",
                    help="key strictly on the accession each paper printed")
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

    # ---- collect one contribution per (dataset row)
    contributions, unknown = [], []
    for path in files:
        key = os.path.basename(path).replace(".uniprot.final.tsv", "")
        if key not in PROVENANCE:
            unknown.append(key)
            continue
        with open(path, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["Entry"]:
                    contributions.append((r["Entry"], key,
                                          r["From"].replace("\t", " | "),
                                          r.get("Notes", "")))
    if unknown:
        print("Datasets with no provenance entry (add them to PROVENANCE):")
        for k in unknown:
            print("  -", k)
        return 1

    accessions = {c[0] for c in contributions}
    print(f"{len(files)} datasets, {len(contributions)} rows, "
          f"{len(accessions)} distinct accessions")

    print("Fetching UniProt metadata...")
    records = fetch_records(accessions)
    missing = accessions - set(records)
    if missing:
        print(f"  WARNING: no metadata for {len(missing)}: {sorted(missing)[:5]}")

    mapping, audit = {}, []
    if not args.no_canonicalize:
        n_unrev = sum(1 for r in records.values() if r.get("entryType") != REVIEWED)
        print(f"Canonicalising {n_unrev} unreviewed entries...")
        mapping, audit, wanted = canonicalise(records)
        extra = wanted - set(records)
        if extra:
            records.update(fetch_records(extra))
        os.makedirs(os.path.dirname(os.path.abspath(args.audit)), exist_ok=True)
        with open(args.audit, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(AUDIT_COLUMNS)
            w.writerows(audit)
        moved = sum(1 for a in audit if a[5] == "canonicalised")
        print(f"  {moved} moved to a reviewed entry, {len(audit) - moved} kept as-is")

    # ---- aggregate on the canonical accession
    proteins = OrderedDict()
    for acc, key, label, notes in contributions:
        target = mapping.get(acc, acc)
        study, line, system, evidence, caveat = PROVENANCE[key]
        p = proteins.setdefault(target, {
            "studies": [], "lines": [], "systems": [], "datasets": [],
            "labels": [], "caveats": [], "notes": [], "accessions": [],
            "evidence": set(),
        })
        p["studies"].append(study)
        p["lines"].append(line)
        p["systems"].append(system)
        p["datasets"].append(key)
        p["evidence"].add(evidence)
        p["accessions"].append(acc)
        p["labels"].append(label)
        if caveat:
            p["caveats"].append(f"{key}: {caveat}")
        if notes:
            p["notes"].append(f"{key}: {notes}")

    rows_out = []
    for acc, p in proteins.items():
        rec = records.get(acc, {})
        pname, genes = describe(rec) if rec else ("", [])
        reviewed = rec.get("entryType") == REVIEWED if rec else False
        gene_sym, nom, sub = rp.get(acc, ("", "", ""))
        # A dataset counts once per protein even where a paper printed it on
        # several rows, or where two of its accessions collapsed onto one entry.
        datasets = dedup(p["datasets"])
        rows_out.append({
            "Entry": acc,
            "Entry Name": rec.get("uniProtkbId", ""),
            "Gene": genes[0] if genes else "",
            "Gene Names": " ".join(genes),
            "Protein names": pname,
            "Reviewed": "yes" if reviewed else "no",
            "Ribosomal": "yes" if acc in rp else "",
            "RP Gene": gene_sym, "RP Nomenclature": nom, "RP Subunit": sub,
            "Evidence": "both" if len(p["evidence"]) > 1 else next(iter(p["evidence"])),
            "N studies": len(set(p["studies"])), "N datasets": len(datasets),
            "Studies": join(p["studies"]), "Cell lines": join(p["lines"]),
            "Systems": join(p["systems"]), "Datasets": "; ".join(datasets),
            "Accessions": join(p["accessions"]),
            "Source labels": join(p["labels"]),
            "Caveats": join(p["caveats"]), "Notes": join(p["notes"]),
        })

    rows_out.sort(key=lambda r: (-r["N studies"], -r["N datasets"],
                                 r["Ribosomal"] != "yes", r["Entry"]))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)

    # ------------------------------------------------------------- report
    print(f"\n{len(accessions)} accessions -> {len(rows_out)} proteins\n")
    by_ev = {}
    for r in rows_out:
        by_ev[r["Evidence"]] = by_ev.get(r["Evidence"], 0) + 1
    print("Evidence:")
    for k in ("in-cell", "cell-free", "both"):
        if k in by_ev:
            print(f"  {k:<10} {by_ev[k]:5d}")
    print(f"\nReviewed (Swiss-Prot): {sum(1 for r in rows_out if r['Reviewed'] == 'yes')}"
          f" / unreviewed: {sum(1 for r in rows_out if r['Reviewed'] == 'no')}")
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
    if audit:
        print(f"Wrote {args.audit}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

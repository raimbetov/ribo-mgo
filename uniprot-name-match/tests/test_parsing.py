#!/usr/bin/env python3
"""Offline tests for input parsing. No network -- these cover the layer that
decides *what gets searched*, which is where silent corruption hides.

    python3 uniprot-name-match/tests/test_parsing.py

Every case below is a real line from 01-data-extracted/, not a synthetic one.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "scripts"))

from fetch_uniprot_candidates import clean_name, infer_record  # noqa: E402

failures = []


def check(label, got, want):
    if got != want:
        failures.append(f"{label}\n     got  {got!r}\n     want {want!r}")


# --- clean_name: list numbering is stripped, meaningful leading digits are not
# A leading number is only numbering when a separator AND whitespace follow it.
for raw, want in [
    ("1. Vimentin", "Vimentin"),
    ("2) Vimentin", "Vimentin"),
    ("3 - Vimentin", "Vimentin"),
    ("12: Vimentin", "Vimentin"),
    ("• Vimentin", "Vimentin"),
    # Regression: these are protein names, not numbered list items. Each one
    # was corrupted before NUMBERING_RE required a separator.
    ("40S ribosomal protein S4", "40S ribosomal protein S4"),
    ("60S ribosomal protein L7", "60S ribosomal protein L7"),
    ("28S ribosomal protein S29, mitochondrial", "28S ribosomal protein S29, mitochondrial"),
    ("39S ribosomal protein L14, mitochondrial", "39S ribosomal protein L14, mitochondrial"),
    ("14-3-3 protein eta", "14-3-3 protein eta"),
    ("14-3-3 protein zeta/delta", "14-3-3 protein zeta/delta"),
    ("78 kDa glucose-regulated protein", "78 kDa glucose-regulated protein"),
    ("10 kDa heat shock protein, mitochondrial", "10 kDa heat shock protein, mitochondrial"),
    ("26S proteasome non-ATPase regulatory subunit 7", "26S proteasome non-ATPase regulatory subunit 7"),
    ("3-ketoacyl-CoA thiolase, mitochondrial", "3-ketoacyl-CoA thiolase, mitochondrial"),
    ("2-amino-3-ketobutyrate coenzyme A ligase, mitochondrial",
     "2-amino-3-ketobutyrate coenzyme A ligase, mitochondrial"),
    ("5'-3' exoribonuclease 1", "5'-3' exoribonuclease 1"),
    ("3'(2'),5'-bisphosphate nucleotidase 1", "3'(2'),5'-bisphosphate nucleotidase 1"),
    # Entry names beginning with digits -- mangling these loses the identifier.
    ("1433Z_HUMAN", "1433Z_HUMAN"),
    ("1433E_HUMAN", "1433E_HUMAN"),
    ("6PGD_HUMAN", "6PGD_HUMAN"),
]:
    check(f"clean_name({raw!r})", clean_name(raw), want)


# --- infer_record: every field on the line has to land somewhere
def fields(line):
    r = infer_record(line)
    return (r.protein_name, r.gene_name, r.identifier)


for line, want in [
    # identifier + name (Donnellan, Zheng S1)
    ("P08670\tVimentin", ("Vimentin", "", "P08670")),
    ("1433Z_HUMAN\t14-3-3 protein zeta/delta",
     ("14-3-3 protein zeta/delta", "", "1433Z_HUMAN")),
    # gene + name (Alhujaily)
    ("VIM\tVimentin", ("Vimentin", "VIM", "")),
    # name + mixed-case hint (Irshad HAEC). Regression: the hint used to be
    # dropped on the floor because it is not upper-case.
    ("Rho GDP-dissociation inhibitor 2\tRhoGDI2",
     ("Rho GDP-dissociation inhibitor 2", "RHOGDI2", "")),
    # name + upper-case hint (Irshad HMEC-1, split from "Name (HINT)")
    ("14-3-3 protein eta\t1433F", ("14-3-3 protein eta", "1433F", "")),
    ("Pyruvate kinase-M\tPKM", ("Pyruvate kinase-M", "PKM", "")),
    # bare gene symbol (Zheng Fig. 4)
    ("VIM", ("VIM", "VIM", "")),
    # bare descriptive name, digits intact
    ("40S ribosomal protein S4", ("40S ribosomal protein S4", "", "")),
    ("60S ribosomal protein L7", ("60S ribosomal protein L7", "", "")),
    # entry name and accession alone
    ("VIME_HUMAN", ("VIME_HUMAN", "", "VIME_HUMAN")),
    ("P08670", ("P08670", "", "P08670")),
]:
    check(f"infer_record({line!r})", fields(line), want)


from fetch_uniprot_candidates import (  # noqa: E402
    exact_name_bonus, entry_name_mnemonic, build_search_strategies, infer_record,
    candidate_rank_key,
)

SWISS = "UniProtKB reviewed (Swiss-Prot)"


def entry(recommended=None, alternatives=(), submissions=(), reviewed=True, uid=""):
    pd = {}
    if recommended:
        pd["recommendedName"] = {"fullName": {"value": recommended}}
    if alternatives:
        pd["alternativeNames"] = [{"fullName": {"value": v}} for v in alternatives]
    if submissions:
        pd["submissionNames"] = [{"fullName": {"value": v}} for v in submissions]
    return {"proteinDescription": pd, "uniProtkbId": uid,
            "entryType": SWISS if reviewed else "UniProtKB unreviewed (TrEMBL)"}


# --- exact_name_bonus: tiered by curation, capped for unreviewed
# "Prosaposin" must beat "Prosaposin receptor GPR37": the first is an exact
# recommended name, the second only shares a token.
check("bonus recommended", exact_name_bonus(entry("Prosaposin"), "Prosaposin"), 110)
check("bonus not-a-match",
      exact_name_bonus(entry("Prosaposin receptor GPR37"), "Prosaposin"), 0)
check("bonus alternative",
      exact_name_bonus(entry("Large ribosomal subunit protein uL30",
                             ["60S ribosomal protein L7"]), "60S ribosomal protein L7"), 80)
check("bonus submission", exact_name_bonus(entry(None, submissions=["Vimentin"]), "Vimentin"), 40)
# An uncurated TrEMBL fragment must not outrank the curated entry.
check("bonus unreviewed capped",
      exact_name_bonus(entry("40S ribosomal protein S3", reviewed=False),
                       "40S ribosomal protein S3"), 80)
check("bonus case-insensitive", exact_name_bonus(entry("Vimentin"), "vimentin"), 110)

check("mnemonic", entry_name_mnemonic({"uniProtkbId": "RL15_HUMAN"}), "RL15")
check("mnemonic no underscore", entry_name_mnemonic({"uniProtkbId": "RL15"}), "RL15")


# --- strategies: a hint is searched as an entry-name mnemonic, and a trailing
# parenthetical is retried without it.
def strategy_names(line):
    return {s.name for s in build_search_strategies(infer_record(line), "9606")}


# --- candidate_rank_key: a supplied identifier outranks a better-scoring entry
def cand(acc, score, reviewed=True, uid=""):
    return {"accession": acc, "entry_name": uid, "heuristic_score": score,
            "reviewed": reviewed}


# Regression: the TrEMBL entry the accession names must win even though the
# reviewed entry scores higher on the description.
trembl = cand("E7EW49", 300, reviewed=False, uid="E7EW49_HUMAN")
swiss = cand("O75122", 380, reviewed=True, uid="CLAP2_HUMAN")
check("identifier outranks score",
      sorted([swiss, trembl], key=lambda i: candidate_rank_key(i, "E7EW49"))[0]["accession"],
      "E7EW49")
# With no identifier, score decides.
check("no identifier -> score wins",
      sorted([trembl, swiss], key=lambda i: candidate_rank_key(i, ""))[0]["accession"],
      "O75122")
# Entry name counts as an identifier match too.
check("entry name matches identifier",
      sorted([swiss, trembl], key=lambda i: candidate_rank_key(i, "E7EW49_HUMAN"))[0]["accession"],
      "E7EW49")
# Equal score: reviewed first, then accession -- deterministic across reruns.
check("tie -> reviewed then accession",
      [c["accession"] for c in sorted(
          [cand("P62701", 230), cand("P15880", 230), cand("B2R491", 230, reviewed=False)],
          key=lambda i: candidate_rank_key(i, ""))],
      ["P15880", "P62701", "B2R491"])

if "hint_entry_name" not in strategy_names("Moesin\tMOES"):
    failures.append("Moesin<TAB>MOES: no hint_entry_name strategy")
if "protein_exact_no_paren" not in strategy_names(
        "Carbohydrate-response element-binding protein (Mondo A)"):
    failures.append("parenthetical name: no protein_exact_no_paren strategy")
if "protein_exact_no_paren" in strategy_names("Vimentin"):
    failures.append("Vimentin: emitted a no-paren strategy with no parenthetical")

if failures:
    print(f"FAIL  {len(failures)} case(s):\n")
    for f in failures:
        print("  -", f)
    sys.exit(1)
print("ok - all parsing cases pass")

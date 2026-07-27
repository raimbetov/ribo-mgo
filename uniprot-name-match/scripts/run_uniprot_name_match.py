#!/usr/bin/env python3
"""
Run an end-to-end UniProt name matching workflow.

This script accepts mixed inputs such as:
- descriptive protein names
- protein name + gene symbol
- current or deprecated-looking UniProt accessions
- UniProt entry names

It writes:
1. a polished final TSV,
2. a review TSV for uncertain cases,
3. optional raw candidate JSONL for auditability.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import requests
import time

from fetch_uniprot_candidates import (
    DEFAULT_ORGANISM_ID,
    InputRecord,
    build_search_strategies,
    candidate_rank_key,
    entry_name,
    extract_gene_names,
    extract_protein_names,
    fetch_results,
    heuristic_score,
    organism_label,
    read_input_records,
)

FINAL_COLUMNS = ["From", "Entry", "Entry Name", "Protein names", "Gene Names", "Organism"]
REVIEW_COLUMNS = [
    "From",
    "Selected Entry",
    "Selected Entry Name",
    "Confidence",
    "Reason",
    "Alternative Entries",
]


def candidate_record(result: dict, record: InputRecord, strategy_name: str) -> dict:
    return {
        "accession": result.get("primaryAccession", ""),
        "entry_name": entry_name(result),
        "protein_names": extract_protein_names(result),
        "gene_names": extract_gene_names(result),
        "organism": organism_label(result),
        "reviewed": result.get("entryType") == "UniProtKB reviewed (Swiss-Prot)",
        "strategy": strategy_name,
        "heuristic_score": heuristic_score(result, record, strategy_name),
    }


def collect_candidates(
    records: list[InputRecord],
    organism_id: str,
    api_size: int,
    candidate_limit: int,
    timeout: float,
    delay: float,
) -> list[dict]:
    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": "codex-uniprot-name-match/1.0",
            "Accept": "application/json",
        }
    )

    cache: dict[tuple[str, str, str], list[dict]] = {}
    payloads: list[dict] = []

    for record in records:
        cache_key = (record.protein_name, record.gene_name, record.identifier)
        if cache_key not in cache:
            candidate_map: dict[str, dict] = {}

            for strategy in build_search_strategies(record, organism_id):
                try:
                    results = fetch_results(
                        session=session,
                        query=strategy.query,
                        size=api_size,
                        timeout=timeout,
                    )
                except requests.RequestException:
                    continue

                for result in results:
                    accession = result.get("primaryAccession", "")
                    if not accession:
                        continue
                    current = candidate_record(result, record, strategy.name)
                    existing = candidate_map.get(accession)
                    if existing is None or current["heuristic_score"] > existing["heuristic_score"]:
                        candidate_map[accession] = current

            ranked = sorted(
                candidate_map.values(),
                key=lambda item: candidate_rank_key(item, record.identifier),
            )
            cache[cache_key] = ranked[:candidate_limit]
            if delay > 0:
                time.sleep(delay)

        payloads.append(
            {
                "input_name": record.display_name,
                "raw_line": record.raw_line,
                "protein_name": record.protein_name,
                "gene_name": record.gene_name,
                "identifier": record.identifier,
                "candidates": cache[cache_key],
            }
        )

    return payloads


def confidence_for(payload: dict) -> tuple[str, str]:
    candidates = payload["candidates"]
    gene_name = payload.get("gene_name", "")
    identifier = payload.get("identifier", "")

    if not candidates:
        return "low", "No UniProt candidates were retrieved."

    best = candidates[0]
    second = candidates[1] if len(candidates) > 1 else None
    best_score = best["heuristic_score"]
    second_score = second["heuristic_score"] if second else -1
    gap = best_score - second_score

    gene_tokens = set(best["gene_names"].upper().split())
    best_mnemonic = best["entry_name"].split("_", 1)[0].upper()
    # A hint agrees if it is one of the entry's gene names OR its entry-name
    # mnemonic -- papers use both ("ARHGDIB" and "MOES" are equally valid ways
    # to point at an entry).
    hint_agrees = bool(gene_name) and (
        gene_name.upper() in gene_tokens or best_mnemonic == gene_name.upper()
    )
    hint_conflicts = bool(gene_name) and not hint_agrees

    if identifier:
        if best["entry_name"].upper() == identifier.upper() or best["accession"].upper() == identifier.upper():
            return "high", "Identifier matched the selected UniProt record directly."
        if best.get("strategy") == "secondary_accession":
            return (
                "high",
                "Input is a secondary/deprecated accession; resolved to the "
                "current UniProt entry that superseded it.",
            )

    if gene_name and hint_agrees and best["reviewed"] and best_score >= 180:
        return "high", "Gene symbol and protein evidence both support the same reviewed entry."

    # The name evidence is deliberately strong, so a hint that points somewhere
    # else has to hold the row back rather than be outvoted. These conflicts are
    # where the papers' own errors surface: Irshad prints "Aspartate
    # aminotransferase, cytoplasmic (AATM)" -- AATM is the *mitochondrial*
    # mnemonic -- and pairs 40S ribosomal protein S7 with the unrelated hint
    # TPD54. Promoting those on name evidence alone would bury a real defect.
    if hint_conflicts:
        return (
            "medium",
            f"Protein name and the hint '{gene_name}' point to different entries "
            f"(selected {best['entry_name']}, genes {best['gene_names'] or 'n/a'}); "
            f"check the source for a naming error.",
        )

    if identifier and best_score >= 180 and gap >= 30:
        return "high", "Identifier-guided search produced a strong leading match."

    if not best["reviewed"]:
        return "medium", "Best candidate is unreviewed, so manual verification is recommended."

    if best_score >= 190 and gap >= 35:
        if identifier:
            return (
                "high",
                "Identifier did not resolve to a live human entry; matched on the "
                "protein name, which selects a reviewed entry decisively.",
            )
        return "high", "Reviewed candidate strongly outranks the alternatives."
    if best_score >= 140 and gap >= 20:
        return "medium", "Best candidate is plausible, but at least one alternative remains credible."

    return "low", "Candidate ranking is ambiguous or only weakly supported."


def select_best_candidate(payload: dict) -> tuple[dict | None, str, str]:
    candidates = payload["candidates"]
    if not candidates:
        return None, "low", "No UniProt candidates were retrieved."

    confidence, reason = confidence_for(payload)
    return candidates[0], confidence, reason


def final_row(input_name: str, candidate: dict | None) -> dict[str, str]:
    if not candidate:
        return {
            "From": input_name,
            "Entry": "",
            "Entry Name": "",
            "Protein names": "",
            "Gene Names": "",
            "Organism": "",
        }

    return {
        "From": input_name,
        "Entry": candidate["accession"],
        "Entry Name": candidate["entry_name"],
        "Protein names": candidate["protein_names"],
        "Gene Names": candidate["gene_names"],
        "Organism": candidate["organism"],
    }


def review_row(
    input_name: str,
    candidate: dict | None,
    confidence: str,
    reason: str,
    candidates: list[dict],
) -> dict[str, str]:
    alternatives = ", ".join(
        f"{item['accession']} ({item['entry_name']})" for item in candidates[1:4]
    )
    return {
        "From": input_name,
        "Selected Entry": "" if not candidate else candidate["accession"],
        "Selected Entry Name": "" if not candidate else candidate["entry_name"],
        "Confidence": confidence,
        "Reason": reason,
        "Alternative Entries": alternatives,
    }


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run end-to-end UniProt name matching.")
    parser.add_argument("input_file", help="Text file with one identifier or clue-set per line")
    parser.add_argument("-o", "--output", required=True, help="Final TSV output path")
    parser.add_argument(
        "--review-output",
        help="Review TSV output path (default: <output stem>.review.tsv)",
    )
    parser.add_argument(
        "--candidates-output",
        help="Optional JSONL path to save raw candidate sets",
    )
    parser.add_argument(
        "--organism-id",
        default=DEFAULT_ORGANISM_ID,
        help=f"NCBI taxonomy ID to search (default: {DEFAULT_ORGANISM_ID})",
    )
    parser.add_argument(
        "--candidate-limit",
        type=int,
        default=5,
        help="Maximum candidates to keep per input (default: 5)",
    )
    parser.add_argument(
        "--api-size",
        type=int,
        default=10,
        help="Maximum results to inspect per strategy (default: 10)",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.2,
        help="Delay in seconds between unique UniProt lookups (default: 0.2)",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="HTTP timeout in seconds (default: 30)",
    )
    parser.add_argument(
        "--accept-confidence",
        choices=["high", "medium"],
        default="high",
        help="Minimum confidence that should be accepted into the final TSV (default: high)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_records = read_input_records(args.input_file)
    if not input_records:
        raise SystemExit(f"No usable identifiers found in {args.input_file}")

    output_path = Path(args.output)
    review_path = Path(args.review_output) if args.review_output else output_path.with_suffix(".review.tsv")

    candidate_payloads = collect_candidates(
        records=input_records,
        organism_id=args.organism_id,
        api_size=args.api_size,
        candidate_limit=args.candidate_limit,
        timeout=args.timeout,
        delay=args.delay,
    )

    accepted_levels = {"high"} if args.accept_confidence == "high" else {"high", "medium"}
    final_rows = []
    review_rows = []

    for payload in candidate_payloads:
        input_name = payload["input_name"]
        candidates = payload["candidates"]
        selected, confidence, reason = select_best_candidate(payload)

        if confidence in accepted_levels and selected is not None:
            final_rows.append(final_row(input_name, selected))
        else:
            final_rows.append(final_row(input_name, None))

        if confidence != "high" or selected is None:
            review_rows.append(review_row(input_name, selected, confidence, reason, candidates))

    write_tsv(output_path, FINAL_COLUMNS, final_rows)
    write_tsv(review_path, REVIEW_COLUMNS, review_rows)

    if args.candidates_output:
        with Path(args.candidates_output).open("w", encoding="utf-8") as handle:
            for payload in candidate_payloads:
                handle.write(json.dumps(payload, ensure_ascii=False) + "\n")

    print(f"Wrote final TSV: {output_path}")
    print(f"Wrote review TSV: {review_path}")
    if args.candidates_output:
        print(f"Wrote candidate JSONL: {args.candidates_output}")

    found = sum(1 for row in final_rows if row["Entry"])
    print(f"Accepted mappings: {found}/{len(final_rows)}")
    print(f"Review cases: {len(review_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

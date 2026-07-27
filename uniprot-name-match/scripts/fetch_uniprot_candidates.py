#!/usr/bin/env python3
"""
Fetch UniProt candidate entries for mixed protein identifiers.

Supported input styles:
- descriptive protein names only
- protein name + gene symbol
- UniProt accessions or entry names
- mixed rows containing several clues separated by tabs, pipes, or semicolons

The output is JSONL with one object per input line.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from dataclasses import dataclass

import requests

BASE_URL = "https://rest.uniprot.org/uniprotkb/search"
DEFAULT_ORGANISM_ID = "9606"

ACCESSION_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$",
    re.IGNORECASE,
)
ENTRY_NAME_RE = re.compile(r"^[A-Z0-9]{1,12}_[A-Z0-9]{2,8}$", re.IGNORECASE)
GENE_SYMBOL_RE = re.compile(r"^[A-Z0-9][A-Z0-9\-]{1,19}$")

# Strips list numbering ("1. Vimentin", "2) Vimentin", "3 - Vimentin") only.
# The separator is REQUIRED and must be followed by whitespace. It used to be
# optional, which meant any leading integer was eaten: "40S ribosomal protein
# S4" -> "S ribosomal protein S4", "14-3-3 protein eta" -> "3-3 protein eta",
# and the entry name "1433Z_HUMAN" -> "Z_HUMAN". Protein names that begin with
# a digit are the norm, not the exception (ribosomal subunits 28S/39S/40S/60S,
# the 14-3-3 family, "78 kDa glucose-regulated protein"), so a bare leading
# number is now left alone -- a missed bullet costs one stray search token, a
# wrongly stripped digit silently destroys the identity of the protein.
NUMBERING_RE = re.compile(r"^\s*\d+\s*(?:[.)]|[-:>]|→)\s+")

# Ribosome/proteasome sedimentation coefficients: 40S, 60S, 28S, 39S, 26S, 5.8S.
SEDIMENTATION_RE = re.compile(r"^\d+(?:\.\d+)?S$")

# Trailing "(synonym)" on a descriptive name.
PARENTHETICAL_TAIL_RE = re.compile(r"\s*\([^()]*\)\s*$")


@dataclass(frozen=True)
class SearchStrategy:
    name: str
    query: str


@dataclass(frozen=True)
class InputRecord:
    raw_line: str
    display_name: str
    protein_name: str
    gene_name: str
    identifier: str


def clean_name(name: str) -> str:
    name = name.strip()
    name = NUMBERING_RE.sub("", name)
    name = re.sub(r"^[*•\-]+\s*", "", name)
    name = re.sub(r"\s+", " ", name)
    return name.strip(" ;,")


def tokenize(text: str) -> list[str]:
    return re.findall(r"[A-Za-z0-9]+", text.lower())


def looks_like_gene_symbol(value: str) -> bool:
    compact = value.replace("_", "")
    if SEDIMENTATION_RE.fullmatch(compact):
        # "40S", "60S", "28S", "5.8S" pass the gene-symbol shape test but are
        # sedimentation coefficients. Treating one as a gene symbol splits
        # "40S ribosomal protein S4" into gene "40S" + protein "ribosomal
        # protein S4", which loses the subunit the name is about.
        return False
    return bool(GENE_SYMBOL_RE.fullmatch(compact)) and any(ch.isalpha() for ch in compact)


def is_explicit_gene_symbol(value: str) -> bool:
    stripped = value.strip()
    return stripped == stripped.upper() and " " not in stripped and looks_like_gene_symbol(stripped)


def looks_like_accession(value: str) -> bool:
    return bool(ACCESSION_RE.fullmatch(value))


def looks_like_entry_name(value: str) -> bool:
    return bool(ENTRY_NAME_RE.fullmatch(value))


def normalize_identifier(value: str) -> str:
    return value.strip().upper()


def split_fields(line: str) -> list[str]:
    if "\t" in line:
        parts = [clean_name(part) for part in line.split("\t")]
    elif " | " in line:
        parts = [clean_name(part) for part in line.split(" | ")]
    elif ";" in line:
        parts = [clean_name(part) for part in line.split(";")]
    else:
        parts = [clean_name(line)]
    return [part for part in parts if part]


def is_inline_gene_symbol(value: str) -> bool:
    """Stricter than is_explicit_gene_symbol, for splitting a space-separated
    line into gene + name.

    Here a false positive silently rewrites the protein name: "40S ribosomal
    protein S4" would become gene "S4" + name "40S ribosomal protein", and the
    search then ranks ribosomal proteins by token overlap with a query that no
    longer says which one. Requiring >=3 characters and >=2 letters keeps real
    labels ("THRA Thyroid hormone receptor alpha") and rejects the subunit
    suffixes that make up the tail of a ribosomal protein name (S4, L7, S28).
    """
    stripped = value.strip()
    if not is_explicit_gene_symbol(stripped):
        return False
    compact = stripped.replace("_", "").replace("-", "")
    return len(compact) >= 3 and sum(ch.isalpha() for ch in compact) >= 2


def infer_space_separated_fields(line: str) -> tuple[str, str]:
    parts = line.split(" ", 1)
    if len(parts) == 2 and is_inline_gene_symbol(parts[0]) and parts[1].strip():
        return parts[1].strip(), normalize_identifier(parts[0])

    parts = line.rsplit(" ", 1)
    if len(parts) == 2 and is_inline_gene_symbol(parts[1]) and parts[0].strip():
        return parts[0].strip(), normalize_identifier(parts[1])

    return "", ""


def infer_record(raw_line: str) -> InputRecord:
    parts = split_fields(raw_line)
    protein_name = ""
    gene_name = ""
    identifier = ""

    if len(parts) == 1:
        inferred_protein, inferred_gene = infer_space_separated_fields(parts[0])
        if inferred_protein:
            protein_name = inferred_protein
            gene_name = inferred_gene
        else:
            only = parts[0]
            only_upper = normalize_identifier(only)
            # A line that is nothing but a gene symbol IS a gene symbol. Without
            # this it falls through to `protein_name` below, leaving `gene_name`
            # empty -- and the gene-agreement rule in confidence_for() can then
            # never fire, so a correctly retrieved entry is gated to `low`.
            if (
                is_explicit_gene_symbol(only)
                and not looks_like_accession(only_upper)
                and not looks_like_entry_name(only_upper)
            ):
                gene_name = only_upper

    if parts and not protein_name and not gene_name:
        first = parts[0]
        first_upper = normalize_identifier(first)
        if looks_like_accession(first_upper) or looks_like_entry_name(first_upper):
            identifier = first_upper
            if len(parts) > 1:
                protein_name = parts[1]
        elif len(parts) > 1 and is_explicit_gene_symbol(first) and not is_explicit_gene_symbol(parts[1]):
            gene_name = first_upper
            protein_name = parts[1]
        else:
            protein_name = first

    for index, part in enumerate(parts):
        upper = normalize_identifier(part)
        # Skip any field already consumed. This used to apply to index 0 only,
        # so the field that had become `protein_name` was reconsidered further
        # down the chain -- "P08670<TAB>Vimentin" ended up with gene "VIMENTIN".
        if protein_name == part or gene_name == upper or identifier == upper:
            continue
        if not identifier and (looks_like_accession(upper) or looks_like_entry_name(upper)):
            identifier = upper
            continue
        if not gene_name and is_explicit_gene_symbol(part) and len(parts) > 1:
            gene_name = upper
            continue
        if not protein_name:
            protein_name = part
        elif not gene_name and is_explicit_gene_symbol(part):
            gene_name = upper
        elif not identifier and (looks_like_accession(upper) or looks_like_entry_name(upper)):
            identifier = upper
        elif not gene_name and " " not in part.strip() and looks_like_gene_symbol(upper):
            # A single-token second field is a hint the source deliberately
            # paired with the name, even when it is not upper-case: Irshad
            # prints "RhoGDI2", Alhujaily "C19orf68". is_explicit_gene_symbol()
            # rejects those on letter-case alone, and before this branch the
            # clue fell off the end of the loop -- never searched, never
            # scored, and nothing said so.
            gene_name = upper

    if not protein_name:
        if identifier:
            protein_name = identifier
        elif gene_name:
            protein_name = gene_name
        else:
            protein_name = clean_name(raw_line)

    return InputRecord(
        raw_line=raw_line,
        display_name=raw_line,
        protein_name=protein_name,
        gene_name=gene_name,
        identifier=identifier,
    )


def build_search_strategies(record: InputRecord, organism_id: str) -> list[SearchStrategy]:
    protein_name = record.protein_name
    escaped_name = protein_name.replace('"', "")
    gene_name = record.gene_name
    identifier = record.identifier
    tokens = [token for token in tokenize(protein_name) if len(token) > 1]
    strategies: list[SearchStrategy] = []

    if identifier:
        strategies.append(
            SearchStrategy(
                "identifier_exact",
                f'("{identifier}") AND organism_id:{organism_id}',
            )
        )
        if looks_like_accession(identifier):
            strategies.append(
                SearchStrategy(
                    "accession_exact",
                    f"accession:{identifier} AND organism_id:{organism_id}",
                )
            )
            # `accession:` does NOT match deprecated or merged accessions --
            # UniProt exposes those only through `sec_acc:`. Without this a stale
            # ID returns no candidates at all, which is indistinguishable from a
            # genuinely unmappable input.
            strategies.append(
                SearchStrategy(
                    "secondary_accession",
                    f"sec_acc:{identifier} AND organism_id:{organism_id}",
                )
            )
        if looks_like_entry_name(identifier):
            strategies.append(
                SearchStrategy(
                    "entry_name_exact",
                    f"id:{identifier} AND organism_id:{organism_id}",
                )
            )

    if gene_name:
        strategies.append(
            SearchStrategy(
                "gene_symbol",
                f'(gene:{gene_name} OR gene_exact:{gene_name}) AND organism_id:{organism_id}',
            )
        )
        # A paper's parenthetical hint is often the UniProt entry-name mnemonic
        # rather than a gene symbol -- Irshad prints "Moesin (MOES)",
        # "Glyceraldehyde-3-phosphate dehydrogenase (G3P)", "Ribosomal protein
        # L15 (RL15)". Those are MOES_HUMAN, G3P_HUMAN, RL15_HUMAN: the exact
        # answers. `gene:` cannot see them, so without this the strongest clue
        # on the line went unused and the name alone chose the entry -- which
        # picked NHRF1 (ezrin-radixin-moesin-binding) for "Moesin".
        strategies.append(
            SearchStrategy(
                "hint_entry_name",
                f"id:{gene_name}_* AND organism_id:{organism_id}",
            )
        )

    if protein_name and is_explicit_gene_symbol(protein_name) and not gene_name:
        compact = normalize_identifier(protein_name)
        strategies.append(
            SearchStrategy(
                "name_as_gene_symbol",
                f'(gene:{compact} OR gene_exact:{compact}) AND organism_id:{organism_id}',
            )
        )

    if protein_name:
        strategies.append(
            SearchStrategy(
                "protein_exact",
                f'protein_name:"{escaped_name}" AND organism_id:{organism_id}',
            )
        )
        strategies.append(
            SearchStrategy(
                "all_fields_exact",
                f'"{escaped_name}" AND organism_id:{organism_id}',
            )
        )
        # Retry without a trailing parenthetical. Supplements often append a
        # synonym the way Alhujaily prints "Carbohydrate-response
        # element-binding protein (Mondo A)", and the full string matches
        # nothing at all -- that row returned zero candidates.
        bare = PARENTHETICAL_TAIL_RE.sub("", escaped_name).strip()
        if bare and bare != escaped_name:
            strategies.append(
                SearchStrategy(
                    "protein_exact_no_paren",
                    f'protein_name:"{bare}" AND organism_id:{organism_id}',
                )
            )

    if protein_name and gene_name:
        strategies.append(
            SearchStrategy(
                "protein_plus_gene",
                f'(protein_name:"{escaped_name}") AND (gene:{gene_name} OR gene_exact:{gene_name}) AND organism_id:{organism_id}',
            )
        )

    if tokens:
        token_query = " AND ".join(f'"{token}"' for token in tokens)
        strategies.append(
            SearchStrategy(
                "token_and",
                f"({token_query}) AND organism_id:{organism_id}",
            )
        )

    deduped: list[SearchStrategy] = []
    seen_queries: set[str] = set()
    for strategy in strategies:
        if strategy.query not in seen_queries:
            deduped.append(strategy)
            seen_queries.add(strategy.query)
    return deduped


def format_name_block(name_data: dict) -> str:
    if not name_data:
        return ""

    parts: list[str] = []
    full_name = name_data.get("fullName", {}).get("value")
    if full_name:
        parts.append(full_name)

    for short_name in name_data.get("shortNames", []):
        value = short_name.get("value")
        if value:
            parts.append(f"({value})")

    for ec_number in name_data.get("ecNumbers", []):
        value = ec_number.get("value")
        if value:
            parts.append(f"(EC {value})")

    return " ".join(parts).strip()


def extract_protein_names(result: dict) -> str:
    protein_desc = result.get("proteinDescription", {})
    names: list[str] = []

    recommended = format_name_block(protein_desc.get("recommendedName", {}))
    if recommended:
        names.append(recommended)

    for alt_name in protein_desc.get("alternativeNames", []):
        formatted = format_name_block(alt_name)
        if formatted:
            if names:
                names.append(f"({formatted})")
            else:
                names.append(formatted)

    # TrEMBL entries carry `submissionNames` and no recommendedName. Without
    # this they come back with an empty name string, so token overlap scores
    # zero and the entry is effectively invisible to name-based matching.
    for sub_name in protein_desc.get("submissionNames", []):
        formatted = format_name_block(sub_name)
        if formatted:
            names.append(f"({formatted})" if names else formatted)

    return " ".join(names)


def name_components(result: dict) -> tuple[list[str], list[str], list[str]]:
    """Raw fullName values split by curation tier: (recommended, alternative,
    submission).

    Kept separate from extract_protein_names() because that function returns
    one concatenated blob ("Recommended (Alt1) (Alt2)") plus short names and EC
    numbers. An exact-name test against the blob almost never fires for an
    entry that has alternative names -- which is how "Prosaposin" tied
    "Prosaposin receptor GPR37" on generic token overlap alone.
    """
    protein_desc = result.get("proteinDescription", {})
    recommended: list[str] = []
    alternative: list[str] = []
    submission: list[str] = []

    value = protein_desc.get("recommendedName", {}).get("fullName", {}).get("value")
    if value:
        recommended.append(value)
    for alt_name in protein_desc.get("alternativeNames", []):
        value = alt_name.get("fullName", {}).get("value")
        if value:
            alternative.append(value)
    for sub_name in protein_desc.get("submissionNames", []):
        value = sub_name.get("fullName", {}).get("value")
        if value:
            submission.append(value)
    return recommended, alternative, submission


def exact_name_bonus(result: dict, protein_name: str) -> int:
    """Score an exact match against one individual protein name."""
    query = protein_name.strip().lower()
    if not query:
        return 0
    recommended, alternative, submission = name_components(result)
    if query in [name.strip().lower() for name in recommended]:
        bonus = 110
    elif query in [name.strip().lower() for name in alternative]:
        bonus = 80
    elif query in [name.strip().lower() for name in submission]:
        bonus = 40
    else:
        return 0
    # An exact hit on an uncurated TrEMBL name is weaker evidence than the same
    # hit on a curated one: UniProt holds many near-duplicate TrEMBL fragments
    # per gene, all auto-named from the same descriptive convention. Uncapped,
    # such a fragment ties the Swiss-Prot entry and drags a correct,
    # unambiguous row into review.
    if result.get("entryType") != "UniProtKB reviewed (Swiss-Prot)":
        bonus = min(bonus, 80)
    return bonus


def entry_name_mnemonic(result: dict) -> str:
    """`RL15_HUMAN` -> `RL15`."""
    return entry_name(result).split("_", 1)[0].upper()


def candidate_rank_key(item: dict, identifier: str):
    """Sort key for candidates: identifier match, score, reviewed, accession.

    An identifier the input supplied outranks score. Score alone is not enough:
    a reviewed entry whose recommended name matches the description can
    out-score the very TrEMBL entry the accession names -- E7EW49
    "CLIP-associating protein 2" loses to CLAP2_HUMAN -- and then the accession
    the paper printed is silently discarded. The descriptive name is the
    fallback for when an identifier matches nothing, never an override for when
    it matches.

    Accession is the final tiebreak so a rerun selects the same entry. Genuine
    ties happen: "40S ribosomal protein S4" is an exact alternative name of both
    P15880 (RPS2, historical RPS4/LLRep3 naming) and P62701 (RPS4X).
    """
    ident = (identifier or "").upper()
    exact_identifier = bool(ident) and ident in (
        item["accession"].upper(),
        item["entry_name"].upper(),
    )
    return (
        not exact_identifier,
        -item["heuristic_score"],
        not item["reviewed"],
        item["accession"],
    )


def extract_gene_names(result: dict) -> str:
    names = []
    for gene in result.get("genes", []):
        primary = gene.get("geneName", {}).get("value")
        if primary and primary not in names:
            names.append(primary)
        for synonym in gene.get("synonyms", []):
            value = synonym.get("value")
            if value and value not in names:
                names.append(value)
    return " ".join(names)


def entry_name(result: dict) -> str:
    return result.get("uniProtkbId", "")


def organism_label(result: dict) -> str:
    organism = result.get("organism", {})
    scientific = organism.get("scientificName", "")
    common = organism.get("commonName", "")
    if scientific and common:
        return f"{scientific} ({common})"
    return scientific


def heuristic_score(result: dict, record: InputRecord, strategy_name: str) -> int:
    score = 0
    candidate_entry_name = entry_name(result)
    candidate_gene_names = extract_gene_names(result)
    candidate_protein_names = extract_protein_names(result)
    candidate_text = " ".join(
        [candidate_entry_name, candidate_gene_names, candidate_protein_names]
    ).lower()
    candidate_tokens = set(tokenize(candidate_text))
    query_tokens = set(tokenize(record.protein_name))

    if result.get("entryType") == "UniProtKB reviewed (Swiss-Prot)":
        score += 100

    if strategy_name in {
        "identifier_exact",
        "accession_exact",
        "entry_name_exact",
        "secondary_accession",
    }:
        score += 70

    if record.identifier:
        identifier = record.identifier.upper()
        if candidate_entry_name.upper() == identifier:
            score += 140
        if result.get("primaryAccession", "").upper() == identifier:
            score += 160
        # A sec_acc hit is exact identifier evidence too: the input accession is
        # this entry's own retired ID. The candidate's primary accession will not
        # equal the input, so it earns the same weight here instead.
        elif strategy_name == "secondary_accession":
            score += 160

    if record.gene_name:
        primary_gene = ""
        genes = result.get("genes", [])
        if genes:
            primary_gene = genes[0].get("geneName", {}).get("value", "")
        if primary_gene.upper() == record.gene_name.upper():
            score += 120
        elif record.gene_name.upper() in candidate_gene_names.upper().split():
            score += 70
        elif entry_name_mnemonic(result) == record.gene_name.upper():
            # The hint is this entry's own mnemonic (MOES -> MOES_HUMAN). That
            # is identifier-grade evidence, not a lexical coincidence, so it
            # earns the same weight as a primary-gene match.
            score += 120

    if candidate_entry_name.upper() == record.protein_name.upper():
        score += 60

    score += exact_name_bonus(result, record.protein_name)

    overlap = len(query_tokens & candidate_tokens)
    score += overlap * 10
    if query_tokens and query_tokens.issubset(candidate_tokens):
        score += 20

    return score


def fetch_results(
    session: requests.Session,
    query: str,
    size: int,
    timeout: float,
) -> list[dict]:
    response = session.get(
        BASE_URL,
        params={
            "query": query,
            "format": "json",
            "fields": "accession,id,protein_name,gene_names,organism_name,reviewed",
            "size": size,
        },
        timeout=timeout,
    )
    response.raise_for_status()
    return response.json().get("results", [])


def read_input_records(path: str) -> list[InputRecord]:
    records: list[InputRecord] = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            record = infer_record(line)
            if record.protein_name or record.gene_name or record.identifier:
                records.append(record)
    return records


def read_names(path: str) -> list[str]:
    return [record.display_name for record in read_input_records(path)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch UniProt candidates for mixed protein inputs.")
    parser.add_argument("input_file", help="Text file with one identifier or clue-set per line")
    parser.add_argument("-o", "--output", required=True, help="Output JSONL path")
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
        help="Maximum results to inspect per search strategy (default: 10)",
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
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    records = read_input_records(args.input_file)
    if not records:
        print(f"No usable identifiers found in {args.input_file}", file=sys.stderr)
        return 1

    cache: dict[tuple[str, str, str], list[dict]] = {}
    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": "codex-uniprot-name-match/1.0",
            "Accept": "application/json",
        }
    )

    with open(args.output, "w", encoding="utf-8") as handle:
        for record in records:
            cache_key = (record.protein_name, record.gene_name, record.identifier)
            if cache_key not in cache:
                candidate_map: dict[str, dict] = {}

                for strategy in build_search_strategies(record, args.organism_id):
                    try:
                        results = fetch_results(
                            session=session,
                            query=strategy.query,
                            size=args.api_size,
                            timeout=args.timeout,
                        )
                    except requests.RequestException as exc:
                        print(
                            f"Error searching for '{record.display_name}' with strategy '{strategy.name}': {exc}",
                            file=sys.stderr,
                        )
                        continue

                    for result in results:
                        accession = result.get("primaryAccession", "")
                        if not accession:
                            continue

                        score = heuristic_score(result, record, strategy.name)
                        existing = candidate_map.get(accession)
                        if existing is None or score > existing["heuristic_score"]:
                            candidate_map[accession] = {
                                "accession": accession,
                                "entry_name": entry_name(result),
                                "protein_names": extract_protein_names(result),
                                "gene_names": extract_gene_names(result),
                                "organism": organism_label(result),
                                "reviewed": result.get("entryType") == "UniProtKB reviewed (Swiss-Prot)",
                                "strategy": strategy.name,
                                "heuristic_score": score,
                            }

                ranked = sorted(
                    candidate_map.values(),
                    key=lambda item: candidate_rank_key(item, record.identifier),
                )
                cache[cache_key] = ranked[: args.candidate_limit]
                if args.delay > 0:
                    time.sleep(args.delay)

            payload = {
                "input_name": record.display_name,
                "raw_line": record.raw_line,
                "protein_name": record.protein_name,
                "gene_name": record.gene_name,
                "identifier": record.identifier,
                "candidates": cache[cache_key],
            }
            handle.write(json.dumps(payload, ensure_ascii=False) + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

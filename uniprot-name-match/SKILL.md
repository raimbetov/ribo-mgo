---
name: uniprot-name-match
description: "Map mixed protein identifiers to up-to-date UniProt entries using deterministic candidate retrieval plus confidence-aware selection. Use when Codex needs to convert descriptive protein names, protein name plus gene symbol rows, UniProt accessions, UniProt entry names, or deprecated-looking identifiers from text, tables, papers, or supplements into reviewable UniProt mappings."
---

# UniProt Name Match

## Overview

Use this skill to turn mixed protein identifiers into reliable UniProt mappings without asking the model to invent IDs. The skill is self-contained: it has a single runner that produces a polished TSV plus a review TSV, and it also exposes lower-level candidate retrieval when a user wants to inspect ambiguity more closely.

## Contract

**In** — any one of these, alone or paired on one line:

| Input | Example |
|---|---|
| UniProt accession | `P08670` |
| UniProt entry name / ID | `VIME_HUMAN` |
| Gene symbol | `VIM` |
| Protein name or long description | `Small ribosomal subunit protein uS11 (40S ribosomal protein S14)` |

**Out** — the *current* UniProt entry. A deprecated or merged accession resolves to the entry that superseded it (via `sec_acc:`), so the answer is up to date rather than a dead ID echoed back.

**Fallback** — when the identifier does not match anything, the descriptive name on the same line is used to find the best match. So pair an identifier with its name whenever the source has one: `P08670<TAB>Vimentin` survives the accession being retired, while `P08670` alone has nothing to fall back to.

An input with no identifier *and* no usable name is reported in the review file with no candidates. That is the one case the skill cannot answer, and it says so rather than guessing.

## Workflow

1. Parse each line into available clues.
2. Retrieve candidate UniProt entries.
3. Score and select the best candidate.
4. Auto-accept only confident matches by default.
5. Emit a polished final TSV and a separate review file for uncertain cases.

## Input Expectations

- Prefer a plain text file with one clue-set per line.
- Accept these common styles:
  `Vimentin`
  `Vimentin<TAB>VIM`
  `VIM<TAB>Vimentin`
  `Thyroid hormone receptor alpha | THRA`
  `THRA Thyroid hormone receptor alpha`
  `P08670`
  `VIME_HUMAN`
- Preserve the original input label in the final output under `From`. `From` is the raw input line, before any cleaning.
- Remove obvious bullets, numbering, and surrounding punctuation, but do not silently rewrite biological meaning. A leading number is stripped **only** when a separator and a space follow it (`1. Vimentin`, `2) Vimentin`, `3 - Vimentin`). A bare leading number is part of the name and is kept — `40S ribosomal protein S4`, `14-3-3 protein eta`, `78 kDa glucose-regulated protein`, and the entry name `1433Z_HUMAN` all begin with digits.
- Do not read a sedimentation coefficient (`40S`, `60S`, `28S`, `5.8S`) as a gene symbol.
- Splitting a space-separated line into gene + name requires the gene-looking token to have at least 3 characters and 2 letters. `THRA Thyroid hormone receptor alpha` splits; `40S ribosomal protein S4` does not, because its trailing `S4` is part of the name, not a gene.
- Handle crude copy-pastes with minimal preprocessing, including both `protein<TAB>gene` and `gene<TAB>protein`.
- Use gene symbols and identifier-like fields as stronger evidence than descriptive name tokens when they are present.
- A line containing nothing but a gene symbol is treated as a gene symbol, not as a descriptive name, so it can reach `high` on gene agreement.
- Strip isoform suffixes before mapping: UniProt does not resolve `P04075-2`, so pass `P04075` and keep the isoform number in your own source column.
- Assume human proteins by default unless the user specifies another organism.

## Candidate Retrieval

For a full run, use the end-to-end runner:

```bash
python3 uniprot-name-match/scripts/run_uniprot_name_match.py input.txt -o output.tsv
```

This runner:
- queries UniProt with several search strategies,
- handles mixed rows containing names, gene symbols, accessions, or entry names,
- gathers multiple candidate records per input,
- scores the candidates deterministically,
- accepts only confident mappings by default,
- writes a final TSV and a review TSV,
- can optionally save raw candidates as JSONL for auditability.

Use the lower-level fetcher only when you need raw candidate sets:

```bash
python3 uniprot-name-match/scripts/fetch_uniprot_candidates.py input.txt -o candidates.jsonl
```

If the user wants a different organism, pass `--organism-id` to either script.

## Selection Rules

Read [references/matching-policy.md](./references/matching-policy.md) before resolving ambiguous cases. Apply these rules:

- Choose only from returned UniProt candidates.
- Treat explicit UniProt accessions, UniProt entry names, and gene symbols as stronger evidence than descriptive-name token overlap.
- A parenthetical hint may be a gene symbol **or** an entry-name mnemonic. `Moesin (MOES)` means `MOES_HUMAN`; `gene:` cannot see it, so the mnemonic is searched and scored as identifier-grade evidence.
- Score an exact match against each protein name individually, not against the concatenated name string, and rank a curated name above an uncurated TrEMBL submission name. `Prosaposin` must select P07602 (whose recommended name is exactly that) over O15354, "Prosaposin **receptor** GPR37", which merely shares the token.
- When a hint disagrees with the protein name, hold the row in review and say so. The disagreement is usually a defect in the source: Irshad prints `Aspartate aminotransferase, cytoplasmic (AATM)` — `AATM` is the *mitochondrial* mnemonic — and pairs `40S ribosomal protein S7` with the unrelated hint `TPD54`.
- Prefer reviewed Swiss-Prot human entries when the biological meaning matches.
- Prefer semantic match over raw token overlap.
- If an old or deprecated-looking UniProt identifier resolves to a current entry via search, prefer the current primary accession in the final TSV.
- Penalize pseudogenes, fragments, uncharacterized proteins, and unrelated family members unless the input strongly indicates them.
- Flag low-confidence or near-tie matches for review instead of forcing certainty.
- Exact ties are ranked by score, then reviewed status, then accession, so a rerun selects the same entry. A tie is a real outcome, not noise: `40S ribosomal protein S4` is an exact alternative name of both P15880 (RPS2, via the historical `RPS4`/LLRep3 naming) and P62701 (RPS4X). Both score identically and the case belongs in review for a human to settle.

## Output Contract

Produce two files when doing a full mapping task:

- Final TSV for accepted mappings with columns:
  `From`, `Entry`, `Entry Name`, `Protein names`, `Gene Names`, `Organism`
- Review TSV or JSON for uncertain mappings with fields such as:
  `From`, `Selected Entry`, `Confidence`, `Reason`, `Alternatives`

Keep unresolved items in the review output. Do not silently drop them.

## LLM Use

The runner works on its own, but an agent may still help with edge cases. Use the model for review and explanation, not for freeform ID generation. A strong pattern is:

1. Fetch candidates.
2. For each input, compare the top candidates.
3. Return a structured decision with:
   `selected_accession`, `confidence`, `reason`, `needs_review`, `alternatives`
4. Auto-accept only high-confidence decisions.

## Quick Checks

- If the chosen accession is not in the candidate list, stop and correct the workflow.
- If two candidates are biologically plausible, mark the case for review.
- If the input carries a gene symbol or UniProt-like identifier, use that clue to boost confidence and narrow the search.
- If the input names come from a paper supplement, preserve their original spelling in `From` even when the chosen UniProt entry uses a cleaner synonym.

## Entry Points

- `scripts/run_uniprot_name_match.py`: End-to-end workflow from raw names to final and review TSVs.
- `scripts/fetch_uniprot_candidates.py`: Candidate retrieval only.

## Resources

- `scripts/run_uniprot_name_match.py`: Produce final TSV, review TSV, and optional raw candidate JSONL.
- `scripts/fetch_uniprot_candidates.py`: Retrieve candidate UniProt entries for each protein name.
- `references/matching-policy.md`: Confidence rubric and ambiguity-handling guidance.

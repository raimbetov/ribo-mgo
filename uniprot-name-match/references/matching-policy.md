# Matching Policy

Use this file when the best UniProt entry is not obvious from the top candidate list.

## Strong Match

Treat a candidate as strong when most of these are true:
- The recommended protein name clearly matches the biological meaning of the input.
- Important synonyms or short names line up with the input wording.
- The candidate is a reviewed Swiss-Prot human entry.
- The gene name and entry name are consistent with the expected protein family.
- Any supplied gene symbol, accession, or UniProt entry name points to the same record.

## Weak Match

Treat a candidate as weak when any of these are true:
- The overlap is only lexical and not biological.
- The result is a pseudogene, fragment, isoform-only record, or uncharacterized protein without strong evidence.
- The input could reasonably refer to several family members.
- The best candidate depends on assumptions not present in the source text.
- A supplied gene symbol or identifier conflicts with the descriptive-name match.

## Confidence Rubric

- `high`: clear semantic match, strong UniProt candidate, no close runner-up
- `medium`: plausible best match, but one alternative remains credible
- `low`: ambiguous, weakly supported, or likely to require manual review

Auto-accept only `high` confidence by default. Put `medium` and `low` into the review output unless the user explicitly asks for aggressive best-effort mapping.

## Tie-Breaking

When two candidates are both plausible:
- Prefer a candidate that satisfies a supplied accession, UniProt entry name, or gene symbol.
- Prefer the reviewed human entry.
- Prefer the candidate whose recommended name matches the descriptive core of the input.
- Prefer broader canonical entries over fragment-like or poorly annotated records.
- Prefer entries whose synonyms explain the exact wording in the source input.

## Hard Rule

Never invent a UniProt accession. Select only from retrieved candidates.

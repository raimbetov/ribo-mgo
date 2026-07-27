# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Research repository for a project on methylglyoxal (MGO) and translation/ribosomes. Two strands of work:

1. **Literature data mining** — collecting MGO-modified protein lists from published proteomics papers (Donnellan 2022, Irshad 2019, Ashour 2020, Zheng 2024, Alhujaily 2021), normalizing their heterogeneous protein identifiers to current UniProt entries, then intersecting against the human cytoplasmic ribosomal protein set.
2. **Wet-lab data analysis** — Jupyter notebooks analyzing dual-luciferase (Rluc/Fluc) translation-fidelity assays, e.g. [dual-luciferase-analysis.ipynb](dual-luciferase-analysis.ipynb) (WT vs. H245R/D357X Fluc reporters ± MGO).

**Sun et al. 2019 is not an MGO paper** and is deliberately excluded, despite sitting in `00-data-raw/`. Its Methods search exactly one glycation delta-mass, `+162.0528 Da` — a hexose, i.e. glucose-derived fructosyl-lysine. No methylglyoxal adduct mass (MG-H1 +54.0106, CEL/CEA +72.0211) appears anywhere in its search space. See [mgo-adduct-extract/references/inclusion-policy.md](mgo-adduct-extract/references/inclusion-policy.md).

Source data is gitignored deliberately: `00-data-raw/`, `00-full-text/`, `rp-script/` and the venv are local-only. Do not "fix" this by adding them to git. The one exception is `01-data-extracted/`, which **is** tracked — it is small, derived, and byte-reproducible from `mgo-adduct-extract/`, so keeping it in history makes the extraction auditable.

## Environment

```bash
source ribo-mgo-env/bin/activate   # python 3.12 venv, gitignored, holds requests + jupyter + pandas/scipy/seaborn
```

There is no package manifest, test suite, or linter. Scripts are standalone `python3` files run directly.

## Directory pipeline

Numbered prefixes encode pipeline stage:

- `00-data-raw/<Author et al. - YEAR>/` — source PDFs, supplements (xlsx/zip) as downloaded
- `00-full-text/` — the papers themselves, one PDF per `<Author et al. - YEAR>`
- `01-data-extracted/<Author et al. YEAR> - <CellLine>.tsv` — **raw, pre-UniProt-mapping** extractions: a verbatim mirror of one published table, its columns exactly as printed, plus `source_file` / `source_table` / `source_row` provenance on every row. One file per paper × cell line × experimental system (in-cell and cell-free are kept apart). Built and checked by `mgo-adduct-extract/`; [01-data-extracted/MANIFEST.md](01-data-extracted/MANIFEST.md) records per-file provenance, coverage gaps and published errata.
- `02-data-processed/` — working stage for UniProt mapping: `uniprot-input/` (clue-sets built from `01-data-extracted/`), `uniprot-mapped/` (machine mapping + per-row review files + candidate JSONL), and `uniprot-mapped/user-reviewed/` (hand-resolved review decisions)
- `03-data-final/` — one `<key>.uniprot.final.tsv` per dataset: the completed, UniProt-mapped protein lists that the ribosomal intersection is computed from. Columns `From`, `Entry`, `Entry Name`, `Protein names`, `Gene Names`, `Organism`, `Source` (`auto`|`user-reviewed`), `Notes`. Built by `apply_user_reviews.py`; do not hand-edit — change the review file and rebuild
- `04-data-merged/mgo-glycated-proteins.tsv` — the canonical list for downstream analysis: one row per protein, keyed on a **canonical (reviewed, where one exists) accession**, carrying the provenance of every dataset that reported it (study, cell line, system, evidence class, and the label each paper printed). `accession-canonicalization.tsv` beside it logs every substitution. Built by `merge_final_datasets.py`
- `rp-script/` — canonical human ribosomal protein reference (85 entries: gene symbol, UniProt accession, Ban et al. 2014 unified nomenclature `uS/uL/eS/eL`). `HUMAN_RIBOSOMAL_PROTEINS_COMPLETE_LIST.md` documents special cases (ubiquitin-fusion RPS27A/FAU/UBA52, RPS4X/Y1/Y2, RACK1) that break naive gene-symbol joins. Regenerate with `python3 rp-script/fetch_ribosomal_proteins.py`.

## mgo-adduct-extract

A self-contained agent skill ([mgo-adduct-extract/SKILL.md](mgo-adduct-extract/SKILL.md)) that is the workhorse for stage 00 → 01: it turns published supplements into the TSVs in `01-data-extracted/`.

```bash
cd mgo-adduct-extract && pip install -r requirements.txt   # pdfplumber, openpyxl, PyYAML
python3 scripts/extract.py --all      # rebuild 01-data-extracted/
python3 scripts/verify.py             # structural + count checks
```

The work splits in two, and the split is the point:

- **Scoping** (judgment, not automatable): decide which table in a paper qualifies. Papers routinely publish adduct lists, abundance-change lists and hotspot *predictions* side by side, and only the first counts. [references/inclusion-policy.md](mgo-adduct-extract/references/inclusion-policy.md) gives the rule, the procedure, and four worked cases — including why Sun 2019 is out and why Ashour's "proteins at risk of MG modification" is in.
- **Extraction** (deterministic): one YAML spec per dataset in `datasets/`, holding the source table, page range, column geometry, expected row count *with the paper's own sentence as evidence*, the adduct masses searched, known publication errata, and acceptance checks. `scripts/extract.py` carries no per-paper knowledge.

Nothing is written for a dataset whose checks fail. Checks prefer evidence the paper itself supplies — Irshad's Table S3 prints `TOTAL: 220` and `TOTAL 411`, and the spec asserts the declared counts and the parsed site tokens each reach 411 independently. Known publication defects are pinned by row (`expect_offenders: [57, 72]`) so a *new* mismatch fails the build rather than hiding behind a blanket allowance.

To add a paper: scope it, run `scripts/measure_columns.py` to read the table geometry off an x-position histogram, write a spec, build, verify.

## uniprot-name-match

A self-contained agent skill ([uniprot-name-match/SKILL.md](uniprot-name-match/SKILL.md)) for stage 01 → 02. It maps mixed identifiers (accessions, entry names, gene symbols, descriptive names, deprecated IDs) to **current** UniProt entries without letting a model invent accessions. A retired accession resolves to the entry that superseded it; when an identifier matches nothing, the name on the same line is used instead.

Build its inputs from `01-data-extracted/` first — each paper identifies proteins differently, and the per-dataset rules (isoform stripping, splitting Irshad's `Name (HINT)`, pairing accession with description) live in one place:

```bash
python3 mgo-adduct-extract/scripts/make_uniprot_input.py   # -> 02-data-processed/uniprot-input/
```

Then map every dataset and QA the result. `summarize_uniprot_mapping.py` checks the mapper's contract (one output row per input, in order, `From` verbatim, well-formed accessions) and prints the cytoplasmic-ribosomal intersection per dataset — the project's actual question, and the earliest smoke test that the right tables were extracted:

```bash
for f in 02-data-processed/uniprot-input/*.txt; do
  python3 uniprot-name-match/scripts/run_uniprot_name_match.py "$f" \
    -o "02-data-processed/uniprot-mapped/$(basename "${f%.txt}").uniprot.tsv"
done
python3 mgo-adduct-extract/scripts/summarize_uniprot_mapping.py
```

The reference set is 85 gene symbols but **84 distinct accessions**: FAU and RPS30 share P62861 (the ubiquitin-fusion precursor whose ribosomal moiety is eS30).

Review cases are resolved by hand into `02-data-processed/uniprot-mapped/user-reviewed/<key>.uniprot.user-reviewed.tsv` (`From`, entry, entry name, `Status` = `reviewed`|`drop`, `Notes`), then folded in:

```bash
python3 mgo-adduct-extract/scripts/apply_user_reviews.py   # -> 03-data-final/
```

This verifies every hand-assigned accession against live UniProt before accepting it, and fills `Protein names`/`Gene Names`/`Organism` from the record so a hand-entered row is as complete as a machine one. `Source` marks each row `auto` or `user-reviewed`. **`03-data-final/` is the only stage where output is not one row per input line**: `drop` removes a row (no UniProtKB entry — the note carries the UniParc ID), and one input may yield several rows (Irshad's `Nucleoside diphosphate kinase (NME1-NME2)` is one label over two gene products).

Finally, merge the per-dataset lists into the canonical one:

```bash
python3 mgo-adduct-extract/scripts/merge_final_datasets.py   # -> 04-data-merged/mgo-glycated-proteins.tsv
```

Per-dataset provenance (study, cell line, system, evidence class) lives in that script's `PROVENANCE` table, mirroring [01-data-extracted/MANIFEST.md](01-data-extracted/MANIFEST.md); a new dataset must be added there or the merge refuses to run.

**Accessions are canonicalized before merging.** The papers key their tables differently — Donnellan's supplement is Proteome Discoverer output carrying TrEMBL and isoform accessions, the others give gene symbols or Swiss-Prot IDs — so the same protein otherwise appears under several accessions, and a protein reported *only* under a TrEMBL accession never matches a Swiss-Prot-keyed reference set at all. That is how RPS14 (A0A2R8Y811) and RPL10 (F8W7C6) went missing from the ribosomal count. Each unreviewed accession is resolved to its reviewed counterpart by exact gene-symbol lookup; anything ambiguous (no gene, no reviewed entry, or a gene with two reviewed entries) is kept as-is rather than guessed at, and every original accession survives in the `Accessions` column. Use `--no-canonicalize` to key strictly on what each paper printed.

**Read the `Evidence` column before drawing conclusions.** Irshad HMEC-1 and Ashour PDLF cell-free are lysates spiked with exogenous MG, not adducts formed in living cells, and are far deeper for that reason — 25 of the 45 ribosomal hits are cell-free only, and just 20 have any in-cell support.

Irshad's Table S3 prints per-row `MG Modification sites` and `Total arg`, which makes assignments falsifiable: an MG-H1 site at R77 requires an arginine at position 77, so a candidate whose residue 77 is serine is excluded regardless of how well its name matches. Use it when resolving that dataset's review cases.

```bash
# end-to-end: final TSV + review TSV (+ optional audit JSONL)
python3 uniprot-name-match/scripts/run_uniprot_name_match.py input.txt \
  -o input.uniprot.tsv --candidates-output input.uniprot.candidates.jsonl

# candidate retrieval only, when you want to inspect ambiguity yourself
python3 uniprot-name-match/scripts/fetch_uniprot_candidates.py input.txt -o candidates.jsonl

# non-human input
python3 uniprot-name-match/scripts/run_uniprot_name_match.py in.txt -o out.tsv --organism-id 10090
```

Output naming convention used throughout: `<stem>.uniprot.tsv`, `<stem>.uniprot.review.tsv` (auto-derived from `-o`), `<stem>.uniprot.candidates.jsonl`.

**How it works** — all scoring logic lives in [fetch_uniprot_candidates.py](uniprot-name-match/scripts/fetch_uniprot_candidates.py); the runner only adds confidence gating and file writing:

1. `infer_record()` splits each line on tabs/` | `/`;`/spaces and classifies fields into `protein_name` / `gene_name` / `identifier` using regexes for accession, entry-name (`VIME_HUMAN`), and gene-symbol shapes.
2. `build_search_strategies()` emits several named UniProt REST queries per record (`identifier_exact`, `accession_exact`, `secondary_accession`, `entry_name_exact`, `gene_symbol`, `hint_entry_name`, `name_as_gene_symbol`, `protein_exact`, `protein_exact_no_paren`, `all_fields_exact`, `protein_plus_gene`, `token_and`), all scoped by `organism_id` (default 9606).
3. `heuristic_score()` ranks every returned record: reviewed Swiss-Prot +100, exact accession match +160, exact entry-name match +140, primary-gene *or* entry-name-mnemonic match +120, exact protein-name match +110/+80/+40 by curation tier (recommended / alternative / TrEMBL submission, capped at +80 when unreviewed), plus token-overlap terms. Candidates are deduped by accession keeping the best-scoring strategy, and ties break on reviewed-then-accession so reruns are reproducible.
4. `confidence_for()` in the runner turns (top score, gap to runner-up, identifier/gene agreement, reviewed flag) into `high`/`medium`/`low` with a human-readable reason. A hint that disagrees with the protein name caps the row at `medium` — that disagreement is usually a defect in the source, and promoting it on name evidence alone would bury it.
5. Only `high` is written into the final TSV by default (`--accept-confidence medium` loosens this). Everything else gets a blank-Entry row in the final TSV *and* a row in the review TSV — inputs are never dropped, so row order and count match the input.

Column contracts:

- final: `From`, `Entry`, `Entry Name`, `Protein names`, `Gene Names`, `Organism`
- review: `From`, `Selected Entry`, `Selected Entry Name`, `Confidence`, `Reason`, `Alternative Entries`

`From` preserves the original input spelling verbatim, which is what makes results traceable back to a paper's supplement.

Gotchas:

- `run_uniprot_name_match.py` does a bare `from fetch_uniprot_candidates import ...`, so it only works when invoked as a script path (Python puts `scripts/` on `sys.path`). It is not importable as a package module.
- Pair an identifier with its protein name (`P08670<TAB>Vimentin`) wherever the source has both. An accession alone has nothing to fall back on if it is ever retired; the paired name makes the fallback work. `mgo-adduct-extract/scripts/make_uniprot_input.py` does this for every dataset.
- Strip isoform suffixes first: UniProt does not resolve `P04075-2`. 53 Donnellan accessions are isoforms; the input generator strips them and the original stays in `01-data-extracted/`.
- When a clue-set has two fields, `From` contains the TAB and the row is CSV-quoted, so `cut -f` mis-splits it — parse with a real CSV reader. Output is CRLF; the extracted TSVs are LF.
- Requests are rate-limited by `--delay` (default 0.2 s) and cached per `(protein_name, gene_name, identifier)` within a run; network failures for one strategy are swallowed and the run continues with fewer candidates.
- Protein names in this domain routinely start with a digit (`40S ribosomal protein S4`, `14-3-3 protein eta`, `78 kDa glucose-regulated protein`, entry name `1433Z_HUMAN`). The input cleaner strips a leading number only when a separator and a space follow it. `uniprot-name-match/tests/test_parsing.py` pins this offline — run it after touching the parser; it needs no network.
- Some review cases are genuine UniProt ambiguity, not weak matching. `40S ribosomal protein S4` is an exact alternative name of **both** P15880 (RPS2, historical `RPS4`/LLRep3 naming) and P62701 (RPS4X); they tie on score and a human has to choose.

When resolving review cases by hand or with a model, follow [references/matching-policy.md](uniprot-name-match/references/matching-policy.md): select only from retrieved candidates, prefer reviewed human entries and identifier/gene evidence over lexical overlap, and leave genuine ties in review rather than forcing a pick.

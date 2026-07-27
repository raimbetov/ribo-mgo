---
name: mgo-adduct-extract
description: "Extract, from a published proteomics paper, the list of proteins shown by mass spectrometry to carry a methylglyoxal-derived adduct — and only those, excluding abundance-change lists and computational hotspot predictions. Use when adding a paper to the MGO literature-mining set, re-deriving 01-data-extracted/, or auditing why a given protein was included or left out."
---

# MGO Adduct Extract

## Overview

Turns published supplements into one TSV per paper × cell line × experimental system,
listing proteins **proteomically shown to carry a methylglyoxal adduct**. Output is a
verbatim mirror of the published table — published columns as printed, no UniProt
mapping, no dedup, no merging — with provenance columns so every row traces back.

The work splits into a judgment stage and a mechanical stage, and the split matters:

- **Scoping** (judgment): decide *which table in the paper qualifies*. This is the hard
  part and the part that goes wrong. It cannot be automated; it can be made auditable.
  See [references/inclusion-policy.md](references/inclusion-policy.md).
- **Extraction** (mechanical): fully deterministic. Each dataset is one YAML spec in
  `datasets/`; `scripts/extract.py` contains no per-paper knowledge. The spec is the
  reviewable artifact.

## Setup

```bash
pip install -r requirements.txt        # pdfplumber, openpyxl, PyYAML
```

`pdftotext` (poppler-utils) must be on PATH. It parses the `pdf_text_rows` model, and
provides the independent second parse used by the `cross_method_identifiers` check.
That check applies only where a table has an *atomic* identifier column (a UniProt
accession or entry name) — Ashour cell-free and Zheng Table S1 declare it. It cannot be
used on a free-text protein-name column, because `pdftotext` word-splitting will not
reconstruct multi-word names; Irshad's Table S3 relies instead on its two independent
routes to the published total of 411.

## Running

```bash
cd scripts
python3 extract.py --all                       # build every dataset
python3 extract.py ../datasets/irshad-2019-hmec1.yaml
python3 verify.py                              # structural + count checks
```

Both default to the repository containing this skill; override with `--repo-root` and
`--out`. `extract.py` writes nothing for a dataset whose checks fail, and exits non-zero.

## Adding a paper

1. **Scope it.** Follow [references/inclusion-policy.md](references/inclusion-policy.md):
   read the Methods for the adduct masses searched, read the Results to find which table
   reports *detected* adducts, then have an independent reviewer try to refute the choice.
   Record the decision and its quoted evidence in the spec's `adduct_evidence` and
   `count_evidence` fields — those quotes are what make the file auditable later.

2. **Measure the table.**
   ```bash
   python3 measure_columns.py supp.pdf --find "Table S3"      # which pages
   python3 measure_columns.py supp.pdf --page 12 --baselines  # how rows lay out
   python3 measure_columns.py supp.pdf --pages 12 18 --anchor-x 55 78
   ```
   Read column boundaries off the gaps in the histogram. The `--baselines` dump tells
   you which `assign` mode to use (see below).

3. **Write the spec.** Copy the closest existing one. Set `expected_rows` from the
   paper's own stated count wherever the paper states one, and record the quote in
   `count_evidence`. Where it states none — 4 of the 11 shipped datasets are in this
   position, including the Donnellan sheets — say so explicitly in `count_evidence`
   ("sheet row count; no count stated in the paper") rather than implying an authority
   the number does not have. Never adjust it to match a parse that disagrees.

4. **Build and verify.** `python3 extract.py ../datasets/new.yaml && python3 verify.py`

## Spec reference

A dataset is a list of `sources`, each parsed by one `model`, concatenated into one TSV.
Nothing is deduplicated — a protein in two source tables appears twice, distinguishable
by provenance.

| model | for |
|---|---|
| `xlsx` | one worksheet; header row found by a marker cell (sheets in one workbook often differ) |
| `pdf_anchor` | PDF table with a row number in the left margin; columns cut by x-position |
| `pdf_text_rows` | regex over `pdftotext -layout`, for pages whose embedded text has no word spacing |
| `literal` | rows transcribed from prose or a figure, carried in the spec |
| `label_union` | union of named label lists (figure panels), with membership columns |

**`assign` is the one setting that silently corrupts data if wrong.** `below` attaches a
word to the nearest anchor at or above it — correct when a row's number is always on its
first baseline. `nearest` attaches to the closest anchor in either direction — required
when a row prints its number *below* its name, as Irshad's rows 1–4 do. Using `below`
there would attach those names to the previous row. Check with `--baselines` first.

### Other spec keys

`lookups` fills a column from another dataset built **in the same run** (the Fig. 4
dataset takes its accessions from Zheng Table S1). Dependencies are ordered
topologically, and a target column must also appear in some source's `emit` or the build
fails — otherwise it would be populated and then silently discarded. Reading a
previously *written* file is deliberately not supported: a stale or failed dependency
could otherwise supply values without anyone noticing.

`caveats` are printed after a successful build and belong in the manifest — coverage
gaps, published errata, anything a downstream reader would misread the file without.

### Checks

Declared per dataset; all must pass before anything is written.

| check | asserts |
|---|---|
| `contiguous_numbering` | row numbers run 1..N with no gaps or repeats |
| `column_sum` | a count column sums to a total the paper prints |
| `token_sum` | the *listed* items independently reach a stated total |
| `declared_vs_listed` | per-row declared count matches the printed list |
| `subset` | one column's tokens are a subset of another's; `on_fail: drop_column` for a misaligned annotation column |
| `regex_column` | identifier syntax |
| `unique_column` | no duplicate identifiers |
| `blank_count` | exactly N rows have an empty cell in this column |
| `no_blank` | required fields present |
| `any_nonempty_matching` | every row carries adduct evidence in at least one column matching a pattern |
| `cross_method_identifiers` | re-parses one source with `pdftotext` and compares the identifier sets — catches both invented identifiers and dropped rows, without depending on row alignment |
| `covers_all_of` | every identifier in a dependency dataset also appears here |

A check naming a column the dataset does not produce is a **spec error**, not a pass.

**Pin known defects; never blanket-allow them.** Where a paper's own numbers disagree,
list the offending rows explicitly — `expect_offenders: [57, 72]` on `declared_vs_listed`,
`expect_duplicates` on `unique_column`, `equals` on `blank_count`. The build then fails if
a *new* mismatch appears **or** a pinned one disappears. A bare "tolerate mismatches" flag
would let a parse silently rot behind a check that was meant to document two typos.

Prefer a check the *paper itself* can settle. Irshad's Table S3 prints `TOTAL: 220` and
`TOTAL 411`; the spec asserts both the declared counts and the parsed site tokens reach
411 independently. Two independent routes to a number the authors published is much
stronger evidence than "the parse looked right".

## Figure-derived data

`zheng-2024-shsy5y-fig4` is the only dataset not from a text layer: Zheng's Fig. 4 is a
raster image and its row labels exist nowhere else. Such data is reproducible as a
checked-in transcription with a consensus protocol, not as a parser. The protocol used:

1. Render the figure page at 400+ dpi and crop each panel's label strip.
2. Transcribe every label, top to bottom, exactly as printed — never "correcting" a
   symbol toward a more likely gene, never padding to reach an expected count.
3. Have **at least two independent readers per panel** re-read blind, each re-rendering
   from the source PDF rather than trusting the supplied crop. Any disagreement is
   resolved against the page, not by majority.
4. Count rows **structurally**, by pixel analysis of the heatmap colour bands, so the
   count does not depend on counting your own text.
5. Reconcile against something external — here, all 77 Table S1 genes had to appear in
   the union.

For this dataset all eight passes returned byte-identical lists and the four pixel row
counts confirmed 33/52/57/57. It still recovers 152 gene symbols against a stated 153
proteins; the spec records that gap rather than papering over it.

## Non-negotiables

- **Transcribe as printed.** Papers contain real errors — Irshad rows 57 and 72 declare
  the wrong site count; Alhujaily gives SLC6A9 and RBM22 the wrong protein names. Record
  them in `caveats` and let the check report them. Never silently fix.
- **`expected_rows` comes from the paper** wherever the paper states a count. Where it
  does not, `count_evidence` must say that plainly. Either way, never edit it to match
  a parse that disagrees — investigate the parse.
- **Never drop an input row** to make a count work. If a count does not reconcile, the
  spec fails and nothing is written — investigate instead.
- **In-cell and cell-free stay in separate files.** A lysate spiked with exogenous MG is
  not the same evidence as adducts formed in living cells, and the cell-free lists are
  far deeper (Irshad 220 vs 2; Ashour 172 vs 5). Pooling them silently overstates the
  in-cell result.
- **Record coverage gaps.** Alhujaily names 142 MG-H1-positive proteins out of the 493 +
  120 + 42 + 26 its TABLE 2 reports per fraction (the total, 681, is a sum of those four
  — the paper never prints it); Zheng tabulates 77 of a stated 153. A downstream zero-overlap result may be a data-availability
  artifact rather than biology.

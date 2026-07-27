# Raw extraction: proteomically-detected methylglyoxal adducts

One TSV per paper × cell line × experimental system, listing proteins **shown by mass
spectrometry to carry a methylglyoxal-derived adduct**. Verbatim mirrors of the published
tables: published columns kept exactly as printed, no UniProt mapping, no dedup, no merging.
Every row carries `source_file` / `source_table` / `source_row` so it traces back to its origin.

Produced by the [`mgo-adduct-extract`](../mgo-adduct-extract/SKILL.md) skill. Every file here is
described by a YAML spec in `../mgo-adduct-extract/datasets/`, which records the source table,
its geometry, the expected row count with the paper's own wording as evidence, the adducts
searched, and the known publication errata. To regenerate and re-check:

```bash
cd ../mgo-adduct-extract && pip install -r requirements.txt
python3 scripts/extract.py --all      # writes this directory
python3 scripts/verify.py
```

Nothing is written for a dataset whose checks fail. The skill's `references/inclusion-policy.md`
documents *why* each table was included or excluded, and how to scope a new paper the same way.

## Inclusion rule

Included only where the adduct itself was detected — a localized MG-H1/MG-H (+54.0106 Da on Arg)
or carboxyethyl (+72.0211 Da; CEA on Arg, CEL on Lys) modification identified by MS/MS.

Deliberately excluded: proteins merely increased/decreased in abundance under MGO or high glucose;
pathway/enrichment analyses; receptor-binding-domain (RBD) predictions; primer and correlation tables.

## Files

| File | Cell line | System | Source | Rows |
|---|---|---|---|---|
| `Donnellan et al. 2022 - WIL2-NS.tsv` | WIL2-NS | in-cell, 500 µM MGO 24 h | S2 xlsx, sheet `WIL2-NS` | 519 |
| `Donnellan et al. 2022 - OV90.tsv` | OV90 parental | in-cell, basal | sheet `OV90 (Parental)` | 119 |
| `Donnellan et al. 2022 - Caov3.tsv` | Caov3 parental | in-cell, basal | sheet `Caov3 (Parental)` | 126 |
| `Donnellan et al. 2022 - PBL.tsv` | PBL, 3 donors pooled | ex vivo, basal | sheet `PBL` | 57 |
| `Irshad et al. 2019 - HMEC-1.tsv` | HMEC-1 | **cell-free** lysate + 500 µM MG, 24 h | supp Table S3 | 220 |
| `Irshad et al. 2019 - HAEC.tsv` | HAEC | in-cell, 20 mM glucose | main text Results (prose) | 2 |
| `Ashour et al. 2020 - PDLF cell-free.tsv` | PDLF (primary) | **cell-free** lysate + MG | supp Table S2 | 172 |
| `Ashour et al. 2020 - PDLF in-cell.tsv` | PDLF (primary) | in-cell, endogenous | main text Results (prose) | 5 |
| `Alhujaily et al. 2021 - HEK293.tsv` | HEK293 | in-cell, 131 µM MG 6 h | S3+S4+S11, Table 4, Table 2 fn a/b | 147 (142 unique) |
| `Zheng et al. 2024 - SH-SY5Y.tsv` | SH-SY5Y | in-cell, 0–1000 µM MGO | supp Table S1 | 77 |
| `Zheng et al. 2024 - SH-SY5Y (Fig4 recovered).tsv` | SH-SY5Y | as above | **Fig. 4 heatmap labels** | 152 |

**In-cell vs cell-free matters.** Irshad's HMEC-1 and Ashour's PDLF cell-free files are lysates
spiked with exogenous MG, not adducts formed in living cells. They are far deeper than the in-cell
lists (220 and 172 vs 2 and 5) precisely because of that. Do not pool the two kinds of evidence
without saying so.

## Adducts searched, per paper

| Paper | Adducts | Evidence |
|---|---|---|
| Donnellan 2022 | MG-H (+54.010565, Arg), carboxyethyl (+72.021129, Arg/Lys) | Methods 4.7, dynamic modifications, 1% FDR |
| Irshad 2019 | MG-H1 (+54 Da, Arg) | Results: "+54 Da mass increment on arginine residues, reflecting MG-H1 formation" |
| Ashour 2020 | MG-H1 (+54.01 Da, Arg) | Results: "proteins with MG-H1 (+54.01 Da mass increment on arginine residues)" |
| Alhujaily 2021 | MG-H1 | Methods: MG-H1 as a Mascot variable modification (no delta stated in paper) |
| Zheng 2024 | Arg MG-H (+54.010565); carboxyethyl/hemiaminal (+72.021129) on Arg, Lys, Cys | Methods 2.6, MaxQuant variable modifications |

## Excluded paper: Sun et al. 2019

**Not a methylglyoxal dataset — excluded in full.** Its Methods specify the only glycation
delta-mass searched: *"variable modifications: oxidation of methionine (+ 15.9949 Da), glycation of
lysine or arginine (+ 162.0528 Da)"*. +162.0528 Da is a hexose (Amadori / fructosyl-lysine), i.e.
glucose-derived. No MGO adduct mass appears anywhere in its search space. "Methylglyoxal" occurs in
the paper only in the Discussion, citing others' work.

Its 5 cell-line datasets (HEK293T, Jurkat, MCF7, plus Expt 1/2 of MCF7) are therefore **not**
extracted. → **`CLAUDE.md` lists Sun among the MGO papers and should be corrected.**

## Coverage gaps — these files are not the complete adduct lists

**Alhujaily: 142 of 681.** Table 2 reports 681 proteins detected with MG-H1 (cytoplasm 493, nucleus
120, mito matrix+IMS 42, mito membrane 26). The article names only those that *also* changed
significantly in abundance, plus the 16 spliceosome proteins of Table 4. The remaining ~539 exist
only in PRIDE **PXD029315** and would require re-searching the raw data.

**Zheng: 77 tabulated, 152 recovered, 153 stated.** Table S1 holds only the 77 proteins passing a
≥3-of-4-replicate filter. The full set appears solely as Fig. 4 heatmap row labels. There is no
repository deposit ("Data will be made available on request").

**Donnellan: CBPR sublines not extracted.** `OV90 (CBPR)` (133 rows) and `Caov3 (CBPR)` (102 rows)
are carboplatin-resistant derivatives, excluded as selection-pressure-derived. Parental lines only.

## The Fig. 4 recovery file

`Zheng et al. 2024 - SH-SY5Y (Fig4 recovered).tsv` is the one file **not** derived from a text layer.
Fig. 4 is a raster image, so its labels were read visually, then re-read independently **eight times**
(two blind passes per panel, each re-rendering from the PDF at 600–1200 dpi). All eight passes
returned byte-identical lists. Each panel's row count was separately confirmed by pixel analysis of
the heatmap colour bands: **A=33, B=52, C=57, D=57**. All 77 Table S1 genes appear in the union.

It nonetheless yields **152 unique gene symbols against the paper's stated 153 proteins**. The most
likely cause is a counting-unit difference — Fig. 4 labels rows by *gene symbol* while the text
counts *proteins* (accessions), so two accessions sharing one symbol collapse to a single label. This
cannot be resolved from the published material. **Treat the file as 152 of 153, and as
figure-derived rather than table-derived.** It is the only source of Zheng's sole ribosomal protein,
RPL8.

## Publication inconsistencies, transcribed as printed and never silently corrected

- **Irshad Table S3** — row 57 (FUBP2/KHSRP) declares 1 site but prints two (R331, R340); row 72
  (RhoGDI2/ARHGD1B) declares 2 but prints one (R131). Both the declared counts and the printed site
  tokens independently sum to the table's own printed **TOTAL 411**, so these offset each other.
- **Irshad Table S3, RBD column dropped.** The "MG Modification sites in the RBD" column is shifted
  one row down in the PDF's own text layer over ~24 rows (e.g. row 86 RPS3A shows row 85 RPS2's
  site). Verified by a subset check: RBD sites must be a subset of that row's detected sites. Since
  RBD is a *predicted* domain annotation and not adduct evidence, the column is omitted rather than
  shipped misaligned.
- **Ashour Table S2** — the declared "No of MG modification sites" column sums to 353 (matching the
  Results text and the table's own total), but the printed site lists contain only 343 tokens; 8 rows
  disagree (rows 1, 2, 8, 11, 17, 43, 51, 58). Vimentin, for instance, declares 20 and prints 19.
  Confirmed against the rendered page — it is the paper's discrepancy, not an extraction artifact.
- **Alhujaily** — two protein-name cells are copy-paste errors in the published PDF: Table S3 row 27
  (SLC6A9) carries FAAP20's name, and Table 4 row 11 (RBM22) carries HSPA1B's name. **Map from the
  Gene column, not the name.** Table S11 has 12 rows while Table 2 states 2 up / 9 down for that
  fraction; the 12 printed rows are taken.
- **Donnellan WIL2-NS** — 519 rows but 514 unique accessions: four histone H2B accessions (P58876 ×2,
  O60814 ×3, Q99879 ×2, Q5QNW6 ×2) recur with different site sets. Preserved, not deduped.

## Verification performed

- Row counts matched against the papers' own stated totals (evidence quoted in `verify_all.py`).
- Site-count reconciliation where a total is printed: Irshad 411 (both declared and parsed), Ashour 353.
- Contiguous row numbering 1..N with no gaps or repeats in every numbered source table.
- Every PDF table extracted **twice by independent methods** — pdfplumber word-coordinates and
  `pdftotext -layout` — and cross-checked on identifiers and counts. (`pdftotext` cannot parse
  Irshad's 4 wrapped rows; those were confirmed against the rendered page instead.)
- UniProt accession / entry-name syntax validated; bidirectional check that no identifier was
  invented and none dropped.
- Structural: no ragged rows, no missing provenance, no stray tabs or newlines.
- Regression against the earlier UniProt-mapped extraction that previously occupied this directory,
  run before those files were deleted as superseded: Caov3 126=126, OV90 119=119, PBL 57=57,
  WIL2-NS 519 raw vs 514 mapped (+5 duplicate histone rows), Irshad 220 raw vs 218 mapped (+2 that
  failed mapping). All differences expected; nothing was lost relative to the earlier files.

Cytoplasmic ribosomal proteins present per file, as a smoke test that the intended tables were
extracted: Irshad HMEC-1 **29**, Donnellan WIL2-NS **16**, Ashour PDLF cell-free **3**
(as entry names: RS17, RS28, RLA2), Donnellan OV90 **2**, Donnellan Caov3 **1** (RPS14), Alhujaily
HEK293 **1** (UBA52), Zheng Fig4 **1** (RPL8), all others **0**.

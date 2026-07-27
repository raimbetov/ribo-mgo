# Inclusion policy: what counts as a methylglyoxal adduct dataset

The decision this document governs is *which table in a paper to extract*. It is the
step that most often goes wrong, because papers on methylglyoxal routinely publish three
different kinds of protein list side by side, and only one of them qualifies.

## The rule

Include a protein **only where the adduct itself was detected** — a localized
methylglyoxal-derived modification identified by MS/MS on a specific residue.

Methylglyoxal-derived adducts and their delta masses:

| Adduct | Residue | Δ mass |
|---|---|---|
| MG-H1 / MG-H2 / MG-H3 (hydroimidazolone) | Arg | +54.0106 |
| Carboxyethyl — CEA on Arg, CEL on Lys | Arg, Lys | +72.0211 |
| Carboxyethyl / hemiaminal | Cys | +72.0211 |
| Argpyrimidine | Arg | +80.0262 |
| THP (tetrahydropyrimidine) | Arg | +126.0317 |

**Not methylglyoxal**, however the paper uses the word "glycation":

| Adduct | Precursor | Δ mass |
|---|---|---|
| Fructosyl-lysine / Amadori | **glucose** | +162.0528 |
| CML, G-H1 | glyoxal | +58.0055 (CML) |
| 3DG-H | 3-deoxyglucosone | +144.0423 |

## Exclude

- **Abundance changes.** Proteins increased or decreased under MGO or high glucose. These
  are the most common trap because they sit in the same supplement, are much longer, and
  are often enriched for interesting complexes.
- **Computational predictions.** Arginine-hotspot or functional-domain analyses with no
  MS evidence of an adduct.
- **Pathway / enrichment analyses.** A DAVID or REACTOME table naming 51 ribosomal genes
  is an enrichment of some *other* list, not an adduct list.
- **Everything else**: primer tables, expression-correlation tables, MS2 spectra figures.

## Procedure

1. **Read the Methods for the modification searched.** Find the actual delta masses or
   named modifications in the database-search parameters. This single step settles most
   cases and is not negotiable — a title saying "glycation" tells you nothing.
2. **Read the Results for how each table is introduced.** The sentence that cites a table
   usually states plainly whether it is detected or predicted, and how many proteins.
3. **Take the count from the paper, not the parse.** Record the sentence verbatim in the
   spec's `count_evidence`.
4. **Separate in-cell from cell-free.** Many papers report a small endogenous list from
   living cells and a much larger list from a lysate spiked with exogenous MG. Both are
   real MS evidence, but they answer different questions and belong in different files.
5. **Have the choice adversarially reviewed.** An independent reviewer should try to
   *refute* it on five axes: inclusion error (abundance or prediction admitted), chemistry
   error (wrong precursor), omission error (a qualifying table missed), count error,
   split error (cell lines merged or wrongly separated). Default to reporting a problem
   when unsure.
6. **Record the decision and its evidence** in the spec, including which tables were
   deliberately excluded and why.

## Worked cases

These four are the reason the policy reads as it does.

### 1. Sun et al. 2019 — excluded; not methylglyoxal at all

Titled "Comprehensive Analysis of Protein Glycation", it sits in a folder of MGO papers
and reports glycated proteins across three human cell lines. Its Methods name exactly one
glycation delta mass:

> variable modifications: oxidation of methionine (+ 15.9949 Da), glycation of lysine or
> arginine (+ 162.0528 Da)

+162.0528 Da is a hexose — glucose-derived fructosyl-lysine. No MGO adduct mass appears
anywhere in its search space. "Methylglyoxal" occurs only in its Discussion, citing
others' work. **Excluded in full**, despite five otherwise usable cell-line datasets.

*Lesson: check the delta mass before anything else. The word "glycation" is not evidence.*

### 2. Ashour et al. 2020 — included; "at risk of" turned out to mean measured

Its Table S2 is headed *"Proteins at risk of MG modification"*, which reads like a
hotspot prediction and would disqualify it. The Results settle it the other way:

> we interrogated proteomics data for evidence of proteins with MG-H1 (+54.01 Da mass
> increment on arginine residues) … We then detected MG-H1 modification on 172 proteins
> in 353 unique modification sites (online supplemental table S2)

A little further in the same paragraph the authors isolate the genuinely predicted layer — *"From RBD analysis, 115 of
the 353 (33%) modifications were in predicted function domains"* — confirming that RBD is
a downstream annotation of already-detected sites, not a protein filter. Irshad et al.
2019 uses the same phrase for the same kind of list. **Both included**; their RBD columns
carry no adduct evidence of their own.

*Lesson: a table title is not the decision. Find the sentence that cites the table.*

### 3. Alhujaily et al. 2021 — included, with a coverage gap that must be recorded

Only tables whose titles say *"containing MG-H1 residues"* (S3, S4, S11), main-text
TABLE 4 (*"Spliceosome proteins detected with MG-H1 modification"*), and the two TABLE 2
footnotes qualify. Tables S1, S2 and S5–S10 are abundance-change lists.

The trap here is TABLE 3, a DAVID enrichment whose "Ribosome" row names 51 ribosomal
genes. In a ribosome-focused project that looks exactly like the answer. It is an
enrichment of the abundance-**decreased** set and carries no adduct evidence; none of
those genes may be admitted. The only cytoplasmic ribosomal protein with real MG-H1
evidence in this paper is UBA52.

TABLE 2 also reports MG-H1-positive protein counts per fraction — 493 + 120 + 42 + 26,
**681** in total, though the paper never prints that sum — while the
article names only the ~142 that also changed in abundance. That gap belongs in the
spec's `caveats`, because a downstream analysis will otherwise read 142 as the whole
picture.

*Lesson: the most attractive-looking table is often the enrichment. And a paper can
detect far more than it names.*

### 4. Zheng et al. 2024 — included; the full list exists only in a figure

Table S1 tabulates 77 proteins, but Results 3.2 states 153 were detected; the other 76
appear only as Fig. 4 heatmap row labels, and there is no repository deposit. Recovering
them needs the vision protocol in SKILL.md, and even then reconciles to 152 of 153.

Note also that Table S1's J-statistic and p-value columns describe the concentration
*trend* of an already-detected adduct. Rows with p = 1.00 are still adduct-positive and
must not be filtered out — a filter here would silently discard real data.

*Lesson: distinguish "was an adduct detected" from "did the adduct respond to dose".*

## When the answer is genuinely unclear

Leave the paper out and say why, rather than guessing. A wrongly included abundance list
contaminates every downstream intersection and is very hard to detect later; a missing
paper is visible and easy to add.

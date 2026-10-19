# Audit of FDR reporting in DIA proteomics articles

A convenience sample of 63 articles from the 2024-2026 literature, scored on
ten criteria and a set of descriptors, from their main text.

| File | Content |
|---|---|
| `scoring.csv` | the scoring sheet: one row per article, with every score and, beside it, the sentence of the article that supports it, quoted verbatim |
| `audit_pipeline.py` | sentence extraction from the PDFs (assigns no score), five verification checks, the numbers quoted in the paper, Tables S3 and S5, and the workbook |
| `Data_S1_audit.xlsx` | Supplementary Data S1, generated from `scoring.csv`: sheets README, Scores, Evidence, Summary, Codebook |
| `out/` | `audit_numbers.txt` (every number of Section 2.2 and Note S7), `tableS3_audit.tex`, `tableS5_corpus.tex` |

## Commands

```
python audit_pipeline.py extract   --pdf-dir pdfs --out-dir work
python audit_pipeline.py check     --sheet scoring.csv --txt-dir work/txt
python audit_pipeline.py summarise --sheet scoring.csv
python audit_pipeline.py build     --sheet scoring.csv --out-dir out
```

`extract` needs `pdftotext` (poppler-utils); `build` needs `openpyxl` for the
workbook. Everything else uses the Python standard library.

The PDFs are not redistributed. Each row of `scoring.csv` gives the DOI and the
file name (`pdf_file`) expected in `pdfs/`. Without the PDFs, `check` still runs
the checks that do not need the article texts, and `summarise` and `build`
run in full.

## The scoring sheet

Each item has a score
column and a quote column prefixed `q_`:

| Column | Item | Values |
|---|---|---|
| `id_fdr` | identification-level threshold stated | Yes / Partial / No |
| `unit` | criterion 1, controlled unit named | Yes / Partial / No |
| `count_any` | an identification count is reported | Yes / No |
| `count_joint` | criterion 3, count and threshold reported together | Yes / No |
| `count_unit` | criterion 4, the count is of the controlled unit | Yes / No |
| `scope` | criterion 5, run-wise or global | Yes / Partial / No |
| `mbr` | criterion 6, match-between-runs setting | Yes / Partial / No |
| `qfield` | criterion 7, thresholded confidence column named | Yes / Partial / No |
| `decoy_count` | criterion 8, accepted decoy count | Yes / No |
| `label_rule` | criterion 9, label-propagation rule | Yes / Partial / No |
| `external` | criterion 10, external validation | Yes / No |
| `da_fdr`, `contradiction`, `pep`, `decoy_strategy`, `defaults`, `version` | descriptors | Yes / No |

Criterion 2 (count and threshold both reported, anywhere) is derived:
`count_any` is Yes and `id_fdr` is Yes. Definitions and the rule for each
value are in the Codebook sheet of the workbook. ` [...] ` in a quote joins
passages that are separate in the article or interrupted by a page or column
break.

To change a score, edit `scoring.csv` (score and quote together) and run
`check`, then `build`.

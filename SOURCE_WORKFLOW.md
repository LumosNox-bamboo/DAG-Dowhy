# Audit of the original DAG/DoWhy workflow

The project's active research directory contains original analysis code. It was missed during the first audit because the first path supplied was an archive directory without the `DAG py` subdirectory. This audit corrects that record. No original workbook, notebook outputs, participant records, or manuscript files are included in this public repository.

## Code found

| Source artifact | Role |
| --- | --- |
| `DAG py/DAG+DoWHy.ipynb` | Original 28-cell Colab notebook with several analysis branches. |
| `DAG py/run_full_reproduction.py` | Entry point for the later local reconstruction of Tables 1–4 and related diagnostics. |
| `DAG py/analysis/repro/dataset_build.py` | Dataset preparation and source-column alias resolution. |
| `DAG py/analysis/repro/dag_model.py` | HillClimbSearch/BDeu graph learning and bootstrap edge summaries. |
| `DAG py/analysis/repro/dowhy_analysis.py` | DoWhy propensity-score-weighted edge estimates. |
| `DAG py/analysis/repro/counterfactual.py` | Logistic g-formula scenarios and bootstrap intervals. |
| `DAG py/analysis/final_outputs/reproduction_report.md` | Prior local reproduction summary. |

The reconstructed pipeline uses 200 bootstrap samples for the full graph, 400 for the two-hop core graph, and 300 for g-formula intervals, with seed 2025 in those branches. The data-preparation functions produced 420 complete rows across 28 graph-learning variables and 421 complete rows across 14 effect-analysis variables from the local 423-row workbook in a read-only check. These counts are audit context, not example data.

## What the existing reconstruction establishes

The local reproduction report says Table 1 counts and age summary matched the manuscript, and Tables 2–4 were reconstructed from selected notebook branches. It also explicitly records unresolved discrepancies: competing notebook branches, source-column/name mismatches, hard-coded edge choices, omitted error/non-significant rows, and Table 3 confidence intervals not directly reproduced by that branch. A separate consistency report notes that Table 2 and Tables 3–4 do not use identical edge sets. This is **partial reconstruction**, not an exact, clean-environment replication of every manuscript number.

The original scripts are not directly reusable on an arbitrary dataset. They hard-code the source workbook path and column aliases; the DoWhy and counterfactual stages read a hard-coded edge table rather than the immediately preceding bootstrap result. Several binary recodes turn missing raw values into 0/1, which requires study-specific review before reuse. The notebook also has Colab installation commands and multiple branches. Copying it unmodified into a public template would risk confusing branch selection and, if outputs were included, disclosing participant-level data.

## How this repository differs

`causal_pipeline.py` is an **independent generic estimator** for one binary exposure and one binary outcome. It requires an expert-specified DAG, validates the supplied backdoor adjustment set, uses cross-fitted AIPW, and reports aggregate diagnostics. It does not perform exploratory graph discovery or call DoWhy, so its estimates should not be expected to equal the source notebook's IPW or g-formula numbers. Those method differences are intentional and documented in the README.

To reproduce the original study, a separate private analysis would need to freeze the exact workbook version, variable recodes, selected notebook branch, graph/edge choices, package versions, and output comparison rules, then rerun the original suite without publishing patient-level files. The source suite and its reports provide the starting point for that work.

# Reusable DAG-guided causal analysis framework

This package is a reusable, data-free causal-analysis template developed after reviewing the source project's DAG/DoWhy workflow. It accepts a new CSV or XLSX file and a JSON configuration. It does **not** contain the original patient data, names, study results, figures, or a dataset-specific DAG. See [SOURCE_WORKFLOW.md](SOURCE_WORKFLOW.md) for the audit of the original code.

**Reproducibility boundary:** The synthetic example and this generic pipeline can be rerun with fixed versions and seed. The original study **does have** a notebook and a separate local reproduction script suite. This repository does not implement their exact `pgmpy`/DoWhy/g-formula branches and therefore does not independently reproduce their numeric results. The local reproduction report documents remaining disagreements with the manuscript. Do not cite this package as independent replication of those results.

The executable analysis estimates the **average treatment effect (ATE)** of a binary exposure on a binary outcome as an absolute **risk difference**. It uses an investigator-specified directed acyclic graph (DAG), validates the proposed backdoor adjustment set, estimates nuisance models with stratified cross-fitting, and reports an augmented inverse-probability-weighted (AIPW) estimate. It also reports standardized outcome risks under the two exposure settings, propensity overlap, weighted effective sample size, and covariate balance.

## Quick start

Python 3.9–3.12 is supported with the pinned dependency versions in `requirements.txt`. From this directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python examples/make_synthetic.py --out data/synthetic.csv
python causal_pipeline.py --config examples/config.json --data data/synthetic.csv --out results/synthetic
python -m unittest discover -s tests -v
```

On Windows, activate with `.venv\Scripts\activate`. The example dataset is entirely synthetic. Its output is ignored by Git.

## Use your own data

1. Place a **copy** of your data in a private directory. Keep the original file unchanged. One row must represent one independent participant or other defined analysis unit. Repeated observations, clustering, survival outcomes, continuous exposures, and continuous outcomes need a different estimator.
2. Define the exposure, outcome, and candidate **pre-exposure** covariates before looking at effect estimates. Draw a DAG based on timing and substantive knowledge. Document what each arrow means. Cross-sectional associations alone cannot orient arrows or resolve feedback loops; use time-indexed nodes if feedback is scientifically plausible.
3. Copy `examples/config.json` and edit `columns`, `binary_coding`, `dag`, and `analysis`. `columns` maps your conceptual node names to exact input headers. `binary_coding` maps raw exposure and outcome values to 0/1; unknown nonmissing values cause an error rather than silent recoding. The DAG must contain a direct exposure-to-outcome arrow and must be acyclic. `adjust_for` must block all backdoor paths in that DAG and cannot include descendants of the exposure. You must still judge whether the DAG itself is credible.
4. Run:

```bash
python causal_pipeline.py --config path/to/config.json --data path/to/data.csv --out path/to/results
```

For XLSX, supply the `.xlsx` path and optionally set top-level `"sheet": "Sheet Name"` in the JSON. The first sheet is used by default. The code reads only mapped analysis columns. Identifiers and free text should never be mapped as covariates. Avoid storing input data or generated results in a public repository.

If the file has one unique participant identifier per row, optionally set top-level `"id_column": "your ID header"`. Missing or duplicate IDs then stop the run. The ID values are never exported.

## Configuration reference

| Key | Meaning |
| --- | --- |
| `seed` | Reproducible stratified fold assignment. |
| `id_column` | Optional unique row identifier used only for validation. |
| `columns` | Conceptual DAG node to exact CSV/XLSX header. Only nodes used in the analysis require mappings. |
| `binary_coding` | String-form raw values mapped to integer 0 and 1 for exposure and outcome. Both labels must be represented. |
| `dag.nodes`, `dag.edges` | Prespecified causal graph. An edge is `["source", "target"]`. |
| `analysis.treatment`, `analysis.outcome` | Binary exposure and outcome nodes. |
| `analysis.adjust_for` | Investigator-chosen, pre-exposure backdoor adjustment variables. |
| `analysis.numeric_covariates` | Subset of `adjust_for` treated as continuous; other covariates are one-hot encoded. |
| `analysis.folds` | Stratified cross-fitting folds, default 5. Each exposure group needs at least `2 × folds` records. |
| `analysis.propensity_clip` | Prediction-only numerical clipping, default 0.01. Clipped count is reported. |

All mapped analysis variables use complete cases. The report records input rows, retained rows, rows excluded for missingness, and missing counts by variable. Missing categorical or numeric values are **not** imputed. If missingness is substantial or related to observed characteristics, define a defensible missing-data strategy before interpreting effects. An imputed source file may be used only if its imputation procedure was developed and validated separately; do not impute outcome or exposure without an explicit analysis plan.

## Outputs and interpretation

`results.json` contains aggregate estimates and diagnostics. `balance.csv` contains before- and after-weighting standardized mean differences (SMDs) for encoded covariate columns. The program writes no participant-level records.

The ATE is `E[Y(1) − Y(0)]`. For example, `0.08` means an estimated 8 percentage-point higher outcome risk under exposure 1 than 0 in the analyzed population; it does **not** mean an 8% relative increase. The 95% interval and Wald p value use the empirical cross-fitted AIPW influence curve and a normal approximation. They do not account for uncertainty from DAG selection, measurement error, unmeasured confounding, or nonindependent observations.

`standardized_risk_if_all_treated` and `standardized_risk_if_none_treated` are model-based marginal predictions under two hypothetical exposure settings. They are **not** observed intervention outcomes or proof that changing an exposure will reproduce the estimate. Their difference may differ slightly from the AIPW ATE because the former use the outcome models alone.

Review `warnings`, propensity range, clipping count, weighted effective sample size, and `balance.csv` before interpreting the ATE. An absolute weighted SMD below 0.1 is a useful diagnostic target, not a guarantee of exchangeability. Poor overlap, large weights, residual imbalance, or sparse outcome cells can invalidate a simple causal interpretation. Sensitivity analysis for unmeasured confounding, alternative reasonable DAGs/adjustment sets, missing-data methods, and external or longitudinal validation should be planned for substantive studies.

## Relationship to the original analysis process

The source code describes graph exploration, effect estimation, hypothetical exposure changes, and robustness checks. This template addresses the causal-estimation portion with several safeguards; it is not a line-by-line port of the source notebook:

| Source-process element | Reusable implementation or required analyst action |
| --- | --- |
| Data-driven DAG learning | Use it only for exploratory hypothesis generation outside this package. Supply an expert-reviewed, time-ordered DAG here. Bootstrap edge frequency does not establish causal direction. |
| DoWhy backdoor/IPW effect estimates | Explicit backdoor validation plus cross-fitted AIPW for the supported binary/binary case. The implementation is self-contained; it does not claim to reproduce the original DoWhy numeric results. |
| Hypothetical removal of an exposure | Report model-standardized risks under exposure 0 and 1, clearly labeled as hypothetical. |
| Subset/placebo refuters | Not presented as causal proof. A refuter p value cannot validate no unmeasured confounding or the direction of a cross-sectional edge. |

The framework is a starting point for transparent observational analysis, **not a turnkey causal conclusion**. Its assumptions include a correct time-ordered graph and adjustment set, consistency/well-defined exposure, no interference, no unmeasured confounding, positivity, and adequate nuisance models. The analyst must state the target population and whether selection into the dataset could bias the estimate. A DAG cannot contain both `A → B` and `B → A` at the same time; feedback requires temporal ordering such as `A_t → B_t+1 → A_t+2` and corresponding longitudinal data.

## Reproducibility and privacy

- Keep input data and analytic outputs outside the tracked package, or under ignored `data/` and `results/` directories.
- Save the exact config, software versions, data provenance, sample selection rules, and variable coding alongside each study's private analysis record.
- Freeze the DAG and primary estimand before inspecting the final effect whenever confirmatory claims are intended.
- Independently review coding, missingness, overlap, model fit, and alternate plausible causal structures.
- Do not treat a fixed random seed, a significant p value, or a successful refuter as evidence that a causal assumption holds.

## Files

| File | Purpose |
| --- | --- |
| `causal_pipeline.py` | Command-line analysis. |
| `examples/config.json` | Generic DAG and column-mapping example. |
| `examples/make_synthetic.py` | Fictional-data smoke test generator. |
| `tests/test_pipeline.py` | DAG validation and synthetic-run checks. |
| `.github/workflows/test.yml` | Runs tests and a synthetic CLI smoke test on GitHub. |
| `SOURCE_WORKFLOW.md` | Audit of the original notebook/reproduction scripts and the scope of this template. |
| `requirements.txt` | Runtime dependencies. |

# DAG-guided causal analysis

Run a binary exposure–outcome analysis on your own CSV or Excel data. The pipeline produces an estimated risk difference and aggregate diagnostics.

## Install

Python 3.9–3.12:

```bash
git clone https://github.com/LumosNox-bamboo/DAG-Dowhy.git
cd DAG-Dowhy
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, activate the environment with `.venv\Scripts\activate`.

## Configure your data

Use one row per person. Put your CSV or XLSX file in `data/`, then copy `examples/config.json` and edit these fields:

| Field | What to enter |
| --- | --- |
| `columns` | Your exact data-column headers for the exposure, outcome, and adjustment variables. |
| `binary_coding` | How the raw exposure and outcome values map to `0` and `1`. |
| `dag.nodes` / `dag.edges` | Your variables and directed edges. |
| `analysis.treatment` / `analysis.outcome` | The exposure and outcome node names. |
| `analysis.adjust_for` | Variables to adjust for. |
| `analysis.numeric_covariates` | Adjustment variables to treat as numeric; the rest are categorical. |

For an Excel file, set `"sheet": "Sheet Name"` at the top level of the config if your data are not on the first sheet. You can also set `"id_column": "ID header"` to check for duplicate or missing IDs.

## Run

```bash
python causal_pipeline.py --config my_config.json --data data/my_data.csv --out results/my_analysis
```

The output folder contains:

- `results.json`: estimated risk difference, 95% interval, sample counts, and diagnostics.
- `balance.csv`: covariate balance before and after weighting.

## Try the included example

```bash
python examples/make_synthetic.py --out data/synthetic.csv
python causal_pipeline.py --config examples/config.json --data data/synthetic.csv --out results/synthetic
python -m unittest discover -s tests -v
```

The example uses generated data. No original study data are included in this repository.

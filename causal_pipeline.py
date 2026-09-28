"""Small, auditable binary-treatment causal analysis from an expert-specified DAG."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import norm


def fail(message):
    raise ValueError(message)


def check_dag(nodes, edges):
    if len(nodes) != len(set(nodes)):
        fail("DAG nodes must be unique")
    parents = {node: set() for node in nodes}
    for source, target in edges:
        if source not in parents or target not in parents or source == target:
            fail(f"Invalid DAG edge: {source} -> {target}")
        parents[target].add(source)
    remaining = set(nodes)
    while remaining:
        ready = {n for n in remaining if not (parents[n] & remaining)}
        if not ready:
            fail("DAG contains a directed cycle")
        remaining -= ready
    return parents


def descendants(node, edges):
    found, frontier = set(), {node}
    while frontier:
        frontier = {b for a, b in edges if a in frontier} - found
        found |= frontier
    return found


def valid_backdoor(treatment, outcome, adjust, nodes, edges):
    """Check d-separation after deleting arrows leaving treatment."""
    parents = {node: set() for node in nodes}
    for source, target in edges:
        if source != treatment:
            parents[target].add(source)
    ancestral = {treatment, outcome, *adjust}
    while True:
        expanded = ancestral | set().union(*(parents[n] for n in ancestral))
        if expanded == ancestral:
            break
        ancestral = expanded
    neighbors = {n: set() for n in ancestral}
    for child in ancestral:
        ps = parents[child] & ancestral
        for parent in ps:
            neighbors[child].add(parent)
            neighbors[parent].add(child)
        for parent in ps:
            neighbors[parent] |= ps - {parent}
    seen, frontier = set(adjust), {treatment}
    while frontier:
        if outcome in frontier:
            return False
        seen |= frontier
        frontier = set().union(*(neighbors[n] for n in frontier)) - seen
    return True


def read_input(path, sheet=None):
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    if path.suffix.lower() == ".xlsx":
        return pd.read_excel(path, sheet_name=sheet or 0)
    fail("Input must be .csv or .xlsx")


def encode_binary(series, mapping, name):
    if not mapping or set(mapping.values()) != {0, 1}:
        fail(f"{name}: coding must map observed values to both 0 and 1")
    converted = series.astype("string").str.strip().map(mapping)
    bad = series.notna() & converted.isna()
    if bad.any():
        fail(f"{name}: {int(bad.sum())} nonmissing values are absent from coding")
    return converted.astype("float64")


def design_matrix(frame, covariates, numeric):
    pieces = []
    for name in covariates:
        if name in numeric:
            col = pd.to_numeric(frame[name], errors="raise").astype(float)
            sd = col.std(ddof=0)
            if sd == 0:
                fail(f"Constant covariate: {name}")
            pieces.append(((col - col.mean()) / sd).to_frame(name))
        else:
            dummies = pd.get_dummies(frame[name].astype("string"), prefix=name, drop_first=True, dtype=float)
            if dummies.empty:
                fail(f"Constant covariate: {name}")
            pieces.append(dummies)
    return pd.concat(pieces, axis=1).to_numpy(float) if pieces else np.empty((len(frame), 0))


def fit_logit(x, y, penalty=1e-4):
    x = np.column_stack([np.ones(len(x)), x])
    if len(set(y)) != 2:
        fail("Both binary classes are required in every fitted group")

    def objective(beta):
        z = np.einsum("ij,j->i", x, beta)
        loss = np.logaddexp(0, (1 - 2 * y) * z).sum() + penalty * np.sum(beta[1:] ** 2)
        grad = np.einsum("ij,i->j", x, expit(z) - y)
        grad[1:] += 2 * penalty * beta[1:]
        return loss, grad

    result = minimize(objective, np.zeros(x.shape[1]), jac=True, method="L-BFGS-B",
                      bounds=[(-30, 30)] * x.shape[1])
    if not result.success:
        fail(f"Logistic model did not converge: {result.message}")
    return result.x


def predict_logit(beta, x):
    return expit(np.einsum("ij,j->i", np.column_stack([np.ones(len(x)), x]), beta))


def folds(treatment, count, seed):
    if count < 2:
        fail("folds must be at least 2")
    rng = np.random.default_rng(seed)
    groups = [np.where(treatment == value)[0] for value in (0, 1)]
    if min(map(len, groups)) < count * 2:
        fail("Each treatment group needs at least twice the number of folds")
    for group in groups:
        rng.shuffle(group)
    parts = [np.array([], dtype=int) for _ in range(count)]
    for group in groups:
        for i, part in enumerate(np.array_split(group, count)):
            parts[i] = np.concatenate([parts[i], part])
    return parts


def estimate(frame, spec, seed):
    treatment = frame[spec["treatment"]].to_numpy(int)
    outcome = frame[spec["outcome"]].to_numpy(int)
    covariates = spec["adjust_for"]
    x = design_matrix(frame, covariates, set(spec.get("numeric_covariates", [])))
    n, k = len(frame), int(spec.get("folds", 5))
    p, m0, m1 = (np.empty(n) for _ in range(3))
    for test in folds(treatment, k, seed):
        train = np.setdiff1d(np.arange(n), test)
        p[test] = predict_logit(fit_logit(x[train], treatment[train]), x[test])
        for value, dest in ((0, m0), (1, m1)):
            selected = train[treatment[train] == value]
            dest[test] = predict_logit(fit_logit(x[selected], outcome[selected]), x[test])
    clip = float(spec.get("propensity_clip", 0.01))
    if not 0 < clip < 0.5:
        fail("propensity_clip must be between 0 and 0.5")
    pc = np.clip(p, clip, 1 - clip)
    influence = m1 - m0 + treatment * (outcome - m1) / pc - (1 - treatment) * (outcome - m0) / (1 - pc)
    ate = float(influence.mean())
    se = float(influence.std(ddof=1) / math.sqrt(n))
    weights = treatment / pc + (1 - treatment) / (1 - pc)
    balance = []
    labels = list(pd.concat([
        ((pd.to_numeric(frame[c]) - pd.to_numeric(frame[c]).mean()) / pd.to_numeric(frame[c]).std(ddof=0)).to_frame(c)
        if c in set(spec.get("numeric_covariates", [])) else
        pd.get_dummies(frame[c].astype("string"), prefix=c, drop_first=True, dtype=float)
        for c in covariates
    ], axis=1).columns) if covariates else []
    for j, name in enumerate(labels):
        a, b = x[treatment == 1, j], x[treatment == 0, j]
        sd = math.sqrt((a.var() + b.var()) / 2)
        raw = float((a.mean() - b.mean()) / sd) if sd else 0.0
        wa, wb = weights[treatment == 1], weights[treatment == 0]
        weighted = float((np.average(a, weights=wa) - np.average(b, weights=wb)) / sd) if sd else 0.0
        balance.append({"covariate": name, "smd_before": raw, "smd_after": weighted})
    return {
        "n": n, "treated": int(treatment.sum()), "controls": int(n - treatment.sum()),
        "ate_risk_difference": ate, "se": se,
        "ci95": [ate - 1.96 * se, ate + 1.96 * se],
        "p_value_wald": float(2 * norm.sf(abs(ate / se))) if se > 0 else None,
        "standardized_risk_if_all_treated": float(m1.mean()),
        "standardized_risk_if_none_treated": float(m0.mean()),
        "propensity_min": float(p.min()), "propensity_max": float(p.max()),
        "propensity_clipped_n": int(((p < clip) | (p > 1 - clip)).sum()),
        "effective_sample_size_weighted": float(weights.sum() ** 2 / np.square(weights).sum()),
        "max_absolute_weighted_smd": max((abs(b["smd_after"]) for b in balance), default=0.0),
        "balance": balance,
    }


def run(config_path, input_path, output_dir):
    config = json.loads(config_path.read_text())
    dag = config["dag"]
    nodes, edges = dag["nodes"], [tuple(e) for e in dag["edges"]]
    check_dag(nodes, edges)
    spec = config["analysis"]
    treatment, outcome = spec["treatment"], spec["outcome"]
    adjust = spec["adjust_for"]
    if treatment not in nodes or outcome not in nodes or treatment == outcome:
        fail("Treatment and outcome must be distinct DAG nodes")
    if [treatment, outcome] not in dag["edges"]:
        fail("The prespecified treatment -> outcome edge is required")
    if len(adjust) != len(set(adjust)) or not set(adjust) <= set(nodes) - {treatment, outcome}:
        fail("Adjustment variables must be distinct DAG nodes other than treatment/outcome")
    if set(adjust) & descendants(treatment, edges):
        fail("Adjustment set contains a treatment descendant")
    if not valid_backdoor(treatment, outcome, adjust, nodes, edges):
        fail("Adjustment set does not block every backdoor path in the supplied DAG")
    if not set(spec.get("numeric_covariates", [])) <= set(adjust):
        fail("numeric_covariates must be a subset of adjust_for")
    raw = read_input(input_path, config.get("sheet"))
    id_column = config.get("id_column")
    if id_column:
        if id_column not in raw.columns:
            fail("id_column is absent from the input")
        if raw[id_column].isna().any() or raw[id_column].duplicated().any():
            fail("id_column contains missing or duplicate values")
    mapping = config["columns"]
    needed = [treatment, outcome] + adjust
    if not set(needed) <= set(mapping):
        fail("Every analysis variable needs a source column mapping")
    sources = [mapping[n] for n in needed]
    if len(sources) != len(set(sources)) or not set(sources) <= set(raw.columns):
        fail("Source column mappings are duplicated or absent")
    frame = raw[sources].rename(columns={mapping[n]: n for n in needed}).copy()
    for name in (treatment, outcome):
        frame[name] = encode_binary(frame[name], config["binary_coding"][name], name)
    rows_input = len(frame)
    missing_by_variable = {name: int(frame[name].isna().sum()) for name in needed}
    frame = frame.dropna(subset=needed)
    if len(frame) < 30:
        fail("Fewer than 30 complete cases; analysis stopped")
    if not set(frame[treatment].unique()) == {0, 1} or not set(frame[outcome].unique()) == {0, 1}:
        fail("Treatment and outcome must each contain both binary classes")
    for name in spec.get("numeric_covariates", []):
        frame[name] = pd.to_numeric(frame[name], errors="raise")
    result = estimate(frame, spec, int(config.get("seed", 2026)))
    warnings = []
    if result["propensity_min"] < 0.05 or result["propensity_max"] > 0.95:
        warnings.append("Weak overlap: estimated propensity extends outside [0.05, 0.95].")
    if result["max_absolute_weighted_smd"] >= 0.1:
        warnings.append("Residual weighted covariate imbalance: absolute SMD >= 0.1.")
    if result["effective_sample_size_weighted"] < 0.5 * result["n"]:
        warnings.append("Weighted effective sample size is below half the complete-case sample.")
    result.update({"input_file": str(input_path.resolve()), "rows_input": rows_input,
                   "rows_excluded_missing": rows_input - len(frame), "estimand": "ATE, binary outcome risk difference",
                   "missing_by_variable": missing_by_variable,
                   "treatment": treatment, "outcome": outcome, "adjust_for": adjust,
                   "warnings": warnings,
                   "assumptions": "User-specified DAG and sufficient pre-treatment adjustment; no unmeasured confounding, consistency, positivity, correct nuisance models"})
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "results.json").write_text(json.dumps(result, indent=2))
    pd.DataFrame(result.pop("balance")).to_csv(output_dir / "balance.csv", index=False)
    # No row-level source or derived data are written.
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.config, args.data, args.out)
    print(f"n={result['n']}; ATE risk difference={result['ate_risk_difference']:.3f}; output={args.out}")

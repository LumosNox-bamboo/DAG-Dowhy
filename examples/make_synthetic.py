"""Generate fictional data for the README smoke test; never uses source patients."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit

parser = argparse.ArgumentParser()
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
rng = np.random.default_rng(19)
n = 1200
baseline = rng.normal(size=n)
exposure = rng.binomial(1, expit(-0.2 + 0.9 * baseline))
outcome = rng.binomial(1, expit(-0.8 + 0.8 * exposure + 0.7 * baseline))
args.out.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame({"baseline_risk": baseline, "exposure": exposure, "outcome": outcome}).to_csv(args.out, index=False)
print(args.out)

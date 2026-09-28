import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit

from causal_pipeline import check_dag, run, valid_backdoor


class PipelineTests(unittest.TestCase):
    def test_backdoor_and_cycle(self):
        nodes = ["c", "t", "y"]
        edges = [("c", "t"), ("c", "y"), ("t", "y")]
        self.assertFalse(valid_backdoor("t", "y", [], nodes, edges))
        self.assertTrue(valid_backdoor("t", "y", ["c"], nodes, edges))
        with self.assertRaises(ValueError):
            check_dag(nodes, edges + [("y", "c")])

    def test_synthetic_run_writes_aggregate_only(self):
        rng = np.random.default_rng(7)
        n = 900
        c = rng.normal(size=n)
        t = rng.binomial(1, expit(0.8 * c))
        y = rng.binomial(1, expit(-0.8 + 0.8 * t + 0.7 * c))
        config = json.loads((Path(__file__).parents[1] / "examples/config.json").read_text())
        config["columns"] = {"baseline_risk": "c", "exposure": "t", "outcome": "y"}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pd.DataFrame({"c": c, "t": t, "y": y}).to_csv(root / "input.csv", index=False)
            (root / "config.json").write_text(json.dumps(config))
            result = run(root / "config.json", root / "input.csv", root / "out")
            self.assertEqual(result["n"], n)
            self.assertGreater(result["ate_risk_difference"], 0.05)
            self.assertLess(result["ate_risk_difference"], 0.4)
            self.assertEqual(sorted(p.name for p in (root / "out").iterdir()), ["balance.csv", "results.json"])


if __name__ == "__main__":
    unittest.main()

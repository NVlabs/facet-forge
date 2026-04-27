#!/usr/bin/env python3
from pathlib import Path
import sys

SCRIPT = Path(__file__).resolve()
sys.path.insert(0, str(SCRIPT.parents[1] / "python"))

from facetforge_test import float_arg, int_arg, run_eval_sample_cli


if __name__ == "__main__":
    raise SystemExit(
        run_eval_sample_cli(
            source=SCRIPT.with_suffix(".cpp"),
            output=SCRIPT.with_suffix(".png"),
            title="Student-T vNDF",
            sample_label="vNDF.sample():",
            eval_label="vNDF.eval():",
            sample_title="vNDF.sample()",
            eval_title="vNDF.eval()",
            specs=[
                float_arg("rough_x", 0.8, "x roughness"),
                float_arg("rough_y", 0.6, "y roughness"),
                float_arg("gamma", 2.9, "Student-T shape parameter"),
                float_arg("theta_i", 0.8, "incident polar angle in radians"),
                float_arg("phi", 1.2, "incident azimuth in radians"),
                int_arg("samples", 1000000, "number of vNDF samples"),
                int_arg("eval_samples", 100, "eval samples per output bin"),
            ],
        )
    )

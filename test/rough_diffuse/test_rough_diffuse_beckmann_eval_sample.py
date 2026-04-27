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
            title="Rough diffuse Beckmann eval/sample",
            specs=[
                float_arg("rough_x", 0.5, "x roughness"),
                float_arg("rough_y", 0.1, "y roughness"),
                float_arg("kd", 0.7, "Lambert diffuse albedo"),
                float_arg("theta_i", 1.4, "incident polar angle in radians"),
                int_arg("samples", 1000000, "number of sample() calls"),
                int_arg("eval_samples", 100, "eval() samples per output bin"),
            ],
        )
    )

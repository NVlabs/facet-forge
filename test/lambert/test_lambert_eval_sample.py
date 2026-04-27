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
            title="Lambert eval/sample",
            specs=[
                float_arg("kd", 0.8, "Lambert diffuse albedo"),
                float_arg("theta_i", 0.0, "incident polar angle in radians"),
                int_arg("samples", 100000000, "number of sample() calls"),
                int_arg("eval_samples", 10, "eval() samples per output bin"),
            ],
        )
    )

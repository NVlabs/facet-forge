#!/usr/bin/env python3
from pathlib import Path
import sys

SCRIPT = Path(__file__).resolve()
sys.path.insert(0, str(SCRIPT.parents[1] / "python"))

from facetforge_test import float_arg, run_series_cli


if __name__ == "__main__":
    raise SystemExit(
        run_series_cli(
            source=SCRIPT.with_suffix(".cpp"),
            output=SCRIPT.with_suffix(".png"),
            title="Null Student-T sigma",
            labels=["quadrature", "monte carlo"],
            starts=(-0.999, -1.0),
            line_indices=(0,),
            specs=[
                float_arg("roughness", 0.75, "isotropic roughness"),
                float_arg("majorant", 3.0, "null-collision majorant"),
                float_arg("gamma", 2.8, "Student-T shape parameter"),
                float_arg("du", 0.01, "u step size"),
            ],
        )
    )

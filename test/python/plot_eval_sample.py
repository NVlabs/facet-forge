#!/usr/bin/env python3
"""Compatibility launcher for the first Python Lambert workflow."""

from __future__ import annotations

from pathlib import Path

from facetforge_test import float_arg, int_arg, run_eval_sample_cli


SCRIPT = Path(__file__).resolve()
REPO_ROOT = SCRIPT.parents[2]


if __name__ == "__main__":
    raise SystemExit(
        run_eval_sample_cli(
            source=REPO_ROOT / "test" / "lambert" / "test_lambert_eval_sample.cpp",
            output=Path("/tmp/facetforge-python-tests/lambert_eval_sample.png"),
            title="lambert",
            specs=[
                float_arg("kd", 0.8, "Lambert diffuse albedo"),
                float_arg("theta_i", 0.0, "incident polar angle in radians"),
                int_arg("samples", 100000000, "number of sample() calls"),
                int_arg("eval_samples", 10, "eval() samples per output bin"),
            ],
        )
    )

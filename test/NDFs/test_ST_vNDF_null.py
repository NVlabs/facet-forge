#!/usr/bin/env python3
import argparse
from pathlib import Path
import sys

import numpy as np

SCRIPT = Path(__file__).resolve()
sys.path.insert(0, str(SCRIPT.parents[1] / "python"))

from facetforge_test import (  # noqa: E402
    DEFAULT_DISPLAY_GAMMA,
    add_driver_args,
    arg_values,
    float_arg,
    int_arg,
    parse_grid,
    plot_grid_panels,
    run_compiled,
    summarize_grid,
    title_args,
)


SPECS = [
    float_arg("roughness", 0.8, "isotropic roughness"),
    float_arg("majorant", 0.6, "null-collision majorant"),
    float_arg("gamma", 2.7, "Student-T shape parameter"),
    float_arg("theta_i", 0.6, "incident polar angle in radians"),
    float_arg("phi", 1.2, "incident azimuth in radians"),
    int_arg("samples", 10000000, "number of vNDF samples"),
    int_arg("eval_samples", 100, "eval samples per output bin"),
]


def diff_stats(lhs, rhs):
    diff = lhs - rhs
    return {
        "mean_abs_error": float(np.mean(np.abs(diff))),
        "rms_error": float(np.sqrt(np.mean(diff * diff))),
        "max_abs_error": float(np.max(np.abs(diff))),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run the null Student-T vNDF test and compare it to the regular Student-T vNDF."
    )
    parser.add_argument("--compiler", default="g++", help="C++ compiler command")
    parser.add_argument("--output", type=Path, default=SCRIPT.with_suffix(".png"), help="path to the output PNG")
    parser.add_argument(
        "--display-gamma",
        type=float,
        default=DEFAULT_DISPLAY_GAMMA,
        help="display gamma for heatmaps; use 1 for linear display",
    )
    add_driver_args(parser, SPECS)
    args = parser.parse_args(argv)

    try:
        null_stdout = run_compiled(SCRIPT.with_suffix(".cpp"), args.compiler, arg_values(args, SPECS))
        null_sample = parse_grid(null_stdout, "vNDF.sample():")
        null_eval = parse_grid(null_stdout, "vNDF.eval():")

        st_source = SCRIPT.with_name("test_ST_vNDF.cpp")
        st_driver_args = [
            args.roughness,
            args.roughness,
            args.gamma,
            args.theta_i,
            args.phi,
            args.samples,
            args.eval_samples,
        ]
        st_stdout = run_compiled(st_source, args.compiler, st_driver_args)
        st_sample = parse_grid(st_stdout, "vNDF.sample():")
        st_eval = parse_grid(st_stdout, "vNDF.eval():")

        null_metrics = summarize_grid(null_sample, null_eval)
        st_metrics = summarize_grid(st_sample, st_eval)
        cross_metrics = diff_stats(null_sample, st_sample)
        panels = [
            ("null eval", null_eval),
            ("null sample", null_sample),
            ("Student-T eval", st_eval),
            ("Student-T sample", st_sample),
            ("abs(null eval - null sample)", np.abs(null_eval - null_sample)),
            ("abs(null sample - Student-T sample)", np.abs(null_sample - st_sample)),
        ]
        vmax = max(float(np.max(panel)) for panel in (null_eval, null_sample, st_eval, st_sample))
        metric_text = (
            f"null MAE: {null_metrics['mean_abs_error']:.3g}    "
            f"ST MAE: {st_metrics['mean_abs_error']:.3g}    "
            f"null/ST sample MAE: {cross_metrics['mean_abs_error']:.3g}    "
            f"display gamma: {args.display_gamma:.4g}"
        )
        plot_grid_panels(
            panels,
            args.output,
            f"Null Student-T vNDF: {title_args(args, SPECS)}",
            args.display_gamma,
            metric_text,
            vmax=max(vmax, 1e-12),
        )
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"wrote {args.output}")
    for prefix, metrics in (
        ("null", null_metrics),
        ("student_t", st_metrics),
        ("null_vs_student_t_sample", cross_metrics),
    ):
        for name, value in metrics.items():
            print(f"{prefix}_{name}: {value:.12g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

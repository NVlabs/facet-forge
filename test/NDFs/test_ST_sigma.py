#!/usr/bin/env python3
import math
from pathlib import Path
import sys

import numpy as np
from scipy import special

SCRIPT = Path(__file__).resolve()
sys.path.insert(0, str(SCRIPT.parents[1] / "python"))

from facetforge_test import float_arg, run_series_cli


def student_t_sigma_reference(us, args):
    u = np.asarray(us, dtype=np.float64)
    alpha = float(args.rough_x)
    gamma = float(args.gamma)
    result = np.empty_like(u)

    lower = u <= -1.0
    upper = u >= 1.0
    interior = ~(lower | upper)
    result[lower] = 0.0
    result[upper] = 1.0

    if np.any(interior):
        ui = u[interior]
        one_minus_u2 = 1.0 - ui * ui
        z = (ui * ui) / ((ui * ui - 1.0) * alpha * alpha * (gamma - 1.0))
        gamma_minus_1 = special.gamma(gamma - 1.0)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            term0 = 0.5 * ui
            term1 = (
                alpha
                * np.power(1.0 - z, 1.5 - gamma)
                * np.sqrt(one_minus_u2 * (gamma - 1.0))
                * special.gamma(gamma - 1.5)
                / (2.0 * math.sqrt(math.pi) * gamma_minus_1)
            )
            term2 = (
                ui
                * ui
                * special.gamma(gamma - 0.5)
                * special.hyp2f1(0.5, gamma - 0.5, 1.5, z)
                / (
                    math.sqrt(math.pi)
                    * alpha
                    * np.sqrt(one_minus_u2 * (gamma - 1.0))
                    * gamma_minus_1
                )
            )
        result[interior] = np.real_if_close(term0 + term1 + term2)

    return result


if __name__ == "__main__":
    raise SystemExit(
        run_series_cli(
            source=SCRIPT.with_suffix(".cpp"),
            output=SCRIPT.with_suffix(".png"),
            title="Student-T sigma",
            labels=["driver sigma"],
            reference=student_t_sigma_reference,
            reference_label="analytic sigma",
            specs=[
                float_arg("rough_x", 0.8, "x roughness"),
                float_arg("rough_y", 0.6, "y roughness"),
                float_arg("gamma", 2.1, "Student-T shape parameter"),
                float_arg("du", 0.01, "u step size"),
                float_arg("phi", 2.2, "azimuth in radians; kept for parity with C++ usage"),
            ],
        )
    )

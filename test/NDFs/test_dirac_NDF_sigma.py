#!/usr/bin/env python3
import math
from pathlib import Path
import sys

import numpy as np

SCRIPT = Path(__file__).resolve()
sys.path.insert(0, str(SCRIPT.parents[1] / "python"))

from facetforge_test import float_arg, run_series_cli


def dirac_sigma_reference(us, args):
    u = np.asarray(us, dtype=np.float64)
    un = float(args.un)
    result = np.zeros_like(u)
    same_hemisphere = u * un > 0.0
    result[same_hemisphere] += 2.0 * math.pi * u[same_hemisphere] * un

    normal_radius_sq = max(0.0, 1.0 - un * un)
    inside = np.abs(u) < math.sqrt(normal_radius_sq)
    if np.any(inside) and normal_radius_sq > 0.0:
        ui = u[inside]
        root = np.sqrt(np.maximum(0.0, 1.0 - un * un - ui * ui))
        denom = np.sqrt(np.maximum(1e-30, (1.0 - ui * ui) * normal_radius_sq))
        acos_arg = np.clip(-ui * un / denom, -1.0, 1.0)
        result[inside] += 2.0 * root + 2.0 * ui * un * np.arccos(acos_arg)
        result[inside & same_hemisphere] -= 2.0 * math.pi * u[inside & same_hemisphere] * un

    return result


if __name__ == "__main__":
    raise SystemExit(
        run_series_cli(
            source=SCRIPT.with_suffix(".cpp"),
            output=SCRIPT.with_suffix(".png"),
            title="Dirac NDF sigma",
            labels=["driver sigma"],
            reference=dirac_sigma_reference,
            reference_label="analytic sigma",
            specs=[
                float_arg("du", 0.01, "u step size"),
                float_arg("un", 0.4, "Dirac microfacet normal z component"),
            ],
        )
    )

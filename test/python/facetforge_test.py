#!/usr/bin/env python3
"""Shared Python helpers for FacetForge C++ reference tests."""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np


NUM_THETA = 100
NUM_PHI = 101
DEFAULT_DISPLAY_GAMMA = 1.0 / 0.3
REPO_ROOT = Path(__file__).resolve().parents[2]


class SignedPowerNorm(colors.Normalize):
    """Symmetric power-law normalization for signed difference images."""

    def __init__(self, gamma: float, vmax: float):
        super().__init__(vmin=-vmax, vmax=vmax, clip=True)
        self.gamma = gamma
        self._vmax = max(vmax, 1e-12)

    def __call__(self, value, clip=None):
        data = np.ma.asarray(value)
        scaled = np.ma.clip(data / self._vmax, -1.0, 1.0)
        signed = np.sign(scaled) * np.ma.power(np.ma.abs(scaled), self.gamma)
        return 0.5 + 0.5 * signed

    def inverse(self, value):
        scaled = 2.0 * np.asarray(value) - 1.0
        signed = np.sign(scaled) * np.power(np.abs(scaled), 1.0 / self.gamma)
        return self._vmax * signed


@dataclass(frozen=True)
class ArgSpec:
    name: str
    default: float | int
    type: Callable[[str], float | int]
    help: str

    @property
    def flag(self) -> str:
        return "--" + self.name.replace("_", "-")

    @property
    def dest(self) -> str:
        return self.name.replace("-", "_")


def float_arg(name: str, default: float, help_text: str) -> ArgSpec:
    return ArgSpec(name=name, default=default, type=float, help=help_text)


def int_arg(name: str, default: int, help_text: str) -> ArgSpec:
    return ArgSpec(name=name, default=default, type=int, help=help_text)


def compile_driver(source: Path, output: Path, compiler: str) -> None:
    command = [
        compiler,
        "-std=c++11",
        "-I",
        str(REPO_ROOT / "include"),
        str(source),
        "-o",
        str(output),
    ]
    subprocess.run(command, cwd=REPO_ROOT, check=True)


def run_driver(binary: Path, driver_args: Sequence[float | int]) -> str:
    command = [str(binary), *(str(value) for value in driver_args)]
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )

    # Several current C++ test drivers return 1 after successful output.
    if completed.returncode not in (0, 1):
        raise RuntimeError(
            f"driver failed with exit code {completed.returncode}\n"
            f"stderr:\n{completed.stderr}"
        )

    return completed.stdout


def parse_grid(output: str, label: str) -> np.ndarray:
    lines = output.splitlines()
    try:
        start = lines.index(label) + 1
    except ValueError as exc:
        raise ValueError(f"could not find output label {label!r}") from exc

    values: list[float] = []
    for line in lines[start:]:
        if line.endswith(":") and values:
            break
        if not line.strip():
            if values:
                break
            continue
        values.extend(float(token) for token in line.split())
        if len(values) >= NUM_THETA * NUM_PHI:
            break

    expected = NUM_THETA * NUM_PHI
    if len(values) != expected:
        raise ValueError(f"{label} had {len(values)} values, expected {expected}")

    return np.asarray(values, dtype=np.float64).reshape(NUM_THETA, NUM_PHI)


def parse_numeric_rows(output: str) -> list[np.ndarray]:
    rows: list[np.ndarray] = []
    for line in output.splitlines():
        tokens = line.split()
        if not tokens:
            continue
        try:
            values = [float(token) for token in tokens]
        except ValueError:
            continue
        rows.append(np.asarray(values, dtype=np.float64))
    return rows


def summarize_grid(sampled: np.ndarray, evaluated: np.ndarray) -> dict[str, float]:
    diff = sampled - evaluated
    return {
        "sample_sum": float(sampled.sum()),
        "eval_sum": float(evaluated.sum()),
        "mean_abs_error": float(np.mean(np.abs(diff))),
        "rms_error": float(math.sqrt(np.mean(diff * diff))),
        "max_abs_error": float(np.max(np.abs(diff))),
    }


def summarize_rows(rows: Sequence[np.ndarray]) -> dict[str, float]:
    if len(rows) < 2 or rows[0].shape != rows[1].shape:
        return {}

    diff = rows[0] - rows[1]
    return {
        "mean_abs_error": float(np.mean(np.abs(diff))),
        "rms_error": float(math.sqrt(np.mean(diff * diff))),
        "max_abs_error": float(np.max(np.abs(diff))),
    }


def plot_grid_comparison(
    sampled: np.ndarray,
    evaluated: np.ndarray,
    output_path: Path,
    title: str,
    metrics: dict[str, float],
    sample_title: str,
    eval_title: str,
    display_gamma: float,
) -> None:
    diff = sampled - evaluated
    vmax = max(float(sampled.max()), float(evaluated.max()), 1e-12)
    diff_vmax = max(float(np.max(np.abs(diff))), 1e-12)
    extent = [-math.pi, math.pi, -0.5 * math.pi, 0.5 * math.pi]
    display_exponent = 1.0 / display_gamma if display_gamma > 0.0 else 1.0

    fig = plt.figure(figsize=(17.0, 6.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 4, width_ratios=[1.0, 1.0, 1.0, 1.05])
    axes = [
        fig.add_subplot(grid[:, 0]),
        fig.add_subplot(grid[:, 1]),
        fig.add_subplot(grid[:, 2]),
    ]
    theta_axis = fig.add_subplot(grid[0, 3])
    phi_axis = fig.add_subplot(grid[1, 3])

    positive_norm = colors.PowerNorm(gamma=display_exponent, vmin=0.0, vmax=vmax)
    diff_norm = SignedPowerNorm(gamma=display_exponent, vmax=diff_vmax)
    images = [
        axes[0].imshow(sampled, origin="lower", aspect="auto", extent=extent, norm=positive_norm),
        axes[1].imshow(evaluated, origin="lower", aspect="auto", extent=extent, norm=positive_norm),
        axes[2].imshow(diff, origin="lower", aspect="auto", extent=extent, norm=diff_norm, cmap="coolwarm"),
    ]

    axes[0].set_title(sample_title)
    axes[1].set_title(eval_title)
    axes[2].set_title("sample - eval")
    for axis in axes:
        axis.set_xlabel("phi")
        axis.set_ylabel("theta")

    fig.colorbar(images[0], ax=axes[:2], shrink=0.85, label="histogram mass")
    fig.colorbar(images[2], ax=axes[2], shrink=0.85, label="difference")

    theta_centers = -0.5 * math.pi + (np.arange(NUM_THETA) + 0.5) * math.pi / NUM_THETA
    phi_centers = -math.pi + (np.arange(NUM_PHI) + 0.5) * 2.0 * math.pi / NUM_PHI

    theta_axis.plot(theta_centers, sampled.sum(axis=1), label="sample", linewidth=1.3)
    theta_axis.plot(theta_centers, evaluated.sum(axis=1), label="eval", linewidth=1.3)
    theta_axis.set_title("theta marginal")
    theta_axis.set_xlabel("theta")
    theta_axis.set_ylabel("mass")
    theta_axis.grid(True, alpha=0.25)
    theta_axis.legend()

    phi_axis.plot(phi_centers, sampled.sum(axis=0), label="sample", linewidth=1.3)
    phi_axis.plot(phi_centers, evaluated.sum(axis=0), label="eval", linewidth=1.3)
    phi_axis.set_title("phi marginal")
    phi_axis.set_xlabel("phi")
    phi_axis.set_ylabel("mass")
    phi_axis.grid(True, alpha=0.25)
    phi_axis.legend()

    metric_text = (
        f"sample sum: {metrics['sample_sum']:.6g}    "
        f"eval sum: {metrics['eval_sum']:.6g}    "
        f"MAE: {metrics['mean_abs_error']:.3g}    "
        f"RMSE: {metrics['rms_error']:.3g}    "
        f"display gamma: {display_gamma:.4g}"
    )
    fig.suptitle(f"{title}\n{metric_text}", fontsize=11)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def plot_grid_panels(
    panels: Sequence[tuple[str, np.ndarray]],
    output_path: Path,
    title: str,
    display_gamma: float,
    metric_text: str = "",
    vmax: float | None = None,
) -> None:
    if not panels:
        raise ValueError("no panels to plot")

    cols = min(3, len(panels))
    rows = int(math.ceil(len(panels) / cols))
    extent = [-math.pi, math.pi, -0.5 * math.pi, 0.5 * math.pi]
    display_exponent = 1.0 / display_gamma if display_gamma > 0.0 else 1.0
    if vmax is None:
        vmax = max(float(np.max(np.abs(panel))) for _, panel in panels)
    positive_norm = colors.PowerNorm(gamma=display_exponent, vmin=0.0, vmax=max(vmax, 1e-12))

    fig, axes = plt.subplots(rows, cols, figsize=(5.2 * cols, 4.2 * rows), constrained_layout=True)
    axes_array = np.asarray(axes).reshape(-1)
    image = None
    for axis, (panel_title, panel) in zip(axes_array, panels):
        image = axis.imshow(np.abs(panel), origin="lower", aspect="auto", extent=extent, norm=positive_norm)
        axis.set_title(panel_title)
        axis.set_xlabel("phi")
        axis.set_ylabel("theta")

    for axis in axes_array[len(panels):]:
        axis.set_visible(False)

    if image is not None:
        fig.colorbar(image, ax=axes_array[: len(panels)], shrink=0.85, label="mass / absolute error")

    if metric_text:
        fig.suptitle(f"{title}\n{metric_text}", fontsize=11)
    else:
        fig.suptitle(title, fontsize=11)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def plot_rows(
    rows: Sequence[np.ndarray],
    output_path: Path,
    title: str,
    du: float,
    labels: Sequence[str],
    starts: Sequence[float],
    metrics: dict[str, float],
    reference_indices: Sequence[int] = (),
    line_indices: Sequence[int] = (),
    plot_offset: float = 0.5,
) -> None:
    fig, axis = plt.subplots(1, 1, figsize=(8.5, 4.8), constrained_layout=True)
    reference_index_set = set(reference_indices)
    line_index_set = set(line_indices)

    for index, row in enumerate(rows):
        start = starts[index] if index < len(starts) else -1.0
        label = labels[index] if index < len(labels) else f"row {index + 1}"
        if index in reference_index_set:
            xs = start + du * np.arange(row.size)
            axis.plot(xs, row, linewidth=1.4, label=label)
        elif index in line_index_set:
            xs = start + du * (np.arange(row.size) + plot_offset)
            axis.plot(xs, row, linewidth=1.2, label=label)
        else:
            xs = start + du * (np.arange(row.size) + plot_offset)
            point_indices = np.arange(0, row.size, 2)
            axis.plot(
                xs[point_indices],
                row[point_indices],
                linestyle="None",
                marker="o",
                markersize=3.6,
                markeredgewidth=0.8,
                label=label,
            )

    axis.set_xlabel("u = cos(theta)")
    axis.set_ylabel("value")
    axis.grid(True, alpha=0.25)
    axis.legend()

    if metrics:
        metric_text = (
            f"MAE: {metrics['mean_abs_error']:.3g}    "
            f"RMSE: {metrics['rms_error']:.3g}    "
            f"max abs: {metrics['max_abs_error']:.3g}"
        )
        fig.suptitle(f"{title}\n{metric_text}", fontsize=11)
    else:
        fig.suptitle(title, fontsize=11)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def add_driver_args(parser: argparse.ArgumentParser, specs: Sequence[ArgSpec]) -> None:
    for spec in specs:
        parser.add_argument(spec.flag, dest=spec.dest, type=spec.type, default=spec.default, help=spec.help)


def arg_values(args: argparse.Namespace, specs: Sequence[ArgSpec]) -> list[float | int]:
    return [getattr(args, spec.dest) for spec in specs]


def title_args(args: argparse.Namespace, specs: Sequence[ArgSpec]) -> str:
    return ", ".join(f"{spec.name}={getattr(args, spec.dest)}" for spec in specs)


def run_compiled(source: Path, compiler: str, driver_args: Sequence[float | int]) -> str:
    with tempfile.TemporaryDirectory(prefix="facetforge-test-") as temp_dir:
        binary = Path(temp_dir) / source.stem
        compile_driver(source, binary, compiler)
        return run_driver(binary, driver_args)


def run_eval_sample_cli(
    *,
    source: Path,
    specs: Sequence[ArgSpec],
    output: Path,
    title: str,
    sample_label: str = "BSDF.sample():",
    eval_label: str = "BSDF.eval():",
    sample_title: str = "sample()",
    eval_title: str = "eval()",
    argv: Sequence[str] | None = None,
) -> int:
    parser = argparse.ArgumentParser(description=f"Run and plot {source.name}.")
    parser.add_argument("--compiler", default="g++", help="C++ compiler command")
    parser.add_argument("--output", type=Path, default=output, help="path to the output PNG")
    parser.add_argument(
        "--display-gamma",
        type=float,
        default=DEFAULT_DISPLAY_GAMMA,
        help="display gamma for heatmaps; use 1 for linear display",
    )
    add_driver_args(parser, specs)
    args = parser.parse_args(argv)

    try:
        stdout = run_compiled(source, args.compiler, arg_values(args, specs))
        sampled = parse_grid(stdout, sample_label)
        evaluated = parse_grid(stdout, eval_label)
        metrics = summarize_grid(sampled, evaluated)
        plot_grid_comparison(
            sampled,
            evaluated,
            args.output,
            f"{title}: {title_args(args, specs)}",
            metrics,
            sample_title,
            eval_title,
            args.display_gamma,
        )
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"wrote {args.output}")
    for name, value in metrics.items():
        print(f"{name}: {value:.12g}")
    return 0


def run_series_cli(
    *,
    source: Path,
    specs: Sequence[ArgSpec],
    output: Path,
    title: str,
    labels: Sequence[str],
    starts: Sequence[float] = (-1.0,),
    reference: Callable[[np.ndarray, argparse.Namespace], np.ndarray] | None = None,
    reference_label: str = "reference",
    line_indices: Sequence[int] = (),
    plot_offset: float = 0.5,
    argv: Sequence[str] | None = None,
) -> int:
    parser = argparse.ArgumentParser(description=f"Run and plot {source.name}.")
    parser.add_argument("--compiler", default="g++", help="C++ compiler command")
    parser.add_argument("--output", type=Path, default=output, help="path to the output PNG")
    add_driver_args(parser, specs)
    args = parser.parse_args(argv)

    try:
        stdout = run_compiled(source, args.compiler, arg_values(args, specs))
        rows = parse_numeric_rows(stdout)
        if not rows:
            raise ValueError("driver produced no numeric rows")
        du = float(getattr(args, "du"))
        if reference is not None:
            start = starts[0] if starts else -1.0
            xs = start + du * np.arange(rows[0].size)
            reference_row = np.asarray(reference(xs, args), dtype=np.float64)
            if reference_row.shape != rows[0].shape:
                raise ValueError(
                    f"reference row had shape {reference_row.shape}, expected {rows[0].shape}"
                )
            rows = [rows[0], reference_row, *rows[1:]]
            labels = [labels[0] if labels else "driver", reference_label, *labels[1:]]
            starts = [start, start, *starts[1:]]
            reference_indices = (1,)
            line_indices = tuple(index + 1 if index >= 1 else index for index in line_indices)
        else:
            reference_indices = ()
        metrics = summarize_rows(rows)
        plot_rows(
            rows,
            args.output,
            f"{title}: {title_args(args, specs)}",
            du,
            labels,
            starts,
            metrics,
            reference_indices=reference_indices,
            line_indices=line_indices,
            plot_offset=plot_offset,
        )
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"wrote {args.output}")
    for index, row in enumerate(rows):
        label = labels[index] if index < len(labels) else f"row {index + 1}"
        print(f"{label}_count: {row.size}")
    for name, value in metrics.items():
        print(f"{name}: {value:.12g}")
    return 0

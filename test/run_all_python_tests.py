#!/usr/bin/env python3
"""Run every Python FacetForge test wrapper and build an HTML image report."""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as _datetime
import fnmatch
import html
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


SCRIPT = Path(__file__).resolve()
TEST_ROOT = SCRIPT.parent
REPO_ROOT = TEST_ROOT.parent
DEFAULT_OUTPUT_DIR = TEST_ROOT / "_python_report" / "latest"
METRIC_RE = re.compile(r"^([A-Za-z0-9_]+):\s*([-+0-9.eE]+)\s*$")


@dataclass(frozen=True)
class TestCase:
    script: Path
    relative: Path
    output: Path


@dataclass(frozen=True)
class TestResult:
    case: TestCase
    command: list[str]
    passed: bool
    timed_out: bool
    elapsed_seconds: float
    returncode: int | None
    stdout: str
    stderr: str
    metrics: dict[str, float]


def discover_tests() -> list[Path]:
    tests: list[Path] = []
    helper_dir = TEST_ROOT / "python"
    for path in sorted(TEST_ROOT.rglob("test_*.py")):
        try:
            path.relative_to(helper_dir)
            in_helper_dir = True
        except ValueError:
            in_helper_dir = False
        if in_helper_dir:
            continue
        tests.append(path)
    return tests


def select_tests(tests: Iterable[Path], patterns: list[str]) -> list[Path]:
    if not patterns:
        return list(tests)

    selected: list[Path] = []
    for path in tests:
        rel = path.relative_to(TEST_ROOT).as_posix()
        name = path.name
        if any(fnmatch.fnmatch(rel, pattern) or fnmatch.fnmatch(name, pattern) or pattern in rel for pattern in patterns):
            selected.append(path)
    return selected


def parse_test_args(raw_args: list[str]) -> list[str]:
    parsed: list[str] = []
    for raw in raw_args:
        if "=" in raw and not raw.startswith("--"):
            name, value = raw.split("=", 1)
            parsed.extend(["--" + name.replace("_", "-"), value])
        else:
            parsed.append(raw)
    return parsed


def safe_clean(output_dir: Path) -> None:
    resolved = output_dir.resolve()
    repo = REPO_ROOT.resolve()
    if resolved == repo or repo in resolved.parents and resolved.name in ("test", "include", "assets"):
        raise ValueError(f"refusing to clean suspicious output directory: {resolved}")
    if resolved.exists():
        shutil.rmtree(resolved)


def metric_lines(stdout: str) -> dict[str, float]:
    metrics: dict[str, float] = {}
    for line in stdout.splitlines():
        match = METRIC_RE.match(line.strip())
        if match:
            metrics[match.group(1)] = float(match.group(2))
    return metrics


def build_case(script: Path, output_dir: Path) -> TestCase:
    relative = script.relative_to(TEST_ROOT)
    output = output_dir / relative.parent / (script.stem + ".png")
    return TestCase(script=script, relative=relative, output=output)


def run_case(case: TestCase, compiler: str, extra_args: list[str], timeout: float | None) -> TestResult:
    case.output.parent.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable,
        str(case.script),
        "--compiler",
        compiler,
        "--output",
        str(case.output),
        *extra_args,
    ]

    start = time.perf_counter()
    try:
        completed = subprocess.run(
            command,
            cwd=REPO_ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )
        elapsed = time.perf_counter() - start
        passed = completed.returncode == 0 and case.output.exists()
        return TestResult(
            case=case,
            command=command,
            passed=passed,
            timed_out=False,
            elapsed_seconds=elapsed,
            returncode=completed.returncode,
            stdout=completed.stdout,
            stderr=completed.stderr,
            metrics=metric_lines(completed.stdout),
        )
    except subprocess.TimeoutExpired as exc:
        elapsed = time.perf_counter() - start
        stdout = exc.stdout if isinstance(exc.stdout, str) else ""
        stderr = exc.stderr if isinstance(exc.stderr, str) else ""
        return TestResult(
            case=case,
            command=command,
            passed=False,
            timed_out=True,
            elapsed_seconds=elapsed,
            returncode=None,
            stdout=stdout,
            stderr=stderr,
            metrics=metric_lines(stdout),
        )


def relative_href(target: Path, base: Path) -> str:
    return os.path.relpath(target, base).replace(os.sep, "/")


def format_seconds(seconds: float) -> str:
    if seconds < 60.0:
        return f"{seconds:.1f}s"
    minutes, sec = divmod(seconds, 60.0)
    return f"{int(minutes)}m {sec:.0f}s"


def metric_summary(metrics: dict[str, float]) -> str:
    preferred = [
        "sample_sum",
        "eval_sum",
        "mean_abs_error",
        "rms_error",
        "max_abs_error",
    ]
    names = [name for name in preferred if name in metrics]
    names.extend(name for name in sorted(metrics) if name not in names and not name.endswith("_count"))
    if not names:
        return ""
    parts = [f"{html.escape(name)}: {metrics[name]:.6g}" for name in names[:8]]
    return " · ".join(parts)


def write_report(results: list[TestResult], report_path: Path, started_at: _datetime.datetime) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    total = len(results)
    passed = sum(1 for result in results if result.passed)
    failed = total - passed
    elapsed = sum(result.elapsed_seconds for result in results)
    generated = _datetime.datetime.now().astimezone()

    nav_items = []
    cards = []
    for index, result in enumerate(results, 1):
        rel = result.case.relative.as_posix()
        status_class = "pass" if result.passed else "fail"
        status_text = "PASS" if result.passed else ("TIMEOUT" if result.timed_out else "FAIL")
        anchor = f"test-{index}"
        metrics = metric_summary(result.metrics)
        nav_items.append(
            f'<a class="pill {status_class}" href="#{anchor}">'
            f"<span>{status_text}</span>{html.escape(rel)}</a>"
        )

        image_html = ""
        if result.case.output.exists():
            href = relative_href(result.case.output, report_path.parent)
            image_html = f'<a href="{html.escape(href)}"><img src="{html.escape(href)}" alt="{html.escape(rel)}"></a>'
        else:
            image_html = '<div class="missing">No PNG was produced.</div>'

        stdout = html.escape(result.stdout.strip())
        stderr = html.escape(result.stderr.strip())
        command = html.escape(" ".join(result.command))
        returncode = "timeout" if result.timed_out else str(result.returncode)
        details = ""
        if stdout or stderr:
            details = (
                "<details><summary>log</summary>"
                f"<p><strong>command</strong><br><code>{command}</code></p>"
                f"<p><strong>return code</strong>: {html.escape(returncode)}</p>"
                f"<pre>{stdout}</pre>"
                + (f"<pre class=\"stderr\">{stderr}</pre>" if stderr else "")
                + "</details>"
            )

        cards.append(
            f'<section class="card {status_class}" id="{anchor}">'
            f"<header><div><h2>{html.escape(rel)}</h2>"
            f'<p>{html.escape(metrics) if metrics else "Completed with no numeric summary."}</p></div>'
            f'<span class="badge {status_class}">{status_text}</span></header>'
            f"{image_html}"
            f'<footer>{format_seconds(result.elapsed_seconds)}</footer>'
            f"{details}"
            "</section>"
        )

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>FacetForge Python Test Report</title>
  <style>
    :root {{
      color-scheme: light;
      --blue: #1559a8;
      --blue-bg: #eaf3ff;
      --blue-border: #8bb9ea;
      --red: #a32121;
      --red-bg: #fff0f0;
      --red-border: #ec9b9b;
      --ink: #172033;
      --muted: #667085;
      --line: #d9dee8;
      --paper: #f7f8fb;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font: 14px/1.45 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      color: var(--ink);
      background: var(--paper);
    }}
    .top {{
      position: sticky;
      top: 0;
      z-index: 10;
      padding: 16px 22px 14px;
      background: rgba(247, 248, 251, 0.96);
      border-bottom: 1px solid var(--line);
      backdrop-filter: blur(10px);
    }}
    h1 {{
      margin: 0 0 6px;
      font-size: 22px;
      font-weight: 700;
      letter-spacing: 0;
    }}
    .meta {{
      margin: 0 0 12px;
      color: var(--muted);
    }}
    .summary {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      max-height: 145px;
      overflow: auto;
      padding-right: 4px;
    }}
    .pill {{
      display: inline-flex;
      align-items: center;
      gap: 7px;
      max-width: 360px;
      padding: 6px 9px;
      border: 1px solid;
      border-radius: 8px;
      color: inherit;
      text-decoration: none;
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
    }}
    .pill span, .badge {{
      font-size: 11px;
      font-weight: 700;
      letter-spacing: 0.04em;
    }}
    .pass {{ border-color: var(--blue-border); background: var(--blue-bg); }}
    .fail {{ border-color: var(--red-border); background: var(--red-bg); }}
    .pill.pass span, .badge.pass {{ color: var(--blue); }}
    .pill.fail span, .badge.fail {{ color: var(--red); }}
    main {{
      padding: 22px;
      display: grid;
      gap: 18px;
    }}
    .card {{
      background: white;
      border: 1px solid var(--line);
      border-left-width: 6px;
      border-radius: 8px;
      padding: 14px;
      box-shadow: 0 1px 2px rgba(20, 25, 40, 0.05);
    }}
    .card.pass {{ border-left-color: var(--blue); background: white; }}
    .card.fail {{ border-left-color: var(--red); background: white; }}
    .card header {{
      display: flex;
      align-items: flex-start;
      justify-content: space-between;
      gap: 12px;
      margin-bottom: 12px;
    }}
    .card h2 {{
      margin: 0 0 4px;
      font-size: 17px;
      letter-spacing: 0;
    }}
    .card p {{
      margin: 0;
      color: var(--muted);
    }}
    .badge {{
      border: 1px solid;
      border-radius: 999px;
      padding: 4px 8px;
      flex: 0 0 auto;
    }}
    img {{
      width: 100%;
      max-height: 760px;
      object-fit: contain;
      border: 1px solid var(--line);
      border-radius: 6px;
      background: #fff;
    }}
    footer {{
      margin-top: 8px;
      color: var(--muted);
    }}
    details {{
      margin-top: 10px;
      border-top: 1px solid var(--line);
      padding-top: 8px;
    }}
    code, pre {{
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 12px;
    }}
    pre {{
      overflow: auto;
      max-height: 260px;
      padding: 10px;
      background: #f4f6fa;
      border-radius: 6px;
    }}
    .stderr {{
      background: #fff4f4;
    }}
    .missing {{
      padding: 48px;
      text-align: center;
      color: var(--red);
      border: 1px dashed var(--red-border);
      border-radius: 6px;
      background: var(--red-bg);
    }}
  </style>
</head>
<body>
  <div class="top">
    <h1>FacetForge Python Test Report</h1>
    <p class="meta">
      {passed}/{total} passed, {failed} failed · generated {html.escape(generated.strftime("%Y-%m-%d %H:%M:%S %Z"))}
      · started {html.escape(started_at.strftime("%Y-%m-%d %H:%M:%S %Z"))}
      · summed runtime {html.escape(format_seconds(elapsed))}
    </p>
    <div class="summary">
      {''.join(nav_items)}
    </div>
  </div>
  <main>
    {''.join(cards)}
  </main>
</body>
</html>
"""
    report_path.write_text(html_text, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run all FacetForge Python test wrappers and write an HTML report.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help="directory for PNGs and index.html")
    parser.add_argument("--compiler", default="g++", help="C++ compiler command passed to each wrapper")
    parser.add_argument("--jobs", type=int, default=1, help="number of tests to run concurrently")
    parser.add_argument("--timeout", type=float, default=None, help="per-test timeout in seconds")
    parser.add_argument("--only", action="append", default=[], help="run tests matching this glob or substring; repeatable")
    parser.add_argument("--test-arg", action="append", default=[], help="extra wrapper argument, e.g. samples=100000")
    parser.add_argument("--clean", action="store_true", help="delete the output directory before running")
    parser.add_argument("--list", action="store_true", help="list selected tests and exit")
    args = parser.parse_args(argv)

    tests = select_tests(discover_tests(), args.only)
    if not tests:
        print("no tests matched", file=sys.stderr)
        return 2

    if args.list:
        for test in tests:
            print(test.relative_to(TEST_ROOT).as_posix())
        return 0

    output_dir = args.output_dir.resolve()
    if args.clean:
        safe_clean(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    extra_args = parse_test_args(args.test_arg)
    cases = [build_case(test, output_dir) for test in tests]
    started_at = _datetime.datetime.now().astimezone()
    jobs = max(1, args.jobs)

    print(f"running {len(cases)} tests with {jobs} job(s)")
    results: list[TestResult] = []
    if jobs == 1:
        for case in cases:
            print(f"[run] {case.relative.as_posix()}")
            result = run_case(case, args.compiler, extra_args, args.timeout)
            print(f"[{'pass' if result.passed else 'fail'}] {case.relative.as_posix()} ({format_seconds(result.elapsed_seconds)})")
            results.append(result)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
            future_to_case = {
                executor.submit(run_case, case, args.compiler, extra_args, args.timeout): case
                for case in cases
            }
            for future in concurrent.futures.as_completed(future_to_case):
                case = future_to_case[future]
                result = future.result()
                print(f"[{'pass' if result.passed else 'fail'}] {case.relative.as_posix()} ({format_seconds(result.elapsed_seconds)})")
                results.append(result)
        results.sort(key=lambda result: result.case.relative.as_posix())

    report_path = output_dir / "index.html"
    write_report(results, report_path, started_at)
    print(f"wrote {report_path}")
    passed = sum(1 for result in results if result.passed)
    print(f"{passed}/{len(results)} passed")
    return 0 if passed == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())

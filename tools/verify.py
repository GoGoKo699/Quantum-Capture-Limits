#!/usr/bin/env python3
"""Run all scientific checks and validate the reader-facing repository.

Use --record to refresh verification/validated.json after a successful run.
Raw reports are outputs, never inputs to the scientific calculations.
"""
from __future__ import annotations
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
SUITES = (("simulation", 4), ("control", 5), ("capacity", 5), ("objectives", 4), ("upper_bound", 4))


def source_hashes():
    files = sorted(p for d in ("capture", "tests", "tools") for p in (ROOT/d).glob("*.py"))
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def structural_checks():
    parsed = 0
    for path in source_hashes():
        ast.parse((ROOT/path).read_text(encoding="utf-8"), filename=path)
        parsed += 1
    links = 0
    for f in sorted(ROOT.rglob("*.md")):
        if any(part.startswith(".") or part in ("artifacts", "venv") for part in f.relative_to(ROOT).parts):
            continue
        text = f.read_text(encoding="utf-8")
        if text.count("```") % 2:
            raise ValueError(f"Unbalanced Markdown fences: {f}")
        if re.search(r"\\[\[\]()]", text):
            raise ValueError(f"Use GitHub-rendered inline/display math: {f}")
        if f.parent.name == "docs" and re.search(r"\b(pilot|round[- ]\d|PRL)\b", text, re.I):
            raise ValueError(f"Stale exploration or journal label: {f}")
        for url in re.findall(r"https?://[^\s<>]+", text):
            host = url.split("/", 3)[2]
            if not re.fullmatch(r"[A-Za-z0-9.-]+(?::[0-9]+)?", host):
                raise ValueError(f"Malformed external URL in {f}: {url}")
        # Only the path is checked; heading anchors can change with renderer rules.
        for target in re.findall(r"\]\(([^)]+)\)", text):
            if re.match(r"[a-zA-Z][a-zA-Z0-9+.-]*:", target) or target.startswith("#"):
                continue
            rel = target.split("#", 1)[0]
            if not (f.parent/rel).exists():
                raise FileNotFoundError(f"Broken local link in {f}: {target}")
            links += 1
    return {"python_files_parsed": parsed, "local_links_checked": links}


def main():
    if not __debug__ or os.environ.get("PYTHONOPTIMIZE"):
        raise RuntimeError("Assertions must be enabled; do not use python -O or PYTHONOPTIMIZE.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT/"artifacts")
    parser.add_argument("--record", action="store_true")
    args = parser.parse_args()
    out = args.output_dir.resolve()
    if out == ROOT or out == ROOT/"verification":
        raise ValueError("Use a separate raw-report output directory.")
    out.mkdir(parents=True, exist_ok=True)
    structure = structural_checks()
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    results = {}
    for name, groups in SUITES:
        report = out/f"{name}.json"
        subprocess.run([sys.executable, "-B", str(ROOT/f"tests/check_{name}.py"), "--output", str(report)],
                       cwd=ROOT, env=env, check=True, timeout=180)
        data = json.loads(report.read_text())
        if data.get("status") != "PASS":
            raise RuntimeError(f"Suite {name} did not report PASS")
        results[name] = {"status": "PASS", "groups": groups, "versions": data["versions"]}
    record = {
        "status": "PASS",
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "scope": f"{sum(n for _, n in SUITES)} scientific check groups plus structural checks supporting the analytical results. See verification/README.md for the evidence covered.",
        "scientific_groups": sum(n for _, n in SUITES),
        "structure": structure, "suites": results, "source_sha256": source_hashes(),
    }
    (out/"summary.json").write_text(json.dumps(record, indent=2)+"\n")
    if args.record:
        (ROOT/"verification/validated.json").write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps({"status": "PASS", "scientific_groups": record["scientific_groups"], "structure": structure, "reports": str(out)}, indent=2))


if __name__ == "__main__":
    main()

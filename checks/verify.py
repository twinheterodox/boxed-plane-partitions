"""Run all six finite checks and compare their results with the bundled records."""

import json
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory


CHECKS = Path(__file__).resolve().parent


def run(script):
    return subprocess.run(
        [sys.executable, "-E", "-B", str(script)],
        cwd=script.parent,
        stdout=subprocess.PIPE,
        text=True,
        check=True,
    ).stdout


def compare(script_name, actual, expected):
    if actual != expected:
        raise SystemExit(f"{script_name} produced a different result from the reference.")
    print(f"PASS  {script_name}", flush=True)


def main():
    # These two programs write JSON beside their source. Run a temporary copy
    # so that a failing check cannot overwrite its own reference result.
    with TemporaryDirectory() as temp:
        for program, record in [
            ("verify_smith.py", "smith-checks.json"),
            ("verify_interval.py", "interval-checks.json"),
        ]:
            script = Path(temp) / program
            shutil.copyfile(CHECKS / program, script)
            run(script)
            compare(program, (script.parent / record).read_bytes(), (CHECKS / record).read_bytes())

    arithmetic = CHECKS / "arithmetic"
    for program, record in [
        ("verify.py", "check-results.json"),
        ("check-addenda.py", "check-addenda-results.json"),
    ]:
        result = json.loads(run(arithmetic / program))
        compare("arithmetic/" + program, result, json.loads((arithmetic / record).read_text()))

    reference = json.loads((CHECKS / "reference-output.json").read_text())
    cofactor_lines = run(CHECKS / "verify_cofactors.py").splitlines()
    compare("verify_cofactors.py", cofactor_lines[-1], reference["cofactor_summary"])
    independent_lines = run(CHECKS / "independent-checks.py").splitlines()
    compare("independent-checks.py", independent_lines, reference["independent_audit"])


if __name__ == "__main__":
    main()

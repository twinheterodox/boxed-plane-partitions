"""Rebuild the manuscript from its TeX source and supplied vector figures."""

from pathlib import Path
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"
STEM = "boxed-plane-partitions"


def main():
    build = ROOT / "build" / "latex"
    build.mkdir(parents=True, exist_ok=True)
    command = [
        "pdflatex", "-no-shell-escape", "-interaction=nonstopmode", "-halt-on-error",
        "-file-line-error", f"-output-directory={build}", STEM + ".tex",
    ]
    previous = None
    for _ in range(5):
        result = subprocess.run(command, cwd=PAPER, capture_output=True, text=True)
        if result.returncode:
            raise SystemExit(result.stdout[-5000:] + result.stderr)
        state = tuple((build / (STEM + ext)).read_bytes() for ext in (".aux", ".out"))
        if state == previous:
            break
        previous = state
    else:
        raise SystemExit("Cross-references did not settle after five runs.")
    log = (build / (STEM + ".log")).read_text(errors="replace")
    # All supplied plates are PDFs. epstopdf reports that it cannot invoke an
    # external converter under -no-shell-escape even when no EPS file is used.
    flags = ("Warning:", "Overfull ", "Underfull ", "Rerun to get")
    issues = [line for line in log.splitlines()
              if any(flag in line for flag in flags)
              and not line.startswith("Package epstopdf Warning: Shell escape feature is not enabled.")]
    if issues:
        raise SystemExit("\n".join(issues))
    output = ROOT / "output" / "pdf" / (STEM + ".pdf")
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(build / output.name, output)
    print(output.relative_to(ROOT))


if __name__ == "__main__":
    main()

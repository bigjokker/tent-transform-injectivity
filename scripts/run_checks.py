"""Run every numerical check in code/ and report pass/fail with timings.

Usage (from the repository root):
    python scripts/run_checks.py          # all checks
    python scripts/run_checks.py --quick  # skip the slow kernel-lemma stress test

Each check script exits with a nonzero status on failure. The checks are numerical
illustrations of the analytic results in papers/ and docs/; they are not proofs.
"""
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECKS = [
    ("real_line_checks.py", "autoconvolution identity; escaping-extremum example"),
    ("periodic_linearization_checks.py", "periodic linearization multipliers and remainder"),
    ("controlled_tail_checks.py", "controlled-tail approximation g ~ pi/(4 f^2)"),
    ("periodic_reconstruction.py", "periodic reconstruction with certified mean search"),
    ("kernel_lemma_checks.py", "kernel lemmas of the uniqueness proof (slow)"),
]


def main() -> int:
    quick = "--quick" in sys.argv[1:]
    failures = 0
    for script, what in CHECKS:
        if quick and script == "kernel_lemma_checks.py":
            print(f"SKIP  {script:36s} {what}")
            continue
        start = time.time()
        proc = subprocess.run([sys.executable, str(ROOT / "code" / script)],
                              cwd=ROOT, capture_output=True, text=True)
        status = "PASS" if proc.returncode == 0 else "FAIL"
        print(f"{status}  {script:36s} {what}  ({time.time() - start:.0f} s)")
        if proc.returncode != 0:
            failures += 1
            print(proc.stdout[-3000:])
            print(proc.stderr[-3000:])
    print("All checks passed." if failures == 0 else f"{failures} check(s) failed.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

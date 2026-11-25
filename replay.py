"""Replay the campaign's exact certificates. Run with Python 3.12 and SymPy 1.14."""

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent


def read_json(path):
    return json.loads(path.read_bytes())


def run(script, *args):
    result = subprocess.run(
        [sys.executable, "-E", "-B", str(HERE / script), *map(str, args)],
        cwd=HERE, stdout=subprocess.PIPE, text=True, check=True,
    )
    return result.stdout


def require_same(actual, expected, label):
    if actual != expected:
        raise ValueError(f"Replayed output differs from the archived {label}")


def verify_hashes(manifest, folder):
    for record in manifest.read_text(encoding="utf-8").splitlines():
        expected, relative = record.split(maxsplit=1)
        actual = sha256((folder / relative).read_bytes()).hexdigest()
        require_same(actual, expected, relative)


def main():
    verify_hashes(HERE / "verification/FROZEN.sha256", HERE)
    witness_record = read_json(HERE / "pivot/witness.json")
    require_same(json.loads(run("pivot/verify.py")), read_json(HERE / "pivot/baseline.json"), "pivot baseline")
    witness = witness_record["witness"]
    parameters = [witness["x"], witness["z"], *witness["y_isolating_interval"]]
    regenerated = json.loads(run("pivot/verify.py", "--witness", *parameters))
    require_same(regenerated, witness_record, "pivot witness")
    require_same(json.loads(run("pivot/region.py")), read_json(HERE / "pivot/region.json"), "pivot box")
    comparator = json.loads(run("pivot/comparator.py"))
    require_same(comparator, read_json(HERE / "pivot/comparator.json"), "pivot comparator")
    print(run("verification/check_lean_inputs.py").strip(), flush=True)
    print(run("verification/pivot_subbox.py").strip(), flush=True)
    print(run("verification/pivot_attainment.py").strip(), flush=True)
    print("Pivot benchmark, exact witness, box exclusion and sharpness checks reproduced.")


if __name__ == "__main__":
    main()

"""Bind the Lean integer inputs to the benchmark and derived eliminant."""

from ast import literal_eval
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read_definition(source, name):
    literal = source.split(f"def {name} :", 1)[1].split(":=", 1)[1].split("\n\n", 1)[0]
    return literal_eval(literal.strip())


def main():
    source = (ROOT / "formalization/BenchmarkSigns.lean").read_text(encoding="utf-8")
    coefficients = (ROOT / "pivot/p5_coefficients.txt").read_text().splitlines()
    if read_definition(source, "benchmarkCoeffs") != list(map(int, coefficients)):
        raise ValueError("Lean benchmark coefficients differ from the published polynomial")
    archive = json.loads((ROOT / "pivot/comparator.json").read_bytes())
    if read_definition(source, "eliminantCoeffs") != archive["coefficients_ascending_z_g"]:
        raise ValueError("Lean eliminant differs from the polynomial reconstructed by replay.py")
    print("Lean benchmark and eliminant coefficients match the exact inputs.")


if __name__ == "__main__":
    main()

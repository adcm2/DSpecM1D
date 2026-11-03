#!/usr/bin/env python3
"""Generate a hash-bound, standalone v4 analytic-pole control receiver gate."""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
DSPEC = HERE.parent
ROOT = HERE.parents[3]
OLD_ROOT = DSPEC / "source-v4-tunya"
CONTROL_ROOT = DSPEC / "source-v4-analytic-pole-control"
V6_ROOT = DSPEC / "source-v6-frozen"
OLD_HEADER = OLD_ROOT / "engine_headers/Windowed3DFields.h"
CONTROL_HEADER = CONTROL_ROOT / "engine_headers/Windowed3DFields.h"
V6_HEADER = V6_ROOT / "engine_headers/Windowed3DFields.h"
OLD_MANIFEST = DSPEC / "source-v4-tunya.sha256"
CONTROL_MANIFEST = DSPEC / "source-v4-analytic-pole-control.sha256"
V6_MANIFEST = DSPEC / "source-v6-frozen.sha256"
TEMPLATE = HERE / "control_receiver_gate_template.cpp.in"
OUTPUT = HERE / "control_receiver_gate.cpp"


def expected_hash(manifest: Path, name: str) -> str:
    for line in manifest.read_text().splitlines():
        digest, separator, entry = line.partition("  ")
        if separator and entry == name:
            return digest
    raise RuntimeError(f"{name} absent from {manifest}")


def verify_hash(path: Path, manifest: Path, entry: str) -> str:
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    expected = expected_hash(manifest, entry)
    if actual != expected:
        raise RuntimeError(f"SHA-256 mismatch for {path}: {actual} != {expected}")
    return actual


def extract_function(source: str, signature: str) -> tuple[str, int]:
    start = source.find(signature)
    if start < 0:
        raise RuntimeError(f"signature not found: {signature}")
    brace = source.find("{", start)
    if brace < 0:
        raise RuntimeError(f"opening brace not found: {signature}")
    depth = 0
    for index in range(brace, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[start:index + 1], source.count("\n", 0, start) + 1
    raise RuntimeError(f"unclosed function: {signature}")


def check_include_lines(source: str) -> None:
    for number, line in enumerate(source.splitlines(), 1):
        if line.startswith("#include") and not re.fullmatch(
            r'#include\s+(?:<[^>\n]+>|"[^"\n]+")', line
        ):
            raise RuntimeError(f"malformed #include at generated line {number}: {line}")


def main() -> None:
    old_sha = verify_hash(OLD_HEADER, OLD_MANIFEST,
                          "source-v4-tunya/engine_headers/Windowed3DFields.h")
    control_sha = verify_hash(CONTROL_HEADER, CONTROL_MANIFEST,
                              "source-v4-analytic-pole-control/engine_headers/Windowed3DFields.h")
    v6_sha = verify_hash(V6_HEADER, V6_MANIFEST,
                         "source-v6-frozen/engine_headers/Windowed3DFields.h")
    old_text, control_text, v6_text = OLD_HEADER.read_text(), CONTROL_HEADER.read_text(), V6_HEADER.read_text()
    harmonic_sig = "static Complex harmonic(int degree, int m, int n, double theta, double phi)"
    receiver_sig = "static Eigen::MatrixXcd receiverAtPole(SensitivitySolveFamily family,int degree)"
    old_harmonic, old_harmonic_line = extract_function(old_text, harmonic_sig)
    control_harmonic, control_harmonic_line = extract_function(control_text, harmonic_sig)
    old_receiver, old_receiver_line = extract_function(old_text, receiver_sig)
    control_receiver, control_receiver_line = extract_function(control_text, receiver_sig)
    v6_receiver, v6_receiver_line = extract_function(v6_text, receiver_sig)

    old_open = old_harmonic.find("{", old_harmonic.find(harmonic_sig))
    branch = (
        "if (degree >= 1 && m >= -degree && m <= degree && n >= -degree && n <= degree &&\n"
        "        theta == 0.0 && phi == 0.0) {\n"
        "      const double norm = std::sqrt((2.0 * degree + 1.0) / (4.0 * EIGEN_PI));\n"
        "      return m == n ? Complex(norm, 0.0) : Complex{};\n"
        "    }"
    )
    original_fallback = old_harmonic[old_open + 1:]
    expected_control = old_harmonic[:old_open + 1] + "\n    " + branch + original_fallback
    if control_harmonic != expected_control:
        raise RuntimeError("control scalar helper is not exactly one pole-only early branch")
    if control_receiver != old_receiver:
        raise RuntimeError("control receiver aggregation differs from original v4")
    radial_call = "const Complex y00=harmonic(0,0,0,0.0,0.0);"
    if old_text.count(radial_call) != 1 or control_text.count(radial_call) != 1:
        raise RuntimeError("radial L=0 scalar-helper callsite changed or is ambiguous")
    if old_receiver.count("harmonic(") != 3:
        raise RuntimeError("unexpected nonradial scalar-harmonic receiver call graph")
    if control_receiver.count("harmonic(") != 3:
        raise RuntimeError("control receiver scalar-harmonic call graph changed")
    receiver_calls = {
        "result(0,col)=harmonic(degree,m,0,0.0,0.0);",
        "const Complex ym=harmonic(degree,m,-1,0.0,0.0);",
        "const Complex yp=harmonic(degree,m,1,0.0,0.0);",
    }
    scalar_calls = [line.strip() for line in old_text.splitlines()
                    if re.search(r"\bharmonic\(\s*(?:degree|0)\s*,", line)]
    if len(scalar_calls) != 4 or sum(line in receiver_calls for line in scalar_calls) != 3:
        raise RuntimeError(f"unexpected scalar harmonic call graph: {scalar_calls}")
    if scalar_calls.count(radial_call) != 1:
        raise RuntimeError("scalar L=0 helper is not uniquely isolated from nonradial receiver")
    control_scalar_calls = [line.strip() for line in control_text.splitlines()
                            if re.search(r"\bharmonic\(\s*(?:degree|0)\s*,", line)]
    if control_scalar_calls != scalar_calls:
        raise RuntimeError("control changed the scalar harmonic call graph")

    replacements = {
        "@OLD_HEADER@": str(OLD_HEADER.relative_to(ROOT)),
        "@OLD_HARMONIC_LINE@": str(old_harmonic_line),
        "@OLD_SHA@": old_sha,
        "@CONTROL_HEADER@": str(CONTROL_HEADER.relative_to(ROOT)),
        "@CONTROL_SHA@": control_sha,
        "@CONTROL_HARMONIC_LINE@": str(control_harmonic_line),
        "@CONTROL_RECEIVER_LINE@": str(control_receiver_line),
        "@V6_HEADER@": str(V6_HEADER.relative_to(ROOT)),
        "@V6_SHA@": v6_sha,
        "@V6_RECEIVER_LINE@": str(v6_receiver_line),
        "@OLD_HARMONIC@": old_harmonic,
        "@CONTROL_HARMONIC@": control_harmonic,
        "@CONTROL_RECEIVER@": control_receiver,
        "@V6_RECEIVER@": v6_receiver,
    }
    generated = TEMPLATE.read_text()
    for token, value in replacements.items():
        if token not in generated:
            raise RuntimeError(f"template token missing: {token}")
        generated = generated.replace(token, value)
    if re.search(r"@[A-Z0-9_]+@", generated):
        raise RuntimeError("unexpanded template token remains")
    check_include_lines(generated)
    OUTPUT.write_text(generated)
    print(f"wrote {OUTPUT.relative_to(DSPEC)}; complete includes and source proofs passed")


if __name__ == "__main__":
    main()

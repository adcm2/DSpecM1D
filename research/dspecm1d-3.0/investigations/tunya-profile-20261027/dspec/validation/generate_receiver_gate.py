#!/usr/bin/env python3
"""Extract frozen v4/v6 receiver bodies into a standalone, reviewable C++ gate."""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
DSPEC = HERE.parent
ROOT = HERE.parents[3]
V4_ROOT = ROOT / "investigations/dspec-sk-profile-20261026/physical-profile/dspec/source-v4-frozen"
V4_HEADER = V4_ROOT / "engine_headers/Windowed3DFields.h"
V4_MANIFEST = V4_ROOT.parent / "source-v4-frozen.sha256"
V6_ROOT = DSPEC / "source-v6-frozen"
V6_HEADER = V6_ROOT / "engine_headers/Windowed3DFields.h"
V6_MANIFEST = DSPEC / "source-v6-frozen.sha256"
TEMPLATE = HERE / "receiver_gate_template.cpp.in"
OUTPUT = HERE / "receiver_gate.cpp"


def expected_hash(manifest: Path, relative_entry: str) -> str:
    for line in manifest.read_text().splitlines():
        digest, separator, name = line.partition("  ")
        if separator and name == relative_entry:
            return digest
    raise RuntimeError(f"{relative_entry} absent from {manifest}")


def verify_hash(path: Path, manifest: Path, relative_entry: str) -> str:
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    expected = expected_hash(manifest, relative_entry)
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
                end = index + 1
                line = source.count("\n", 0, start) + 1
                return source[start:end], line
    raise RuntimeError(f"unclosed function body: {signature}")


def check_include_directives(source: str) -> None:
    for number, line in enumerate(source.splitlines(), 1):
        if line.startswith("#include") and not re.fullmatch(
            r'#include\s+(?:<[^>\n]+>|"[^"\n]+")', line
        ):
            raise RuntimeError(f"malformed preprocessor include on generated line {number}: {line}")


def main() -> None:
    v4_sha = verify_hash(V4_HEADER, V4_MANIFEST, "source-v4-frozen/engine_headers/Windowed3DFields.h")
    v6_sha = verify_hash(V6_HEADER, V6_MANIFEST, "source-v6-frozen/engine_headers/Windowed3DFields.h")
    v4_text = V4_HEADER.read_text()
    v6_text = V6_HEADER.read_text()
    v4_harmonic, v4_harmonic_line = extract_function(
        v4_text, "static Complex harmonic(int degree, int m, int n, double theta, double phi)"
    )
    v4_receiver, v4_receiver_line = extract_function(
        v4_text, "static Eigen::MatrixXcd receiverAtPole(SensitivitySolveFamily family,int degree)"
    )
    v6_receiver, v6_receiver_line = extract_function(
        v6_text, "static Eigen::MatrixXcd receiverAtPole(SensitivitySolveFamily family,int degree)"
    )
    v6_estimator, v6_estimator_line = extract_function(
        v6_text, "static std::size_t angularWignerBytesPerAngle(int maxDegree,int mMax)"
    )
    replacements = {
        "@V4_HEADER@": str(V4_HEADER.relative_to(ROOT)),
        "@V4_RECEIVER_LINE@": str(v4_receiver_line),
        "@V4_HARMONIC_LINE@": str(v4_harmonic_line),
        "@V4_SHA256@": v4_sha,
        "@V6_HEADER@": str(V6_HEADER.relative_to(ROOT)),
        "@V6_RECEIVER_LINE@": str(v6_receiver_line),
        "@V6_ESTIMATOR_LINE@": str(v6_estimator_line),
        "@V6_SHA256@": v6_sha,
        "@V4_HARMONIC@": v4_harmonic,
        "@V4_RECEIVER@": v4_receiver,
        "@V6_RECEIVER@": v6_receiver,
        "@V6_ESTIMATOR@": v6_estimator,
    }
    generated = TEMPLATE.read_text()
    for token, value in replacements.items():
        if token not in generated:
            raise RuntimeError(f"template token missing: {token}")
        generated = generated.replace(token, value)
    if re.search(r"@[A-Z0-9_]+@", generated):
        raise RuntimeError("unexpanded token remains in generated receiver gate")
    check_include_directives(generated)
    OUTPUT.write_text(generated)
    print(f"wrote {OUTPUT.relative_to(DSPEC)}; all generated #include directives are complete")


if __name__ == "__main__":
    main()

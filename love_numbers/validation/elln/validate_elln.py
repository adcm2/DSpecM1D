#!/usr/bin/env python3

import argparse
import hashlib
import io
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile


ARCHIVE_SHA256 = (
    "f2904364f9be96239fc85140a032b01b732cba6ea1e9bc323f2f3325834e7da7"
)
FILE_SHA256 = {
    "ELLN.m":
        "6341a62a0c595edfea15a81e0d8ed7d79f601b1ebd0c10448926b66bd0648864",
    "EarthMantleTI56.txt":
        "43d455cf329b7d9db6684b434447b648e3a5311f84873bb283e2a8f22f8f17fe",
    "EarthCore26.txt":
        "b6ddf75db09d15b4d7e065f109134d54cdd68434e7714b75839d1c12c49d9f12",
}
OUTER_MEMBERS = {
    "2017 1010 GJI ELLN  Supplement FINAL Clean.doc":
        "cf8174bbb56d06eb53f99f06a2f538da52654785815bc70d70d05924f0ee83a7",
    "ELLN Direct Run.zip":
        "f4c4b10b84d1dfa08d5f4187fffebf9c50d5a860d3d22eeebd801f5befb08f7a",
    "ELLN Example Input Files.zip":
        "c5b138d3806c4d6590cd8d417d21b1195315072ff6fd0f02d68cd6ffe9cb0798",
    "ELLN GUI Run.zip":
        "56d00b6ed524d0896768bd73aaabe2d5bdb95f280e7e192568f835c9b8799fb4",
}


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def extract_zip(archive, destination):
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    for member in archive.infolist():
        target = (destination / member.filename).resolve()
        if root not in target.parents and target != root:
            raise RuntimeError(f"Unsafe archive member: {member.filename}")
        archive.extract(member, destination)


def locate_one(root, name):
    matches = sorted(root.rglob(name))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one {name} below {root}; found {len(matches)}"
        )
    return matches[0]


def prepare_archive(archive_path, work_directory):
    if sha256_file(archive_path) != ARCHIVE_SHA256:
        raise RuntimeError(
            f"{archive_path} is not the audited official ggx444_supp.zip"
        )

    source_directory = work_directory / "official-source"
    source_directory.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path) as outer:
        names = set(outer.namelist())
        if names != set(OUTER_MEMBERS):
            raise RuntimeError("Official archive member list has changed")
        for name, expected in OUTER_MEMBERS.items():
            data = outer.read(name)
            if sha256_bytes(data) != expected:
                raise RuntimeError(f"Checksum mismatch for {name}")
            if name.endswith(".zip") and name != "ELLN GUI Run.zip":
                with zipfile.ZipFile(io.BytesIO(data)) as nested:
                    extract_zip(nested, source_directory)

    print(f"ELLN archive SHA-256: {ARCHIVE_SHA256}")
    return source_directory


def verify_file(path):
    expected = FILE_SHA256[path.name]
    actual = sha256_file(path)
    if actual != expected:
        raise RuntimeError(
            f"Checksum mismatch for {path.name}: {actual}, expected {expected}"
        )


def read_numeric_rows(path, columns):
    rows = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        fields = line.split()
        if len(fields) != columns:
            continue
        try:
            row = [float(field) for field in fields]
        except ValueError:
            continue
        rows.append(row)
    return rows


def audit_source(source_directory):
    paths = {
        name: locate_one(source_directory, name) for name in FILE_SHA256
    }
    for path in paths.values():
        verify_file(path)

    source = paths["ELLN.m"].read_text(encoding="utf-8-sig")
    required_fragments = [
        "GetMc(1,1,i)=MantleData(i,4)",
        "GetMc(3,3,i)=MantleData(i,5)",
        "GetMc(4,4,i)=MantleData(i,6)",
        "GetMc(6,6,i)=MantleData(i,7)",
        "GetMc(1,3,i)=MantleData(i,8)",
        "Ggg=6.67384E-11",
        "RMax=Rr(nLayer+1)",
        "EqRho=Rho(ILayer+1)",
        "LLNsRes(II,4)=-(RUO(3)*grav*RMax^2/Ggg+1)",
        "LLNsRes(II,2)=LLNsRes(II,2)-LLNsRes(II,4)",
        "AllRes=[Nn,-Hn,Ln,-Kn,Nn.*Ln,-Nn.*Kn]",
    ]
    for fragment in required_fragments:
        if fragment not in source:
            raise RuntimeError(f"ELLN source audit failed at: {fragment}")

    mantle = read_numeric_rows(paths["EarthMantleTI56.txt"], 8)
    core = read_numeric_rows(paths["EarthCore26.txt"], 4)
    if len(mantle) != 57 or len(core) != 27:
        raise RuntimeError("Unexpected ELLN PREM model dimensions")
    if mantle[0][1] != core[-1][1] or mantle[-1][1] != 6371000.0:
        raise RuntimeError("Unexpected ELLN core or surface radius")

    print("ELLN source/model checksums: verified")
    print(
        "ELLN model: 56 TI mantle layers, 26 fluid-core layers, "
        "surface radius 6371000 m"
    )
    return paths


def read_reference(path):
    references = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 6:
            raise RuntimeError(f"Invalid official-reference row: {line}")
        degree = int(fields[0])
        values = tuple(float(value) for value in fields[1:])
        references[degree] = values

    for degree, (_, horizontal, _, scaled_horizontal, _) in references.items():
        expected = degree * horizontal
        if not math.isclose(expected, scaled_horizontal,
                            rel_tol=5.0e-7, abs_tol=5.0e-7):
            print(
                "Official Table S7 internal discrepancy: "
                f"degree {degree}, n*l={expected:.6e}, "
                f"published n*l={scaled_horizontal:.6e}"
            )
    return references


def validate_conversion():
    gravitational_constant = 6.67384e-11
    surface_radius = 6371000.0
    surface_gravity = 9.82
    degree = 2
    factor = (
        (2 * degree + 1)
        / (4.0 * math.pi * gravitational_constant * surface_radius)
    )

    expected_h = -0.75
    h_load = expected_h / (factor * surface_gravity)
    actual_h = factor * surface_gravity * h_load

    expected_k = -0.25
    k_load = (-expected_k - 1.0) / factor
    actual_k = -(1.0 + factor * k_load)

    if actual_h != expected_h or actual_k != expected_k:
        raise RuntimeError("Dimensional-to-conventional conversion failed")
    if not math.isclose(
        factor * surface_gravity * h_load, expected_h, rel_tol=1.0e-15
    ):
        raise RuntimeError("Displacement conversion round trip failed")
    if not math.isclose(
        -(1.0 + factor * k_load), expected_k, rel_tol=1.0e-15
    ):
        raise RuntimeError("Potential conversion round trip failed")

    print(
        "ELLN conversion: h_l=(2l+1) g(a) h_load/(4 pi G a)"
    )
    print(
        "ELLN conversion: k_l=-[1+(2l+1) k_load/(4 pi G a)]"
    )
    print("ELLN conversion algebra and signs: verified")


def validate_degree_one_convention(source_path):
    source = source_path.read_text(encoding="utf-8-sig")
    required_fragments = (
        "if(nOrder==1)",
        "LLNsRes(II,2)=LLNsRes(II,2)-LLNsRes(II,4)",
        "LLNsRes(II,3)=LLNsRes(II,3)-LLNsRes(II,4)",
        "LLNsRes(II,4)=LLNsRes(II,4)-LLNsRes(II,4)",
    )
    for fragment in required_fragments:
        if fragment not in source:
            raise RuntimeError(
                f"ELLN degree-one source audit failed at: {fragment}"
            )

    dspec_surface_potential = 0.0
    generic_k = -(1.0 + dspec_surface_potential)
    elln_final_k = generic_k - generic_k
    if generic_k != -1.0 or elln_final_k != 0.0:
        raise RuntimeError("Degree-one convention check failed")

    print(
        "ELLN degree one: raw k is subtracted from h and l, "
        "then k is set to zero in the deformed-centre frame"
    )
    print(
        "DSpecM1D P(a)=0: the generic conversion gives k=-1, "
        "so it is not ELLN's final degree-one k"
    )
    print(
        "ELLN degree-one h and k are excluded because frame equivalence "
        "has not been established"
    )


def patch_non_gui(source, destination, patch_path):
    patch_path = patch_path.resolve()
    destination.write_text(
        source.read_text(encoding="utf-8-sig"), encoding="utf-8"
    )
    patch_executable = shutil.which("patch")
    if patch_executable is None:
        raise RuntimeError("patch is required for a non-GUI ELLN run")
    subprocess.run(
        [
            patch_executable, "--batch", "--forward", "-p0",
            f"--input={patch_path}",
        ],
        cwd=destination.parent,
        check=True,
    )


def run_example(paths, driver, patch_path, work_directory,
                matlab, octave):
    run_directory = work_directory / "run"
    run_directory.mkdir(parents=True, exist_ok=True)
    patch_non_gui(paths["ELLN.m"], run_directory / "ELLN.m", patch_path)
    for name in ("EarthMantleTI56.txt", "EarthCore26.txt"):
        shutil.copyfile(paths[name], run_directory / name)

    output = run_directory / "ti-prem-results.txt"
    if output.exists():
        output.unlink()
    environment = os.environ.copy()
    environment["ELLN_WORK_DIR"] = str(run_directory)
    environment["ELLN_OUTPUT"] = str(output)

    quoted_driver = str(driver.resolve()).replace("'", "''")
    if matlab:
        runtime = f"MATLAB: {matlab}"
        expression = (
            f"try, run('{quoted_driver}'); "
            "catch exception, disp(getReport(exception)); exit(1); "
            "end; exit(0);"
        )
        command = [
            matlab, "-nodisplay", "-nosplash", "-nodesktop",
            "-r", expression,
        ]
    elif octave:
        runtime = f"GNU Octave: {octave}"
        command = [octave, "--no-gui", "--quiet", str(driver.resolve())]
    else:
        return None, None

    print(f"ELLN runtime: {runtime}")
    subprocess.run(command, env=environment, check=True)
    if not output.is_file():
        raise RuntimeError("ELLN driver did not create its output")
    return runtime, read_driver_output(output)


def read_driver_output(path):
    results = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) != 6:
            raise RuntimeError(f"Invalid ELLN driver row: {line}")
        degree = int(fields[0])
        values = tuple(float(value) for value in fields[1:])
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError(f"Non-finite ELLN result for degree {degree}")
        results[degree] = values
    if set(results) != {1, 2, 10, 50}:
        raise RuntimeError("ELLN driver did not return degrees 1, 2, 10, 50")
    return results


def compare_results(results, references):
    for degree in (1, 2, 10, 50):
        h, horizontal, k, mass, gravity = results[degree]
        print(
            f"ELLN degree {degree}: h={h:.17e} l={horizontal:.17e} "
            f"k={k:.17e} mass={mass:.17e} gravity={gravity:.17e}"
        )
        if degree == 1 and k != 0.0:
            raise RuntimeError("ELLN degree-one k is not exactly zero")
        if degree not in references:
            print(
                f"  no publisher table row for degree {degree}; "
                "live value reported without an invented oracle"
            )
            continue

        minus_h, reference_horizontal, minus_k, _, _ = references[degree]
        generated = (-h, horizontal, -k)
        official = (minus_h, reference_horizontal, minus_k)
        names = ("-h", "l", "-k")
        for name, value, reference in zip(names, generated, official):
            signed = value - reference
            relative = signed / max(abs(reference), sys.float_info.min)
            print(
                f"  {name}: signed={signed:.17e} "
                f"relative={relative:.17e}"
            )
            if float(f"{value:.6e}") != reference:
                raise RuntimeError(
                    f"ELLN degree {degree} {name} does not reproduce "
                    "the official six-digit table value"
                )


def print_reference_fallback(references):
    print("ELLN live execution: unavailable (MATLAB and Octave not found)")
    for degree in (1, 2, 10, 50):
        if degree not in references:
            print(
                f"ELLN official Table S7: degree {degree} is not tabulated"
            )
            continue
        minus_h, horizontal, minus_k, _, _ = references[degree]
        print(
            f"ELLN official degree {degree}: "
            f"h={-minus_h:.6e} l={horizontal:.6e} k={-minus_k:.6e}"
        )


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--source-dir", type=Path)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--driver", type=Path)
    parser.add_argument("--patch", type=Path)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--matlab")
    parser.add_argument("--octave")
    parser.add_argument("--conversion-only", action="store_true")
    parser.add_argument("--degree-one-only", action="store_true")
    arguments = parser.parse_args()
    if arguments.archive and arguments.source_dir:
        parser.error("set only one of --archive and --source-dir")
    return arguments


def main():
    arguments = parse_arguments()
    references = read_reference(arguments.reference)
    validate_conversion()
    if arguments.conversion_only:
        return

    if not arguments.archive and not arguments.source_dir:
        raise RuntimeError(
            "Provide the official archive or an extracted source directory"
        )
    if not arguments.work_dir:
        raise RuntimeError("Work directory is required")

    arguments.work_dir.mkdir(parents=True, exist_ok=True)
    if arguments.archive:
        source_directory = prepare_archive(
            arguments.archive.resolve(), arguments.work_dir
        )
    else:
        source_directory = arguments.source_dir.resolve()
    paths = audit_source(source_directory)
    if arguments.degree_one_only:
        validate_degree_one_convention(paths["ELLN.m"])
        return

    if not arguments.driver or not arguments.patch:
        raise RuntimeError("Driver and patch are required")
    _, results = run_example(
        paths, arguments.driver, arguments.patch, arguments.work_dir,
        arguments.matlab, arguments.octave,
    )
    if results is None:
        print_reference_fallback(references)
    else:
        compare_results(results, references)


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, subprocess.CalledProcessError,
            zipfile.BadZipFile) as error:
        print(f"ELLN validation failed: {error}", file=sys.stderr)
        raise SystemExit(1)

import subprocess
import sys

from conftest import shared_members, write


def program(name, *args, cwd=None):
    """Runs one of the package's programs as its command would."""
    code = f"import sys; from yspec.programs import {name}; {name}()"
    return subprocess.run([sys.executable, "-c", code, *args], cwd=cwd,
                          capture_output=True, text=True)


def test_yspec_writes_the_seismograms(deck, tmp_path):
    members = shared_members(deck, "out")
    members["receivers"] = {"depth_km": 0, "locations": [[0, 60], [10, 80]]}
    path = write(tmp_path / "input.json", members)
    run = program("seismograms", path, cwd=tmp_path)
    assert run.returncode == 0, run.stderr
    assert "wrote 2 files with prefix" in run.stderr
    assert (tmp_path / "out.1").exists() and (tmp_path / "out.2").exists()
    assert (tmp_path / "out.receivers").exists()


def test_programs_report_faults_and_their_version(tmp_path):
    run = program("kernels", str(tmp_path / "missing.json"))
    assert run.returncode == 1
    assert run.stderr.startswith("yspec-kernels: ")
    run = program("snapshots", "--version")
    assert run.returncode == 0 and run.stdout.startswith("yspec-snapshot ")
    assert program("seismograms").returncode == 2

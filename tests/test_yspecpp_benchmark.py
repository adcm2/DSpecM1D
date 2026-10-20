import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "benchmarks" / "yspecpp" / "run_benchmark.py"


class OutputDirectoryGuardTests(unittest.TestCase):
    def test_populated_run_directory_stops_before_subprocesses_or_writes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            hook = root / "hooks"
            hook.mkdir()
            (hook / "sitecustomize.py").write_text(
                "import os\n"
                "from pathlib import Path\n"
                "import subprocess\n"
                "def blocked(*args, **kwargs):\n"
                "    Path(os.environ['PROCESS_SENTINEL']).write_text('called')\n"
                "    raise RuntimeError('subprocess sentinel called')\n"
                "subprocess.run = blocked\n"
                "subprocess.check_output = blocked\n"
                "subprocess.Popen = blocked\n",
                encoding="utf-8",
            )
            output = root / "saved-run"
            (output / "nested").mkdir(parents=True)
            (output / ".provenance").write_text("keep", encoding="utf-8")
            (output / "nested" / "input.json").write_text(
                '{"frozen": true}\n', encoding="utf-8"
            )
            before = {
                path.relative_to(output).as_posix():
                ("directory", None) if path.is_dir() else
                ("file", path.read_bytes())
                for path in sorted(output.rglob("*"))
            }
            sentinel = root / "subprocess-called"
            env = os.environ.copy()
            env.update({
                "PYTHONDONTWRITEBYTECODE": "1",
                "PYTHONPATH": str(hook),
                "PROCESS_SENTINEL": str(sentinel),
            })
            result = subprocess.run(
                [sys.executable, "-B", str(RUNNER),
                 "--dspec-source", "missing-dspec",
                 "--yspec-source", "missing-yspec",
                 "--dependency-root", "missing-dependencies",
                 "--output", str(output)],
                cwd=root, env=env, capture_output=True, text=True, check=False,
            )
            after = {
                path.relative_to(output).as_posix():
                ("directory", None) if path.is_dir() else
                ("file", path.read_bytes())
                for path in sorted(output.rglob("*"))
            }
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn("output directory must be empty", result.stderr)
            self.assertFalse(sentinel.exists())
            self.assertEqual(after, before)


if __name__ == "__main__":
    unittest.main()

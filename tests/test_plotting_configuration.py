import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "plotting" / "functions"))

from helpers.mode_annotations import get_modes_in_range
from helpers.plot_paths import build_output_path, default_build_dir


class PlotPathTests(unittest.TestCase):
    def test_standard_and_nested_build_trees_resolve_to_output_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            ordinary = Path(tmp) / "build" / "release"
            nested = Path(tmp) / "build" / "ci" / "release"
            with patch.dict(os.environ, {"DSPECM1D_BUILD_DIR": str(ordinary)}):
                self.assertEqual(default_build_dir(), ordinary.resolve())
            self.assertEqual(
                build_output_path("ex1_w.out", nested),
                str((nested / "output" / "ex1_w.out").resolve()),
            )


class ModeCatalogueTests(unittest.TestCase):
    def test_explicit_catalogues_are_sorted_and_keep_mode_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "modes"
            path.write_text(
                " mode phs vel w(mhz)\n"
                " n type l x w\n"
                "MINEOS mode catalogue\n"
                "01 S 02 x 0.40\n"
                "2 R 3 x 0.30\n"
                "End of catalogue: two mode rows were listed above.\n",
                encoding="utf-8",
            )

            names, frequencies = get_modes_in_range(0.2, 0.5, [path])

        self.assertEqual(frequencies, [0.3, 0.4])
        self.assertEqual(names, ["${_2}$R${_{3}}$", "${_01}$S${_{02}}$"])

    def test_absent_catalogues_skip_optional_annotations(self):
        self.assertEqual(get_modes_in_range(0.2, 0.5), ([], []))

    def test_missing_and_malformed_explicit_catalogues_fail_clearly(self):
        with tempfile.TemporaryDirectory() as tmp:
            missing = Path(tmp) / "missing.modes"
            malformed = Path(tmp) / "malformed.modes"
            malformed_rows = {
                "one-token-truncated": "1",
                "signed-negative-incomplete": "-1 S 2",
                "negative-radial-order": "-1 S 2 x 0.40",
                "negative-angular-order": "1 S -2 x 0.40",
                "invalid-mode-class": "1 SS 2 x 0.40",
            }
            malformed.write_text("not a MINEOS mode table\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "Cannot read mode catalogue"):
                get_modes_in_range(0.2, 0.5, [missing])
            with self.assertRaisesRegex(ValueError, "no MINEOS table header"):
                get_modes_in_range(0.2, 0.5, [malformed])
            for name, row in malformed_rows.items():
                path = Path(tmp) / f"{name}.modes"
                path.write_text(
                    f" mode phs vel w(mhz)\n n type l x w\n{row}\n",
                    encoding="utf-8",
                )
                with self.subTest(row=name), self.assertRaisesRegex(
                    ValueError, f"{path.name}' at line 3"
                ):
                    get_modes_in_range(0.2, 0.5, [path])


class Ex1HeadlessTests(unittest.TestCase):
    def test_default_data_path_and_headless_save_work_from_another_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            temp = Path(tmp)
            build = temp / "nested" / "release"
            output = build / "output"
            output.mkdir(parents=True)
            frequencies = np.linspace(0.0, 5.0, 501)
            data = np.zeros((frequencies.size, 22))
            data[:, 0] = frequencies
            data[:, 3] = 0.5 + 0.4 * np.sin(frequencies)
            data[:, 12] = 0.6 + 0.3 * np.cos(frequencies)
            data[:, 21] = 0.55 + 0.2 * np.sin(2.0 * frequencies)
            np.savetxt(output / "ex1_w.out", data, delimiter=";")
            save_path = temp / "plots" / "ex1.png"
            env = os.environ.copy()
            env["DSPECM1D_BUILD_DIR"] = str(build)
            env["PYTHONDONTWRITEBYTECODE"] = "1"

            result = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(REPO / "plotting" / "functions" / "ex1.py"),
                    "--no-show",
                    "--save",
                    str(save_path),
                ],
                cwd=temp,
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn(str((output / "ex1_w.out").resolve()), result.stdout)
            self.assertTrue(save_path.is_file())
            self.assertGreater(save_path.stat().st_size, 0)

            missing = output / "missing.modes"
            invalid = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(REPO / "plotting" / "functions" / "ex1.py"),
                    "--no-show",
                    "--mode-file",
                    str(missing),
                ],
                cwd=temp,
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(invalid.returncode, 2)
            self.assertIn("Cannot read mode catalogue", invalid.stderr)

            truncated = output / "truncated.modes"
            truncated.write_text(
                " mode phs vel w(mhz)\n n type l x w\n1\n",
                encoding="utf-8",
            )
            invalid = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(REPO / "plotting" / "functions" / "ex1.py"),
                    "--no-show",
                    "--mode-file",
                    str(truncated),
                ],
                cwd=temp,
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(invalid.returncode, 2)
            self.assertIn(str(truncated), invalid.stderr)
            self.assertIn("line 3: expected at least 5 columns", invalid.stderr)


if __name__ == "__main__":
    unittest.main()

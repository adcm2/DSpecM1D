import os
from pathlib import Path


_REPO_ROOT = Path(__file__).resolve().parents[3]


def default_build_dir():
    """Return the configured build tree or the repository's standard Release tree."""
    configured = os.environ.get("DSPECM1D_BUILD_DIR")
    if configured:
        return Path(configured).expanduser().resolve()
    return _REPO_ROOT / "build" / "release"


def build_output_path(filename, build_dir=None):
    """Resolve an example product beneath the selected build tree's output folder."""
    root = Path(build_dir) if build_dir is not None else default_build_dir()
    return str((root / "output" / Path(filename).name).resolve())


def project_file_path(*parts):
    """Resolve a source-controlled plotting input independently of cwd."""
    return str((_REPO_ROOT.joinpath(*parts)).resolve())

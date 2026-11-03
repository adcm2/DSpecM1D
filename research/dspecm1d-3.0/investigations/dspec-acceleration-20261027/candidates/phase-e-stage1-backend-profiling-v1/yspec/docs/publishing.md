# Publishing the Python package

The package is `yspec-seismo` on PyPI, imported as `yspec`. GitHub Actions
builds it and uploads it (`.github/workflows/wheels.yml`); nothing is built or
uploaded from a laptop, and no token is stored anywhere. PyPI trusts the
workflow itself ("trusted publishing").

What a release puts on PyPI:

| file | for |
|------|-----|
| `yspec_seismo-X.Y.Z-cp312-abi3-manylinux_2_28_x86_64.whl` | Linux, Intel/AMD, glibc 2.28 or later |
| `yspec_seismo-X.Y.Z-cp312-abi3-manylinux_2_28_aarch64.whl` | Linux, ARM |
| `yspec_seismo-X.Y.Z-cp312-abi3-macosx_14_0_arm64.whl` | macOS 14 or later on Apple silicon |
| `yspec_seismo-X.Y.Z.tar.gz` | everything else: pip compiles it (C++23 compiler, CMake, network) |

Each wheel serves Python 3.12 and every later Python.

## Once

1. **PyPI.** Log in at <https://pypi.org>, go to *Account settings →
   Publishing*, and add a *pending publisher* (the project does not exist
   yet; the first upload creates it):
   - PyPI project name: `yspec-seismo`
   - Owner: `da380`
   - Repository name: `YSpecpp`
   - Workflow name: `wheels.yml`
   - Environment name: `pypi`
2. **GitHub.** In the repository, *Settings → Environments → New
   environment*, named `pypi`. Adding yourself as a *required reviewer* makes
   every upload wait for your approval, a last look before anything is
   public; worth having.
3. **TestPyPI, for rehearsals** (optional, recommended before the first
   release). The same at <https://test.pypi.org> (a separate account),
   with environment name `testpypi`, and a GitHub environment `testpypi`.

## Each release

1. **Set the version.** It is set in one place, the `project(YSpec VERSION
   ...)` line of `CMakeLists.txt`; the package takes it from there.
2. **Merge to `main` and push.** The wheels are built and tested on every
   push to `main` (*Actions → Wheels*). Wait for them, and for CI, to pass.
3. **Rehearse** (optional). *Actions → Wheels → Run workflow* on `main`
   builds everything and uploads it to TestPyPI. Try it in a fresh
   environment:
   ```sh
   python -m venv /tmp/try && . /tmp/try/bin/activate
   pip install --index-url https://test.pypi.org/simple/ \
     --extra-index-url https://pypi.org/simple/ "yspec-seismo[plot]"
   python -c "import yspec; print(yspec.__version__, yspec.Model.prem())"
   ```
   TestPyPI takes each file once, so a second rehearsal of the same version
   uploads nothing new.
4. **Tag and push the tag**, which is what publishes:
   ```sh
   git tag v0.1.0
   git push origin v0.1.0
   ```
   The workflow builds and tests the wheels again, checks that the tag is
   the version in `CMakeLists.txt`, waits for your approval if you asked for
   it, and uploads. The release appears at
   <https://pypi.org/project/yspec-seismo/>.
5. **Check it.** In a fresh environment, `pip install "yspec-seismo[plot]"`
   and run an example.

## When something goes wrong

- **A build or test fails:** nothing is uploaded; the upload needs every
  wheel and the source. Fix it, delete the tag (`git tag -d v0.1.0`,
  `git push origin :refs/tags/v0.1.0`) and tag again.
- **A release is faulty:** PyPI never takes the same version twice, even
  after a deletion. *Yank* it (on PyPI, *Manage → Releases*), which hides it
  from new installs without breaking pinned ones, and release the next
  patch version.
- **The tag check fails:** the tag and `CMakeLists.txt` disagree; delete the
  tag and tag the right version.

## Building by hand, to look

```sh
pip install build && python -m build --sdist       # the source, in dist/
pip wheel . -w dist                                # a wheel for this machine
```

A wheel built so is tagged `linux_x86_64` and PyPI refuses it: Linux wheels
must be built in the manylinux containers, which the workflow does. With
Docker, `pipx run cibuildwheel --platform linux` makes them locally too.

The page on PyPI is `python/README.md`; its metadata, extras and build
settings are in `pyproject.toml`.

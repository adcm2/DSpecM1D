"""Copy an auditable backup into a separate checkout; never edit source trees."""
from pathlib import Path
import hashlib
import io
import json
import re
import shutil
import subprocess
import tarfile

WORKSPACE = Path(__file__).resolve().parents[2]
ORIGINAL = Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
DEST = WORKSPACE / 'backups/dspecm1d-3.0-20261103'
RESEARCH = DEST / 'research/dspecm1d-3.0'
BRANCH = 'backup/dspecm1d-3.0-20261103'
ROOTS = {
    'dspec-acceleration-20261027': 'complete',
    'dspec-scientific-validation-20261103': 'complete',
    'dspec-yspec-wedges-20261027': 'complete',
    'dspec-sk-profile-20261026': 'support',
    'tunya-profile-20261027': 'support',
}

def git(*args):
    return subprocess.check_output(['git', '-C', str(ORIGINAL), *args])

def sha(data):
    return hashlib.sha256(data).hexdigest()

def reason(path, root, mode):
    relative = path.relative_to(root)
    if any(s in ('.git', '__pycache__', '.pytest_cache', 'tmp', 'CMakeFiles') for s in relative.parts):
        return 'cache, temporary render, or Git/build metadata'
    if path.suffix in ('.pyc', '.pyo', '.o', '.a', '.so', '.aux', '.out', '.tmp'):
        return 'compiled/cache/temporary output'
    with path.open('rb') as f:
        magic = f.read(4)
    if magic == b'\x7fELF':
        return 'compiled executable/shared object'
    if mode == 'support':
        if path.name.endswith(('.tar.gz', '.tar', '.zip')) or re.search(r'\.part\d+$', path.name):
            return 'earlier profiling transfer/archive; retained at original location'
        if path.suffix == '.bin' or path.stat().st_size > 1_000_000:
            return 'earlier profiling raw/large output; retained at original location'
    return None

SECRET = re.compile(rb'-----BEGIN (?:OPENSSH|RSA|EC|DSA) PRIVATE KEY-----|gh[pousr]_[A-Za-z0-9]{30,}|AKIA[A-Z0-9]{16}')

def check_sensitive(name, data):
    assert not SECRET.search(data), f'Possible credential in {name}; do not publish'
    if name.endswith('.tar.gz') or name.endswith('.tar'):
        with tarfile.open(fileobj=io.BytesIO(data), mode='r:*') as archive:
            for member in archive:
                assert not any(s in ('.ssh', '.git', '.env') for s in Path(member.name).parts), f'Unexpected private member in {name}'
                if member.isfile() and member.size < 2_000_000:
                    assert not SECRET.search(archive.extractfile(member).read()), f'Possible credential in archive {name}'

assert subprocess.check_output(['git', '-C', str(DEST), 'branch', '--show-current']).decode().strip() == BRANCH
assert not RESEARCH.exists(), 'Backup must be assembled exclusively once'
RESEARCH.mkdir(parents=True)
status_before = git('status', '--porcelain=v1', '-uall').decode()
base = git('rev-parse', 'HEAD').decode().strip()
baseline = []
changed = git('diff', '--name-only', '-z').decode().split('\0')
untracked = git('ls-files', '--others', '--exclude-standard', '-z').decode().split('\0')
for name in sorted(set(filter(None, changed + untracked))):
    assert name in changed or name.startswith('validation/phase0/'), ('Unexpected untracked scope', name)
    source = ORIGINAL / name
    assert source.is_file() and not source.is_symlink()
    data = source.read_bytes()
    check_sensitive(name, data)
    target = DEST / name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    baseline.append({'path': name, 'bytes': len(data), 'sha256': sha(data)})

(RESEARCH / 'baseline-worktree.patch').write_bytes(git('diff', '--binary', 'HEAD'))
included, excluded = [], []
for name, mode in ROOTS.items():
    root = WORKSPACE / 'investigations' / name
    for source in sorted(root.rglob('*')):
        assert not source.is_symlink(), ('Unexpected link', source)
        if not source.is_file():
            continue
        rel = source.relative_to(WORKSPACE)
        why = reason(source, root, mode)
        if why:
            excluded.append({'source': str(rel), 'bytes': source.stat().st_size, 'reason': why})
            continue
        data = source.read_bytes()
        assert len(data) < 95_000_000, ('GitHub individual-file limit', source)
        check_sensitive(str(rel), data)
        target = RESEARCH / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        assert not target.exists()
        shutil.copy2(source, target)
        assert sha(target.read_bytes()) == sha(data)
        included.append({'source': str(rel), 'destination': str(target.relative_to(DEST)), 'bytes': len(data), 'sha256': sha(data)})

assert git('status', '--porcelain=v1', '-uall').decode() == status_before
for entry in baseline:
    assert sha((ORIGINAL / entry['path']).read_bytes()) == entry['sha256']
manifest = {
    'schema': 'dspecm1d-3-backup-manifest-v1', 'date': '2026-11-03',
    'branch': BRANCH, 'remote': 'https://github.com/adcm2/DSpecM1D.git',
    'base_commit': base, 'original_repository': str(ORIGINAL),
    'workspace': str(WORKSPACE), 'original_status_before_and_after': status_before,
    'baseline_worktree_files': baseline, 'included': included, 'excluded': excluded,
    'included_bytes': sum(e['bytes'] for e in included),
    'note': 'Recent acceleration/scientific source, reviews, evidence and archives preserved. Earlier profiling support is selective; excluded raw dumps/transfer archives remain local/Tunya and are not remotely backed up. Existing provenance files retain their available hashes. This is not a standalone dependency/binary distribution.',
}
(RESEARCH / 'backup-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(json.dumps({'base': base, 'baseline_files': len(baseline), 'included_files': len(included), 'included_MB': manifest['included_bytes']/1e6, 'excluded_files': len(excluded)}, indent=2))

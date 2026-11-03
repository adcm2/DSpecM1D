"""Verify and extract the immutable campaign into fresh paths (Python 3.11)."""
from pathlib import Path, PurePosixPath
import hashlib
import json
import os
import stat
import tarfile
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def safe(name):
    p = PurePosixPath(name)
    if p.is_absolute() or '..' in p.parts or not p.parts:
        raise RuntimeError('unsafe archive member: ' + name)
    if p.parts[0] not in ('original', 'dependencies'):
        raise RuntimeError('unexpected archive root: ' + name)
    target = ROOT.joinpath(*p.parts)
    for parent in target.parents:
        if parent == ROOT:
            break
        if parent.is_symlink():
            raise RuntimeError('symlink parent: ' + str(parent))
    return target


def main():
    manifest = json.loads((HERE / 'manifest.json').read_text())
    directories = json.loads((HERE / 'directories.json').read_text())
    for name in ('original', 'dependencies'):
        if os.path.lexists(ROOT / name):
            raise RuntimeError('refusing existing extraction root: ' + name)
    archive = HERE / 'verified-campaign.tar.gz'
    aggregate = hashlib.sha256()
    count = 0
    with archive.open('xb') as out:
        for part in manifest['parts']:
            h = hashlib.sha256()
            n = 0
            with (HERE / part['name']).open('rb') as inp:
                for block in iter(lambda: inp.read(8 * 1024 * 1024), b''):
                    out.write(block)
                    h.update(block)
                    aggregate.update(block)
                    n += len(block)
            if n != part['bytes'] or h.hexdigest() != part['sha256']:
                raise RuntimeError('part verification failed: ' + part['name'])
            count += n
            print('verified', part['name'], flush=True)
    if count != manifest['archive']['bytes'] or aggregate.hexdigest() != manifest['archive']['sha256']:
        raise RuntimeError('aggregate archive verification failed')
    for name in sorted(directories, key=lambda n: len(PurePosixPath(n).parts)):
        safe(name).mkdir(parents=True, exist_ok=True)
    seen = set()
    expected = manifest['files']
    with tarfile.open(archive, mode='r|gz') as tf:
        for member in tf:
            name = member.name
            if name not in expected or name in seen:
                raise RuntimeError('unexpected/duplicate archive member: ' + name)
            target = safe(name)
            if os.path.lexists(target):
                raise RuntimeError('refusing existing member: ' + name)
            record = expected[name]
            target.parent.mkdir(parents=True, exist_ok=True)
            if member.isfile() and record['kind'] == 'file':
                if member.size != record['bytes']:
                    raise RuntimeError('archive size mismatch: ' + name)
                h = hashlib.sha256()
                with tf.extractfile(member) as inp, target.open('xb') as out:
                    for block in iter(lambda: inp.read(8 * 1024 * 1024), b''):
                        h.update(block)
                        out.write(block)
                if h.hexdigest() != record['sha256'] or digest(target) != record['sha256']:
                    raise RuntimeError('payload hash mismatch: ' + name)
                os.chmod(target, record['mode'])
                os.utime(target, ns=(record['mtime_ns'], record['mtime_ns']))
                if stat.S_IMODE(target.stat().st_mode) != record['mode']:
                    raise RuntimeError('mode mismatch: ' + name)
            elif member.issym() and record['kind'] == 'symlink':
                link = record.get('target', record.get('link_target'))
                if member.linkname != link:
                    raise RuntimeError('symlink mismatch: ' + name)
                os.symlink(link, target)
                os.utime(target, ns=(record['mtime_ns'], record['mtime_ns']), follow_symlinks=False)
            else:
                raise RuntimeError('unsupported or mismatched archive type: ' + name)
            seen.add(name)
    if seen != set(expected):
        raise RuntimeError('missing archive members')
    for name in sorted(directories, key=lambda n: len(PurePosixPath(n).parts), reverse=True):
        target = safe(name)
        os.chmod(target, directories[name]['mode'])
        os.utime(target, ns=(directories[name]['mtime_ns'], directories[name]['mtime_ns']))
    report = {'status': 'verified', 'files_and_symlinks': len(seen),
              'original_directories': len(directories), 'archive': manifest['archive'],
              'dependencies': manifest['dependencies'], 'finished_unix': time.time()}
    with (HERE / 'transfer-verification.json').open('x') as f:
        json.dump(report, f, indent=2)
        f.write('\n')
    with (ROOT / 'implementation_status.md').open('a') as f:
        f.write('\nTransfer verified: all %d file/symlink payloads and %d original directories preserved; exact dependency revisions retained separately. No pre-existing files removed or replaced.\n' % (len(seen), len(directories)))
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()

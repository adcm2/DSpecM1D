"""Archive the complete original campaign and exact dependency sources."""
from pathlib import Path
import datetime, hashlib, json, os, stat, subprocess, tarfile

here = Path(__file__).resolve().parent
source = here.parents[1] / 'dspec-sk-profile-20261026'
archive = here / 'campaign-and-dependencies.tar.gz'
assert not archive.exists()
files = []
for base, dirs, names in os.walk(source, followlinks=False):
    for name in list(dirs):
        path = Path(base) / name
        if path.is_symlink():
            files.append((path, 'original/dspec-sk-profile-20261026/' + str(path.relative_to(source))))
            dirs.remove(name)
    for name in names:
        path = Path(base) / name
        files.append((path, 'original/dspec-sk-profile-20261026/' + str(path.relative_to(source))))
deps = Path('/home/adcm2/Documents/Research/Projects/DSpecM1D/build/phase-d-build-cache')
dep_sources = {name: deps / 'build-dspec/_deps' / (name + '-src') for name in
               ['numericconcepts', 'fftwpp', 'gaussquad', 'interpolation', 'gshtrans']}
dep_sources['planetmodel-dspec'] = deps / 'build-dspec/_deps/planetmodel-src'
dep_sources['planetmodel-yspec'] = deps / 'build-yspecpp/_deps/planetmodel-src'
dependency_records = {}
env = dict(os.environ, GIT_OPTIONAL_LOCKS='0')
for name, path in dep_sources.items():
    def git(*args):
        return subprocess.check_output(['git', '-C', str(path), *args], env=env)
    status = git('status', '--porcelain').decode()
    assert not status, (name, status)
    dependency_records[name] = {'source': str(path), 'git_head': git('rev-parse', 'HEAD').decode().strip()}
    for rel in git('ls-files', '-z').split(b'\0'):
        if rel:
            rel = os.fsdecode(rel)
            files.append((path / rel, f'dependencies/{name}/{rel}'))

class HashReader:
    def __init__(self, stream):
        self.stream, self.hash = stream, hashlib.sha256()
    def read(self, size=-1):
        value = self.stream.read(size)
        self.hash.update(value)
        return value

records = {}
with tarfile.open(archive, 'w:gz', compresslevel=1) as output:
    for index, (path, relative) in enumerate(sorted(files, key=lambda item: item[1])):
        before = path.lstat()
        info = output.gettarinfo(str(path), arcname=relative)
        record = {'source': str(path), 'mode': stat.S_IMODE(before.st_mode),
                  'mtime_ns': before.st_mtime_ns}
        if path.is_symlink():
            record.update(kind='symlink', target=os.readlink(path))
            output.addfile(info)
        else:
            assert path.is_file(), path
            with path.open('rb') as stream:
                reader = HashReader(stream)
                output.addfile(info, reader)
            record.update(kind='file', bytes=before.st_size, sha256=reader.hash.hexdigest())
        after = path.lstat()
        assert (before.st_size, before.st_mtime_ns, before.st_ino) == (after.st_size, after.st_mtime_ns, after.st_ino), path
        records[relative] = record
        if index % 500 == 0:
            print(index, 'of', len(files), relative, flush=True)

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024**2), b''):
            h.update(chunk)
    return h.hexdigest()

parts = []
with archive.open('rb') as stream:
    index = 0
    while True:
        chunk = stream.read(4 * 1024**2)
        if not chunk:
            break
        part = here / f'campaign.part{index:03}'
        h = hashlib.sha256()
        size = 0
        with part.open('xb') as output:
            while chunk and size < 768 * 1024**2:
                output.write(chunk); h.update(chunk); size += len(chunk)
                if size >= 768 * 1024**2:
                    break
                chunk = stream.read(4 * 1024**2)
        parts.append({'name': part.name, 'bytes': size, 'sha256': h.hexdigest()})
        index += 1
manifest = {'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'source_campaign': str(source), 'files': records,
            'dependencies': dependency_records, 'parts': parts,
            'archive': {'bytes': archive.stat().st_size, 'sha256': sha256(archive)}}
(here / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
(here / 'STAGED').write_text('Complete archive and manifest; original files preserved.\n')
print('STAGED', len(records), 'files', archive.stat().st_size, 'compressed bytes', len(parts), 'parts', flush=True)

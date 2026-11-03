"""Detached, bounded wait for all immutable parts, then verified extraction."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import time

p = Path(__file__).resolve().parent
extractor = p / 'extract_on_tunya.py'
expected_hash = sys.argv[1]
deadline = time.monotonic() + 7200
manifest = json.loads((p / 'manifest.json').read_text())
while True:
    ready = all((p / part['name']).is_file() and
                (p / part['name']).stat().st_size == part['bytes']
                for part in manifest['parts'])
    if ready:
        break
    if time.monotonic() > deadline:
        raise RuntimeError('transfer did not finish within two hours; all partial files preserved')
    time.sleep(5)
if hashlib.sha256(extractor.read_bytes()).hexdigest() != expected_hash:
    raise RuntimeError('extractor hash changed')
subprocess.run(['/usr/bin/python3', str(extractor)], check=True)

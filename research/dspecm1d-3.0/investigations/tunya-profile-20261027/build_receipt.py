"""Run the reviewed sequential builder and bind its outputs to a task receipt."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

tool, config, receipt = map(Path, sys.argv[1:])
subprocess.run([sys.executable, str(tool), str(config)], check=True)
cfg = json.loads(config.read_text())
root = Path(cfg['build_root'])
report = json.loads((root / 'build_campaign.json').read_text())
if report['status'] != 'completed':
    raise RuntimeError('build did not complete')
for variant in report['variants']:
    v = next(x for x in cfg['variants'] if x['id'] == variant['id'])
    binary = root / variant['id'] / v.get('output_name', 'program')
    if hashlib.sha256(binary.read_bytes()).hexdigest() != variant['binary_sha256']:
        raise RuntimeError('binary changed after build: ' + str(binary))
with receipt.open('x') as out:
    json.dump(report, out, indent=2)
    out.write('\n')

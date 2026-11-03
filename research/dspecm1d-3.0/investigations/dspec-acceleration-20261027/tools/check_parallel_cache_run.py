#!/usr/bin/env python3
"""Run bounded cache tests and distinguish unsupported TSan from a clean run."""
from pathlib import Path
import argparse,json,subprocess,sys
p=argparse.ArgumentParser();p.add_argument('binary');p.add_argument('output');p.add_argument('--tsan',action='store_true');a=p.parse_args()
r=subprocess.run([a.binary],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=120)
report={'binary':a.binary,'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr}
if a.tsan and r.returncode!=0 and 'ThreadSanitizer: unexpected memory mapping' in r.stderr:
 report.update(status='UNSUPPORTED_RUNTIME',sanitizer_passed=False)
elif r.returncode==0:
 report.update(status='PASS',test=json.loads(r.stdout),sanitizer_passed=a.tsan)
else:report.update(status='FAIL',sanitizer_passed=False)
with Path(a.output).open('x') as f:json.dump(report,f,indent=2)
sys.exit(1 if report['status']=='FAIL' else 0)

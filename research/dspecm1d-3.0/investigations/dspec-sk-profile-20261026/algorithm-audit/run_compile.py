"""Shared two-slot compiler semaphore for all diagnostic implementation scopes."""
from pathlib import Path
import fcntl,os,subprocess,sys,time
root=Path(__file__).resolve().parent/'compile-locks';root.mkdir(exist_ok=True)
command=sys.argv[1:]
if command and command[0]=='--':command=command[1:]
if not command:raise SystemExit('usage: run_compile.py -- COMPILER ARG...')
held=None
while held is None:
 for slot in range(2):
  f=(root/f'slot-{slot}.lock').open('a+')
  try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
  except BlockingIOError:f.close();continue
  held=f;f.seek(0);f.truncate();f.write(f'pid={os.getpid()} start={time.time()} command={command!r}\n');f.flush();break
 if held is None:time.sleep(.2)
try:raise SystemExit(subprocess.run(command).returncode)
finally:fcntl.flock(held,fcntl.LOCK_UN);held.close()

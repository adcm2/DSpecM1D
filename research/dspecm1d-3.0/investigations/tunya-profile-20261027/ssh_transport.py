"""Use the established Tunya route without modifying local SSH state."""
import subprocess
import sys
import shlex

OPTIONS = ['-F', '/tmp/ccdsm-tunya-20261027/ssh_config',
           '-o', 'ControlMaster=no', '-o', 'ControlPath=none',
           '-o', 'UpdateHostKeys=no', '-o', 'AddKeysToAgent=no',
           '-o', 'ProxyCommand=ssh -F /tmp/ccdsm-tunya-20261027/ssh_config -o ControlMaster=no -o ControlPath=none -o UpdateHostKeys=no -o AddKeysToAgent=no -W %h:%p ssh-west']


def python(source):
    result = subprocess.run(['ssh', *OPTIONS, '-T', 'tunya', 'python3 -'],
                            input=source, text=True, capture_output=True)
    if result.returncode:
        sys.stdout.write(result.stdout)
        sys.stderr.write(result.stderr)
        result.check_returncode()
    return result.stdout


def send(paths, destination):
    subprocess.run(['scp', '-B', *OPTIONS, *map(str, paths), 'tunya:' + destination], check=True)


if __name__ == '__main__':
    if sys.argv[1] == 'python':
        print(python(sys.stdin.read()), end='')
    elif sys.argv[1] == 'send':
        send(sys.argv[3:], sys.argv[2])
    else:
        raise SystemExit('usage: ssh_transport.py python | send DEST FILE...')

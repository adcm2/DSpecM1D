# Supervisor v3 snapshot

This directory contains a separate supervisor snapshot that fixes final log-tail loss from nonblocking pipes. `EAGAIN` is treated as “no bytes available yet”; only a zero-length read marks EOF. A normally completed job log is closed and hashed only after EOF, leader completion, and the existing owned-process/thread check. After a watchdog stop, identity cleanup runs first, then the pipe is drained through EOF. The drain is bounded; if cleanup cannot confirm gone writers or EOF remains unavailable, `incomplete_pipe_drain` is recorded in the job state and the log is still hashed without an unbounded wait. The deployed `tools/supervise_campaign.py` and prior snapshots are unchanged.

The snapshot passed Python syntax compilation and 100 synthetic fast-failure cases that alternated short/long final writes to stdout/stderr. Each case verified the exact `job.log` bytes, SHA-256, failed status, and exit status. A synthetic watchdog-timeout case also verified exact final stderr from a SIGTERM handler, a verified-gone cleanup result, and the log hash. These are tooling regression checks; they do not validate scientific results.


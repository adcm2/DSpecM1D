# Supervisor v3 snapshot

This directory contains a separate supervisor snapshot that fixes final log-tail loss from nonblocking pipes. `EAGAIN` is treated as “no bytes available yet”; only a zero-length read marks EOF. A job log is closed and hashed only after EOF, leader completion, and the existing owned-process/thread check. The deployed `tools/supervise_campaign.py` and prior snapshots are unchanged.

The snapshot passed Python syntax compilation and 100 synthetic fast-failure cases that alternated short/long final writes to stdout/stderr. Each case verified the exact `job.log` bytes, SHA-256, failed status, and exit status. This is a tooling regression check; it does not validate scientific results.

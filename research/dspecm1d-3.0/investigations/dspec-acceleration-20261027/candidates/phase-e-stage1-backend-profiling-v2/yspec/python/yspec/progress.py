"""Progress of a calculation, which the library reports frequency by frequency.

Each calculation takes `progress`: None (the default) to report nothing, True
to print a line at every tenth of the frequencies as the programs do, or a
callable taking the number of frequencies done and the number in all.
Whichever it is, Ctrl-C stops the calculation at the next frequency to
finish.
"""

import sys


class Reporter:
    """Prints progress at every tenth of the frequencies, as the programs do."""

    def __init__(self, name, stream=None):
        self.name = name
        self.stream = stream
        self.reported = 0

    def __call__(self, done, total):
        percent = 100 * done // total
        if percent >= self.reported + 10 or done == total:
            self.reported = percent - percent % 10
            print(f"{self.name}: {percent:3d}% of {total} frequencies",
                  file=self.stream or sys.stderr, flush=True)


class StageReporter:
    """The same for the stages of the kernels, a Reporter for each."""

    def __init__(self, name, stream=None):
        self.name = name
        self.stream = stream
        self.stages = {}

    def __call__(self, stage, done, total):
        if stage not in self.stages:
            self.stages[stage] = Reporter(f"{self.name}: {stage}", self.stream)
        self.stages[stage](done, total)


def reporter(progress, name):
    """What the library is given for a `progress` argument."""
    if progress is None or progress is False:
        return None
    if progress is True:
        return Reporter(name)
    if not callable(progress):
        raise TypeError("progress should be None, True or a callable taking "
                        "(done, total)")
    return progress

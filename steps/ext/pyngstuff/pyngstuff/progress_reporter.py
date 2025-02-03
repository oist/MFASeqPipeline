import sys
import progress.bar

class progress_reporter:
    class bar(progress.bar.Bar):
        def __init__(self, reporter, *args, **kwargs):
            super().__init__(suffix='%(statusline)s', *args, **kwargs)
            self.reporter = reporter

        @property
        def statusline(self):
            return self.reporter.current_statusline

    def __init__(self, iterator, every=1000, total=None, show_progressbar=True, target=sys.stderr, statusline=lambda n: "processed %d" % n):
        self.iterator = iterator
        self.total = total
        self.every = min(every, self.total / 100) if self.total is not None else every
        self.target = target
        self.statusline = statusline
        self.n = 0;
        self.reporting = False
        self.iteration_stopped = False
        if self.total is not None and show_progressbar and self.target.isatty():
            self.progressbar = progress_reporter.bar(self, max=self.total)
        else:
            self.progressbar = None

    def __iter__(self):
        return self

    def __next__(self):
        try:
            read = next(self.iterator)
            self.n += 1
            if self.reporting and self.n % self.every == 0:
                if self.progressbar is not None:
                    self.progressbar.next(n=self.every)
                else:
                    print(self.current_statusline, file=self.target, end="\r", flush=True)
            return read
        except StopIteration:
            if self.reporting and not self.iteration_stopped:
                if self.progressbar is not None:
                    self.progressbar.finish()
                else:
                    print(self.current_statusline, file=self.target, end="\n", flush=True)
                self.iteration_stopped = True
            raise

    def __enter__(self):
        if self.progressbar is not None:
            self.progressbar.__enter__()
        self.reporting = True
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.reporting = False
        if self.progressbar is not None:
            self.progressbar.__exit__(exc_type, exc_val, exc_tb)

    @property
    def current_statusline(self):
        return self.statusline(self.n, self.total)

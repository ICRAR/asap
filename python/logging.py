"""This module presents a logging abstraction layer on top of casa.
"""
__all__ = ["asaplog", "print_log", "print_log_dec", "AsapLogger"]

from asap.env import is_casapy
from asap.parameters import rcParams
from asap._asap import LogSink, set_global_sink
try:
    from functools import wraps as wraps_dec
except ImportError:
    from asap.compatibility import wraps as wraps_dec


class AsapLogger(object):
    """Wrapper object to allow for both casapy and asap logging.

    Inside casapy this will connect to `taskinit.casalog`. Otherwise it will
    create its own casa log sink.

    .. note:: Do instantiate a new one - use the :obj:`asaplog` instead.
    
    """
    def __init__(self):
        self._enabled = False
        self._log = ""
        if is_casapy():
            from taskinit import casalog
            self.logger = casalog
        else:
            self.logger = LogSink()
            set_global_sink(self.logger)

    def post(self, level):
        """Post the messages to the logger. This will clear the buffered
        logs.

        Parameters:

            level:  The log level (severity). One of INFO, WARN, ERROR.

        """
        if not self._enabled:
            return
        if not rcParams['verbose']:
            return

        logs = self._log.strip()
        if len(logs) > 0:
           self.logger.post(logs, priority=level)
        if isinstance(self.logger, LogSink):
            logs = self.logger.pop().strip()
            if len(logs) > 0:
                print logs
        self._log = ""

    def push(self, msg, newline=True):
        """Push logs into the buffer. post needs to be called to send them.

        Parameters:

            msg:        the log message (string)

            newline:    should we terminate with a newline (default yes)

        """
        from asap import rcParams
        if self._enabled:
            if rcParams["verbose"]:
                sep = ""
                self._log = sep.join([self._log, msg])
                if newline:
                    self._log += "\n"

    def enable(self, flag=True):
        """Enable (or disable) logging."""
        self._enabled = flag

    def disable(self, flag=False):
        """Disable (or enable) logging"""
        self._enabled = flag

asaplog = AsapLogger()
"""Default asap logger"""

def print_log_dec(f, level='INFO'):
    """Decorator which posts log at completion of the wrapped method.
    Example::

        @print_log_dec
        def test(self):
            do_stuff()
            asaplog.push('testing...', False)
            do_more_stuff()
            asaplog.push('finished')
    """
    @wraps_dec(f)
    def wrap_it(*args, **kw):
        val = f(*args, **kw)
        print_log(level)
        return val
    return wrap_it

def print_log(level='INFO'):
    """Alias for asaplog.post(level)"""
    asaplog.post(level)

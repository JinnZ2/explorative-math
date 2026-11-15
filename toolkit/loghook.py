"""loghook -- every run writes the schema. CC0. Co-authored. Stdlib only.

log_run() wraps emlib.write_run for automatic runs: both interest slots and
the divergence note are void (null), because no one has looked yet.
logged() is a decorator for a function returning (params, outputs, angles,
symbols, definition); calling the wrapped function writes the run file.
"""
import functools

import emlib


def log_run(name, seed, params, outputs, angles=None, symbols=None,
            definition=None, nature_match=None, out_dir=None):
    comp = emlib.compressibility(symbols, definition) if symbols is not None else None
    return emlib.write_run(name, seed, params, outputs, angles, None, comp,
                           nature_match, out_dir=out_dir)


def logged(name, out_dir=None):
    def wrap(fn):
        @functools.wraps(fn)
        def inner(seed, out_dir_override=None):
            params, outputs, angles, symbols, definition = fn(seed)
            return log_run(name, seed, params, outputs, angles, symbols, definition,
                           out_dir=out_dir_override or out_dir)
        return inner
    return wrap

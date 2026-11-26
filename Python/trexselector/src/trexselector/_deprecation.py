"""Deprecation warnings of the package."""
from __future__ import annotations

import sys
import warnings

_PACKAGE = __name__.partition(".")[0]


def warn_deprecated(old: str, new: str) -> None:
    """Raise the ``FutureWarning`` that ``old`` is deprecated in favour of ``new``.

    The warning is attributed to the first caller outside the package, so it
    points at the user's line also when a wrapper of the package sits between.
    """
    level = 2
    frame = sys._getframe(1)
    while frame is not None and frame.f_globals.get("__name__", "").partition(".")[0] == _PACKAGE:
        frame = frame.f_back
        level += 1
    warnings.warn(
        f"{old} is deprecated. Use {new}.",
        FutureWarning,
        stacklevel=level,
    )

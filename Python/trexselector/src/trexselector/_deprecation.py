"""Deprecation warnings of the package."""
from __future__ import annotations

import sys
import warnings

_PACKAGE = __name__.partition(".")[0]


def warn_deprecated(old: str, new: str) -> None:
    """Raise the ``FutureWarning`` that ``old`` is deprecated in favour of ``new``.

    The warning is attributed to the first caller outside the package, so it
    points at the user's line also when a wrapper of the package sits between.
    The compiled bindings call this function where a solver class of the
    coordinate-descent family is constructed and where a solver type of the
    family reaches a selector.
    """
    level = 2
    frame = sys._getframe(1)
    while frame is not None and frame.f_globals.get("__name__", "").partition(".")[0] == _PACKAGE:
        frame = frame.f_back
        level += 1
    warnings.warn(
        f"{old} is deprecated and will be removed in the next minor release of trexselector. "
        f"Use {new}, the homotopy solver of the same problem.",
        FutureWarning,
        stacklevel=level,
    )

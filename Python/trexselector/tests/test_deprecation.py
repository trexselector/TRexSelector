"""
The deprecation warning of the package.

``trexselector._deprecation.warn_deprecated`` raises one FutureWarning that names the old name and
its replacement, and it attributes the warning to the first caller outside the package.
"""
import warnings

import pytest

from trexselector._deprecation import warn_deprecated


def _deprecations(fn):
    """Run fn and return the FutureWarnings it raised."""
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        fn()
    return [w for w in record if issubclass(w.category, FutureWarning)]


def test_warning_names_the_old_name_and_the_replacement():
    record = _deprecations(lambda: warn_deprecated("old_name", "new_name"))
    assert len(record) == 1
    assert str(record[0].message) == "old_name is deprecated. Use new_name."


def test_warning_is_attributed_to_the_caller_outside_the_package():
    record = _deprecations(lambda: warn_deprecated("old_name", "new_name"))
    assert record[0].filename == __file__


def test_warning_turned_into_an_error_raises():
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        with pytest.raises(FutureWarning, match="old_name is deprecated"):
            warn_deprecated("old_name", "new_name")

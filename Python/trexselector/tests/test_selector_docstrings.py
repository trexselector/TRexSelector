"""
Tests for the documentation of the selectors of the package root.
"""
import inspect

import pytest
import trexselector
from trexselector.trex_selector_methods import ad as core

SELECTORS = [
    "TRexSelector",
    "TRexDASelector",
    "TRexGVSSelector",
    "TRexTikhonovSelector",
    "TRexScreeningSelector",
    "TRexBiobankScreeningSelector",
    "TRexSPCASelector",
]


@pytest.mark.parametrize("name", SELECTORS)
def test_root_selector_documents_itself(name):
    """A selector of the package root says what it is and lists its parameters."""
    cls = getattr(trexselector, name)
    doc = inspect.getdoc(cls)
    assert not doc.startswith("Python wrapper")
    assert "\nParameters\n----------\n" in doc
    for parameter in inspect.signature(cls.__init__).parameters:
        if parameter != "self":
            assert f"\n{parameter} : " in doc


@pytest.mark.parametrize("name", SELECTORS)
def test_root_selector_carries_the_core_sections(name):
    """The notes and the reference list are those of the core class."""
    doc = inspect.getdoc(getattr(trexselector, name))
    core_doc = inspect.getdoc(getattr(core, name))
    heading = "Notes\n-----" if "Notes\n-----" in core_doc else "References\n----------"
    sections = core_doc[core_doc.index(heading):]
    assert "References\n----------" in sections
    assert doc.rstrip().endswith(sections.rstrip())


def test_da_selector_carries_its_notes():
    """The dependency-aware selector shows the notes on its penalty and its grid."""
    doc = inspect.getdoc(trexselector.TRexDASelector)
    for phrase in ("The penalty of the last step", "The reference threshold", "FDR limitation"):
        assert phrase in doc

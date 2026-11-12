"""
Tests for TRexBiobankScreeningSelector.
"""
import numpy as np
import pytest
from trexselector import (
    TRexBiobankScreeningSelector,
    BiobankScreenTRexControl,
    BiobankScreenTRexResult,
)


# ---------------------------------------------------------------------------
# BiobankScreenTRexControl defaults
# ---------------------------------------------------------------------------

def test_biobank_control_defaults():
    ctrl = BiobankScreenTRexControl()
    assert isinstance(ctrl.target_FDR_trex, float)
    assert isinstance(ctrl.lower_bound_FDR, float)
    assert isinstance(ctrl.upper_bound_FDR, float)
    assert ctrl.lower_bound_FDR <= ctrl.upper_bound_FDR


def test_biobank_control_assignment():
    ctrl = BiobankScreenTRexControl()
    ctrl.target_FDR_trex = 0.05
    ctrl.lower_bound_FDR = 0.01
    ctrl.upper_bound_FDR = 0.2
    assert ctrl.target_FDR_trex == pytest.approx(0.05)
    assert ctrl.lower_bound_FDR == pytest.approx(0.01)
    assert ctrl.upper_bound_FDR == pytest.approx(0.2)


def test_biobank_nested_screen_control_settable():
    """The nested Screen-TRex control (and its nested base T-Rex control) must
    be reachable and mutable, matching the R package's screen_control/control
    knobs."""
    ctrl = BiobankScreenTRexControl()
    ctrl.trex_screen_ctrl.trex_ctrl.K = 7
    ctrl.trex_screen_ctrl.use_bootstrap_CI = True
    assert ctrl.trex_screen_ctrl.trex_ctrl.K == 7
    assert ctrl.trex_screen_ctrl.use_bootstrap_CI is True


def test_biobank_bio_ctrl_kwarg(signal_data):
    """The control argument is named bio_ctrl (matches the C++ ctor)."""
    X, y, n, p = signal_data
    ctrl = BiobankScreenTRexControl()
    ctrl.trex_screen_ctrl.trex_ctrl.K = 5
    sel = TRexBiobankScreeningSelector(X, y, bio_ctrl=ctrl, verbose=False)
    assert sel is not None


# ---------------------------------------------------------------------------
# Single phenotype (1-D y) — screenPhenotype()
# ---------------------------------------------------------------------------

def test_instantiation_1d(signal_data):
    X, y, n, p = signal_data
    sel = TRexBiobankScreeningSelector(X, y, verbose=False)
    assert sel is not None


def test_select_returns_result(signal_data):
    X, y, n, p = signal_data
    sel = TRexBiobankScreeningSelector(X, y, verbose=False)
    res = sel.select()
    assert res is not None
    assert isinstance(res, BiobankScreenTRexResult)


def test_result_phenotype_index(signal_data):
    X, y, n, p = signal_data
    sel = TRexBiobankScreeningSelector(X, y, verbose=False)
    res = sel.select()
    assert isinstance(res, BiobankScreenTRexResult)
    assert res.phenotype_index == 0


def test_result_selected_indices(signal_data):
    X, y, n, p = signal_data
    sel = TRexBiobankScreeningSelector(X, y, verbose=False)
    res = sel.select()
    assert isinstance(res, BiobankScreenTRexResult)
    assert isinstance(res.selected_indices, list)
    if res.selected_indices:
        assert all(0 <= i < p for i in res.selected_indices)


def test_result_estimated_fdr(signal_data):
    X, y, n, p = signal_data
    sel = TRexBiobankScreeningSelector(X, y, verbose=False)
    res = sel.select()
    assert isinstance(res, BiobankScreenTRexResult)
    assert isinstance(res.estimated_FDR, float)
    assert 0.0 <= res.estimated_FDR <= 1.0


def test_result_method_used(signal_data):
    X, y, n, p = signal_data
    sel = TRexBiobankScreeningSelector(X, y, verbose=False)
    res = sel.select()
    assert isinstance(res, BiobankScreenTRexResult)
    assert isinstance(res.method_used, str)
    assert len(res.method_used) > 0


def test_result_used_fallback_trex(signal_data):
    X, y, n, p = signal_data
    sel = TRexBiobankScreeningSelector(X, y, verbose=False)
    res = sel.select()
    assert isinstance(res, BiobankScreenTRexResult)
    assert isinstance(res.used_fallback_trex, bool)


def test_result_auxiliary_index_lists(signal_data):
    X, y, n, p = signal_data
    sel = TRexBiobankScreeningSelector(X, y, verbose=False)
    res = sel.select()
    assert isinstance(res, BiobankScreenTRexResult)
    assert isinstance(res.selected_indices_screen_ordinary, list)
    assert isinstance(res.selected_indices_screen_bootstrap, list)


# ---------------------------------------------------------------------------
# Multiple phenotypes (2-D Y) — screenPhenotypes()
# ---------------------------------------------------------------------------

def test_instantiation_2d(signal_data):
    X, y, n, p = signal_data
    Y = np.column_stack([y, y])
    sel = TRexBiobankScreeningSelector(X, Y, verbose=False)
    assert sel is not None


def test_select_2d_returns_list(signal_data):
    X, y, n, p = signal_data
    Y = np.column_stack([y, y])
    sel = TRexBiobankScreeningSelector(X, Y, verbose=False)
    results = sel.select()
    assert isinstance(results, list)
    assert len(results) == 2


def test_select_2d_phenotype_indices(signal_data):
    X, y, n, p = signal_data
    Y = np.column_stack([y, y])
    sel = TRexBiobankScreeningSelector(X, Y, verbose=False)
    results = sel.select()
    assert isinstance(results, list)
    indices = [r.phenotype_index for r in results]
    assert sorted(indices) == [0, 1]


def test_select_2d_result_types(signal_data):
    X, y, n, p = signal_data
    Y = np.column_stack([y, y])
    sel = TRexBiobankScreeningSelector(X, Y, verbose=False)
    results = sel.select()
    assert isinstance(results, list)
    for res in results:
        assert isinstance(res, BiobankScreenTRexResult)
        assert isinstance(res.selected_indices, list)
        assert isinstance(res.estimated_FDR, float)
        assert 0.0 <= res.estimated_FDR <= 1.0


def test_phenotypes_use_independent_seeds():
    """Phenotype i runs with seed + i, and a single-phenotype run with that seed reproduces it."""
    rng = np.random.default_rng(31)
    n, p, copies = 80, 40, 12
    X = np.asfortranarray(rng.standard_normal((n, p)))
    y_noise = rng.standard_normal(n)

    # Twelve copies of one noise response make the dummy draw visible: the screening of such
    # a response selects something in about half of the draws
    Y = np.asfortranarray(np.tile(y_noise[:, None], (1, copies)))
    ctrl = BiobankScreenTRexControl()
    ctrl.trex_screen_ctrl.R_boot = 100
    ctrl.trex_screen_ctrl.trex_ctrl.K = 10

    multi = TRexBiobankScreeningSelector(X.copy(order="F"), Y, bio_ctrl=ctrl, seed=11, verbose=False).select()
    assert len(multi) == copies

    fields = ("selected_indices", "selected_indices_screen_ordinary", "selected_indices_screen_bootstrap",
              "estimated_FDR_screen_ordinary", "estimated_FDR_screen_bootstrap", "method_used")
    for i, res in enumerate(multi):
        single = TRexBiobankScreeningSelector(X.copy(order="F"), y_noise.copy(), bio_ctrl=ctrl, seed=11 + i,
                                              verbose=False).select()
        for field in fields:
            assert getattr(res, field) == getattr(single, field), (i, field)

    # With one shared dummy draw the twelve copies would give twelve identical results
    first = multi[0].selected_indices_screen_ordinary
    assert any(res.selected_indices_screen_ordinary != first for res in multi)

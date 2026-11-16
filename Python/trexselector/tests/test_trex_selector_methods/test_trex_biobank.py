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


def test_fallback_result_reports_the_target(signal_data):
    """A fallback result reports the target level, the flag marks it, and the screening
    estimates stay in their own fields."""
    X, y, n, p = signal_data
    ctrl = BiobankScreenTRexControl()
    ctrl.lower_bound_FDR = 0.0001   # an acceptance window no estimate 1 / R can meet
    ctrl.upper_bound_FDR = 0.0002
    ctrl.target_FDR_trex = 0.2
    res = TRexBiobankScreeningSelector(X.copy(order="F"), y.copy(), bio_ctrl=ctrl, seed=3, verbose=False).select()

    assert res.used_fallback_trex is True
    assert res.method_used == "T-Rex (fallback)"
    assert res.estimated_FDR == pytest.approx(0.2)
    assert 0.0 < res.estimated_FDR_screen_ordinary <= 1.0
    assert 0.0 < res.estimated_FDR_screen_bootstrap <= 1.0


def test_phenotypes_use_independent_seeds():
    """Phenotype i runs with seed + i, and a single-phenotype run with that seed reproduces it."""
    rng = np.random.default_rng(31)
    n, p, seed = 80, 40, 11
    X = np.asfortranarray(rng.standard_normal((n, p)))
    y_noise = rng.standard_normal(n)
    ctrl = BiobankScreenTRexControl()
    ctrl.trex_screen_ctrl.R_boot = 100
    ctrl.trex_screen_ctrl.trex_ctrl.K = 10
    fields = ("selected_indices", "selected_indices_screen_ordinary", "selected_indices_screen_bootstrap",
              "estimated_FDR_screen_ordinary", "estimated_FDR_screen_bootstrap", "method_used")

    def single(run_seed):
        res = TRexBiobankScreeningSelector(X.copy(order="F"), y_noise.copy(), bio_ctrl=ctrl, seed=run_seed,
                                           verbose=False).select()
        return tuple(getattr(res, field) for field in fields)

    # The dummy draw has to show in the result. The screening of this noise response has few
    # outcomes, and most seeds give the same one, so no fixed number of seeds is sure to show a
    # second one. The single-phenotype runs go on until a seed changes the result.
    limit = 1000
    singles = [single(seed)]
    while len(singles) < limit and singles[-1] == singles[0]:
        singles.append(single(seed + len(singles)))
    copies = len(singles)
    assert singles[-1] != singles[0], "1,000 seeds in a row gave one result"

    # Copies of the response: phenotype i runs with seed + i, and the single-phenotype run with
    # that seed reproduces it. With one shared dummy draw the last copy would equal the first.
    Y = np.asfortranarray(np.tile(y_noise[:, None], (1, copies)))
    multi = TRexBiobankScreeningSelector(X.copy(order="F"), Y, bio_ctrl=ctrl, seed=seed, verbose=False).select()
    assert len(multi) == copies
    for i, res in enumerate(multi):
        for field, expected in zip(fields, singles[i]):
            assert getattr(res, field) == expected, (i, field)


def _biobank_design():
    rng = np.random.default_rng(11)
    n, p = 100, 60
    X = np.asfortranarray(rng.standard_normal((n, p)))
    y = X[:, :5] @ np.full(5, 3.0) + rng.standard_normal(n)
    return X, y


def _biobank_run(X, y, lower, upper, use_bootstrap_ci=False):
    ctrl = BiobankScreenTRexControl()
    ctrl.lower_bound_FDR = lower
    ctrl.upper_bound_FDR = upper
    ctrl.trex_screen_ctrl.use_bootstrap_CI = use_bootstrap_ci
    ctrl.trex_screen_ctrl.trex_ctrl.K = 10
    return TRexBiobankScreeningSelector(X, y, bio_ctrl=ctrl, seed=3).select()


@pytest.mark.parametrize("lower, upper", [(0.01, 0.05), (0.05, 0.15), (0.1, 0.25), (0.15, 0.3), (0.2, 1.0),
                                           (0.3, 0.5)])
def test_biobank_follows_the_documented_decision_rule(lower, upper):
    """A variant is accepted inside the bounds, the larger accepted estimate wins, the bootstrap on a tie."""
    X, y = _biobank_design()
    res = _biobank_run(X, y, lower, upper)
    ordinary, bootstrap = res.estimated_FDR_screen_ordinary, res.estimated_FDR_screen_bootstrap
    ok_ordinary = lower <= ordinary <= upper
    ok_bootstrap = lower <= bootstrap <= upper
    if ok_bootstrap and (not ok_ordinary or bootstrap >= ordinary):
        expected = "Screen-TRex (bootstrap-CI)"
    elif ok_ordinary:
        expected = "Screen-TRex (ordinary)"
    else:
        expected = "T-Rex (fallback)"
    assert res.method_used == expected
    assert res.used_fallback_trex == (expected == "T-Rex (fallback)")


def test_biobank_computes_both_variants_whatever_use_bootstrap_ci_says():
    """use_bootstrap_CI of the nested screening control has no effect in the biobank screening."""
    X, y = _biobank_design()
    off = _biobank_run(X, y, 0.05, 0.15, use_bootstrap_ci=False)
    on = _biobank_run(X, y, 0.05, 0.15, use_bootstrap_ci=True)
    for field in ("selected_indices", "estimated_FDR", "method_used", "estimated_FDR_screen_ordinary",
                  "estimated_FDR_screen_bootstrap", "selected_indices_screen_ordinary",
                  "selected_indices_screen_bootstrap", "used_fallback_trex"):
        assert getattr(off, field) == getattr(on, field)


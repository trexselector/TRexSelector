"""
Tests for TRexSPCASelector (trexselector.TRexSPCASelector).
"""
import warnings

import numpy as np
import pytest
from trexselector import (
    ENSolverType,
    GVSType,
    LLoopStrategy,
    SPCAMode,
    SolverTypeForTRex,
    TRexSPCAControlParameter,
    TRexSPCASelector,
)
from trexselector.trex_selector_methods.ad import LambdaSelectionMethod


# ---------------------------------------------------------------------------
# Shared fixture — small dataset to keep tests fast
# ---------------------------------------------------------------------------

@pytest.fixture
def spca_fixture():
    rng = np.random.default_rng(42)
    n, p = 30, 10
    X = rng.standard_normal((n, p)).astype(np.float64)
    return X, n, p


M    = 2
TFDR = 0.2


def make_ctrl(**kwargs):
    ctrl = TRexSPCAControlParameter()
    ctrl.gvs_ctrl.trex_ctrl.K = kwargs.get("K", 3)
    ctrl.gvs_ctrl.trex_ctrl.max_dummy_multiplier = kwargs.get("max_dummy_multiplier", 3)
    ctrl.mode = kwargs.get("mode", SPCAMode.ActiveSet)
    if "lambda2_ridge_loadings" in kwargs:
        ctrl.lambda2_ridge_loadings = kwargs["lambda2_ridge_loadings"]
    if "gvs_type" in kwargs:
        ctrl.gvs_ctrl.gvs_type = kwargs["gvs_type"]
    if "lambda_2" in kwargs:
        ctrl.gvs_ctrl.lambda_2 = kwargs["lambda_2"]
    return ctrl


# ---------------------------------------------------------------------------
# Group 1: TRexSPCAControlParameter defaults
# ---------------------------------------------------------------------------

def test_spca_control_default_mode():
    ctrl = TRexSPCAControlParameter()
    assert ctrl.mode == SPCAMode.ActiveSet


def test_spca_control_default_lambda2_ridge_loadings():
    ctrl = TRexSPCAControlParameter()
    assert np.isclose(ctrl.lambda2_ridge_loadings, 1e-6)


def test_spca_control_default_gvs_type():
    ctrl = TRexSPCAControlParameter()
    assert ctrl.gvs_ctrl.gvs_type == GVSType.EN


def test_spca_control_default_lambda_2():
    ctrl = TRexSPCAControlParameter()
    # -1.0 is the "not supplied" sentinel that triggers auto-selection of lambda_2.
    assert ctrl.gvs_ctrl.lambda_2 == -1.0


def test_spca_control_thresholded_mode():
    ctrl = TRexSPCAControlParameter()
    ctrl.mode = SPCAMode.Thresholded
    assert ctrl.mode == SPCAMode.Thresholded


# ---------------------------------------------------------------------------
# Group 2: Basic execution
# ---------------------------------------------------------------------------

def test_spca_select_runs_activeset_en(spca_fixture):
    X, _, _ = spca_fixture
    ctrl = make_ctrl()
    sel = TRexSPCASelector(X, ctrl, seed=1)
    sel.select(M, TFDR)  # must not raise


def test_spca_select_runs_thresholded_en(spca_fixture):
    X, _, _ = spca_fixture
    ctrl = make_ctrl(mode=SPCAMode.Thresholded)
    sel = TRexSPCASelector(X, ctrl, seed=1)
    sel.select(M, TFDR)  # must not raise


def test_spca_rejects_ien(spca_fixture):
    # T-Rex+SPCA is defined on the elastic net: IEN raises instead of running EN
    X, _, _ = spca_fixture
    ctrl = make_ctrl(gvs_type=GVSType.IEN)
    sel = TRexSPCASelector(X, ctrl, seed=1)
    with pytest.raises(ValueError, match="gvs_type must be EN"):
        sel.select(M, TFDR)


def test_spca_select_runs_with_lambda_2_bypass(spca_fixture):
    X, _, _ = spca_fixture
    ctrl = make_ctrl(lambda_2=0.01)
    sel = TRexSPCASelector(X, ctrl, seed=1)
    sel.select(M, TFDR)  # must not raise


# ---------------------------------------------------------------------------
# Group 3: Output dimensions
# ---------------------------------------------------------------------------

def test_spca_Z_shape(spca_fixture):
    X, n, _ = spca_fixture
    sel = TRexSPCASelector(X, make_ctrl(), seed=1)
    sel.select(M, TFDR)
    assert sel.Z.shape == (n, M)


def test_spca_V_shape(spca_fixture):
    X, _, p = spca_fixture
    sel = TRexSPCASelector(X, make_ctrl(), seed=1)
    sel.select(M, TFDR)
    assert sel.V.shape == (p, M)


def test_spca_active_sets_length(spca_fixture):
    X, _, _ = spca_fixture
    sel = TRexSPCASelector(X, make_ctrl(), seed=1)
    sel.select(M, TFDR)
    assert len(sel.active_sets) == M


def test_spca_explained_variance_length(spca_fixture):
    X, _, _ = spca_fixture
    sel = TRexSPCASelector(X, make_ctrl(), seed=1)
    sel.select(M, TFDR)
    assert sel.adjusted_ev.shape == (M,)
    assert sel.cumulative_ev.shape == (M,)


def test_spca_M1_shapes(spca_fixture):
    X, n, p = spca_fixture
    sel = TRexSPCASelector(X, make_ctrl(), seed=1)
    sel.select(1, TFDR)
    assert sel.Z.shape == (n, 1)
    assert sel.V.shape == (p, 1)
    assert len(sel.active_sets) == 1


# ---------------------------------------------------------------------------
# Group 4: Active set bounds (0-based indices)
# ---------------------------------------------------------------------------

def test_spca_active_set_indices_in_range(spca_fixture):
    X, _, p = spca_fixture
    sel = TRexSPCASelector(X, make_ctrl(), seed=1)
    sel.select(M, TFDR)
    for as_m in sel.active_sets:
        for idx in as_m:
            assert 0 <= idx < p


# ---------------------------------------------------------------------------
# Group 5: Explained variance properties
# ---------------------------------------------------------------------------

def test_spca_adjusted_ev_non_negative(spca_fixture):
    X, _, _ = spca_fixture
    sel = TRexSPCASelector(X, make_ctrl(), seed=1)
    sel.select(M, TFDR)
    assert np.all(sel.adjusted_ev >= 0)


def test_spca_cumulative_ev_non_decreasing(spca_fixture):
    X, _, _ = spca_fixture
    sel = TRexSPCASelector(X, make_ctrl(), seed=1)
    sel.select(M, TFDR)
    assert np.all(np.diff(sel.cumulative_ev) >= -1e-10)


def test_spca_cumulative_ev_in_unit_interval(spca_fixture):
    X, _, _ = spca_fixture
    sel = TRexSPCASelector(X, make_ctrl(), seed=1)
    sel.select(M, TFDR)
    assert np.all(sel.cumulative_ev >= 0)
    assert np.all(sel.cumulative_ev <= 1.0 + 1e-10)


# ---------------------------------------------------------------------------
# Group 6: Data integrity — X restored after select()
# ---------------------------------------------------------------------------

def test_spca_X_restored_after_select(spca_fixture):
    X, _, _ = spca_fixture
    X_orig = X.copy()  # genuine deep copy
    sel = TRexSPCASelector(X, make_ctrl(), seed=1)
    sel.select(M, TFDR)
    assert np.allclose(X, X_orig, atol=1e-12)


# ---------------------------------------------------------------------------
# Group 7: Input validation
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad_M", [0, -1, 11])
def test_spca_rejects_M_outside_its_range(spca_fixture, bad_M):
    """The core refuses M below 1 and above min(n, p) = 10."""
    X, _, _ = spca_fixture
    sel = TRexSPCASelector(X, make_ctrl())
    with pytest.raises(ValueError, match="M must"):
        sel.select(bad_M, TFDR)


@pytest.mark.parametrize("bad_tFDR", [-0.1, 1.5, float("nan")])
def test_spca_rejects_tFDR_outside_the_unit_interval(spca_fixture, bad_tFDR):
    X, _, _ = spca_fixture
    sel = TRexSPCASelector(X, make_ctrl())
    with pytest.raises(ValueError, match=r"tFDR must be in \[0, 1\]"):
        sel.select(M, bad_tFDR)


@pytest.mark.parametrize("end", [0.0, 1.0])
def test_spca_takes_the_two_ends_of_tFDR(spca_fixture, end):
    """0 and 1 are valid targets, as for the other selectors."""
    X, _, _ = spca_fixture
    sel = TRexSPCASelector(X, make_ctrl(), seed=1)
    sel.select(M, end)
    assert len(sel.active_sets) == M


def test_spca_properties_raise_before_select(spca_fixture):
    X, _, _ = spca_fixture
    sel = TRexSPCASelector(X, make_ctrl())
    with pytest.raises(RuntimeError):
        _ = sel.Z


# ---------------------------------------------------------------------------
# Group 8: Documented behaviour of the control and the result
# ---------------------------------------------------------------------------

def test_spca_verbose_prints_one_line_per_component(spca_fixture, capfd):
    """The root class is quiet by default and passes verbose on."""
    import ctypes
    import sys
    X, _, _ = spca_fixture

    def c_stdout():
        # The lines are written by the C++ core, so the C stdio buffer is flushed before
        # reading. CDLL(None) is the running process on POSIX. Windows has no such handle,
        # and the UCRT shared by CPython and the extension is loaded by name.
        libc = ctypes.CDLL("ucrtbase" if sys.platform == "win32" else None)
        libc.fflush(None)
        return capfd.readouterr().out

    TRexSPCASelector(X, make_ctrl(), seed=1).select(M, TFDR)
    assert "Processing PC" not in c_stdout()
    TRexSPCASelector(X, make_ctrl(), seed=1, verbose=True).select(M, TFDR)
    assert c_stdout().count("Processing PC") == M


def test_spca_nested_solver_fields_have_no_effect(spca_fixture):
    """The selector replaces gvs_ctrl.en_solver and gvs_ctrl.trex_ctrl.solver_type by en_solver."""
    X, _, _ = spca_fixture
    plain = make_ctrl(lambda_2=1.0)
    nested = make_ctrl(lambda_2=1.0)
    nested.gvs_ctrl.en_solver = ENSolverType.TCDENET
    nested.gvs_ctrl.trex_ctrl.solver_type = SolverTypeForTRex.TLASSO
    # One sweep would end a coordinate-descent run with an error.
    nested.gvs_ctrl.trex_ctrl.solver_params.cd_max_sweeps = 1

    a = TRexSPCASelector(X, plain, seed=1).select(M, TFDR)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        b = TRexSPCASelector(X, nested, seed=1).select(M, TFDR)
    assert a.active_sets == b.active_sets
    assert np.array_equal(a.V, b.V)


def test_spca_outer_tcdenet_is_deprecated(spca_fixture):
    X, _, _ = spca_fixture
    ctrl = make_ctrl(lambda_2=1.0)
    ctrl.en_solver = ENSolverType.TCDENET
    ctrl.gvs_ctrl.trex_ctrl.solver_params.cd_max_sweeps = 200000
    with pytest.warns(FutureWarning, match="TCDENET"):
        TRexSPCASelector(X, ctrl, seed=1).select(M, TFDR)


@pytest.mark.parametrize("rule", ["CV_1SE_IEN_CCD", "CV_MIN_IEN_CCD", "CV_1SE_TIK_SVD", "CV_MIN_TIK_SVD"])
def test_spca_rejects_the_rules_of_the_informed_elastic_net(spca_fixture, rule):
    X, _, _ = spca_fixture
    ctrl = make_ctrl()
    ctrl.gvs_ctrl.lambda2_method = getattr(LambdaSelectionMethod, rule)
    with pytest.raises(ValueError, match="rule of the informed elastic net"):
        TRexSPCASelector(X, ctrl, seed=1).select(M, TFDR)


@pytest.mark.parametrize("strategy", ["PERMUTATION", "PERMUTATION_SEEDED"])
def test_spca_rejects_the_permutation_strategies(spca_fixture, strategy):
    X, _, _ = spca_fixture
    ctrl = make_ctrl()
    ctrl.gvs_ctrl.trex_ctrl.lloop_strategy = getattr(LLoopStrategy, strategy)
    with pytest.raises(ValueError, match="PERMUTATION"):
        TRexSPCASelector(X, ctrl, seed=1).select(M, TFDR)


@pytest.mark.parametrize("strategy", ["STANDARD", "HCONCAT", "SKIPL", "SEEDED"])
def test_spca_runs_the_strategies_of_the_gvs_selector(spca_fixture, strategy):
    X, _, _ = spca_fixture
    ctrl = make_ctrl()
    ctrl.gvs_ctrl.trex_ctrl.lloop_strategy = getattr(LLoopStrategy, strategy)
    sel = TRexSPCASelector(X, ctrl, seed=1).select(M, TFDR)
    assert len(sel.active_sets) == M


def test_spca_cumulative_ev_is_a_fraction_of_the_total_variance(spca_fixture):
    X, n, _ = spca_fixture
    sel = TRexSPCASelector(X, make_ctrl(), seed=1).select(M, TFDR)
    Xc = X - X.mean(axis=0)
    total_variance = np.sum(Xc ** 2) / (n - 1)
    assert np.allclose(sel.cumulative_ev, np.cumsum(sel.adjusted_ev) / total_variance, rtol=1e-12)
    assert np.allclose(sel.Z, Xc @ sel.V, atol=1e-10)


def test_spca_component_without_selection_is_zero(spca_fixture):
    """The rule holds for every component, whatever is selected."""
    X, _, _ = spca_fixture
    for target in (0.0, 0.05, 0.2):
        sel = TRexSPCASelector(X, make_ctrl(), seed=1).select(3, target)
        previous = 0.0
        for m in range(3):
            if len(sel.active_sets[m]) == 0:
                assert not sel.V[:, m].any()
                assert not sel.Z[:, m].any()
                assert sel.adjusted_ev[m] == 0.0
                assert sel.cumulative_ev[m] == previous
            else:
                assert np.isclose(np.linalg.norm(sel.V[:, m]), 1.0, atol=1e-12)
            previous = sel.cumulative_ev[m]


def test_spca_thresholded_keeps_as_many_loadings_as_selected(spca_fixture):
    X, _, _ = spca_fixture
    sel = TRexSPCASelector(X, make_ctrl(mode=SPCAMode.Thresholded), seed=1).select(M, TFDR)
    assert [int(np.count_nonzero(sel.V[:, m])) for m in range(M)] == [len(a) for a in sel.active_sets]


def test_spca_copies_an_array_that_is_not_fortran_ordered(spca_fixture):
    """Only a Fortran-ordered float64 array is used without a copy."""
    X, _, _ = spca_fixture
    X_c = np.ascontiguousarray(X)
    before = X_c.copy()
    sel = TRexSPCASelector(X_c, make_ctrl(), seed=1)
    assert not np.shares_memory(sel.X, X_c)
    sel.select(M, TFDR)
    assert np.array_equal(X_c, before)

    X_f = np.asfortranarray(X)
    assert np.shares_memory(TRexSPCASelector(X_f, make_ctrl(), seed=1).X, X_f)

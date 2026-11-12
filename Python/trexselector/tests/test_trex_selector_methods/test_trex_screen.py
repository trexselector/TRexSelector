"""
Tests for TRexScreeningSelector.
"""
import numpy as np
import pytest
import trexselector


# ---------------------------------------------------------------------------
# Instantiation
# ---------------------------------------------------------------------------

def test_instantiation(signal_data, fast_control):
    X, y, n, p = signal_data
    sel = trexselector.TRexScreeningSelector(X, y, trex_control=fast_control, verbose=False)
    assert sel is not None


def test_dimension_mismatch_raises(signal_data, fast_control):
    X, y, n, p = signal_data
    y_bad = np.random.randn(n + 10)
    with pytest.raises(Exception):
        trexselector.TRexScreeningSelector(X, y_bad, trex_control=fast_control, verbose=False)


# ---------------------------------------------------------------------------
# Execution and result fields
# ---------------------------------------------------------------------------

def test_execution_and_result(signal_data, fast_control):
    X, y, n, p = signal_data
    sel = trexselector.TRexScreeningSelector(X, y, trex_control=fast_control, verbose=False)
    res = sel.select()

    assert res is not None
    assert isinstance(sel.selected_var, np.ndarray)
    assert isinstance(sel.selected_indices, list)

    if sel.selected_indices:
        assert all(0 <= i < p for i in sel.selected_indices)


def test_screen_result_fields(signal_data, fast_control):
    X, y, n, p = signal_data
    sel = trexselector.TRexScreeningSelector(X, y, trex_control=fast_control, verbose=False)
    screen_res = sel.select()
    assert screen_res is not None
    assert isinstance(screen_res.estimated_FDR, float)
    assert isinstance(screen_res.used_bootstrap, bool)
    assert isinstance(screen_res.confidence_level, float)
    assert isinstance(screen_res.beta_mat, np.ndarray)
    assert isinstance(screen_res.dummy_betas, np.ndarray)
    assert isinstance(screen_res.estimated_correlation, float)


# ---------------------------------------------------------------------------
# Screen control parameter fields
# ---------------------------------------------------------------------------

def test_screen_control_n_blocks():
    ctrl = trexselector.ScreenTRexControlParameter()
    ctrl.n_blocks = 5
    assert ctrl.n_blocks == 5


def test_screen_control_fields():
    ctrl = trexselector.ScreenTRexControlParameter()
    ctrl.use_bootstrap_CI = False
    ctrl.R_boot = 50
    ctrl.cor_coef = 0.5
    assert ctrl.use_bootstrap_CI is False
    assert ctrl.R_boot == 50
    assert ctrl.cor_coef == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# Properties inherited from TRexSelector base
# ---------------------------------------------------------------------------

def test_base_properties(signal_data, fast_control):
    X, y, n, p = signal_data
    sel = trexselector.TRexScreeningSelector(X, y, trex_control=fast_control, verbose=False)
    sel.select()

    assert isinstance(sel.tFDR, float)
    assert isinstance(sel.T_stop, int)
    assert isinstance(sel.fdp_hat_mat, np.ndarray)


# ---------------------------------------------------------------------------
# ScreenTRexMethod variants
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("method", [
    trexselector.ScreenTRexMethod.TREX,
    trexselector.ScreenTRexMethod.TREX_DA_AR1,
    trexselector.ScreenTRexMethod.TREX_DA_EQUI,
    trexselector.ScreenTRexMethod.TREX_DA_BLOCK_EQUI,
])
def test_screen_method_variants(signal_data, method):
    X, y, n, p = signal_data
    ctrl = trexselector.TRexControlParameter()
    ctrl.K = 3
    screen_ctrl = trexselector.ScreenTRexControlParameter()
    screen_ctrl.trex_method = method
    if method == trexselector.ScreenTRexMethod.TREX_DA_BLOCK_EQUI:
        screen_ctrl.n_blocks = 2
    sel = trexselector.TRexScreeningSelector(X, y, screen_control=screen_ctrl,
                                               trex_control=ctrl, verbose=False)
    res = sel.select()
    assert res is not None
    assert isinstance(res.estimated_FDR, float)
    assert 0.0 <= res.estimated_FDR <= 1.0


@pytest.mark.parametrize("method, rho_thr", [
    (trexselector.ScreenTRexMethod.TREX, 1.5),
    (trexselector.ScreenTRexMethod.TREX_DA_EQUI, -0.5),
    (trexselector.ScreenTRexMethod.TREX_DA_EQUI, float("nan")),
    (trexselector.ScreenTRexMethod.TREX_DA_AR1, 1.0),
])
def test_screen_rejects_invalid_rho_threshold(signal_data, method, rho_thr):
    """rho_thr_DA is checked at construction: [0, 1], and [0, 1) for TREX_DA_AR1."""
    X, y, n, p = signal_data
    screen_ctrl = trexselector.ScreenTRexControlParameter()
    screen_ctrl.trex_method = method
    screen_ctrl.rho_thr_DA = rho_thr
    with pytest.raises(ValueError, match="rho_thr_DA"):
        trexselector.TRexScreeningSelector(X, y, screen_control=screen_ctrl, verbose=False)


@pytest.mark.parametrize("n_blocks", [0, 21])
def test_screen_rejects_invalid_block_count(signal_data, n_blocks):
    """TREX_DA_BLOCK_EQUI takes between 1 and p blocks (p = 20)."""
    X, y, n, p = signal_data
    screen_ctrl = trexselector.ScreenTRexControlParameter()
    screen_ctrl.trex_method = trexselector.ScreenTRexMethod.TREX_DA_BLOCK_EQUI
    screen_ctrl.n_blocks = n_blocks
    with pytest.raises(ValueError, match="n_blocks"):
        trexselector.TRexScreeningSelector(X, y, screen_control=screen_ctrl, verbose=False)


def test_screen_ar1_window_edges(signal_data):
    """TREX_DA_AR1: the window follows abs(cor_coef), and a zero threshold is the full window."""
    X, y, n, p = signal_data
    ctrl = trexselector.TRexControlParameter()
    ctrl.K = 10

    def run(cor_coef, rho_thr):
        screen_ctrl = trexselector.ScreenTRexControlParameter()
        screen_ctrl.trex_method = trexselector.ScreenTRexMethod.TREX_DA_AR1
        screen_ctrl.cor_coef = cor_coef
        screen_ctrl.rho_thr_DA = rho_thr
        sel = trexselector.TRexScreeningSelector(X.copy(order="F"), y.copy(), screen_control=screen_ctrl,
                                                 trex_control=ctrl, seed=7, verbose=False)
        return sel.select()

    neg, pos = run(-0.5, 0.02), run(0.5, 0.02)
    assert np.array_equal(neg.selected_var, pos.selected_var)
    np.testing.assert_allclose(neg.phi_mat, pos.phi_mat)

    full, unit = run(0.5, 0.0), run(1.0, 0.02)
    assert np.all(np.isfinite(full.phi_mat))
    np.testing.assert_allclose(full.phi_mat, unit.phi_mat)


@pytest.mark.parametrize("method", [
    trexselector.ScreenTRexMethod.TREX_DA_AR1,
    trexselector.ScreenTRexMethod.TREX_DA_EQUI,
    trexselector.ScreenTRexMethod.TREX_DA_BLOCK_EQUI,
])
def test_screen_unusable_coefficient_is_estimated(signal_data, method, capfd):
    """A cor_coef outside [-1, 1] is replaced by the estimate, with a warning on stderr."""
    X, y, n, p = signal_data
    ctrl = trexselector.TRexControlParameter()
    ctrl.K = 5

    def run(cor_coef):
        screen_ctrl = trexselector.ScreenTRexControlParameter()
        screen_ctrl.trex_method = method
        screen_ctrl.cor_coef = cor_coef
        sel = trexselector.TRexScreeningSelector(X.copy(order="F"), y.copy(), screen_control=screen_ctrl,
                                                 trex_control=ctrl, seed=7, verbose=False)
        capfd.readouterr()
        coefficient = sel.select().estimated_correlation
        return coefficient, capfd.readouterr().err

    estimated, warning = run(-2.0)
    assert "cor_coef" not in warning
    for unusable in (-1.5, 1.5, -2.5, float("nan")):
        coefficient, warning = run(unusable)
        assert coefficient == pytest.approx(estimated)
        assert "cor_coef" in warning

    coefficient, warning = run(-0.5)
    assert coefficient == pytest.approx(-0.5)
    assert "cor_coef" not in warning


def test_screen_auto_estimate_cor(signal_data):
    """cor_coef = -2.0 triggers auto-estimation of the correlation coefficient."""
    X, y, n, p = signal_data
    ctrl = trexselector.TRexControlParameter()
    ctrl.K = 3
    screen_ctrl = trexselector.ScreenTRexControlParameter()
    screen_ctrl.trex_method = trexselector.ScreenTRexMethod.TREX_DA_AR1
    screen_ctrl.cor_coef = -2.0  # sentinel: auto-estimate
    sel = trexselector.TRexScreeningSelector(X, y, screen_control=screen_ctrl,
                                               trex_control=ctrl, verbose=False)
    res = sel.select()
    assert res is not None
    # estimated_correlation should have been filled in
    assert isinstance(res.estimated_correlation, float)

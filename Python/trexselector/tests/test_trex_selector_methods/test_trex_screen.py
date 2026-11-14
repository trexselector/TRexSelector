"""
Tests for TRexScreeningSelector.
"""
import warnings

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
    trexselector.ScreenTRexMethod.TREX_DA_NN,
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
    (trexselector.ScreenTRexMethod.TREX_DA_NN, 1.0),
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
    trexselector.ScreenTRexMethod.TREX_DA_NN,
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


@pytest.mark.parametrize("field, value", [
    ("R_boot", 0),
    ("R_boot", 1),
    ("ci_grid_step", 0.0),
    ("ci_grid_step", -0.01),
    ("ci_grid_step", 1.5),
    ("ci_grid_step", float("nan")),
])
@pytest.mark.parametrize("use_bootstrap", [False, True])
def test_screen_rejects_invalid_bootstrap_parameters(signal_data, field, value, use_bootstrap):
    """R_boot is at least 2 and ci_grid_step lies in (0, 1], whether or not the bootstrap is on."""
    X, y, n, p = signal_data
    screen_ctrl = trexselector.ScreenTRexControlParameter()
    screen_ctrl.use_bootstrap_CI = use_bootstrap
    setattr(screen_ctrl, field, value)
    with pytest.raises(ValueError, match=field):
        trexselector.TRexScreeningSelector(X, y, screen_control=screen_ctrl, verbose=False)


@pytest.mark.parametrize("use_bootstrap", [False, True])
def test_screen_empty_selection_reports_one(use_bootstrap):
    """The estimated FDR is one over the number of selected variables, and 1 for an empty selection."""
    rng = np.random.default_rng(21)
    n, p = 80, 40
    X = np.asfortranarray(rng.standard_normal((n, p)))
    y_noise = rng.standard_normal(n)

    # A response orthogonal to the design and to the intercept leaves every variable behind
    # the first dummy, so nothing is selected
    basis = np.column_stack([np.ones(n), X])
    y_orth = y_noise - basis @ np.linalg.lstsq(basis, y_noise, rcond=None)[0]

    screen_ctrl = trexselector.ScreenTRexControlParameter()
    screen_ctrl.use_bootstrap_CI = use_bootstrap
    screen_ctrl.R_boot = 100
    res = trexselector.TRexScreeningSelector(X.copy(order="F"), y_orth, screen_control=screen_ctrl,
                                             seed=5, verbose=False).select()
    assert np.asarray(res.selected_var).sum() == 0
    assert res.estimated_FDR == pytest.approx(1.0)

    y_signal = 5 * X[:, 1] - 5 * X[:, 9] + 5 * X[:, 20] + 0.1 * y_noise
    res = trexselector.TRexScreeningSelector(X.copy(order="F"), y_signal, seed=5, verbose=False).select()
    selected = np.asarray(res.selected_var).sum()
    assert selected > 0
    assert res.estimated_FDR == pytest.approx(1.0 / selected)


def test_screen_nn_groups_and_empty_group(signal_data):
    """TREX_DA_NN: a zero threshold keeps every neighbor (equals TREX_DA_AR1), and a variable
    without a correlated neighbor keeps its relative occurrence (equals plain screening)."""
    X, y, n, p = signal_data
    ctrl = trexselector.TRexControlParameter()
    ctrl.K = 10

    def run(X_run, y_run, method, rho_thr):
        screen_ctrl = trexselector.ScreenTRexControlParameter()
        screen_ctrl.trex_method = method
        screen_ctrl.cor_coef = 0.5
        screen_ctrl.rho_thr_DA = rho_thr
        sel = trexselector.TRexScreeningSelector(X_run.copy(order="F"), y_run.copy(), screen_control=screen_ctrl,
                                                 trex_control=ctrl, seed=7, verbose=False)
        return sel.select()

    method = trexselector.ScreenTRexMethod
    nn_all, ar1_all = run(X, y, method.TREX_DA_NN, 0.0), run(X, y, method.TREX_DA_AR1, 0.0)
    assert np.array_equal(nn_all.selected_var, ar1_all.selected_var)
    np.testing.assert_allclose(nn_all.phi_mat, ar1_all.phi_mat)

    # Centered orthonormal columns have no correlated neighbor
    rng = np.random.default_rng(5)
    centered = rng.standard_normal((n, p))
    centered -= centered.mean(axis=0)
    X_orth = np.linalg.qr(centered)[0]
    y_orth = X_orth[:, 0] * 5 - X_orth[:, 7] * 4 + X_orth[:, 15] * 3 + 0.01 * rng.standard_normal(n)
    nn_none, plain = run(X_orth, y_orth, method.TREX_DA_NN, 0.5), run(X_orth, y_orth, method.TREX, 0.5)
    assert np.array_equal(nn_none.selected_var, plain.selected_var)
    np.testing.assert_allclose(nn_none.phi_mat, plain.phi_mat)
    # A penalty of one half for an empty group would select nothing
    assert np.asarray(nn_none.selected_var).sum() > 0


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


def _screen_design():
    rng = np.random.default_rng(11)
    n, p = 100, 60
    X = np.asfortranarray(rng.standard_normal((n, p)))
    y = X[:, :5] @ np.full(5, 3.0) + rng.standard_normal(n)
    return X, y, n, p


def _screen_phi(X, y, strategy):
    ctrl = trexselector.TRexControlParameter()
    ctrl.K = 10
    ctrl.lloop_strategy = getattr(trexselector.LLoopStrategy, strategy)
    sel = trexselector.TRexScreeningSelector(X, y, trex_control=ctrl, seed=3, verbose=False)
    sel.select()
    return np.asarray(sel.phi_mat).copy(), ctrl


def test_screen_takes_the_four_strategies_of_the_core():
    """The strategy of the control reaches the selector, and each strategy equals its seeded twin."""
    X, y, n, p = _screen_design()
    standard, _ = _screen_phi(X, y, "STANDARD")
    permutation, _ = _screen_phi(X, y, "PERMUTATION")
    assert np.array_equal(_screen_phi(X, y, "SEEDED")[0], standard)
    assert np.array_equal(_screen_phi(X, y, "PERMUTATION_SEEDED")[0], permutation)
    assert not np.array_equal(permutation, standard)


@pytest.mark.parametrize("strategy", ["HCONCAT", "SKIPL"])
def test_screen_refuses_the_growth_strategies(strategy):
    """The screening is one pass at L = p. The root class neither remaps nor changes the control."""
    X, y, n, p = _screen_design()
    ctrl = trexselector.TRexControlParameter()
    ctrl.lloop_strategy = getattr(trexselector.LLoopStrategy, strategy)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(ValueError, match="only the STANDARD, SEEDED, PERMUTATION, and PERMUTATION_SEEDED"):
            trexselector.TRexScreeningSelector(X, y, trex_control=ctrl, verbose=False)
    assert ctrl.lloop_strategy == getattr(trexselector.LLoopStrategy, strategy)


def test_screen_base_attributes_hold_the_fixed_values_of_a_screening_run():
    """One pass at T = 1 and L = p, without a target and without a calibration."""
    X, y, n, p = _screen_design()
    ctrl = trexselector.TRexControlParameter()
    ctrl.K = 10
    sel = trexselector.TRexScreeningSelector(X, y, trex_control=ctrl, seed=3, verbose=False)
    res = sel.select()
    n_selected = len(sel.selected_indices)

    assert sel.T_stop == 1
    assert sel.L == 1
    assert sel.num_dummies == p
    assert sel.v_thresh == 0.5
    assert sel.tFDR == 0.0
    assert np.asarray(sel.r_mat).shape == (1, 1)
    assert np.asarray(sel.r_mat)[0, 0] == n_selected
    assert np.asarray(sel.fdp_hat_mat).shape == (1, 1)
    assert np.asarray(sel.fdp_hat_mat)[0, 0] == 0.0
    assert np.array_equal(np.asarray(sel.phi_prime).ravel(), np.asarray(sel.phi_mat).ravel())
    assert res.estimated_FDR == 1.0 / max(1, n_selected)

    # The fields of the screening without bootstrap and without a filter
    assert not res.used_bootstrap
    assert res.confidence_level == 0.0
    assert res.estimated_correlation == 0.0
    assert np.asarray(res.beta_mat).shape == (p,)
    assert np.asarray(res.dummy_betas).size > 0
    assert np.all(np.asarray(res.dummy_betas) != 0.0)

    # A supplied coefficient of a dependency-aware variant is reported as it is
    screen_ctrl = trexselector.ScreenTRexControlParameter()
    screen_ctrl.trex_method = trexselector.ScreenTRexMethod.TREX_DA_AR1
    screen_ctrl.cor_coef = 0.5
    filtered = trexselector.TRexScreeningSelector(X, y, screen_control=screen_ctrl, trex_control=ctrl, seed=3,
                                                  verbose=False)
    assert filtered.select().estimated_correlation == 0.5


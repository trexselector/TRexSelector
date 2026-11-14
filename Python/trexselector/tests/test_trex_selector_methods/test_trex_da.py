"""
Tests for TRexDASelector.
"""
import numpy as np
import pytest
import trexselector


# ---------------------------------------------------------------------------
# Instantiation
# ---------------------------------------------------------------------------

def test_instantiation(signal_data, fast_control):
    X, y, n, p = signal_data
    sel = trexselector.TRexDASelector(X, y, trex_control=fast_control, verbose=False)
    assert sel is not None


def test_dimension_mismatch_raises(signal_data, fast_control):
    X, y, n, p = signal_data
    y_bad = np.random.randn(n + 10)
    with pytest.raises(Exception):
        trexselector.TRexDASelector(X, y_bad, trex_control=fast_control, verbose=False)


# ---------------------------------------------------------------------------
# Execution and result fields
# ---------------------------------------------------------------------------

def test_execution_and_result(signal_data, fast_control):
    X, y, n, p = signal_data
    sel = trexselector.TRexDASelector(X, y, trex_control=fast_control, verbose=False)
    res = sel.select()

    assert res is not None
    assert isinstance(sel.selected_var, np.ndarray)
    assert isinstance(sel.selected_indices, list)

    # 0-based indices in range
    if sel.selected_indices:
        assert all(0 <= i < p for i in sel.selected_indices)


def test_da_result_fields(signal_data, fast_control):
    X, y, n, p = signal_data
    sel = trexselector.TRexDASelector(X, y, trex_control=fast_control, verbose=False)
    da_res = sel.select()
    assert da_res is not None
    assert isinstance(da_res.rho_thresh, float)
    assert isinstance(da_res.rho_grid, np.ndarray)
    assert isinstance(da_res.method, trexselector.DAMethod)
    assert isinstance(da_res.cor_coef, float)

    # BT arrays may be empty (depends on method) but must be list-like
    assert isinstance(da_res.fdp_hat_array_bt, list)
    assert isinstance(da_res.phi_array_bt, list)


def test_da_control_parameter_fields():
    ctrl = trexselector.TRexDAControlParameter()
    ctrl.rho_thr_DA = 0.5
    assert ctrl.rho_thr_DA == pytest.approx(0.5)

    # bt_selection_mode
    ctrl.bt_selection_mode = trexselector.BTSelectionMode.RFaithful
    assert ctrl.bt_selection_mode == trexselector.BTSelectionMode.RFaithful


# ---------------------------------------------------------------------------
# Properties inherited from TRexSelector base
# ---------------------------------------------------------------------------

def test_base_properties(signal_data, fast_control):
    X, y, n, p = signal_data
    sel = trexselector.TRexDASelector(X, y, trex_control=fast_control, verbose=False)
    sel.select()

    assert isinstance(sel.tFDR, float)
    assert isinstance(sel.T_stop, int)
    assert isinstance(sel.L, int)
    assert isinstance(sel.fdp_hat_mat, np.ndarray)
    assert isinstance(sel.phi_mat, np.ndarray)


# ---------------------------------------------------------------------------
# DAMethod variants
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("method", [
    trexselector.DAMethod.AR1,
    trexselector.DAMethod.EQUI,
    trexselector.DAMethod.BT,
])
def test_da_method_variants(signal_data, method):
    X, y, n, p = signal_data
    ctrl = trexselector.TRexControlParameter()
    ctrl.K = 3
    da_ctrl = trexselector.TRexDAControlParameter()
    da_ctrl.method = method
    sel = trexselector.TRexDASelector(X, y, da_control=da_ctrl,
                                        trex_control=ctrl, verbose=False)
    res = sel.select()
    assert res is not None
    assert res.method == method


def test_da_method_bt_r_array(signal_data):
    """Bootstrap method populates r_array_bt."""
    X, y, n, p = signal_data
    ctrl = trexselector.TRexControlParameter()
    ctrl.K = 3
    da_ctrl = trexselector.TRexDAControlParameter()
    da_ctrl.method = trexselector.DAMethod.BT
    sel = trexselector.TRexDASelector(X, y, da_control=da_ctrl,
                                        trex_control=ctrl, verbose=False)
    res = sel.select()
    assert res is not None
    # r_array_bt should be a non-None array when BT method is used
    assert res.r_array_bt is not None


def test_da_prior_groups(signal_data):
    """PRIOR_GROUPS runs constrained sub-clustering within the given groups."""
    X, y, n, p = signal_data
    ctrl = trexselector.TRexControlParameter()
    ctrl.K = 3
    da_ctrl = trexselector.TRexDAControlParameter()
    da_ctrl.method = trexselector.DAMethod.PRIOR_GROUPS
    # p=20 predictors; assign them to 4 groups
    da_ctrl.prior_groups = [[i % 4 for i in range(p)]]  # one constraint level, 4 groups
    sel = trexselector.TRexDASelector(X, y, da_control=da_ctrl,
                                        trex_control=ctrl, verbose=False)
    res = sel.select()
    assert res is not None
    # Data-driven rho grid: hc_grid_length points (default min(20, p) = 20),
    # ascending and terminated by the conservative rho = 1 singleton anchor.
    rho_grid = np.asarray(res.rho_grid).ravel()
    assert rho_grid.size == min(20, p)
    assert np.all(np.diff(rho_grid) >= 0)
    assert rho_grid[-1] == pytest.approx(1.0)


def test_da_ar1_uses_coefficient_magnitude(signal_data, fast_control):
    """The AR1 window follows abs(cor_coef): -0.5 and 0.5 give the same result."""
    X, y, n, p = signal_data

    def run(cor_coef):
        da_ctrl = trexselector.TRexDAControlParameter()
        da_ctrl.method = trexselector.DAMethod.AR1
        da_ctrl.cor_coef = cor_coef
        sel = trexselector.TRexDASelector(X.copy(order="F"), y.copy(), da_control=da_ctrl,
                                          trex_control=fast_control, seed=7, verbose=False)
        return sel.select()

    neg, pos = run(-0.5), run(0.5)
    assert np.array_equal(neg.selected_var, pos.selected_var)
    np.testing.assert_allclose(neg.phi_prime, pos.phi_prime)
    np.testing.assert_allclose(neg.phi_mat, pos.phi_mat)
    # The field reports the coefficient as supplied
    assert neg.cor_coef == pytest.approx(-0.5)


def test_da_rejects_ward_linkage(signal_data, fast_control):
    """Ward needs squared Euclidean distances, and the correlation distance is not Euclidean."""
    from trexselector.ml_methods.clustering import LinkageMethod

    X, y, n, p = signal_data
    da_ctrl = trexselector.TRexDAControlParameter()
    da_ctrl.hc_linkage = LinkageMethod.Ward
    with pytest.raises(Exception, match="Unsupported linkage method"):
        sel = trexselector.TRexDASelector(X, y, da_control=da_ctrl, trex_control=fast_control,
                                          verbose=False)
        sel.select()


def test_da_ar1_penalty_of_the_last_step_keeps_the_history_monotone():
    """AR1 deflates every step with the penalty of the last step.

    The deflated history is then non-decreasing, so the calibrated occurrences Phi' of the FDP
    estimator never exceed the deflated occurrences of the last step. A penalty per step breaks
    this in about every second run of this design. The design is tall (p < n / 2 - 1), so no
    random experiment exhausts its path within the T-loop.
    """
    n, p, rho = 150, 60, 0.7
    for seed in range(1, 7):
        rng = np.random.default_rng(seed)
        X = np.zeros((n, p), order="F")
        X[:, 0] = rng.standard_normal(n)
        for j in range(1, p):
            X[:, j] = rho * X[:, j - 1] + np.sqrt(1 - rho**2) * rng.standard_normal(n)
        beta = np.zeros(p)
        beta[[9, 29, 49]] = 2.0
        y = X @ beta + rng.standard_normal(n)

        ctrl = trexselector.TRexControlParameter()
        ctrl.K = 10
        da_ctrl = trexselector.TRexDAControlParameter()
        da_ctrl.method = trexselector.DAMethod.AR1
        sel = trexselector.TRexDASelector(X, y, tFDR=0.2, da_control=da_ctrl, trex_control=ctrl,
                                          seed=seed, verbose=False)
        sel.select()

        phi_da = np.asarray(sel.phi_mat)[-1]
        phi_prime = np.asarray(sel.phi_prime).ravel()
        assert np.all(phi_prime <= phi_da + 1e-12), f"Phi' above the deflated Phi, seed {seed}"


@pytest.mark.parametrize("method", [trexselector.DAMethod.AR1, trexselector.DAMethod.EQUI])
def test_da_unusable_coefficient_is_estimated(signal_data, fast_control, method, capfd):
    """A cor_coef outside [-1, 1] is replaced by the estimate, with a warning on stderr."""
    X, y, n, p = signal_data

    def run(cor_coef):
        da_ctrl = trexselector.TRexDAControlParameter()
        da_ctrl.method = method
        da_ctrl.cor_coef = cor_coef
        sel = trexselector.TRexDASelector(X.copy(order="F"), y.copy(), da_control=da_ctrl,
                                          trex_control=fast_control, seed=7, verbose=False)
        capfd.readouterr()
        coefficient = sel.select().cor_coef
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


def test_da_ar1_rejects_unit_threshold(signal_data):
    """AR1 takes thresholds in [0, 1): at 1 the window holds no neighbor. EQUI accepts 1."""
    X, y, n, p = signal_data
    da_ctrl = trexselector.TRexDAControlParameter()
    da_ctrl.method = trexselector.DAMethod.AR1
    da_ctrl.rho_thr_DA = 1.0
    with pytest.raises(ValueError, match="rho_thr_DA"):
        trexselector.TRexDASelector(X, y, da_control=da_ctrl, verbose=False)

    da_ctrl.method = trexselector.DAMethod.EQUI
    trexselector.TRexDASelector(X, y, da_control=da_ctrl, verbose=False)


def test_da_prior_groups_rejects_negative_labels(signal_data):
    """Negative prior-group labels are rejected at validation."""
    X, y, n, p = signal_data
    da_ctrl = trexselector.TRexDAControlParameter()
    da_ctrl.method = trexselector.DAMethod.PRIOR_GROUPS
    da_ctrl.prior_groups = [[-1] + [i % 4 for i in range(p - 1)]]
    with pytest.raises(ValueError, match="non-negative"):
        trexselector.TRexDASelector(X, y, da_control=da_ctrl,
                                         verbose=False).select()


def test_da_prior_groups_method_needs_groups(signal_data):
    """PRIOR_GROUPS without prior_groups is rejected at construction."""
    X, y, n, p = signal_data
    da_ctrl = trexselector.TRexDAControlParameter()
    da_ctrl.method = trexselector.DAMethod.PRIOR_GROUPS
    with pytest.raises(ValueError, match="requires prior_groups"):
        trexselector.TRexDASelector(X, y, da_control=da_ctrl, verbose=False)


@pytest.mark.parametrize("method", [trexselector.DAMethod.BT, trexselector.DAMethod.NN])
def test_da_grid_methods_report_the_reference_threshold(method):
    """The matrices of the selector belong to the reference threshold of the grid."""
    rng = np.random.default_rng(42)
    n, p = 100, 40
    X = np.asfortranarray(rng.standard_normal((n, p)))
    y = 5.0 * X[:, 0] - 4.0 * X[:, 1] + 3.0 * X[:, 2] + rng.standard_normal(n)

    da_ctrl = trexselector.TRexDAControlParameter()
    da_ctrl.method = method
    ctrl = trexselector.TRexControlParameter()
    ctrl.K = 5
    sel = trexselector.TRexDASelector(X, y, tFDR=0.2, da_control=da_ctrl, trex_control=ctrl, seed=42,
                                      verbose=False)
    res = sel.select()

    # Grid of 20 thresholds, reference element 15 (index 14)
    assert res.rho_grid.size == 20
    assert np.array_equal(sel.fdp_hat_mat, res.fdp_hat_array_bt[14])
    assert np.array_equal(sel.phi_mat, res.phi_array_bt[14])

    # The result holds the arrays, its two base matrices stay empty
    assert res.fdp_hat_mat.size == 0
    assert res.phi_mat.size == 0

    # One row per stopping time, one less for the selection counts
    assert res.fdp_hat_array_bt[14].shape[0] == sel.T_stop
    assert res.r_array_bt[14].shape[0] == max(sel.T_stop - 1, 1)

    # The selection counts of the selector belong to the selected threshold
    if res.rho_thresh >= 0:
        grid = res.rho_grid.ravel()
        assert any(np.array_equal(sel.r_mat, res.r_array_bt[r])
                   for r in range(grid.size) if grid[r] == res.rho_thresh)


@pytest.mark.parametrize("method", [trexselector.DAMethod.AR1, trexselector.DAMethod.EQUI])
def test_da_window_methods_report_no_grid(signal_data, fast_control, method):
    """AR1 and EQUI report -2 as rho_thresh, an empty grid and empty arrays."""
    X, y, n, p = signal_data
    da_ctrl = trexselector.TRexDAControlParameter()
    da_ctrl.method = method
    sel = trexselector.TRexDASelector(X, y, da_control=da_ctrl, trex_control=fast_control, verbose=False)
    res = sel.select()
    assert res.method == method
    assert res.rho_thresh == -2.0
    assert res.rho_grid.size == 0
    assert len(res.fdp_hat_array_bt) == 0


def test_da_prior_groups_set_the_method(signal_data, fast_control):
    """Supplied prior groups make PRIOR_GROUPS the method of the run."""
    X, y, n, p = signal_data
    da_ctrl = trexselector.TRexDAControlParameter()
    assert da_ctrl.method == trexselector.DAMethod.BT
    da_ctrl.prior_groups = [[i % 4 for i in range(p)]]
    sel = trexselector.TRexDASelector(X, y, da_control=da_ctrl, trex_control=fast_control, verbose=False)
    res = sel.select()
    assert res.method == trexselector.DAMethod.PRIOR_GROUPS

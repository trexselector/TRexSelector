"""
Tests for TRexTikhonovSelector (general-Tikhonov T-Rex).

Mirrors the R suite (test-trex_tikhonov.R): control round trip, construction
gates, the solver-family warning, a fixed-lambda_2 run with its getters, the
group-mean-K equivalence with GVS-IEN TCDIENET, the K = I collapse onto the
TCDENET run, the two K builders, and the lambda_2 tuners (CV and pilot).
"""
import warnings

import numpy as np
import pytest
import trexselector as tsn


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _grouped_data(seed=7, n=60, p=12, gsize=4, rho=0.7):
    """Correlated groups of size gsize; the signal sits in the first group."""
    rng = np.random.default_rng(seed)
    X = np.empty((n, p))
    for g in range(0, p, gsize):
        z = rng.standard_normal(n)
        for j in range(g, min(g + gsize, p)):
            X[:, j] = rho * z + np.sqrt(1.0 - rho * rho) * rng.standard_normal(n)
    y = 3.0 * X[:, 0] + 2.0 * X[:, 1] - 2.0 * X[:, 2] + rng.standard_normal(n)
    groups = [j // gsize for j in range(p)]
    return X, y, groups


def _group_mean_k(groups):
    """W = sum_m 1_m 1_m^T / p_m (the IEN group-mean matrix)."""
    groups = np.asarray(groups)
    p = groups.size
    K = np.zeros((p, p))
    for g in np.unique(groups):
        idx = np.where(groups == g)[0]
        K[np.ix_(idx, idx)] = 1.0 / idx.size
    return K


def _tik_control(K=None, lambda_2=1.0, **fields):
    ctrl = tsn.TRexTikhonovControlParameter()
    if K is not None:
        ctrl.tikhonov_K = K
    ctrl.lambda_2 = lambda_2
    for name, value in fields.items():
        setattr(ctrl, name, value)
    return ctrl


def _base_control(K=5, mult=2, solver=None):
    ctrl = tsn.TRexControlParameter()
    ctrl.K = K
    ctrl.max_dummy_multiplier = mult
    # A base control that is passed has to name the solver of the selector.
    ctrl.solver_type = solver if solver is not None else tsn.SolverTypeForTRex.TIKNET
    return ctrl


# ---------------------------------------------------------------------------
# Control parameter round trip
# ---------------------------------------------------------------------------

def test_control_fields_roundtrip():
    ctrl = tsn.TRexTikhonovControlParameter()
    # defaults mirror the C++ struct
    assert ctrl.penalty_structure == tsn.TikPenaltyStructure.SPARSE_GENERAL
    assert ctrl.lambda_2 == pytest.approx(-1.0)
    assert ctrl.lambda2_method == tsn.TikLambda2Method.PILOT_KNEE
    assert not hasattr(tsn.TikLambda2Method, "PILOT_MAX")   # the argmax variant is gone
    assert ctrl.cluster_dummies is True
    assert ctrl.fold_dummy_coupling is True
    assert ctrl.trex_ctrl.solver_type == tsn.SolverTypeForTRex.TIKNET   # the homotopy solver
    assert ctrl.tikhonov_K.shape == (0, 0)
    assert ctrl.cv_dfmax_ratio == pytest.approx(0.8)

    K = _group_mean_k([0, 0, 1, 1])
    ctrl.tikhonov_K = K
    np.testing.assert_allclose(ctrl.tikhonov_K.toarray(), K)   # a dense array is held sparse

    ctrl.penalty_structure = tsn.TikPenaltyStructure.GROUP_CENTERING
    ctrl.lambda2_method = tsn.TikLambda2Method.PILOT_KNEE
    ctrl.pilot_statistic = tsn.PilotStatistic.GROUP_MAJORITY
    ctrl.cv_n_folds = 5
    ctrl.cv_n_lambda1 = 15
    ctrl.cv_n_lambda2 = 25
    ctrl.cv_seed = 3
    ctrl.cv_dfmax_ratio = 0.95
    ctrl.pilot_grid = [0.0, 1.0, 10.0]
    ctrl.pilot_K = 4
    ctrl.pilot_T = 2
    ctrl.pilot_L_multiplier = 1
    ctrl.pilot_knee_tolerance = 2.0
    ctrl.pilot_group_fraction = 0.6
    ctrl.pilot_shared_dummies = True
    ctrl.fold_dummy_coupling = False
    ctrl.cluster_dummies = False
    ctrl.prior_groups = [0, 0, 1, 1]
    ctrl.corr_max = 0.3
    assert ctrl.penalty_structure == tsn.TikPenaltyStructure.GROUP_CENTERING
    assert ctrl.lambda2_method == tsn.TikLambda2Method.PILOT_KNEE
    assert ctrl.pilot_statistic == tsn.PilotStatistic.GROUP_MAJORITY
    assert (ctrl.cv_n_folds, ctrl.cv_n_lambda1, ctrl.cv_n_lambda2, ctrl.cv_seed) == (5, 15, 25, 3)
    assert ctrl.cv_dfmax_ratio == pytest.approx(0.95)
    assert list(ctrl.pilot_grid) == [0.0, 1.0, 10.0]
    assert (ctrl.pilot_K, ctrl.pilot_T, ctrl.pilot_L_multiplier) == (4, 2, 1)
    assert ctrl.pilot_knee_tolerance == pytest.approx(2.0)
    assert ctrl.pilot_group_fraction == pytest.approx(0.6)
    assert ctrl.pilot_shared_dummies is True
    assert ctrl.fold_dummy_coupling is False
    assert ctrl.cluster_dummies is False
    assert list(ctrl.prior_groups) == [0, 0, 1, 1]
    assert ctrl.corr_max == pytest.approx(0.3)


# ---------------------------------------------------------------------------
# Construction gates
# ---------------------------------------------------------------------------

def test_tik_control_is_required():
    X, y, _ = _grouped_data()
    with pytest.raises(ValueError, match="tik_control is required"):
        tsn.TRexTikhonovSelector(X, y, trex_control=_base_control(), verbose=False)


def test_empty_or_misdimensioned_k_raises():
    X, y, groups = _grouped_data()
    # SPARSE_GENERAL without a K
    with pytest.raises(Exception):
        tsn.TRexTikhonovSelector(X, y, tik_control=_tik_control(None, 1.0),
                                 trex_control=_base_control(), verbose=False)
    # K of the wrong size
    with pytest.raises(Exception):
        tsn.TRexTikhonovSelector(X, y, tik_control=_tik_control(np.eye(X.shape[1] - 1), 1.0),
                                 trex_control=_base_control(), verbose=False)


@pytest.mark.parametrize("solver", ["TLARS", "TLASSO", "TENET", "TOMP", "TIENET"])
def test_a_solver_outside_the_family_is_refused(solver):
    """The root class, the core class and the power rule refuse it, as C++ does."""
    from trexselector.trex_selector_methods import ad as core

    X, y, groups = _grouped_data()
    X = np.asfortranarray(X)
    X0 = X.copy(order="F")
    tik = _tik_control(_group_mean_k(groups), 1.0, prior_groups=groups)
    base = _base_control(solver=getattr(tsn.SolverTypeForTRex, solver))
    with pytest.raises(ValueError, match="solver_type must be TIKNET"):
        tsn.TRexTikhonovSelector(X, y, tik_control=tik, trex_control=base, seed=1, verbose=False)
    with pytest.raises(ValueError, match="solver_type must be TIKNET"):
        core.TRexTikhonovSelector(X, y, 0.1, tik, base, 1, False)
    with pytest.raises(ValueError, match="solver_type must be TIKNET"):
        tsn.tikhonov_lambda2_power_rule(X, y, tik_control=tik, trex_control=base, seed=1, grid=[0.0, 1.0])
    np.testing.assert_allclose(X, X0, atol=1e-12)     # a refused control leaves X alone


def test_the_coordinate_descent_solver_is_refused():
    """TCDTIKNET, which the selector accepted up to 1.0.0, is an error on every entry point, without a
    deprecation warning before it."""
    from trexselector.trex_selector_methods import ad as core

    X, y, groups = _grouped_data()
    X = np.asfortranarray(X)
    X0 = X.copy(order="F")
    tik = _tik_control(_group_mean_k(groups), 1.0, prior_groups=groups)
    base = _base_control(solver=tsn.SolverTypeForTRex.TCDTIKNET)
    refused = "TCDTIKNET is no longer accepted by the selector"
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(ValueError, match=refused):
            tsn.TRexTikhonovSelector(X, y, tik_control=tik, trex_control=base, seed=1, verbose=False)
        with pytest.raises(ValueError, match=refused):
            core.TRexTikhonovSelector(X, y, 0.1, tik, base, 1, False)
        with pytest.raises(ValueError, match=refused):
            tsn.tikhonov_lambda2_power_rule(X, y, tik_control=tik, trex_control=base, seed=1, grid=[0.0, 1.0])
        # Also when it is named on the nested control
        tik.trex_ctrl.solver_type = tsn.SolverTypeForTRex.TCDTIKNET
        with pytest.raises(ValueError, match=refused):
            tsn.TRexTikhonovSelector(X, y, tik_control=tik, seed=1, verbose=False)
    np.testing.assert_allclose(X, X0, atol=1e-12)     # a refused control leaves X alone


def test_the_decoupled_dummy_block_is_refused():
    """fold_dummy_coupling = False existed with TCDTIKNET only, so it has no solver left."""
    X, y, _ = _grouped_data(seed=5)
    tik = _tik_control(np.eye(X.shape[1]), 1.5, fold_dummy_coupling=False)
    with pytest.raises(ValueError, match="fold_dummy_coupling = false is no longer available"):
        tsn.TRexTikhonovSelector(X.copy(), y.copy(), tik_control=tik, seed=42, verbose=False)
    with pytest.raises(ValueError, match="fold_dummy_coupling = false is no longer available"):
        tsn.TRexTikhonovSelector(X.copy(), y.copy(), tik_control=tik, trex_control=_base_control(),
                                 seed=42, verbose=False)


def test_the_nested_base_control_applies_without_trex_control():
    """Without trex_control the base control of the Tikhonov control is the one of the run, as in C++."""
    X, y, groups = _grouped_data()
    tik = _tik_control(_group_mean_k(groups), 1.0, prior_groups=groups)
    assert tik.trex_ctrl.solver_type == tsn.SolverTypeForTRex.TIKNET
    tik.trex_ctrl.K = 7
    tik.trex_ctrl.max_dummy_multiplier = 2

    nested = tsn.TRexTikhonovSelector(np.asfortranarray(X.copy()), y, tik_control=tik, seed=1, verbose=False)
    assert nested.K == 7
    selected = np.array(nested.select().selected_var)

    # The same control passed as trex_control gives the same run and replaces the nested one.
    tik.trex_ctrl.K = 3
    passed = tsn.TRexTikhonovSelector(np.asfortranarray(X.copy()), y, tik_control=tik,
                                      trex_control=_base_control(K=7, mult=2), seed=1, verbose=False)
    assert passed.K == 7
    assert np.array_equal(np.array(passed.select().selected_var), selected)
    # Neither control is changed by the selector.
    assert tik.trex_ctrl.K == 3


def test_the_default_of_the_base_control_is_refused():
    """TRexControlParameter() names TLARS, so a base control that is passed has to set the solver."""
    X, y, groups = _grouped_data()
    tik = _tik_control(_group_mean_k(groups), 1.0, prior_groups=groups)
    base = tsn.TRexControlParameter()
    assert base.solver_type == tsn.SolverTypeForTRex.TLARS
    with pytest.raises(ValueError, match="solver_type must be TIKNET"):
        tsn.TRexTikhonovSelector(np.asfortranarray(X.copy()), y, tik_control=tik, trex_control=base,
                                 verbose=False)
    base.solver_type = tsn.SolverTypeForTRex.TIKNET
    tsn.TRexTikhonovSelector(np.asfortranarray(X.copy()), y, tik_control=tik, trex_control=base,
                             verbose=False).restore()


# ---------------------------------------------------------------------------
# A fixed-lambda_2 run and the getters
# ---------------------------------------------------------------------------

def test_fixed_lambda2_run_and_getters():
    X, y, groups = _grouped_data()
    tik = _tik_control(_group_mean_k(groups), 2.5, prior_groups=groups)
    sel = tsn.TRexTikhonovSelector(X, y, tik_control=tik,
                                   trex_control=_base_control(), seed=1, verbose=False)
    assert sel.lambda2_used == pytest.approx(-1.0)          # before select()
    res = sel.select()
    assert isinstance(res, tsn.SelectionResult)
    assert sel.lambda2_used == pytest.approx(2.5)
    assert sel.lambda2_cv_diagnostics.ran is False
    assert sel.lambda2_pilot.ran is False
    assert all(0 <= i < X.shape[1] for i in sel.selected_indices)
    assert sel.T_stop >= 1
    assert 1 <= sel.L <= 2


def test_group_centering_structure_runs():
    X, y, groups = _grouped_data()
    tik = _tik_control(None, 1.0, prior_groups=groups,
                       penalty_structure=tsn.TikPenaltyStructure.GROUP_CENTERING)
    sel = tsn.TRexTikhonovSelector(X, y, tik_control=tik,
                                   trex_control=_base_control(), seed=1, verbose=False)
    sel.select()
    assert sel.selected_var.shape == (X.shape[1],)


# ---------------------------------------------------------------------------
# K builders
# ---------------------------------------------------------------------------

def test_gamma_to_k_builds_the_chain_laplacian():
    p = 6
    gamma = np.zeros((p - 1, p))
    for i in range(p - 1):
        gamma[i, i] = -1.0
        gamma[i, i + 1] = 1.0
    K = tsn.gamma_to_k(gamma)
    expected = gamma.T @ gamma
    np.testing.assert_allclose(K.toarray(), expected)
    assert K[0, 0] == pytest.approx(1.0) and K[1, 1] == pytest.approx(2.0)
    assert K[0, 1] == pytest.approx(-1.0)


def test_build_group_centering_k():
    K = tsn.build_group_centering_k([0, 0, 1, 1, 1])
    expected = np.zeros((5, 5))
    expected[:2, :2] = np.eye(2) - np.ones((2, 2)) / 2.0
    expected[2:, 2:] = np.eye(3) - np.ones((3, 3)) / 3.0
    np.testing.assert_allclose(K.toarray(), expected)


# ---------------------------------------------------------------------------
# lambda_2 resolution: CV tuner and pilot rule
# ---------------------------------------------------------------------------

def test_lambda2_resolves_via_the_ien_tuner():
    X, y, groups = _grouped_data(seed=3)
    tik = _tik_control(_group_mean_k(groups), -1.0, prior_groups=groups,
                       lambda2_method=tsn.TikLambda2Method.CV_1SE_IEN_CCD,
                       cv_n_folds=5, cv_n_lambda1=15, cv_n_lambda2=12, cv_seed=1,
                       cv_dfmax_ratio=0.9)
    sel = tsn.TRexTikhonovSelector(X, y, tik_control=tik,
                                   trex_control=_base_control(), seed=1, verbose=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")   # a boundary pick may warn on tiny data
        sel.select()
    diag = sel.lambda2_cv_diagnostics
    assert diag.ran is True
    assert diag.grid_size == 12
    assert 0 <= diag.grid_index < 12
    assert sel.lambda2_used >= 0.0
    assert sel.lambda2_pilot.ran is False
    # The support-cap ratio reached the tuner and is recorded as passed.
    assert diag.dfmax_ratio == pytest.approx(0.9)
    assert diag.profile_at_edge or not diag.profile_at_support_cap


def test_non_finite_cv_dfmax_ratio_raises():
    X, y, groups = _grouped_data()
    tik = _tik_control(_group_mean_k(groups), -1.0, prior_groups=groups,
                       cv_dfmax_ratio=float("nan"))
    with pytest.raises(ValueError, match="cv_dfmax_ratio must be finite"):
        tsn.TRexTikhonovSelector(X, y, tik_control=tik,
                                 trex_control=_base_control(), seed=1, verbose=False)


def test_lambda2_resolves_via_the_pilot_rule():
    X, y, groups = _grouped_data(seed=3)
    grid = [0.0, 1.0, 10.0, 100.0]
    tik = _tik_control(_group_mean_k(groups), -1.0, prior_groups=groups,
                       lambda2_method=tsn.TikLambda2Method.PILOT_KNEE,
                       pilot_grid=grid, pilot_K=3, pilot_T=2, pilot_L_multiplier=1)
    sel = tsn.TRexTikhonovSelector(X, y, tik_control=tik,
                                   trex_control=_base_control(), seed=1, verbose=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sel.select()
    pilot = sel.lambda2_pilot
    assert pilot.ran is True
    assert list(pilot.grid) == grid
    assert len(pilot.counts) == len(grid)
    assert pilot.lambda2 in grid
    assert pilot.index == grid.index(pilot.lambda2)
    assert pilot.pick.max_selected >= 0
    # L_0 is the pilot's dummy count: pilot_L_multiplier * p
    assert (pilot.K_pilot, pilot.T_pilot, pilot.L_0) == (3, 2, X.shape[1])
    assert sel.lambda2_used == pytest.approx(pilot.lambda2)
    # The CV record notes that a tuner ran (and its seconds); the grid lives
    # in the pilot record, so the CV grid fields stay unset.
    diag = sel.lambda2_cv_diagnostics
    assert diag.ran is True
    assert diag.grid_size == 0 and diag.grid_index == -1


# ---------------------------------------------------------------------------
# The full-profile knee (one selection per grid value)
# ---------------------------------------------------------------------------

def test_lambda2_power_rule_profiles_the_grid_and_keeps_the_contract():
    X, y, groups = _grouped_data(seed=3)
    X_before = X.copy()
    tik = _tik_control(_group_mean_k(groups), -1.0, prior_groups=groups)
    base = _base_control()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        prof = tsn.tikhonov_lambda2_power_rule(X, y, tFDR=0.2, tik_control=tik, trex_control=base,
                                               seed=1, grid=[10.0, 0.0, 1.0, 10.0])
    grid = [p.lambda2 for p in prof.probes]
    assert grid == [0.0, 1.0, 10.0]                       # sorted, de-duplicated
    counts = list(prof.counts())
    assert counts == [len(p.selected) for p in prof.probes]
    assert prof.index == next(k for k, c in enumerate(counts) if c >= max(counts) - 1)   # the knee
    assert prof.lambda2 == grid[prof.index]
    assert list(prof.selection()) == list(prof.probes[prof.index].selected)
    assert prof.pick.max_selected == max(counts)
    assert prof.seconds > 0
    np.testing.assert_array_equal(X, X_before)            # every probe hands X back

    # Tolerance 0 reads the largest selection, ties toward the smallest lambda_2.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        argmax = tsn.tikhonov_lambda2_power_rule(X, y, tFDR=0.2, tik_control=tik, trex_control=base,
                                                 seed=1, grid=[0.0, 1.0, 10.0], knee_tolerance=0.0)
    assert argmax.index == counts.index(max(counts))

    # The contract: a fixed-lambda_2 run at the pick with the same seed selects the same variables.
    tik_fixed = _tik_control(_group_mean_k(groups), prof.lambda2, prior_groups=groups)
    sel = tsn.TRexTikhonovSelector(X, y, tFDR=0.2, tik_control=tik_fixed, trex_control=base,
                                   seed=1, verbose=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sel.select()
    assert sorted(sel.selected_indices) == list(prof.selection())

    with pytest.raises(ValueError, match="seed"):
        tsn.tikhonov_lambda2_power_rule(X, y, tik_control=tik, trex_control=base, seed=-1)
    with pytest.raises(ValueError):
        tsn.tikhonov_lambda2_power_rule(X, y, tik_control=tik, trex_control=base, knee_tolerance=-1.0)
    with pytest.raises(ValueError, match="tik_control"):
        tsn.tikhonov_lambda2_power_rule(X, y)


# ---------------------------------------------------------------------------
# The default run on a chain Laplacian
# ---------------------------------------------------------------------------

def test_default_run_on_a_chain_laplacian_selects_the_planted_block():
    """The default control (homotopy solver, pilot rule over its default grid up to 3000).

    The coordinate-descent solver stops unconverged at the large grid values, the pilot then
    takes its count there for the knee, picks 3000 and the selection is dominated by nulls.
    """
    rng = np.random.default_rng(2718)
    n, p = 120, 60
    X = np.asfortranarray(rng.standard_normal((n, p)))
    support = np.arange(20, 30)
    y = X[:, support].sum(axis=1) + rng.standard_normal(n)
    gamma = np.diff(np.eye(p), axis=0)

    tik = tsn.TRexTikhonovControlParameter()   # lambda_2 < 0: the pilot rule
    tik.tikhonov_K = gamma.T @ gamma
    sel = tsn.TRexTikhonovSelector(X, y, tFDR=0.1, tik_control=tik, trex_control=_base_control(K=10, mult=10),
                                   seed=4, verbose=False)
    sel.select()

    selected = np.flatnonzero(sel.selected_var)
    assert np.isin(support, selected).sum() >= 8
    assert np.setdiff1d(selected, support).size <= 4
    assert sel.lambda2_used < 1000.0   # not the top of the pilot grid


# ---------------------------------------------------------------------------
# Checks at construction
# ---------------------------------------------------------------------------

def _small():
    rng = np.random.default_rng(31)
    n, p = 40, 10
    X = np.asfortranarray(rng.standard_normal((n, p)))
    y = X[:, 0] + rng.standard_normal(n)
    return X, y, n, p


def _build(X, y, **fields):
    tik = tsn.TRexTikhonovControlParameter()
    for name, value in fields.items():
        setattr(tik, name, value)
    return tsn.TRexTikhonovSelector(X, y, tFDR=0.3, tik_control=tik, trex_control=_base_control(K=5, mult=10),
                                    seed=2, verbose=False)


def test_the_content_of_k_is_checked_at_construction():
    X, y, _, p = _small()
    X0 = X.copy(order="F")
    K_asym = np.eye(p)
    K_asym[0, 1] = 0.5
    K_neg = np.eye(p)
    K_neg[p - 1, p - 1] = -1.0
    with pytest.raises(ValueError, match="TRexTikhonovSelector.*not symmetric"):
        _build(X, y, tikhonov_K=K_asym, lambda_2=1.0)
    with pytest.raises(ValueError, match="TRexTikhonovSelector.*positive semi-definite"):
        _build(X, y, tikhonov_K=K_neg, lambda_2=1.0)
    np.testing.assert_allclose(X, X0, atol=1e-12)     # a refused control leaves X alone


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_lambda_2_must_be_finite(bad):
    X, y, _, p = _small()
    with pytest.raises(ValueError, match="lambda_2 must be finite"):
        _build(X, y, tikhonov_K=np.eye(p), lambda_2=bad)


def test_the_folds_of_a_cv_rule_are_checked_at_construction():
    X, y, n, p = _small()
    for folds in (1, n + 1):
        with pytest.raises(ValueError, match=r"cv_n_folds must be in \[2, n\]"):
            _build(X, y, tikhonov_K=np.eye(p), lambda2_method=tsn.TikLambda2Method.CV_1SE_IEN_CCD,
                   cv_n_folds=folds)
    # A fixed lambda_2 does not use the folds.
    _build(X, y, tikhonov_K=np.eye(p), lambda_2=1.0, cv_n_folds=1).restore()


def test_group_centering_refuses_another_matrix():
    X, y, _, p = _small()
    groups = [0] * 5 + [1] * 5

    def run(K):
        fields = dict(penalty_structure=tsn.TikPenaltyStructure.GROUP_CENTERING, prior_groups=groups,
                      lambda_2=1.0)
        if K is not None:
            fields["tikhonov_K"] = K
        return np.array(_build(X, y, **fields).select().selected_var)

    with pytest.raises(ValueError, match="must be the centering matrix"):
        run(np.eye(p))
    # The right matrix and no matrix give the same run.
    assert np.array_equal(run(tsn.build_group_centering_k(groups)), run(None))


def test_the_power_rule_and_the_pilot_rule_have_one_default_grid():
    import inspect

    grid = list(tsn.TRexTikhonovControlParameter().pilot_grid)
    assert grid == [0.0, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0]
    # The grid stands once, in C++. The root functions take it from there when none is passed.
    assert list(tsn.default_lambda2_grid()) == grid
    for rule in (tsn.tikhonov_lambda2_power_rule, tsn.gvs_lambda2_power_rule):
        assert inspect.signature(rule).parameters["grid"].default is None


def test_the_core_power_rule_refuses_a_negative_seed():
    """The probes share their dummies through the seed, so the core refuses a negative one."""
    from trexselector.trex_selector_methods import ad as core

    X, y, _, p = _small()
    tik = tsn.TRexTikhonovControlParameter()
    tik.tikhonov_K = np.eye(p)
    base = _base_control(K=5, mult=10)
    with pytest.raises(ValueError, match="seed must be >= 0"):
        core.tikhonov_lambda2_power_rule(X, y, 0.3, tik, base, -1, [0.0, 1.0], 1.0, False)
    # Without a grid the core takes the default grid of the pilot rule.
    res = core.tikhonov_lambda2_power_rule(X, y, 0.3, tik, base, 3)
    assert [pr.lambda2 for pr in res.probes] == list(tik.pilot_grid)


def test_a_lambda_2_grid_means_the_same_shrinkage_under_zscore():
    """Under ZSCORE the pilot and the power rule run at the grid times n - 1."""
    X, y, groups = _grouped_data()
    n = X.shape[0]
    grid = [0.0, 1.0, 10.0]

    def controls(mode):
        tik = tsn.TRexTikhonovControlParameter()
        tik.penalty_structure = tsn.TikPenaltyStructure.GROUP_CENTERING
        tik.prior_groups = list(groups)
        tik.pilot_grid = grid
        tik.pilot_K = 5
        base = _base_control(K=10, mult=10)
        base.scaling_mode = mode
        return tik, base

    def pilot(mode):
        tik, base = controls(mode)
        sel = tsn.TRexTikhonovSelector(np.asfortranarray(X.copy()), y, tFDR=0.2, tik_control=tik,
                                       trex_control=base, seed=5, verbose=False)
        selected = np.array(sel.select().selected_var)
        return sel, selected

    l2, sel_l2 = pilot(tsn.ScalingMode.L2)
    z, sel_z = pilot(tsn.ScalingMode.ZSCORE)
    assert l2.lambda2_pilot.grid_scale == 1.0
    assert z.lambda2_pilot.grid_scale == n - 1
    assert list(l2.lambda2_pilot.grid) == grid
    assert np.allclose(z.lambda2_pilot.grid, np.array(grid) * (n - 1))
    assert list(z.lambda2_pilot.counts) == list(l2.lambda2_pilot.counts)
    assert z.lambda2_pilot.index == l2.lambda2_pilot.index
    assert np.isclose(z.lambda2_used, l2.lambda2_used * (n - 1))
    assert np.array_equal(sel_z, sel_l2)

    def rule(mode):
        tik, base = controls(mode)
        return tsn.tikhonov_lambda2_power_rule(np.asfortranarray(X.copy()), y, tFDR=0.2, tik_control=tik,
                                               trex_control=base, seed=5, grid=grid)

    a = rule(tsn.ScalingMode.L2)
    b = rule(tsn.ScalingMode.ZSCORE)
    assert a.grid_scale == 1.0
    assert b.grid_scale == n - 1
    assert np.allclose([pr.lambda2 for pr in b.probes], np.array(grid) * (n - 1))
    assert list(b.counts()) == list(a.counts())
    assert list(b.selection()) == list(a.selection())


# ---------------------------------------------------------------------------
# The penalty matrix is sparse
# ---------------------------------------------------------------------------

def test_the_penalty_matrix_is_held_and_returned_sparse():
    import scipy.sparse as sp

    K = _group_mean_k([0, 0, 1, 1, 1])
    for given in (K, sp.csc_matrix(K), sp.csr_matrix(K), sp.coo_matrix(K)):
        ctrl = tsn.TRexTikhonovControlParameter()
        ctrl.tikhonov_K = given
        back = ctrl.tikhonov_K
        assert isinstance(back, sp.csc_matrix)
        np.testing.assert_allclose(back.toarray(), K)
    assert isinstance(tsn.TRexTikhonovControlParameter().tikhonov_K, sp.csc_matrix)

    gamma = np.diff(np.eye(6), axis=0)
    for given in (gamma, sp.csr_matrix(gamma)):
        built = tsn.gamma_to_k(given)
        assert isinstance(built, sp.csc_matrix)
        np.testing.assert_allclose(built.toarray(), gamma.T @ gamma)
    assert isinstance(tsn.build_group_centering_k([0, 0, 1, 1, 1]), sp.csc_matrix)


def test_a_sparse_and_a_dense_penalty_matrix_give_the_same_run():
    import scipy.sparse as sp

    X, y, groups = _grouped_data()
    K = _group_mean_k(groups)

    def run(matrix):
        sel = tsn.TRexTikhonovSelector(np.asfortranarray(X.copy()), y, tFDR=0.2,
                                       tik_control=_tik_control(matrix, 1.0, prior_groups=groups),
                                       trex_control=_base_control(), seed=3, verbose=False)
        return np.array(sel.select().selected_var)

    assert np.array_equal(run(K), run(sp.csr_matrix(K)))


def test_a_large_penalty_matrix_stays_sparse():
    """The centering matrix of 20,000 variables has a million entries, where a dense array has 3.2 GB."""
    p, size = 20_000, 50
    groups = [j // size for j in range(p)]
    K = tsn.build_group_centering_k(groups)
    assert K.shape == (p, p)
    assert K.nnz == p * size
    ctrl = tsn.TRexTikhonovControlParameter()
    ctrl.tikhonov_K = K
    assert ctrl.tikhonov_K.nnz == p * size
    chain = tsn.gamma_to_k(K[:1000, :1000])
    assert chain.shape == (1000, 1000)


# ---------------------------------------------------------------------------
# What the C++ class has: the pilot alone, the clusters of the dummies, the GVS power rule
# ---------------------------------------------------------------------------

def _pilot_selector(X, y, groups, **fields):
    tik = tsn.TRexTikhonovControlParameter()
    tik.penalty_structure = tsn.TikPenaltyStructure.GROUP_CENTERING
    tik.prior_groups = list(groups)
    tik.pilot_grid = [0.0, 1.0, 10.0]
    tik.pilot_K = 5
    for name, value in fields.items():
        setattr(tik, name, value)
    return tsn.TRexTikhonovSelector(np.asfortranarray(X.copy()), y, tFDR=0.2, tik_control=tik,
                                    trex_control=_base_control(K=10, mult=10), seed=5, verbose=False)


def test_run_pilot_only_gives_the_pilot_of_a_selection():
    X, y, groups = _grouped_data()
    full = _pilot_selector(X, y, groups)
    full.select()

    Xf = np.asfortranarray(X.copy())
    X0 = Xf.copy(order="F")
    tik = tsn.TRexTikhonovControlParameter()
    tik.penalty_structure = tsn.TikPenaltyStructure.GROUP_CENTERING
    tik.prior_groups = list(groups)
    tik.pilot_grid = [0.0, 1.0, 10.0]
    tik.pilot_K = 5
    only = tsn.TRexTikhonovSelector(Xf, y, tFDR=0.2, tik_control=tik, trex_control=_base_control(K=10, mult=10),
                                    seed=5, verbose=False)
    record = only.run_pilot_only()
    assert record.ran
    assert list(record.counts) == list(full.lambda2_pilot.counts)
    assert record.lambda2 == full.lambda2_pilot.lambda2
    assert only.lambda2_used == full.lambda2_used
    assert list(only.lambda2_pilot.counts) == list(record.counts)
    np.testing.assert_allclose(Xf, X0, atol=1e-12)        # X is given back

    # A selector runs once.
    with pytest.raises(RuntimeError, match="has run already"):
        only.select()
    with pytest.raises(RuntimeError, match="has run already"):
        only.run_pilot_only()
    with pytest.raises(RuntimeError, match="has run already"):
        full.run_pilot_only()


def test_run_pilot_only_needs_the_pilot_rule():
    X, y, groups = _grouped_data()
    fixed = _pilot_selector(X, y, groups, lambda_2=1.0)
    with pytest.raises(ValueError, match="PILOT_KNEE"):
        fixed.run_pilot_only()
    fixed.select()                                           # the refused call did not spend the selector


def test_dummy_clusters_are_reported():
    X, y, groups = _grouped_data()
    p = X.shape[1]

    before = _pilot_selector(X, y, groups, lambda_2=1.0)
    record = before.dummy_clusters
    assert len(record.groups_vec) == 0 and record.max_clusters == 0 and record.hc_method_used == ""
    before.restore()

    prior = _pilot_selector(X, y, groups, lambda_2=1.0)
    prior.select()
    assert list(prior.dummy_clusters.groups_vec) == list(groups)
    assert prior.dummy_clusters.max_clusters == len(set(groups))
    assert prior.dummy_clusters.hc_method_used == "prior"

    def found(**fields):
        tik = _tik_control(np.eye(p), 1.0)
        for name, value in fields.items():
            setattr(tik, name, value)
        sel = tsn.TRexTikhonovSelector(np.asfortranarray(X.copy()), y, tFDR=0.2, tik_control=tik,
                                       trex_control=_base_control(), seed=5, verbose=False)
        sel.select()
        return sel.dummy_clusters

    single = found()
    assert single.hc_method_used == "single"
    assert len(single.groups_vec) == p
    assert sorted(set(int(g) for g in single.groups_vec)) == list(range(single.max_clusters))

    none = found(cluster_dummies=False)
    assert len(none.groups_vec) == 0 and none.max_clusters == 0 and none.hc_method_used == ""


def test_the_gvs_power_rule_runs_the_grid():
    X, y, groups = _grouped_data()
    gvs = tsn.TRexGVSControlParameter()
    gvs.prior_groups = list(groups)
    base = tsn.TRexControlParameter()
    base.solver_type = tsn.SolverTypeForTRex.TENET      # a base control names the solver of the track
    base.K = 5
    base.max_dummy_multiplier = 2
    grid = [0.0, 1.0, 10.0]

    def rule(seed=3, **kw):
        return tsn.gvs_lambda2_power_rule(np.asfortranarray(X.copy()), y, tFDR=0.2, gvs_control=gvs,
                                          trex_control=base, seed=seed, grid=grid, **kw)

    a = rule()
    assert [pr.lambda2 for pr in a.probes] == grid
    counts = list(a.counts())
    assert a.pick.max_selected == max(counts)
    assert len(a.selection()) == counts[a.index]
    assert a.lambda2 == grid[a.index]
    assert a.grid_scale == 1.0

    b = rule()                                               # the same seed reproduces the profile
    assert [list(pr.selected) for pr in b.probes] == [list(pr.selected) for pr in a.probes]

    with pytest.raises(ValueError, match="non-negative"):
        rule(seed=-1)
    assert gvs.lambda_2 == -1.0                              # the control of the caller is not changed


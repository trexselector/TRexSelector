"""
Tests for TRexGVSSelector.
"""
import warnings

import numpy as np
import pytest
import trexselector

TENET = trexselector.SolverTypeForTRex.TENET
TIENET = trexselector.SolverTypeForTRex.TIENET


def _base(solver=TENET, K=3):
    """A base control for the GVS selector: a control that is passed names a solver of the track."""
    ctrl = trexselector.TRexControlParameter()
    ctrl.solver_type = solver
    ctrl.K = K
    return ctrl


@pytest.fixture
def gvs_base():
    return _base()

# ---------------------------------------------------------------------------
# Instantiation
# ---------------------------------------------------------------------------

def test_instantiation(signal_data, gvs_base):
    X, y, n, p = signal_data
    sel = trexselector.TRexGVSSelector(X, y, trex_control=gvs_base, verbose=False)
    assert sel is not None


def test_dimension_mismatch_raises(signal_data, gvs_base):
    X, y, n, p = signal_data
    y_bad = np.random.randn(n + 10)
    with pytest.raises(Exception):
        trexselector.TRexGVSSelector(X, y_bad, trex_control=gvs_base, verbose=False)


# ---------------------------------------------------------------------------
# Execution and result fields
# ---------------------------------------------------------------------------

def test_execution_and_result(signal_data, gvs_base):
    X, y, n, p = signal_data
    sel = trexselector.TRexGVSSelector(X, y, trex_control=gvs_base, verbose=False)
    res = sel.select()

    assert res is not None
    assert isinstance(sel.selected_var, np.ndarray)
    assert isinstance(sel.selected_indices, list)

    if sel.selected_indices:
        assert all(0 <= i < p for i in sel.selected_indices)


# ---------------------------------------------------------------------------
# GVS control parameter fields
# ---------------------------------------------------------------------------

def test_gvs_control_new_fields():
    gvs_ctrl = trexselector.TRexGVSControlParameter()

    # corr_max
    gvs_ctrl.corr_max = 0.7
    assert gvs_ctrl.corr_max == pytest.approx(0.7)

    # hc_linkage — just check assignment doesn't raise
    from trexselector._core.ml_methods.clustering import LinkageMethod
    gvs_ctrl.hc_linkage = LinkageMethod.Complete
    assert gvs_ctrl.hc_linkage == LinkageMethod.Complete
    # lambda_2
    gvs_ctrl.lambda_2 = 0.01
    assert gvs_ctrl.lambda_2 == pytest.approx(0.01)


def test_gvs_control_en_solver_and_cv_fields():
    gvs_ctrl = trexselector.TRexGVSControlParameter()

    # en_solver
    gvs_ctrl.en_solver = trexselector.ENSolverType.TENET_AUG
    assert gvs_ctrl.en_solver == trexselector.ENSolverType.TENET_AUG

    # lambda2_method
    gvs_ctrl.lambda2_method = trexselector.LambdaSelectionMethod.CV_MIN_CCD
    assert gvs_ctrl.lambda2_method == trexselector.LambdaSelectionMethod.CV_MIN_CCD

    # CV knobs
    gvs_ctrl.cv_n_folds = 5
    gvs_ctrl.cv_n_lambda = 250
    gvs_ctrl.cv_seed = 7
    assert gvs_ctrl.cv_n_folds == 5
    assert gvs_ctrl.cv_n_lambda == 250
    assert gvs_ctrl.cv_seed == 7

    # tenet_aug_use_lars
    gvs_ctrl.tenet_aug_use_lars = True
    assert gvs_ctrl.tenet_aug_use_lars is True


def test_gvs_en_tenet_aug_runs(signal_data, gvs_base):
    X, y, n, p = signal_data
    gvs_ctrl = trexselector.TRexGVSControlParameter()
    gvs_ctrl.gvs_type = trexselector.GVSType.EN
    gvs_ctrl.en_solver = trexselector.ENSolverType.TENET_AUG
    sel = trexselector.TRexGVSSelector(X, y, gvs_control=gvs_ctrl,
                                        trex_control=gvs_base, verbose=False)
    res = sel.select()
    assert res is not None
    assert res.gvs_type == trexselector.GVSType.EN


def test_gvs_prior_groups(signal_data, gvs_base):
    X, y, n, p = signal_data

    gvs_ctrl = trexselector.TRexGVSControlParameter()
    gvs_ctrl.prior_groups = [i % 4 for i in range(p)]  # 5 groups of 4

    sel = trexselector.TRexGVSSelector(X, y, gvs_control=gvs_ctrl,
                                        trex_control=gvs_base, verbose=False)
    res = sel.select()
    assert res is not None


def test_gvs_result_fields(signal_data, gvs_base):
    X, y, n, p = signal_data
    sel = trexselector.TRexGVSSelector(X, y, trex_control=gvs_base, verbose=False)
    gvs_res = sel.select()
    assert gvs_res is not None
    assert isinstance(gvs_res.gvs_type, trexselector.GVSType)


# ---------------------------------------------------------------------------
# Properties inherited from TRexSelector base
# ---------------------------------------------------------------------------

def test_base_properties(signal_data, gvs_base):
    X, y, n, p = signal_data
    sel = trexselector.TRexGVSSelector(X, y, trex_control=gvs_base, verbose=False)
    sel.select()

    assert isinstance(sel.tFDR, float)
    assert isinstance(sel.T_stop, int)
    assert isinstance(sel.L, int)
    assert isinstance(sel.fdp_hat_mat, np.ndarray)


# ---------------------------------------------------------------------------
# GVSType variants
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("gvs_type, solver", [
    (trexselector.GVSType.EN, TENET),
    (trexselector.GVSType.IEN, TIENET),
])
def test_gvs_type_variants(signal_data, gvs_type, solver):
    X, y, n, p = signal_data
    gvs_ctrl = trexselector.TRexGVSControlParameter()
    gvs_ctrl.gvs_type = gvs_type
    sel = trexselector.TRexGVSSelector(X, y, gvs_control=gvs_ctrl,
                                         trex_control=_base(solver), verbose=False)
    res = sel.select()
    assert res is not None
    assert res.gvs_type == gvs_type


def test_gvs_lambda2_used(signal_data, gvs_base):
    X, y, n, p = signal_data
    sel = trexselector.TRexGVSSelector(X, y, trex_control=gvs_base, verbose=False)
    res = sel.select()
    assert isinstance(res.lambda2_used, float)
    assert res.lambda2_used > 0.0


# ---------------------------------------------------------------------------
# The solver of the base control
# ---------------------------------------------------------------------------

def _gvs(gvs_type, en_solver=None):
    gvs_ctrl = trexselector.TRexGVSControlParameter()
    gvs_ctrl.gvs_type = gvs_type
    gvs_ctrl.lambda_2 = 1.0
    if en_solver is not None:
        gvs_ctrl.en_solver = en_solver
    return gvs_ctrl


@pytest.mark.parametrize("solver", ["TLARS", "TOMP", "TLASSO", "TIENET", "TIKNET"])
def test_gvs_en_refuses_a_solver_outside_the_elastic_net(signal_data, solver):
    """EN: C++ refuses every solver outside the elastic net, TLARS among them. The binding changes nothing."""
    X, y, n, p = signal_data
    X0 = X.copy(order="F")
    ctrl = _base(getattr(trexselector.SolverTypeForTRex, solver))
    with pytest.raises(ValueError, match="gvs_type = EN drives the ENET solver family only"):
        trexselector.TRexGVSSelector(X, y, gvs_control=_gvs(trexselector.GVSType.EN), trex_control=ctrl,
                                     verbose=False)
    with pytest.raises(ValueError, match="gvs_type = EN drives the ENET solver family only"):
        trexselector.gvs_lambda2_power_rule(X, y, gvs_control=_gvs(trexselector.GVSType.EN), trex_control=ctrl,
                                            seed=1, grid=[0.0, 1.0])
    assert ctrl.solver_type == getattr(trexselector.SolverTypeForTRex, solver)   # the control is not changed
    np.testing.assert_allclose(X, X0, atol=1e-12)                                # a refused control leaves X alone


def test_gvs_en_takes_its_solver_from_en_solver(signal_data):
    """EN: en_solver chooses the solver. The base control keeps TENET or names the same solver."""
    X, y, n, p = signal_data
    en = trexselector.GVSType.EN
    aug = trexselector.ENSolverType.TENET_AUG
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        trexselector.TRexGVSSelector(X.copy(order="F"), y, gvs_control=_gvs(en, aug), trex_control=_base(),
                                     verbose=False)
        trexselector.TRexGVSSelector(X.copy(order="F"), y, gvs_control=_gvs(en, aug),
                                     trex_control=_base(trexselector.SolverTypeForTRex.TENET_AUG), verbose=False)
    # A base control that names another elastic-net solver than en_solver is refused
    with pytest.raises(ValueError, match="selects its solver via en_solver"):
        trexselector.TRexGVSSelector(X, y, gvs_control=_gvs(en),
                                     trex_control=_base(trexselector.SolverTypeForTRex.TENET_AUG), verbose=False)


@pytest.mark.parametrize("solver", ["TLARS", "TOMP", "TENET", "TENET_AUG", "TIKNET"])
def test_gvs_ien_refuses_a_solver_outside_its_family(signal_data, solver):
    """IEN: the base control has to name TIENET or TIENET_AUG, as in C++."""
    X, y, n, p = signal_data
    ctrl = _base(getattr(trexselector.SolverTypeForTRex, solver))
    with pytest.raises(ValueError, match="gvs_type = IEN requires solver_type"):
        trexselector.TRexGVSSelector(X, y, gvs_control=_gvs(trexselector.GVSType.IEN), trex_control=ctrl,
                                     verbose=False)
    with pytest.raises(ValueError, match="gvs_type = IEN requires solver_type"):
        trexselector.gvs_lambda2_power_rule(X, y, gvs_control=_gvs(trexselector.GVSType.IEN), trex_control=ctrl,
                                            seed=1, grid=[0.0, 1.0])
    assert ctrl.solver_type == getattr(trexselector.SolverTypeForTRex, solver)


def test_gvs_ien_needs_its_solver_named(signal_data):
    """IEN: the nested base control names TENET, so IEN alone is refused, and a named solver runs."""
    X, y, n, p = signal_data
    ien = trexselector.GVSType.IEN
    with pytest.raises(ValueError, match="gvs_type = IEN requires solver_type"):
        trexselector.TRexGVSSelector(X, y, gvs_control=_gvs(ien), verbose=False)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        for solver in (TIENET, trexselector.SolverTypeForTRex.TIENET_AUG):
            sel = trexselector.TRexGVSSelector(X.copy(order="F"), y, gvs_control=_gvs(ien),
                                               trex_control=_base(solver), seed=3, verbose=False)
            assert sel.select().gvs_type == ien
        # The solver named on the nested control runs without a base control
        nested = _gvs(ien)
        nested.trex_ctrl.solver_type = TIENET
        nested.trex_ctrl.K = 3
        sel = trexselector.TRexGVSSelector(X.copy(order="F"), y, gvs_control=nested, seed=3, verbose=False)
        assert sel.K == 3
        assert sel.select().gvs_type == ien


def test_gvs_nested_base_control_applies_without_trex_control(signal_data):
    """Without trex_control the base control of the GVS control is the one of the run, as in C++."""
    X, y, n, p = signal_data
    gvs_ctrl = _gvs(trexselector.GVSType.EN)
    assert gvs_ctrl.trex_ctrl.solver_type == TENET
    gvs_ctrl.trex_ctrl.K = 7

    nested = trexselector.TRexGVSSelector(X.copy(order="F"), y, gvs_control=gvs_ctrl, seed=1, verbose=False)
    assert nested.K == 7
    selected = np.array(nested.select().selected_var)

    # The same control passed as trex_control gives the same run and replaces the nested one
    gvs_ctrl.trex_ctrl.K = 3
    passed = trexselector.TRexGVSSelector(X.copy(order="F"), y, gvs_control=gvs_ctrl, trex_control=_base(K=7),
                                          seed=1, verbose=False)
    assert passed.K == 7
    assert np.array_equal(np.array(passed.select().selected_var), selected)
    assert gvs_ctrl.trex_ctrl.K == 3                    # neither control is changed by the selector


def test_gvs_warns_on_lambda2_zero(signal_data):
    """lambda_2 == 0 is the degenerate no-ridge case, not automatic CV; the
    wrapper must warn so users don't silently get pure T-LASSO."""
    X, y, n, p = signal_data
    gvs_ctrl = trexselector.TRexGVSControlParameter()
    gvs_ctrl.lambda_2 = 0.0
    with pytest.warns(UserWarning, match="degenerate no-ridge"):
        trexselector.TRexGVSSelector(X, y, gvs_control=gvs_ctrl, verbose=False)


def test_gvs_no_lambda2_warning_for_auto_or_positive(signal_data):
    """lambda_2 < 0 (auto/CV) and lambda_2 > 0 (fixed) must not warn."""
    X, y, n, p = signal_data
    for lam in (-1.0, 0.05):
        gvs_ctrl = trexselector.TRexGVSControlParameter()
        gvs_ctrl.lambda_2 = lam
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            trexselector.TRexGVSSelector(X, y, gvs_control=gvs_ctrl, verbose=False)


# ---------------------------------------------------------------------------
# LinkageMethod variants (non-Ward only — Ward requires Euclidean)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("linkage", [
    pytest.param("Average", id="Average"),
    pytest.param("Complete", id="Complete"),
    pytest.param("Single", id="Single"),
])
def test_gvs_linkage_methods(signal_data, linkage):
    from trexselector._core.ml_methods.clustering import LinkageMethod
    X, y, n, p = signal_data
    gvs_ctrl = trexselector.TRexGVSControlParameter()
    gvs_ctrl.hc_linkage = getattr(LinkageMethod, linkage)
    sel = trexselector.TRexGVSSelector(X, y, gvs_control=gvs_ctrl,
                                         trex_control=_base(), verbose=False)
    res = sel.select()
    assert res is not None


# ---------------------------------------------------------------------------
# CCD family + IEN-geometry lambda_2 tuners (bindings-parity additions)
# ---------------------------------------------------------------------------

def test_lambda2_method_ien_entries_exist():
    m = trexselector.LambdaSelectionMethod
    for name in ("CV_1SE_IEN_CCD", "CV_MIN_IEN_CCD",
                 "CV_1SE_TIK_SVD", "CV_MIN_TIK_SVD"):
        assert hasattr(m, name)


def test_lambda2_method_track_resolution_and_restore(signal_data):
    X, y, n, p = signal_data
    gvs_ctrl = trexselector.TRexGVSControlParameter()
    gvs_ctrl.gvs_type = trexselector.GVSType.IEN
    gvs_ctrl.lambda_2 = 1.0
    assert gvs_ctrl.lambda2_method == trexselector.LambdaSelectionMethod.CV_1SE_CCD
    trexselector.TRexGVSSelector(X, y, gvs_control=gvs_ctrl,
                                      trex_control=_base(TIENET), verbose=False)
    # The C++ constructor resolves the default to the IEN-geometry tuner on
    # its own copy; the caller's object stays as set.
    assert gvs_ctrl.lambda2_method == trexselector.LambdaSelectionMethod.CV_1SE_CCD


def test_lambda2_method_ien_tuner_rejected_on_en_track(signal_data):
    X, y, n, p = signal_data
    gvs_ctrl = trexselector.TRexGVSControlParameter()
    gvs_ctrl.gvs_type = trexselector.GVSType.EN
    gvs_ctrl.lambda2_method = trexselector.LambdaSelectionMethod.CV_1SE_IEN_CCD
    with pytest.raises(ValueError, match="requires gvs_type = IEN"):
        trexselector.TRexGVSSelector(X, y, gvs_control=gvs_ctrl, verbose=False)


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_gvs_non_finite_cv_dfmax_ratio_raises(signal_data, bad):
    X, y, n, p = signal_data
    gvs_ctrl = trexselector.TRexGVSControlParameter()
    gvs_ctrl.gvs_type = trexselector.GVSType.IEN
    gvs_ctrl.cv_dfmax_ratio = bad
    with pytest.raises(ValueError, match="cv_dfmax_ratio must be finite"):
        trexselector.TRexGVSSelector(X, y, gvs_control=gvs_ctrl, trex_control=_base(TIENET), verbose=False)


def test_gvs_en_tcdenet_runs(signal_data):
    X, y, n, p = signal_data
    ctrl = _base()
    ctrl.solver_params.cd_max_sweeps = 200000   # converges on every draw of the dummies
    gvs_ctrl = trexselector.TRexGVSControlParameter()
    gvs_ctrl.en_solver = trexselector.ENSolverType.TCDENET
    gvs_ctrl.lambda_2 = 1.0
    sel = trexselector.TRexGVSSelector(X, y, gvs_control=gvs_ctrl,
                                            trex_control=ctrl, seed=42, verbose=False)
    res = sel.select()
    assert res is not None


def test_gvs_ien_solver_choice_tcdienet(signal_data):
    X, y, n, p = signal_data
    ctrl = trexselector.TRexControlParameter()
    ctrl.K = 3
    ctrl.solver_type = trexselector.SolverTypeForTRex.TCDIENET
    ctrl.solver_params.cd_max_sweeps = 200000   # converges on every draw of the dummies
    gvs_ctrl = trexselector.TRexGVSControlParameter()
    gvs_ctrl.gvs_type = trexselector.GVSType.IEN
    gvs_ctrl.lambda_2 = 1.0
    sel = trexselector.TRexGVSSelector(X, y, gvs_control=gvs_ctrl,
                                            trex_control=ctrl, seed=42, verbose=False)
    res = sel.select()
    assert res is not None


def test_gvs_ien_tuner_grids_and_diagnostics(signal_data):
    X, y, n, p = signal_data
    ctrl = _base(TIENET)

    # Fixed lambda_2: no tuner ran.
    fixed = trexselector.TRexGVSControlParameter()
    fixed.gvs_type = trexselector.GVSType.IEN
    fixed.lambda_2 = 1.0
    sel = trexselector.TRexGVSSelector(X, y, gvs_control=fixed,
                                            trex_control=ctrl, verbose=False)
    res = sel.select()
    assert res.lambda2_cv.ran is False
    assert sel.getLambda2CVDiagnostics().ran is False
    assert res.lambda2_cv.dfmax_ratio == 0.0     # no IEN-geometry tuner ran

    # IEN-geometry CV with reduced grids: the tuner ran and is reported.
    gvs_ctrl = trexselector.TRexGVSControlParameter()
    gvs_ctrl.gvs_type = trexselector.GVSType.IEN
    gvs_ctrl.lambda_2 = -1.0
    gvs_ctrl.cv_n_lambda1 = 20
    gvs_ctrl.cv_n_lambda2 = 10
    assert gvs_ctrl.cv_n_lambda1 == 20 and gvs_ctrl.cv_n_lambda2 == 10
    assert gvs_ctrl.cv_dfmax_ratio == pytest.approx(0.8)   # default
    gvs_ctrl.cv_dfmax_ratio = 0.9
    assert gvs_ctrl.cv_dfmax_ratio == pytest.approx(0.9)
    sel = trexselector.TRexGVSSelector(X, y, gvs_control=gvs_ctrl,
                                            trex_control=ctrl, verbose=False)
    res = sel.select()
    diag = res.lambda2_cv
    assert diag.ran is True
    assert diag.seconds >= 0.0
    assert diag.n_profile_at_edge >= 0
    assert res.lambda2_used > 0

    # The support-cap ratio reached the tuner and is recorded as passed.
    assert diag.dfmax_ratio == pytest.approx(0.9)
    assert isinstance(diag.profile_at_support_cap, bool)
    assert diag.profile_at_edge or not diag.profile_at_support_cap

    # Boundary-pick guard fields: the chosen pick on its grid.
    assert diag.grid_size == 10
    assert 0 <= diag.grid_index < diag.grid_size
    assert isinstance(diag.at_grid_edge, bool) and isinstance(diag.profile_at_edge, bool)
    assert diag.boundary_pick == (diag.at_grid_edge or diag.profile_at_edge)
    assert diag.at_grid_edge == (diag.grid_index in (0, diag.grid_size - 1))


def test_gvs_result_names_the_origin_of_the_groups(signal_data, gvs_base):
    """Supplied groups give "prior", their number and their names, and the supplied lambda_2."""
    X, y, n, p = signal_data
    gvs_ctrl = trexselector.TRexGVSControlParameter()
    gvs_ctrl.prior_groups = [i % 4 for i in range(p)]
    gvs_ctrl.group_labels = ["a", "b", "c", "d"]
    gvs_ctrl.lambda_2 = 1.0
    sel = trexselector.TRexGVSSelector(X, y, gvs_control=gvs_ctrl, trex_control=gvs_base, verbose=False)
    res = sel.select()
    assert res.hc_method_used == "prior"
    assert res.max_clusters == 4
    assert res.group_labels == ["a", "b", "c", "d"]
    assert res.lambda2_used == 1.0
    assert not res.lambda2_cv.ran


@pytest.mark.parametrize("linkage", ["Single", "Complete", "Average", "WPGMA"])
def test_gvs_result_names_the_linkage_in_lower_case(signal_data, gvs_base, linkage):
    """Without groups hc_method_used holds the linkage of the clustering in lower case."""
    from trexselector.ml_methods.clustering import LinkageMethod

    X, y, n, p = signal_data
    gvs_ctrl = trexselector.TRexGVSControlParameter()
    gvs_ctrl.hc_linkage = getattr(LinkageMethod, linkage)
    gvs_ctrl.lambda_2 = 1.0
    sel = trexselector.TRexGVSSelector(X, y, gvs_control=gvs_ctrl, trex_control=gvs_base, verbose=False)
    res = sel.select()
    assert res.hc_method_used == linkage.lower()
    assert res.group_labels == []


@pytest.mark.parametrize("strategy", ["PERMUTATION", "PERMUTATION_SEEDED"])
def test_gvs_refuses_the_permutation_strategies(signal_data, strategy):
    """A row permutation destroys the covariance of the dummies within a cluster."""
    X, y, n, p = signal_data
    ctrl = _base()
    ctrl.lloop_strategy = getattr(trexselector.LLoopStrategy, strategy)
    with pytest.raises(ValueError, match="PERMUTATION and PERMUTATION_SEEDED"):
        sel = trexselector.TRexGVSSelector(X, y, trex_control=ctrl, verbose=False)
        sel.select()


def test_unconverged_coordinate_descent_raises_and_gives_the_design_back():
    """TCDIENET at a lambda_2 beyond its sweep budget, one group per variable.

    select() ends with an error instead of a selection from an unconverged path, the design is
    given back at once and the selector is retired.
    """
    rng = np.random.default_rng(2718)
    n, p = 120, 60
    X = np.asfortranarray(rng.standard_normal((n, p)))
    y = X[:, 20:30].sum(axis=1) + rng.standard_normal(n)
    X0 = X.copy(order="F")

    gvs_ctrl = trexselector.TRexGVSControlParameter()
    gvs_ctrl.gvs_type = trexselector.GVSType.IEN
    gvs_ctrl.prior_groups = list(range(p))
    gvs_ctrl.lambda_2 = 3000.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)      # the solver is deprecated
        sel = trexselector.TRexGVSSelector(X, y, tFDR=0.1, gvs_control=gvs_ctrl,
                                           trex_control=_base(trexselector.SolverTypeForTRex.TCDIENET, K=10),
                                           seed=4, verbose=False)
    assert not np.allclose(X, X0)                     # standardized in place
    with pytest.raises(RuntimeError, match="sweep budget.*before it converged.*cd_max_sweeps"):
        sel.select()
    np.testing.assert_allclose(X, X0, atol=1e-10)     # given back
    with pytest.raises(RuntimeError, match="can no longer run"):
        sel.select()                                  # retired


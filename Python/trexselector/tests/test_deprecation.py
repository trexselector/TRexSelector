"""
Deprecation of the coordinate-descent solver family.

A solver class of the family warns when it is constructed, and a selector warns
when a solver type of the family reaches it: one FutureWarning per construction,
attributed to the caller, and none from select().
"""
import warnings

import numpy as np
import pytest

import trexselector as tsn
from trexselector.trex_selector_methods import ad as raw
from trexselector.tsolvers.ad import SolverOptions
from trexselector.tsolvers.ad.cd_based import (
    TCDENET_Solver, TCDIENET_Solver, TCDLASSO_Solver, TCDTIKNET_Solver,
)
from trexselector.tsolvers.ad.lars_based import (
    TENET_Solver, TIENET_Solver, TIKNET_Solver, TLASSO_Solver,
)

N, P = 60, 12
GROUPS = np.repeat(np.arange(P // 2), 2).astype(np.int32)


def _message(old, new):
    return (f"{old} is deprecated and will be removed in the next minor release of trexselector. "
            f"Use {new}, the homotopy solver of the same problem.")


def _deprecations(fn):
    """Run fn and return its result with the FutureWarnings it raised."""
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        out = fn()
    return out, [w for w in record if issubclass(w.category, FutureWarning)]


def _assert_one(record, old, new):
    assert [str(w.message) for w in record] == [_message(old, new)]
    assert record[0].filename == __file__   # attributed to the caller, not to the package


@pytest.fixture
def data():
    rng = np.random.default_rng(7)
    X = np.asfortranarray(rng.standard_normal((N, P)))
    D = np.asfortranarray(rng.standard_normal((N, P)))
    y = X[:, :3].sum(axis=1) * 3.0 + rng.standard_normal(N)
    return X, D, y


def _control(solver=None):
    ctrl = tsn.TRexControlParameter()
    ctrl.K = 3
    ctrl.solver_params.cd_max_sweeps = 200000   # the runs converge on every draw of the dummies
    if solver is not None:
        ctrl.solver_type = getattr(tsn.SolverTypeForTRex, solver)
    return ctrl


def _tik_control():
    tik = tsn.TRexTikhonovControlParameter()
    tik.tikhonov_K = np.eye(P)
    tik.lambda_2 = 0.5
    return tik


# ---------------------------------------------------------------------------
# Solver classes
# ---------------------------------------------------------------------------

_RIDGE = dict(lambda2=0.1)

CD_CLASSES = [
    ("TCDLASSO_Solver", "TLASSO_Solver", lambda X, D, y: TCDLASSO_Solver(X, D, y)),
    ("TCDENET_Solver", "TENET_Solver", lambda X, D, y: TCDENET_Solver(X, D, y, SolverOptions(**_RIDGE))),
    ("TCDIENET_Solver", "TIENET_Solver",
     lambda X, D, y: TCDIENET_Solver(X, D, y, GROUPS, SolverOptions(**_RIDGE))),
    ("TCDTIKNET_Solver", "TIKNET_Solver",
     lambda X, D, y: TCDTIKNET_Solver(X, D, y, np.eye(P), SolverOptions(**_RIDGE))),
]

HOMOTOPY_CLASSES = [
    lambda X, D, y: TLASSO_Solver(X, D, y),
    lambda X, D, y: TENET_Solver(X, D, y, SolverOptions(**_RIDGE)),
    lambda X, D, y: TIENET_Solver(X, D, y, GROUPS, SolverOptions(**_RIDGE)),
    lambda X, D, y: TIKNET_Solver(X, D, y, np.eye(P), SolverOptions(**_RIDGE)),
]


@pytest.mark.parametrize("old, new, make", CD_CLASSES, ids=[c[0] for c in CD_CLASSES])
def test_cd_solver_class_warns_once_and_names_the_replacement(data, old, new, make):
    X, D, y = data
    solver, record = _deprecations(lambda: make(X, D, y))
    _assert_one(record, old, new)

    # The solver works as before
    solver.executeStep(2)
    assert solver.getNumSteps() > 0


@pytest.mark.parametrize("make", HOMOTOPY_CLASSES, ids=["TLASSO", "TENET", "TIENET", "TIKNET"])
def test_homotopy_solver_classes_do_not_warn(data, make):
    X, D, y = data
    _, record = _deprecations(lambda: make(X, D, y))
    assert record == []


# ---------------------------------------------------------------------------
# Selectors
# ---------------------------------------------------------------------------

def _biobank(X, y, ctrl):
    bio = tsn.BiobankScreenTRexControl()
    bio.trex_screen_ctrl.trex_ctrl = ctrl
    return tsn.TRexBiobankScreeningSelector(X, y, bio_ctrl=bio, verbose=False)


BASE_SELECTORS = [
    ("plain", lambda X, y, ctrl: tsn.TRexSelector(X, y, trex_control=ctrl, verbose=False)),
    ("da", lambda X, y, ctrl: tsn.TRexDASelector(X, y, trex_control=ctrl, verbose=False)),
    ("screening", lambda X, y, ctrl: tsn.TRexScreeningSelector(X, y, trex_control=ctrl, verbose=False)),
    ("biobank", _biobank),
    ("plain, subpackage class", lambda X, y, ctrl: raw.TRexSelector(X, y, 0.1, ctrl, 42, False)),
]


@pytest.mark.parametrize("make", [s[1] for s in BASE_SELECTORS], ids=[s[0] for s in BASE_SELECTORS])
def test_selector_warns_once_for_a_cd_solver_type(data, make):
    X, _, y = data
    _, record = _deprecations(lambda: make(X.copy(order="F"), y, _control("TCDENET")))
    _assert_one(record, "SolverTypeForTRex.TCDENET", "SolverTypeForTRex.TENET")

    # The homotopy solver of the same problem and the default raise none
    for ctrl in (_control("TENET"), _control()):
        _, record = _deprecations(lambda: make(X.copy(order="F"), y, ctrl))
        assert record == []


def test_select_raises_no_further_warning(data):
    X, _, y = data
    sel, record = _deprecations(
        lambda: tsn.TRexSelector(X, y, trex_control=_control("TCDLASSO"), seed=42, verbose=False))
    _assert_one(record, "SolverTypeForTRex.TCDLASSO", "SolverTypeForTRex.TLASSO")

    # The K experiments run the deprecated solver without a second warning
    result, record = _deprecations(sel.select)
    assert record == []
    assert result.selected_var.shape[0] == P


def test_grouped_selector_warns_on_both_tracks(data):
    X, _, y = data

    # Elastic-net track: the solver comes from en_solver
    gvs = tsn.TRexGVSControlParameter()
    gvs.en_solver = tsn.ENSolverType.TCDENET
    gvs.lambda_2 = 1.0
    _, record = _deprecations(
        lambda: tsn.TRexGVSSelector(X.copy(order="F"), y, gvs_control=gvs, trex_control=_control(), verbose=False))
    _assert_one(record, "ENSolverType.TCDENET", "ENSolverType.TENET")

    # Informed elastic-net track: the solver comes from the base control
    gvs = tsn.TRexGVSControlParameter()
    gvs.gvs_type = tsn.GVSType.IEN
    gvs.lambda_2 = 1.0
    _, record = _deprecations(
        lambda: tsn.TRexGVSSelector(X.copy(order="F"), y, gvs_control=gvs, trex_control=_control("TCDIENET"),
                                    verbose=False))
    _assert_one(record, "SolverTypeForTRex.TCDIENET", "SolverTypeForTRex.TIENET")

    # The defaults of both tracks raise none
    for gvs_type in (tsn.GVSType.EN, tsn.GVSType.IEN):
        gvs = tsn.TRexGVSControlParameter()
        gvs.gvs_type = gvs_type
        gvs.lambda_2 = 1.0
        _, record = _deprecations(
            lambda: tsn.TRexGVSSelector(X.copy(order="F"), y, gvs_control=gvs, trex_control=_control(),
                                        verbose=False))
        assert record == []


def test_tikhonov_selector_and_power_rule_warn(data):
    X, _, y = data

    _, record = _deprecations(
        lambda: tsn.TRexTikhonovSelector(X.copy(order="F"), y, tik_control=_tik_control(),
                                         trex_control=_control("TCDTIKNET"), verbose=False))
    _assert_one(record, "SolverTypeForTRex.TCDTIKNET", "SolverTypeForTRex.TIKNET")

    result, record = _deprecations(
        lambda: tsn.tikhonov_lambda2_power_rule(X.copy(order="F"), y, tik_control=_tik_control(),
                                                trex_control=_control("TCDTIKNET"), seed=1, grid=[0.0, 1.0]))
    _assert_one(record, "SolverTypeForTRex.TCDTIKNET", "SolverTypeForTRex.TIKNET")
    assert len(result.probes) == 2

    # The default solver, the homotopy TIKNET, raises none
    _, record = _deprecations(
        lambda: tsn.TRexTikhonovSelector(X.copy(order="F"), y, tik_control=_tik_control(), verbose=False))
    assert record == []


def test_spca_selector_warns_for_the_cd_elastic_net(data):
    X, _, _ = data
    ctrl = tsn.TRexSPCAControlParameter()
    ctrl.en_solver = tsn.ENSolverType.TCDENET
    _, record = _deprecations(lambda: raw.TRexSPCASelector(X.copy(order="F"), 1, 0.2, ctrl, 1, False))
    _assert_one(record, "ENSolverType.TCDENET", "ENSolverType.TENET")

    _, record = _deprecations(
        lambda: raw.TRexSPCASelector(X.copy(order="F"), 1, 0.2, tsn.TRexSPCAControlParameter(), 1, False))
    assert record == []


def test_warning_turned_into_an_error_leaves_the_design_untouched(data):
    X, _, y = data
    X0 = X.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        with pytest.raises(FutureWarning, match="SolverTypeForTRex.TCDLASSO is deprecated"):
            tsn.TRexSelector(X, y, trex_control=_control("TCDLASSO"), verbose=False)

    # Nothing was standardized or claimed: the next selector takes the same array
    np.testing.assert_array_equal(X, X0)
    result = tsn.TRexSelector(X, y, trex_control=_control("TLASSO"), seed=42, verbose=False).select()
    assert result.selected_var.shape[0] == P

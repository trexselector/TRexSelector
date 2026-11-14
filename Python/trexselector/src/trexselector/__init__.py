from __future__ import annotations

import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import inspect
import warnings

from importlib.metadata import PackageNotFoundError, version as _pkg_version

import numpy as np

try:
    __version__ = _pkg_version("trexselector")
except PackageNotFoundError:  # Running from a source tree without install
    __version__ = "0+unknown"

from . import ml_methods
from . import tsolvers
from . import utils

# Import the compiled C++ extension inside the nested modules. The root
# re-exports the augmented-dummy (ad) selectors as the everyday spelling;
# the vd / sd families, once bound, stay under trex_selector_methods.vd / .sd.
from . import trex_selector_methods
from .trex_selector_methods.ad import (
    LLoopStrategy, ScalingMode, TRexControlParameter, SelectionResult, TRexSelector as PyTRexSelector,
    TRexDASelector as PyTRexDASelector, TRexDAControlParameter, DAMethod, BTSelectionMode, DASelectionResult,
    TRexGVSSelector as PyTRexGVSSelector, TRexGVSControlParameter, GVSType, ENSolverType,
                        LambdaSelectionMethod, GVSSelectionResult,
    TRexScreeningSelector as PyTRexScreeningSelector, ScreenTRexControlParameter, ScreenTRexMethod,
                             ScreenTRexSelectionResult,
    TRexBiobankScreeningSelector as PyBiobankScreeningSelector, BiobankScreenTRexControl,
                                    BiobankScreenTRexResult,
    SolverTypeForTRex, SolverHyperparameters, DummyDistribution,
    SPCAMode, TRexSPCAControlParameter, TRexSPCAResult,
    TRexSPCASelector as PyTRexSPCASelector,
    Lambda2CVDiagnostics,
    TikLambda2Method, PilotStatistic, TikPenaltyStructure,
    Lambda2PickInfo, Lambda2PilotResult, Lambda2Probe, Lambda2PowerRuleResult, DummyClusters,
    TRexTikhonovControlParameter,
    TRexTikhonovSelector as PyTRexTikhonovSelector,
    gamma_to_k, build_group_centering_k, default_lambda2_grid,
    tikhonov_lambda2_power_rule as _tikhonov_lambda2_power_rule,
    gvs_lambda2_power_rule as _gvs_lambda2_power_rule,
    normalize_design,
)


def _validate_core_solver_type(trex_control: "TRexControlParameter") -> None:
    """Reject solver types the core dispatch cannot serve, at construction.

    The shared C++ dispatch throws the same guidance at SELECT time, but from
    inside the OpenMP experiment loop, where the exception cannot cross the
    region boundary and terminates the process instead of raising. Used by
    the plain/DA/screening/biobank wrappers; the GVS wrapper deliberately
    bypasses it (the IEN family is exactly its vocabulary).
    """
    st = trex_control.solver_type
    if st == SolverTypeForTRex.TIENET:
        raise ValueError(
            'solver_type = TIENET is a GVS-only solver (needs a group '
            'assignment); use TRexGVSSelector with gvs_type = IEN.')
    if st == SolverTypeForTRex.TIENET_AUG:
        raise ValueError(
            'solver_type = TIENET_AUG is a GVS-only solver; use '
            'TRexGVSSelector with gvs_type = IEN.')
    if st == SolverTypeForTRex.TCDIENET:
        raise ValueError(
            'solver_type = TCDIENET needs a group assignment; use '
            'TRexGVSSelector with gvs_type = IEN.')
    if st in (SolverTypeForTRex.TCDTIKNET, SolverTypeForTRex.TIKNET):
        raise ValueError(
            f'solver_type = {st.name} needs a Tikhonov penalty matrix; use '
            'TRexTikhonovSelector.')


def _with_core_notes(core_cls):
    """Class decorator: append the notes and the references of the core class.

    The binding docstring of the core class holds the notes of a selector and
    its rendered reference list. The selector of the package root takes both
    from there, so the text has one source.
    """
    def decorate(cls):
        own = inspect.cleandoc(cls.__doc__ or "")
        lines = inspect.cleandoc(core_cls.__doc__ or "").splitlines()
        # First section heading: a title line above a line of dashes
        start = next((i for i in range(len(lines) - 1)
                      if lines[i].strip() and lines[i + 1].strip()
                      and set(lines[i + 1].strip()) == {"-"}), None)
        if own and start is not None:
            cls.__doc__ = own + "\n\n" + "\n".join(lines[start:]) + "\n"
        return cls
    return decorate


@_with_core_notes(PyTRexSelector)
class TRexSelector:
    """
    T-Rex Selector: FDR-controlled variable selection via dummy comparison.

    X is used without a copy when it is a Fortran-contiguous float64 array, for
    example the ``to_numpy()`` view of a ``MemoryMappedMatrix``. The selector
    standardizes X in place and restores it, so X must be writable. Under
    ``ScalingMode.NONE`` it writes nothing into X and accepts an array that is not
    writable, for example the view of a file opened with ``AccessMode.ReadOnly``.
    The selectors derived from this class follow the same rule.

    The coordinate-descent solver types ``TCDLASSO``, ``TCDENET``, ``TCDIENET``
    and ``TCDTIKNET`` are deprecated and will be removed in the next minor
    release: a selector that receives one of them raises ``FutureWarning``. Use
    ``TLASSO``, ``TENET``, ``TIENET`` and ``TIKNET``, the homotopy solvers of the
    same problems.

    Parameters
    ----------
    X : numpy.ndarray of shape (n, p)
        Design matrix.
    y : numpy.ndarray of shape (n,)
        Response vector.
    tFDR : float, default 0.1
        Target FDR level.
    trex_control : TRexControlParameter, optional
        Control parameters of the T-Rex selector. The defaults apply when it
        is omitted.
    seed : int, default -1
        Random seed. A negative value gives a run that is not reproducible.
    verbose : bool, default True
        Whether to print progress.
    """

    def __init__(self, X: np.ndarray, y: np.ndarray, tFDR: float = 0.1,
                 trex_control: TRexControlParameter | None = None,
                 seed: int = -1, verbose: bool = True):

        # Validate data types and ensure Fortran contiguity for Eigen
        # Eigen Map defaults to ColumnMajor (Fortran order)
        self.X = np.asfortranarray(X, dtype=np.float64)
        self.y = np.asfortranarray(y, dtype=np.float64)

        if trex_control is None:
            trex_control = TRexControlParameter()
        _validate_core_solver_type(trex_control)

        self._selector = PyTRexSelector(self.X, self.y, tFDR, trex_control, seed, verbose)

    def select(self) -> SelectionResult:
        """
        Executes the T-Rex selector algorithm.
        Releases the GIL to allow underlying C++ OpenMP multithreading to run in parallel.

        A selector runs once: ``select()`` gives ``X`` back when it ends, so a
        second call raises ``RuntimeError`` instead of running on a design that
        is no longer standardized. Construct a new selector to run again. A
        ``select()`` that ends with an error gives ``X`` back at once (original
        scale) and retires the selector as well.

        Returns
        -------
        SelectionResult
            A data class containing the selected variables and extensive metadata.
        """
        return self._selector.select()

    def restore(self) -> None:
        """
        Give X back before select(): undo the in-place normalization applied at
        construction and release the buffer claim.

        Idempotent, and a no-op after select() (which restores X itself). A
        restored selector cannot run: select() then raises RuntimeError,
        construct a new one. Without it, an unselected selector gives X back
        only when it is garbage-collected.
        """
        self._selector.restore()

    @property
    def tFDR(self) -> float:
        return self._selector.getTFDR()

    @property
    def K(self) -> int:
        return self._selector.getK()

    @property
    def num_dummies(self) -> int:
        """Number of dummies used by the last select() (the L-loop's calibration)."""
        return self._selector.getNumDummies()

    @property
    def num_selected(self) -> int:
        return self._selector.getNumSelected()

    @property
    def T_stop(self) -> int:
        return self._selector.getTStop()

    @property
    def L(self) -> int:
        return self._selector.getDummyMultiplierL()

    @property
    def v_thresh(self) -> float:
        return self._selector.getVotingThreshold()

    @property
    def voting_grid(self) -> np.ndarray:
        return self._selector.getVotingGrid()

    @property
    def fdp_hat_mat(self) -> np.ndarray:
        return self._selector.getFDPHatMat()

    @property
    def phi_mat(self) -> np.ndarray:
        return self._selector.getPhiMat()

    @property
    def r_mat(self) -> np.ndarray:
        return self._selector.getRMat()

    @property
    def phi_prime(self) -> np.ndarray:
        return self._selector.getPhiPrime()

    @property
    def x_means(self) -> np.ndarray:
        return self._selector.getXMeans()

    @property
    def x_l2_norms(self) -> np.ndarray:
        return self._selector.getXL2Norms()

    @property
    def y_mean(self) -> float:
        return self._selector.getYMean()

    @property
    def selected_var(self) -> np.ndarray:
        return self._selector.getSelectedVar()

    @property
    def selected_indices(self) -> np.ndarray:
        return self._selector.getSelectedIndices()


@_with_core_notes(PyTRexDASelector)
class TRexDASelector(TRexSelector):
    """
    Dependency-aware T-Rex selector (T-Rex+DA): FDR-controlled variable
    selection for dependent variables.

    Parameters
    ----------
    X : numpy.ndarray of shape (n, p)
        Design matrix, used as described for ``TRexSelector``.
    y : numpy.ndarray of shape (n,)
        Response vector.
    tFDR : float, default 0.1
        Target FDR level.
    da_control : TRexDAControlParameter, optional
        Control parameters of the dependency-aware correction. The defaults
        apply when it is omitted.
    trex_control : TRexControlParameter, optional
        Control parameters of the T-Rex selector. The defaults apply when it
        is omitted.
    seed : int, default -1
        Random seed. A negative value gives a run that is not reproducible.
    verbose : bool, default True
        Whether to print progress.
    """
    def __init__(self, X: np.ndarray, y: np.ndarray, tFDR: float = 0.1,
                 da_control: TRexDAControlParameter | None = None,
                 trex_control: TRexControlParameter | None = None,
                 seed: int = -1, verbose: bool = True):

        self.X = np.asfortranarray(X, dtype=np.float64)
        self.y = np.asfortranarray(y, dtype=np.float64)

        if da_control is None:
            da_control = TRexDAControlParameter()
        if trex_control is None:
            trex_control = TRexControlParameter()
        _validate_core_solver_type(trex_control)

        self._selector: PyTRexDASelector = PyTRexDASelector(self.X, self.y, tFDR,
                                                            da_control,
                                                            trex_control,
                                                            seed, verbose)

    def select(self) -> DASelectionResult:
        """
        Run the selector. A selector runs once: a second call raises
        ``RuntimeError``, see ``TRexSelector.select``.

        Returns
        -------
        DASelectionResult
            The selection with the record of the dependency estimation.
        """
        self._selector.select()
        return self._selector.getDAResult()


@_with_core_notes(PyTRexGVSSelector)
class TRexGVSSelector(TRexSelector):
    """
    T-Rex selector for grouped variable selection (T-Rex+GVS): FDR-controlled
    variable selection for groups of highly correlated variables.

    The coordinate-descent solvers of the two tracks, ``ENSolverType.TCDENET``
    and ``SolverTypeForTRex.TCDIENET``, are deprecated and will be removed in
    the next minor release: the selector raises ``FutureWarning`` for them.
    Use ``ENSolverType.TENET`` and ``SolverTypeForTRex.TIENET``.

    Parameters
    ----------
    X : numpy.ndarray of shape (n, p)
        Design matrix, used as described for ``TRexSelector``.
    y : numpy.ndarray of shape (n,)
        Response vector.
    tFDR : float, default 0.1
        Target FDR level.
    gvs_control : TRexGVSControlParameter, optional
        Control parameters of the grouped variable selection. The defaults
        apply when it is omitted.
    trex_control : TRexControlParameter, optional
        Control parameters of the T-Rex selector. The defaults apply when it
        is omitted.
    seed : int, default -1
        Random seed. A negative value gives a run that is not reproducible.
    verbose : bool, default True
        Whether to print progress.
    """
    def __init__(self, X: np.ndarray, y: np.ndarray, tFDR: float = 0.1,
                 gvs_control: TRexGVSControlParameter | None = None,
                 trex_control: TRexControlParameter | None = None,
                 seed: int = -1, verbose: bool = True):

        self.X = np.asfortranarray(X, dtype=np.float64)
        self.y = np.asfortranarray(y, dtype=np.float64)

        if gvs_control is None:
            gvs_control = TRexGVSControlParameter()
        if trex_control is None:
            trex_control = TRexControlParameter()

        # Solver selection mirrors the R surface:
        #   EN  -> derived from en_solver (TENET / TENET_AUG / TCDENET);
        #          trex_control.solver_type is ignored (warn on an explicit
        #          non-matching value; TLARS is the default and cannot be
        #          distinguished from an explicit choice, so it never warns).
        #   IEN -> trex_control.solver_type IS the IEN-solver choice:
        #          TIENET (native pathwise, default), TIENET_AUG, TCDIENET.
        #          The default TLARS maps to TIENET; any other value warns
        #          and falls back to TIENET (in the C++ glue).
        _ien_family = (SolverTypeForTRex.TIENET, SolverTypeForTRex.TIENET_AUG,
                       SolverTypeForTRex.TCDIENET)
        if gvs_control.gvs_type == GVSType.IEN:
            if trex_control.solver_type not in (SolverTypeForTRex.TLARS,) + _ien_family:
                warnings.warn(
                    f"trex_control.solver_type = {trex_control.solver_type.name} "
                    f"is not an IEN-family solver; gvs_type = IEN accepts "
                    f"TIENET, TIENET_AUG, or TCDIENET. Using TIENET.",
                    UserWarning,
                    stacklevel=2,
                )
        else:
            if gvs_control.en_solver == ENSolverType.TENET_AUG:
                derived_solver = SolverTypeForTRex.TENET_AUG
            elif gvs_control.en_solver == ENSolverType.TCDENET:
                derived_solver = SolverTypeForTRex.TCDENET
            else:
                derived_solver = SolverTypeForTRex.TENET
            if trex_control.solver_type not in (SolverTypeForTRex.TLARS, derived_solver):
                warnings.warn(
                    f"trex_control.solver_type = {trex_control.solver_type.name} is "
                    f"ignored by TRexGVSSelector: gvs_type = EN "
                    f"derives solver {derived_solver.name} from en_solver.",
                    UserWarning,
                    stacklevel=2,
                )

        # lambda2_method resolves by track in the C++ constructor: the struct
        # default CV_1SE_CCD becomes CV_1SE_IEN_CCD on the IEN track, and the
        # IEN-geometry tuners are refused on the EN track (ValueError).

        # lambda_2 == 0 is the degenerate no-ridge case (pure T-LASSO), NOT
        # automatic CV (the sentinel for that is lambda_2 < 0). Warn so a user
        # who passed exactly 0 expecting CV is not silently given no ridge.
        if gvs_control.lambda_2 == 0.0:
            warnings.warn(
                "gvs_control.lambda_2 = 0 is the degenerate no-ridge case "
                "(pure T-LASSO), not automatic cross-validation. Use "
                "lambda_2 < 0 (e.g. -1) to auto-determine lambda_2 via CV.",
                UserWarning,
                stacklevel=2,
            )
        self._selector: PyTRexGVSSelector = PyTRexGVSSelector(self.X, self.y, tFDR,
                                                              gvs_control,
                                                              trex_control,
                                                              seed, verbose)

    def select(self) -> GVSSelectionResult:
        """
        Run the selector. A selector runs once: a second call raises
        ``RuntimeError``, see ``TRexSelector.select``.

        Returns
        -------
        GVSSelectionResult
            The selection with the GVS-specific fields.
        """
        self._selector.select()
        return self._selector.getGVSResult()

    def getLambda2CVDiagnostics(self):
        """Health record of the lambda_2 tuner of the last select().

        Returns a Lambda2CVDiagnostics: ``ran`` (a tuner ran; False when
        lambda_2 came from the caller), ``seconds``, ``n_profile_at_edge``
        (lambda_2 values whose profiled lambda_1 optimum sat on a
        truncation edge; 0 is healthy), ``support_cap`` (-1 = disabled),
        ``dfmax_ratio`` (IEN tuner: the ``cv_dfmax_ratio`` it was called
        with; 0 when no IEN-geometry tuner ran), ``low_dimensional``, and
        the boundary-pick guard fields ``grid_index`` / ``grid_size``
        (0-based position of the returned lambda_2 on the tuner's grid),
        ``at_grid_edge`` (first or last grid point: the grid boundary, not a
        CV optimum), ``profile_at_edge`` (IEN tuner: the chosen lambda_2's
        profiled lambda_1 optimum sat on the truncation edge),
        ``profile_at_support_cap`` (that edge is the support cap; raise
        ``cv_dfmax_ratio`` to move it) and ``boundary_pick`` (at_grid_edge
        or profile_at_edge; the selector warned once and returned the pick
        unchanged). The same record is carried in
        ``GVSSelectionResult.lambda2_cv``.
        """
        return self._selector.getLambda2CVDiagnostics()


@_with_core_notes(PyTRexTikhonovSelector)
class TRexTikhonovSelector(TRexSelector):
    """
    T-Rex selector with a general Tikhonov penalty (TIKNET, or the deprecated
    TCDTIKNET) and cluster-aware dummies.

    The penalty comes from ``tik_control``: a p x p PSD matrix ``tikhonov_K``
    (``penalty_structure = SPARSE_GENERAL``), a scipy.sparse matrix or a dense
    array, or ``GROUP_CENTERING`` built from the groups. lambda_2 < 0
    resolves by ``tik_control.lambda2_method`` (CV tuners or the pilot power
    rules); the dummies are cluster-aware by default. Solver: TIKNET, the
    pathwise homotopy solver (default), or TCDTIKNET, the coordinate-descent
    solver, chosen through the base control. Any other solver raises
    ``ValueError``, also ``TLARS``, the default of ``TRexControlParameter``:
    a base control that is passed has to name the solver. TCDTIKNET is not
    finally developed: once lambda_2 times K is large it can stop at its sweep
    budget (``cd_max_sweeps``) before it has converged. ``select()`` then
    raises ``RuntimeError`` instead of returning a selection from an
    unconverged path. TCDTIKNET is deprecated and will be removed in the next
    minor release: the selector raises ``FutureWarning`` for it.

    Parameters
    ----------
    X : numpy.ndarray of shape (n, p)
        Design matrix, used as described for ``TRexSelector``.
    y : numpy.ndarray of shape (n,)
        Response vector.
    tFDR : float, default 0.1
        Target FDR level.
    tik_control : TRexTikhonovControlParameter
        Control parameters of the Tikhonov selector. Required, because it
        carries the penalty matrix ``tikhonov_K`` or the penalty structure.
    trex_control : TRexControlParameter, optional
        Control parameters of the T-Rex selector. When it is omitted, the
        base control nested in ``tik_control`` applies
        (``tik_control.trex_ctrl``, whose solver is TIKNET). A control that
        is passed replaces the nested one as it is, so its ``solver_type``
        must be ``SolverTypeForTRex.TIKNET``.
    seed : int, default -1
        Random seed. A negative value gives a run that is not reproducible.
    verbose : bool, default True
        Whether to print progress.
    """
    def __init__(self, X: np.ndarray, y: np.ndarray, tFDR: float = 0.1,
                 tik_control: TRexTikhonovControlParameter | None = None,
                 trex_control: TRexControlParameter | None = None,
                 seed: int = -1, verbose: bool = True):

        self.X = np.asfortranarray(X, dtype=np.float64)
        self.y = np.asfortranarray(y, dtype=np.float64)

        if tik_control is None:
            raise ValueError(
                "tik_control is required (a TRexTikhonovControlParameter; it "
                "carries the penalty matrix tikhonov_K or the penalty structure).")
        # The core builds the control of the run: without trex_control the base
        # control nested in tik_control applies, and the core validates it.
        self._selector: PyTRexTikhonovSelector = PyTRexTikhonovSelector(
            self.X, self.y, tFDR, tik_control, trex_control, seed, verbose)

    @property
    def lambda2_used(self) -> float:
        """lambda_2 used by the last select() (solver scale); -1 before select()."""
        return self._selector.getLambda2Used()

    @property
    def lambda2_cv_diagnostics(self) -> Lambda2CVDiagnostics:
        """Health record of the lambda_2 tuner of the last select() (Lambda2CVDiagnostics)."""
        return self._selector.getLambda2CVDiagnostics()

    @property
    def lambda2_pilot(self) -> Lambda2PilotResult:
        """Record of the lambda_2 pilot of the last select() (Lambda2PilotResult; ``ran`` is False unless the
        pilot rule resolved lambda_2)."""
        return self._selector.getLambda2PilotResult()

    @property
    def dummy_clusters(self) -> DummyClusters:
        """The clusters the dummies of the last run were drawn for (DummyClusters): ``groups_vec``, the
        cluster of each variable, ``max_clusters`` and ``hc_method_used``. The record is empty before a run
        and with ``cluster_dummies = False``."""
        return self._selector.getDummyClusters()

    def run_pilot_only(self) -> Lambda2PilotResult:
        """
        Run only the lambda_2 pilot, without the selection.

        The pilot builds the cluster model, runs the pilot experiments and
        picks lambda_2, and X is given back. It is for a caller who wants the
        pick alone. It needs ``lambda_2 < 0`` and
        ``lambda2_method = TikLambda2Method.PILOT_KNEE`` in ``tik_control``
        (``ValueError`` otherwise).

        A selector runs once: after ``run_pilot_only()`` a call of ``select()``
        or of ``run_pilot_only()`` raises ``RuntimeError``, and so does
        ``run_pilot_only()`` after ``select()``.

        Returns
        -------
        Lambda2PilotResult
            The record of the pilot, as ``lambda2_pilot`` holds it afterwards.
        """
        return self._selector.runPilotOnly()


def tikhonov_lambda2_power_rule(X: np.ndarray, y: np.ndarray, tFDR: float = 0.1,
                                tik_control: TRexTikhonovControlParameter | None = None,
                                trex_control: TRexControlParameter | None = None,
                                seed: int = 0,
                                grid: "list[float] | np.ndarray | None" = None,
                                knee_tolerance: float = 1.0,
                                verbose: bool = False) -> Lambda2PowerRuleResult:
    """
    The knee of the full lambda_2 profile of the Tikhonov T-Rex.

    Runs one FDR-controlled Tikhonov T-Rex selection per value of ``grid`` with
    the same seed, so the probes share their dummies and differ by the penalty
    weight alone, and picks the smallest lambda_2 whose selection is within
    ``knee_tolerance`` selections of the largest (the knee of the profile; 0
    reads the largest selection, ties toward the smallest lambda_2). This is
    the reference rule behind the default ``TikLambda2Method.PILOT_KNEE``,
    which reads the same knee off a few random experiments per grid value at a
    fraction of the cost; the full profile costs one selection per grid value
    and returns every probe's selection (``Lambda2PowerRuleResult``). The pick
    depends on y and lies outside the fixed-lambda_2 FDR guarantee that each
    probe carries.

    Without ``grid`` the default grid of the pilot rule is used,
    ``default_lambda2_grid()``: 0, 1, 3, 10, 30, 100, 300, 1000, 3000.
    ``grid`` is stated for columns of unit norm. Under ``ScalingMode.ZSCORE``
    the probes run at its values times n - 1, as the pilot rule does, and the
    pick is reported on that working scale (``grid_scale`` of the result holds
    the factor).

    ``X`` and ``y`` are converted to Fortran-ordered float64 arrays once (no
    copy when they already are) and accessed zero-copy: every probe normalizes
    X in place and hands it back, so X is unchanged on return. ``tik_control``
    is required; its ``lambda_2`` is overwritten per probe. ``trex_control``
    is taken as for ``TRexTikhonovSelector``: without it the base control
    nested in ``tik_control`` applies, and a control that is passed must name
    the solver TIKNET. A negative seed would draw fresh dummies per probe and
    is refused.
    """
    if tik_control is None:
        raise ValueError(
            "tik_control is required (a TRexTikhonovControlParameter; it "
            "carries the penalty matrix tikhonov_K or the penalty structure).")
    if seed < 0:
        raise ValueError("seed must be non-negative: the probes share their dummies.")
    X = np.asfortranarray(X, dtype=np.float64)
    y = np.asfortranarray(y, dtype=np.float64)
    grid = default_lambda2_grid() if grid is None else [float(v) for v in grid]
    return _tikhonov_lambda2_power_rule(X, y, tFDR, tik_control, trex_control, int(seed),
                                        grid, float(knee_tolerance), verbose)


def gvs_lambda2_power_rule(X: np.ndarray, y: np.ndarray, tFDR: float = 0.1,
                           gvs_control: TRexGVSControlParameter | None = None,
                           trex_control: TRexControlParameter | None = None,
                           seed: int = 0,
                           grid: "list[float] | np.ndarray | None" = None,
                           knee_tolerance: float = 1.0,
                           verbose: bool = False) -> Lambda2PowerRuleResult:
    """
    The knee of the full lambda_2 profile of the GVS selector.

    Runs one FDR-controlled selection of ``TRexGVSSelector`` per value of
    ``grid`` for ``lambda_2``, the ridge part of the elastic net or of the
    informed elastic net, with the same seed, so the probes share their
    dummies and differ by the penalty weight alone. The pick is the smallest
    lambda_2 whose selection is within ``knee_tolerance`` selections of the
    largest (0 reads the largest selection, ties toward the smallest
    lambda_2). The result holds every probe's selection
    (``Lambda2PowerRuleResult``). The pick depends on y and lies outside the
    fixed-lambda_2 FDR guarantee that each probe carries.

    Without ``grid`` the default grid is used, ``default_lambda2_grid()``.
    ``grid`` is stated for columns of unit norm. Under ``ScalingMode.ZSCORE``
    the probes run at its values times n - 1, and the pick is reported on that
    working scale (``grid_scale`` of the result holds the factor).

    ``gvs_control`` and ``trex_control`` are taken as by ``TRexGVSSelector``,
    and ``lambda_2`` of ``gvs_control`` is replaced per probe. ``X`` and ``y``
    are converted to Fortran-ordered float64 arrays once (no copy when they
    already are): every probe normalizes X in place and hands it back, so X is
    unchanged on return. A negative seed would draw fresh dummies per probe
    and is refused.
    """
    if gvs_control is None:
        gvs_control = TRexGVSControlParameter()
    if trex_control is None:
        trex_control = TRexControlParameter()
    if seed < 0:
        raise ValueError("seed must be non-negative: the probes share their dummies.")
    X = np.asfortranarray(X, dtype=np.float64)
    y = np.asfortranarray(y, dtype=np.float64)
    grid = default_lambda2_grid() if grid is None else [float(v) for v in grid]
    return _gvs_lambda2_power_rule(X, y, tFDR, gvs_control, trex_control, int(seed),
                                   grid, float(knee_tolerance), verbose)


@_with_core_notes(PyTRexScreeningSelector)
class TRexScreeningSelector(TRexSelector):
    """
    T-Rex screening selector (T-Rex+Screening): fast FDR-controlled variable
    screening in one pass of the random experiments, until one dummy is
    included (T = 1), with as many dummies as variables (L = p).

    Parameters
    ----------
    X : numpy.ndarray of shape (n, p)
        Design matrix, used as described for ``TRexSelector``.
    y : numpy.ndarray of shape (n,)
        Response vector.
    screen_control : ScreenTRexControlParameter, optional
        Control parameters of the screening. The defaults apply when it is
        omitted.
    trex_control : TRexControlParameter, optional
        Control parameters of the T-Rex selector. The defaults apply when it
        is omitted.
    seed : int, default -1
        Random seed. A negative value gives a run that is not reproducible.
    verbose : bool, default True
        Whether to print progress.
    """
    def __init__(self, X: np.ndarray, y: np.ndarray,
                 screen_control: ScreenTRexControlParameter | None = None,
                 trex_control: TRexControlParameter | None = None,
                 seed: int = -1, verbose: bool = True):

        self.X = np.asfortranarray(X, dtype=np.float64)
        self.y = np.asfortranarray(y, dtype=np.float64)

        if screen_control is None:
            screen_control = ScreenTRexControlParameter()
        if trex_control is None:
            trex_control = TRexControlParameter()

        # The core validates the L-loop strategy. The control is passed on as it is.
        _validate_core_solver_type(trex_control)

        self._selector: PyTRexScreeningSelector = PyTRexScreeningSelector(self.X, self.y,
                                                                          screen_control,
                                                                          trex_control,
                                                                          seed, verbose)

    def select(self) -> ScreenTRexSelectionResult:
        """
        Run the selector. A selector runs once: a second call raises
        ``RuntimeError``, see ``TRexSelector.select``.

        Returns
        -------
        ScreenTRexSelectionResult
            The selection with the screening-specific fields.
        """
        self._selector.select()
        return self._selector.getScreenResult()


@_with_core_notes(PyBiobankScreeningSelector)
class TRexBiobankScreeningSelector:
    """
    Biobank screening with the T-Rex selector: per phenotype the fast screening
    (T-Rex+Screening), and the T-Rex selector as fallback when no FDR estimate
    of the screening lies in the accepted range.

    ``select()`` returns one result for a response vector and a list of
    results, one per column, for a response matrix.

    The inner selectors standardize X in place and restore it, so X must be
    writable. Under ``ScalingMode.NONE`` they write nothing into X and accept an
    array that is not writable, for example the ``to_numpy()`` view of a
    ``MemoryMappedMatrix`` opened with ``AccessMode.ReadOnly``. A response
    matrix Y is only read, so it may be such an array under every scaling mode.

    Parameters
    ----------
    X : numpy.ndarray of shape (n, p)
        Design matrix.
    Y : numpy.ndarray of shape (n,) or (n, m)
        Response vector of one phenotype, or a matrix with one column per
        phenotype.
    bio_ctrl : BiobankScreenTRexControl, optional
        Control parameters of the biobank screening. The defaults apply when
        it is omitted.
    seed : int, default -1
        Random seed. A negative value gives a run that is not reproducible.
        With several phenotypes, phenotype i (0-based) runs with seed + i.
    verbose : bool, default False
        Whether to print progress.
    """

    def __init__(self, X: np.ndarray, Y,
                 bio_ctrl: BiobankScreenTRexControl | None = None,
                 seed: int = -1, verbose: bool = False):

        self.X = np.asfortranarray(X, dtype=np.float64)
        if bio_ctrl is None:
            bio_ctrl = BiobankScreenTRexControl()

        _validate_core_solver_type(bio_ctrl.trex_screen_ctrl.trex_ctrl)

        Y = np.asarray(Y, dtype=np.float64)
        self._is_mD = Y.ndim > 1
        self._Y = np.asfortranarray(Y)
        self._selector = PyBiobankScreeningSelector(self.X, self._Y,
                                                    bio_ctrl,
                                                    seed, verbose)

    def select(self) -> BiobankScreenTRexResult | list[BiobankScreenTRexResult]:
        """
        Run Biobank screening.

        Returns
        -------
        BiobankScreenTRexResult
            For a 1-D response (single phenotype).
        list[BiobankScreenTRexResult]
            For a m-D response (multiple phenotypes).
        """
        if self._is_mD:
            return self._selector.screenPhenotypes()
        return self._selector.screenPhenotype()


@_with_core_notes(PyTRexSPCASelector)
class TRexSPCASelector:
    """
    T-Rex sparse PCA selector (T-Rex+SPCA): sparse principal components whose
    non-zero loadings are selected with FDR control. Each ordinary principal
    component is the response of a T-Rex+GVS run with the elastic net.

    X is used without a copy when it is a Fortran-contiguous float64 array, for
    example the ``to_numpy()`` view of a ``MemoryMappedMatrix``. The selector
    centers X in place during ``select()`` and restores it on return. The
    centering happens under every scaling mode, so an array that is not
    writable, for example the view of a file opened with
    ``AccessMode.ReadOnly``, raises ``ValueError`` in ``select()``. Any other
    array is copied at construction, and the array of the caller is not
    touched.

    Parameters
    ----------
    X : numpy.ndarray of shape (n, p)
        Design matrix.
    spca_ctrl : TRexSPCAControlParameter, optional
        Control parameters: the loading assembly, the elastic net solver and
        the control of the T-Rex+GVS runs (``gvs_ctrl``). The defaults apply
        when it is omitted. A control that the selector refuses
        (``GVSType.IEN``, a rule of the informed elastic net, a permutation
        strategy) raises ``ValueError`` in ``select()``.
        ``gvs_ctrl.lambda_2 = 0`` gives a ``UserWarning``: it is the T-LASSO
        without the ridge part, not the automatic choice.
    seed : int, default -1
        Random seed. A negative value gives a run that is not reproducible.
        Each component runs with its own seed derived from it.
    verbose : bool, default False
        Whether to print one line per component. The T-Rex+GVS runs stay
        quiet.
    """

    def __init__(self,
                 X: np.ndarray,
                 spca_ctrl: TRexSPCAControlParameter | None = None,
                 seed: int = -1,
                 verbose: bool = False):
        self.X = np.asfortranarray(X, dtype=np.float64)
        self._ctrl = spca_ctrl if spca_ctrl is not None else TRexSPCAControlParameter()
        self._seed = seed
        self._verbose = verbose
        self._result: TRexSPCAResult | None = None

        # lambda_2 == 0 is the T-LASSO without the ridge part, not the
        # automatic choice (a negative lambda_2). The warning is for a user who
        # set 0 and expected the cross-validation.
        if self._ctrl.gvs_ctrl.lambda_2 == 0.0:
            warnings.warn(
                "spca_ctrl.gvs_ctrl.lambda_2 = 0 is the degenerate no-ridge "
                "case (pure T-LASSO), not automatic cross-validation. Use "
                "lambda_2 < 0 (e.g. -1) to auto-determine lambda_2 via CV.",
                UserWarning,
                stacklevel=2,
            )

    def select(self, M: int, tFDR: float) -> "TRexSPCASelector":
        """
        Run the T-Rex sparse PCA selection.

        Parameters
        ----------
        M : int
            Number of sparse principal components, 1 to min(n, p).
        tFDR : float
            Target FDR level of each component, in [0, 1].

        Returns
        -------
        TRexSPCASelector
            The selector itself. The results are in its properties ``Z``,
            ``V``, ``active_sets``, ``adjusted_ev`` and ``cumulative_ev``.
        """
        # The core validates M, tFDR and the control.
        selector = PyTRexSPCASelector(self.X, M, tFDR, self._ctrl, self._seed, self._verbose)
        selector.select()
        self._result = selector.getResult()
        return self

    @property
    def Z(self) -> np.ndarray:
        """Scores of the sparse components (n x M): column m is the centered X
        times column m of ``V``."""
        if self._result is None:
            raise RuntimeError("Call select() first.")
        return self._result.Z

    @property
    def V(self) -> np.ndarray:
        """Sparse loadings (p x M), each column of unit L2 norm. In ActiveSet
        mode the non-zero entries of column m are the variables of
        ``active_sets[m]``. In Thresholded mode they are the largest ordinary
        PCA loadings, as many as ``active_sets[m]`` has entries, which can be
        other variables."""
        if self._result is None:
            raise RuntimeError("Call select() first.")
        return self._result.V

    @property
    def active_sets(self) -> list[list[int]]:
        """Selected variables per component (list of M index lists, 0-based),
        the result of the T-Rex+GVS runs. The FDR control holds for these sets.
        In Thresholded mode they can differ from the non-zero entries of
        ``V``."""
        if self._result is None:
            raise RuntimeError("Call select() first.")
        return self._result.active_sets

    @property
    def adjusted_ev(self) -> np.ndarray:
        """Adjusted explained variance per component (length M):
        r_mm^2 / (n - 1), with r_mm the m-th diagonal element of R in the QR
        decomposition Z = QR. The adjustment accounts for sparse components
        that are not orthogonal."""
        if self._result is None:
            raise RuntimeError("Call select() first.")
        return self._result.adjusted_ev

    @property
    def cumulative_ev(self) -> np.ndarray:
        """Cumulative adjusted explained variance as a fraction of the total
        variance of the centered X (length M), between 0 and 1. It is not the
        PEV of Definition 1 of the paper, whose denominator (signal EV plus
        mixed EV) needs the true support."""
        if self._result is None:
            raise RuntimeError("Call select() first.")
        return self._result.cumulative_ev


__all__ = [
    "LLoopStrategy",
    "ScalingMode",
    "normalize_design",
    "SolverTypeForTRex",
    "SolverHyperparameters",
    "DummyDistribution",
    "TRexControlParameter",
    "SelectionResult",
    "TRexSelector",
    "TRexDASelector",
    "TRexDAControlParameter",
    "DAMethod",
    "BTSelectionMode",
    "DASelectionResult",
    "TRexGVSSelector",
    "TRexGVSControlParameter",
    "GVSType",
    "ENSolverType",
    "LambdaSelectionMethod",
    "GVSSelectionResult",
    "TRexScreeningSelector",
    "ScreenTRexControlParameter",
    "ScreenTRexMethod",
    "ScreenTRexSelectionResult",
    "TRexBiobankScreeningSelector",
    "BiobankScreenTRexControl",
    "BiobankScreenTRexResult",
    "SPCAMode",
    "TRexSPCAControlParameter",
    "TRexSPCAResult",
    "TRexSPCASelector",
    "Lambda2CVDiagnostics",
    "TikLambda2Method",
    "PilotStatistic",
    "TikPenaltyStructure",
    "Lambda2PickInfo",
    "Lambda2PilotResult",
    "TRexTikhonovControlParameter",
    "TRexTikhonovSelector",
    "gamma_to_k",
    "build_group_centering_k",
    "Lambda2Probe",
    "Lambda2PowerRuleResult",
    "tikhonov_lambda2_power_rule",
    "gvs_lambda2_power_rule",
    "default_lambda2_grid",
    "DummyClusters",
    # Submodules
    "ml_methods",
    "tsolvers",
    "trex_selector_methods",
    "utils",
]

"""
Augmented-dummy (AD) T-solvers: the classical family running on [X, D]
"""
from __future__ import annotations
import typing
from . import afs_based
from . import lars_based
from . import omp_based
__all__: list[str] = ['L2', 'ScalingMode', 'SolverOptions', 'ZSCORE', 'afs_based', 'lars_based', 'omp_based']
class ScalingMode:
    """
    Column scaling convention applied by the solver after centering.
    
    Members:
    
      L2 : Divide each centered column by its L2 norm (unit L2 norm).
    
      ZSCORE : Divide each centered column by its sample SD (n-1).
    """
    L2: typing.ClassVar[ScalingMode]  # value = <ScalingMode.L2: 0>
    ZSCORE: typing.ClassVar[ScalingMode]  # value = <ScalingMode.ZSCORE: 1>
    __members__: typing.ClassVar[dict[str, ScalingMode]]  # value = {'L2': <ScalingMode.L2: 0>, 'ZSCORE': <ScalingMode.ZSCORE: 1>}
    @typing.overload
    def __eq__(self, other: ScalingMode) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: ScalingMode) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class SolverOptions:
    """
    Every T-solver input except the data X, D and y. Each solver reads the fields it needs and ignores the rest; the defaults are the solver-level defaults. Construct with keywords (`SolverOptions(lambda2=0.5, normalize=False)`) or set attributes.
    """
    @typing.overload
    def __init__(self) -> None:
        """
        All defaults.
        """
    @typing.overload
    def __init__(self, **kwargs) -> None:
        """
        Keyword constructor: every keyword must name a SolverOptions field.
        """
    def __repr__(self) -> str:
        ...
    @property
    def exch_tie_alpha(self) -> float:
        """
        Exchangeable-tie band width of the greedy solvers (TOMP, TOOLS, TGP, TACGP, TMP, TNCGMP and TAFS); 0 (default) = off. The path solvers ignore it.
        """
    @exch_tie_alpha.setter
    def exch_tie_alpha(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def exch_tie_floor(self) -> float:
        """
        Minimum |correlation| for exchangeable-tie candidates in [0, 1). Default 0.5.
        """
    @exch_tie_floor.setter
    def exch_tie_floor(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def intercept(self) -> bool:
        """
        Center X and D in place (a writer, as for normalize) and the solver's internal copy of y (an intercept is fitted). False when the caller pre-centered. Default True.
        """
    @intercept.setter
    def intercept(self, arg0: bool) -> None:
        ...
    @property
    def lambda2(self) -> float:
        """
        L2 (ridge / group-ridge / Tikhonov) penalty weight, >= 0; read by the elastic-net, IEN and Tikhonov solvers. Default 0.1, as inside the selectors; 0 means no quadratic penalty.
        """
    @lambda2.setter
    def lambda2(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def n_threads(self) -> int:
        """
        Thread budget of the solver's own passes over the design (correlation refreshes, step-length search, standardization). Default 1 (serial); a whole number of at least 1. A pass splits only above about a million entries and never inside a parallel region of the caller. The T-Rex selectors set it from max_inner_threads.
        """
    @n_threads.setter
    def n_threads(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def ncgmp_variant(self) -> omp_based.NCGMPVariant:
        """
        TNCGMP variant (NCGMPVariant). Default LineSearch.
        """
    @ncgmp_variant.setter
    def ncgmp_variant(self, arg0: omp_based.NCGMPVariant) -> None:
        ...
    @property
    def normalize(self) -> bool:
        """
        Scale each column of X (and D) to unit L2 norm (or unit SD, see scaling_mode), IN PLACE through the views (the solver becomes a writer and claims the buffers until restore() or destruction). False when the caller pre-scaled. Default True.
        """
    @normalize.setter
    def normalize(self, arg0: bool) -> None:
        ...
    @property
    def precentered(self) -> bool:
        """
        The caller centered X, D and y (intercept=False on centered data, as in the T-Rex selectors). The centered columns span at most n - 1 dimensions, so the solver counts the intercept as with intercept=True: the active set is capped at n - 1 and the degrees of freedom and Cp include the intercept. The data are not touched. Redundant when intercept is True. Default False.
        """
    @precentered.setter
    def precentered(self, arg0: bool) -> None:
        ...
    @property
    def rho_afs(self) -> float:
        """
        TAFS shrinkage step in (0, 1]. Default 0.3, as inside the selectors; 1 is forward stepwise.
        """
    @rho_afs.setter
    def rho_afs(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def scaling_mode(self) -> ScalingMode:
        """
        Column scaling applied when normalize is True: ScalingMode.L2 (default) or ZSCORE.
        """
    @scaling_mode.setter
    def scaling_mode(self, arg0: ScalingMode) -> None:
        ...
    @property
    def tenet_aug_use_lars(self) -> bool:
        """
        TENETAug / TIENETAug: pure-LARS inner path that never drops variables (True) instead of the LASSO inner path (False, default).
        """
    @tenet_aug_use_lars.setter
    def tenet_aug_use_lars(self, arg0: bool) -> None:
        ...
    @property
    def tie_seed(self) -> int:
        """
        Seed of the tie-break engine. Default 0: a fixed seed, so a solver run is reproducible also where exact ties occur; < 0 seeds from std::random_device (nondeterministic).
        """
    @tie_seed.setter
    def tie_seed(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def tikhonov_validated(self) -> bool:
        """
        TIKNET: True says that the Tikhonov matrix was checked by the caller (finite entries, symmetric, positive semi-definite), and the solver skips its own check. For a caller that builds many solvers on one matrix, as the T-Rex selector does. Default False.
        """
    @tikhonov_validated.setter
    def tikhonov_validated(self, arg0: bool) -> None:
        ...
    @property
    def tol(self) -> float:
        """
        Numerical tolerance of the solver's internal comparisons. Default: machine epsilon. A
        selector sets its own tolerance on the solvers it runs (1e-6 by default).
        """
    @tol.setter
    def tol(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def verbose(self) -> bool:
        """
        Print status and diagnostics. Default False.
        """
    @verbose.setter
    def verbose(self, arg0: bool) -> None:
        ...
L2: ScalingMode  # value = <ScalingMode.L2: 0>
ZSCORE: ScalingMode  # value = <ScalingMode.ZSCORE: 1>

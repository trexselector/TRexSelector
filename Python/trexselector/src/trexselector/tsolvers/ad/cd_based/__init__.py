"""Coordinate-descent T-solvers.

Deprecated: the family will be removed in the next minor release of
trexselector, and constructing one of its classes raises ``FutureWarning``.
Use the homotopy solvers of ``trexselector.tsolvers.ad.lars_based``:
``TLASSO_Solver``, ``TENET_Solver``, ``TIENET_Solver`` and ``TIKNET_Solver``.
"""
from ...._core.tsolvers.ad.cd_based import (
    TCDLASSO_Solver,
    TCDENET_Solver,
    TCDIENET_Solver,
    TCDTIKNET_Solver,
)

__all__ = [
    "TCDLASSO_Solver",
    "TCDENET_Solver",
    "TCDIENET_Solver",
    "TCDTIKNET_Solver",
]

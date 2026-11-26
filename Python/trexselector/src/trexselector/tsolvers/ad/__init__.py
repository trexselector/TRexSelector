"""Augmented-dummy (AD) T-solvers: the classical family running on ``[X, D]``.

``SolverOptions`` and ``ScalingMode`` live here rather than one level up
because every user of them is an AD solver (the same placement rule as
``cpp/src/tsolvers/ad``).
"""
from ..._core.tsolvers.ad import ScalingMode, SolverOptions
from . import lars_based
from . import omp_based
from . import afs_based

__all__ = [
    "ScalingMode",
    "SolverOptions",
    "lars_based",
    "omp_based",
    "afs_based",
]

"""
Tests for a design that is not writable: the view of a memory-mapped file opened
read-only is accepted where nothing is written into it and refused elsewhere.
"""
import hashlib

import numpy as np
import pytest
import trexselector as tsn
from trexselector.trex_selector_methods import ad
from trexselector.utils import AccessMode, MemoryMappedMatrix, numpy_to_memmap


def _control(mode, K=10):
    ctrl = tsn.TRexControlParameter()
    ctrl.K = K
    ctrl.scaling_mode = mode
    return ctrl


def _md5(path):
    with open(path, "rb") as handle:
        return hashlib.md5(handle.read()).hexdigest()


@pytest.fixture
def readonly_design(tmp_path):
    """A pre-normalized design, written to a file and reopened read-only."""
    rng = np.random.default_rng(5)
    n, p = 120, 40
    X = rng.standard_normal((n, p))
    y = 5.0 * X[:, 1] - 3.0 * X[:, 2] + 2.0 * X[:, 3] + rng.standard_normal(n)
    d = tsn.normalize_design(X, y)
    X_pre = np.asfortranarray(d["X"])
    filename = str(tmp_path / "design.bin")
    numpy_to_memmap(filename, X_pre)
    mapped = MemoryMappedMatrix(filename, n, p, AccessMode.ReadOnly)
    return {"X": X_pre, "y": np.asarray(d["y"]), "view": mapped.to_numpy(), "file": filename, "p": p}


def test_readonly_view_serves_every_selector_under_none(readonly_design):
    d = readonly_design
    view, y = d["view"], d["y"]
    md5 = _md5(d["file"])
    none = _control(tsn.ScalingMode.NONE)
    assert not view.flags.writeable

    # The run on the file equals the run in memory.
    ref = tsn.TRexSelector(d["X"].copy(order="F"), y, tFDR=0.2, trex_control=none, seed=11,
                           verbose=False).select()
    sel = tsn.TRexSelector(view, y, tFDR=0.2, trex_control=none, seed=11, verbose=False)
    assert np.shares_memory(sel.X, view)  # no copy on the way in
    res = sel.select()
    np.testing.assert_array_equal(res.selected_var, ref.selected_var)
    np.testing.assert_array_equal(res.phi_prime, ref.phi_prime)
    assert np.flatnonzero(ref.selected_var).size > 0

    tsn.TRexDASelector(view, y, tFDR=0.2, trex_control=none, seed=11, verbose=False).select()
    tsn.TRexGVSSelector(view, y, tFDR=0.2, trex_control=none, seed=11, verbose=False).select()
    tik_ctrl = tsn.TRexTikhonovControlParameter()
    tik_ctrl.tikhonov_K = np.eye(d["p"])
    tik_ctrl.lambda_2 = 0.5
    tik_ctrl.trex_ctrl.scaling_mode = tsn.ScalingMode.NONE   # the nested base control of the run
    tsn.TRexTikhonovSelector(view, y, tFDR=0.2, tik_control=tik_ctrl, seed=11, verbose=False).select()
    tsn.TRexScreeningSelector(view, y, trex_control=none, seed=11, verbose=False).select()

    bio_ctrl = tsn.BiobankScreenTRexControl()
    bio_ctrl.trex_screen_ctrl.trex_ctrl.K = 10
    bio_ctrl.trex_screen_ctrl.trex_ctrl.scaling_mode = tsn.ScalingMode.NONE
    tsn.TRexBiobankScreeningSelector(view, y, bio_ctrl=bio_ctrl, seed=11).select()
    Y = np.column_stack([y, y[::-1]])
    results = tsn.TRexBiobankScreeningSelector(view, Y, bio_ctrl=bio_ctrl, seed=11).select()
    assert len(results) == 2

    # The file is as it was.
    np.testing.assert_array_equal(view, d["X"])
    assert _md5(d["file"]) == md5


@pytest.mark.parametrize("mode", [tsn.ScalingMode.L2, tsn.ScalingMode.ZSCORE])
def test_readonly_view_is_refused_where_the_design_is_written(readonly_design, mode):
    d = readonly_design
    view, y = d["view"], d["y"]
    md5 = _md5(d["file"])
    ctrl = _control(mode)

    # The message names the way out.
    with pytest.raises(ValueError, match="not writable.*ScalingMode.NONE"):
        tsn.TRexSelector(view, y, trex_control=ctrl, verbose=False)
    with pytest.raises(ValueError, match="not writable"):
        tsn.TRexDASelector(view, y, trex_control=ctrl, verbose=False)
    with pytest.raises(ValueError, match="not writable"):
        tsn.TRexGVSSelector(view, y, trex_control=ctrl, verbose=False)
    tik_ctrl = tsn.TRexTikhonovControlParameter()
    tik_ctrl.tikhonov_K = np.eye(d["p"])
    tik_ctrl.trex_ctrl.scaling_mode = mode                   # the nested base control of the run
    with pytest.raises(ValueError, match="not writable"):
        tsn.TRexTikhonovSelector(view, y, tik_control=tik_ctrl, verbose=False)
    with pytest.raises(ValueError, match="not writable"):
        tsn.TRexScreeningSelector(view, y, trex_control=ctrl, verbose=False)

    bio_ctrl = tsn.BiobankScreenTRexControl()
    bio_ctrl.trex_screen_ctrl.trex_ctrl.scaling_mode = mode
    with pytest.raises(ValueError, match="not writable"):
        tsn.TRexBiobankScreeningSelector(view, y, bio_ctrl=bio_ctrl)
    with pytest.raises(ValueError, match="not writable"):
        tsn.TRexBiobankScreeningSelector(view, np.column_stack([y, y]), bio_ctrl=bio_ctrl)

    assert _md5(d["file"]) == md5


def test_sparse_pca_refuses_a_readonly_view(readonly_design):
    d = readonly_design
    md5 = _md5(d["file"])
    # The sparse PCA selector centers the design in every run.
    with pytest.raises(ValueError, match="not writable.*centers the design"):
        tsn.TRexSPCASelector(d["view"], seed=5).select(M=2, tFDR=0.2)
    assert _md5(d["file"]) == md5


def test_compiled_selectors_take_the_readonly_view(readonly_design):
    """The compiled classes under trex_selector_methods.ad follow the same rule."""
    d = readonly_design
    view, y = d["view"], d["y"]
    sel = ad.TRexSelector(view, y, 0.2, _control(tsn.ScalingMode.NONE), 11, False)
    ref = ad.TRexSelector(d["X"].copy(order="F"), y, 0.2, _control(tsn.ScalingMode.NONE), 11, False)
    np.testing.assert_array_equal(sel.select().selected_var, ref.select().selected_var)
    with pytest.raises(ValueError, match="not writable"):
        ad.TRexSelector(view, y, 0.2, _control(tsn.ScalingMode.L2), 11, False)


def test_any_array_that_is_not_writable_follows_the_rule():
    """The rule reads the array's write flag, whatever the source of the array."""
    rng = np.random.default_rng(5)
    X = rng.standard_normal((120, 40))
    y = 5.0 * X[:, 1] - 3.0 * X[:, 2] + rng.standard_normal(120)
    d = tsn.normalize_design(X, y)
    X_pre = np.asfortranarray(d["X"])
    X0 = X_pre.copy(order="F")
    y_pre = np.asarray(d["y"])
    X_pre.setflags(write=False)

    ref = tsn.TRexSelector(X0.copy(order="F"), y_pre, tFDR=0.2,
                           trex_control=_control(tsn.ScalingMode.NONE), seed=11, verbose=False).select()
    res = tsn.TRexSelector(X_pre, y_pre, tFDR=0.2, trex_control=_control(tsn.ScalingMode.NONE), seed=11,
                           verbose=False).select()
    np.testing.assert_array_equal(res.selected_var, ref.selected_var)
    np.testing.assert_array_equal(X_pre, X0)
    with pytest.raises(ValueError, match="not writable"):
        tsn.TRexSelector(X_pre, y_pre, trex_control=_control(tsn.ScalingMode.L2), verbose=False)


def test_biobank_reads_a_response_matrix_that_is_not_writable(tmp_path):
    """The biobank screening only reads Y, so Y may be the view of a file opened read-only."""
    rng = np.random.default_rng(73)
    n, p, q = 80, 40, 4
    X = np.asfortranarray(rng.standard_normal((n, p)))
    Y = np.asfortranarray(np.column_stack(
        [4.0 * X[:, k] - 3.0 * X[:, k + 10] + rng.standard_normal(n) for k in range(q)]))
    filename = str(tmp_path / "responses.bin")
    numpy_to_memmap(filename, Y)
    Y_view = MemoryMappedMatrix(filename, n, q, AccessMode.ReadOnly).to_numpy()
    md5 = _md5(filename)
    assert not Y_view.flags.writeable

    ctrl = tsn.BiobankScreenTRexControl()
    ctrl.trex_screen_ctrl.R_boot = 100
    ctrl.trex_screen_ctrl.trex_ctrl.K = 10

    def run(design, response, control=ctrl):
        return tsn.TRexBiobankScreeningSelector(design, response, bio_ctrl=control, seed=9).select()

    expected = run(X.copy(order="F"), Y.copy(order="F"))
    results = run(X.copy(order="F"), Y_view)   # the default scaling, X writable
    assert len(results) == q
    for res, ref in zip(results, expected):
        assert res.selected_indices == ref.selected_indices
        assert res.method_used == ref.method_used

    # X and Y both not writable, under NONE
    d = tsn.normalize_design(X, Y[:, 0])
    X_pre = np.asfortranarray(d["X"])
    X_pre.setflags(write=False)
    none = tsn.BiobankScreenTRexControl()
    none.trex_screen_ctrl.R_boot = 100
    none.trex_screen_ctrl.trex_ctrl.K = 10
    none.trex_screen_ctrl.trex_ctrl.scaling_mode = tsn.ScalingMode.NONE
    assert len(run(X_pre, Y_view, none)) == q

    # A response vector that is not writable stays refused, as before.
    y_view = Y_view[:, 0]
    with pytest.raises(TypeError):
        run(X.copy(order="F"), y_view)

    np.testing.assert_array_equal(Y_view, Y)
    assert _md5(filename) == md5

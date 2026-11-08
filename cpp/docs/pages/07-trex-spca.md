# T-Rex SPCA: sparse principal components under FDR control {#trex_spca}

[TOC]

## Sparse loadings with a guarantee

Principal component analysis (PCA) finds the directions of largest variance, and every variable
contributes to every loading vector. Sparse PCA asks for loadings with few non-zero entries, so that
a component can be read as a small set of variables. The T-Rex sparse PCA (T-Rex SPCA)
\cite machkour2024spca decides which variables enter a loading vector with the T-Rex selector. For
each of the \f$M\f$ components, ordinary PCA gives the component as a response, the T-Rex selector
selects the variables for it at the target false discovery rate (FDR), and a sparse loading vector
is assembled from the selected set. The selected set of every component is FDR-controlled.

## Data with two sparse factors

The headers and names of this page:

\snippet{trimleft} trex_spca.cpp includes

A small helper prints the selected set of every component:

\snippet{trimleft} trex_spca.cpp helpers

Two latent factors, the first loading on variables 0 to 9 and the second on variables 30 to 39
(0-based), plus noise.

\snippet{trimleft} trex_spca.cpp data

## Running the selector

`TRexSPCA(X, M, tFDR, spca_ctrl, seed, verbose)` takes the design, the number of components `M`, the
target FDR and the control, and `select()` extracts the components. The elastic-net penalty of the
inner selector is chosen by cross-validation by default. `X` is centered in place during `select()`
and restored on return. The base control of the inner selector, with the number of random
experiments `K`, is nested in the group control of the TRexSPCAControlParameter.

\snippet{trimleft} trex_spca.cpp run
\output{trex_spca.run.txt}

The results are the scores `Z` (\f$n \times M\f$), the sparse loadings `V` (\f$p \times M\f$), the
selected sets, and the explained variance.

\snippet{trimleft} trex_spca.cpp fields
\output{trex_spca.fields.txt}

The non-zero entries of the two loading vectors and the share of each loading's squared norm that
falls on its true support. The true supports are variables 0 to 9 and 30 to 39.

\snippet{trimleft} trex_spca.cpp loadings
\output{trex_spca.loadings.txt}

\image html trex_spca_loadings.png "Figure 1: The two sparse loading vectors. The true supports are variables 0 to 9 and 30 to 39." width=90%

Sparse components are not orthogonal in general, so the explained variance is adjusted for the
non-orthogonality by a QR decomposition of the scores. `adjusted_ev` is the marginal adjusted
variance per component and `cumulative_ev` the cumulative fraction of the total variance of the
centered design.

## The control

TRexSPCAControlParameter sets the inner machinery.

- `mode`: `SPCAMode::ActiveSet` (default) assembles the loading by a ridge regression of the
  component on the selected variables with penalty `lambda2_ridge_loadings`, `SPCAMode::Thresholded`
  keeps the top \f$|A|\f$ ordinary loadings and normalizes them.
- `gvs_ctrl.lambda_2` and `gvs_ctrl.lambda2_method`: the elastic-net penalty of the inner group
  selector, fixed or chosen by ridge cross-validation (`CV_1SE_CCD` by default when `lambda_2 < 0`).
  A value of 0 is the degenerate lasso case and gives a warning.
- `en_solver`: `ENSolverType::TENET_AUG` (default) or `TENET`, two forms of the same elastic net.
- `gvs_ctrl.trex_ctrl.K` and `gvs_ctrl.trex_ctrl.max_dummy_multiplier`: the random experiments and
  the dummy cap of the inner selector.

T-Rex SPCA is defined on the elastic net, and `gvs_ctrl.gvs_type = GVSType::IEN` is rejected.

\snippet{trimleft} trex_spca.cpp thresholded
\output{trex_spca.thresholded.txt}

## Memory-mapped data

`X` may be the `getMap()` view of a MemoryMappedMatrix opened in `AccessMode::ReadWrite` (see
\ref memory_mapping "the memory-mapping page"). It is centered in the file during `select()` and
restored. Assigning a matrix to the map writes it into the file.

\snippet{trimleft} trex_spca.cpp mmap
\output{trex_spca.mmap.txt}

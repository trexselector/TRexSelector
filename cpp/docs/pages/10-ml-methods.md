# The ML methods {#ml_methods}

[TOC]

## The building blocks

The machine learning (ML) methods are the numerical building blocks the selectors use: scalers,
ridge and elastic-net regression with cross-validation (the \f$\lambda_2\f$ tuners), the singular
value decomposition (SVD) and principal component analysis (PCA) behind the sparse PCA selector, and
the hierarchical clustering behind the dependency-aware and the grouped selectors. They are exposed
as classes and functions of their own in the namespace `trex::ml_methods` under `src/ml_methods/`,
built on Eigen \cite guennebaud2010eigen. The headers and names of this page:

\snippet{trimleft} ml_methods.cpp includes

The data are drawn with a fixed seed into an Eigen matrix. The in-place methods take an
`Eigen::Map` of it and share its memory without a copy.

\snippet{trimleft} ml_methods.cpp data

## Scalers

\ref trex::ml_methods::scaler_methods::ZScoreScaler "ZScoreScaler" centers and divides by the sample
standard deviation, \ref trex::ml_methods::scaler_methods::LpNormScaler "LpNormScaler" centers and
divides by the L1 or L2 norm. Both work in place on a mapped double matrix (a memory-mapped one
included) and can undo their transformation.

\snippet{trimleft} ml_methods.cpp zscore
\output{ml_methods.zscore.txt}

`fit_transform_inplace()` does both in one call, `get_dropped_indices()` names constant columns, and
`save()` / `load()` persist the statistics.

\snippet{trimleft} ml_methods.cpp lpnorm
\output{ml_methods.lpnorm.txt}

## Ridge regression

\ref trex::ml_methods::ridge::RidgeSolver "RidgeSolver" solves
\f$\min_{\boldsymbol\beta} \|\boldsymbol y - \boldsymbol X \boldsymbol\beta\|_2^2
+ \lambda \|\boldsymbol\beta\|_2^2\f$ for one penalty \cite hoerl2000ridge, by the primal or the
dual Cholesky path depending on the shape of `X`. The class
\ref trex::ml_methods::model_selection::ridge_cv_svd "ridge_cv_svd" chooses the penalty
by \f$K\f$-fold cross-validation over a geometric grid and reports the minimum and the
one-standard-error rule.

\snippet{trimleft} ml_methods.cpp ridge
\output{ml_methods.ridge.txt}

The cross-validation curve is available through `lambdas()`, `cv_mse()` and `cv_sem()`. The
one-standard-error choice is the largest penalty whose error stays within one standard error of the
minimum, so its error lies between the minimum and the largest error on the grid.

\snippet{trimleft} ml_methods.cpp ridge-curve
\output{ml_methods.ridge-curve.txt}

\image html ml_methods_ridge.png "Figure 1: Cross-validation error of the ridge penalty with the minimum and the one-standard-error choice." width=80%

`ridge_cv_svd(true)` runs the dual (kernel) form on column blocks of `X` read in place, with memory
\f$O(n^2)\f$ instead of \f$O(np)\f$, the route for \f$p \gg n\f$.

## The elastic net

The class \ref trex::ml_methods::model_selection::enet_gaussian "enet_gaussian" fits the
elastic-net path \cite zou2005elasticnet by coordinate descent over a glmnet-style grid
\cite friedman2010glmnet. `alpha = 1` is the lasso, `alpha = 0` the ridge. The coefficients are
reported on the original scale.

\snippet{trimleft} ml_methods.cpp enet
\output{ml_methods.enet.txt}

The class \ref trex::ml_methods::model_selection::enet_cv_ccd "enet_cv_ccd" selects the penalty by
cross-validation on a shared grid. The strong rules \cite tibshirani2012strong speed up the sweeps.

\snippet{trimleft} ml_methods.cpp enet-cv
\output{ml_methods.enet-cv.txt}

Both classes take `low_memory = true` in their constructor and then run on `X` itself, a
memory-mapped matrix included, without an allocation of the size of `X`.

## SVD and PCA

The `compute(X, M)` method of \ref trex::ml_methods::svd::SVDSolver "SVDSolver" returns the
top-\f$M\f$ truncated SVD as an \ref trex::ml_methods::svd::SVDResult "SVDResult" with `U`, `S` and
`V`, by a direct, a Gram or a randomized path depending on the shape of `X`.

\snippet{trimleft} ml_methods.cpp svd
\output{ml_methods.svd.txt}

\ref trex::ml_methods::pca::PCA "PCA" fits the top-\f$M\f$ components. It centers and normalizes
`X` in place at construction, the same rule as the T-solvers, `fit()` returns the scores, and
`restore()` undoes the preprocessing.

\snippet{trimleft} ml_methods.cpp pca
\output{ml_methods.pca.txt}

## Hierarchical clustering

The static `cluster()` function of `AgglomerativeClustering` in the namespace
\ref trex::ml_methods::clustering::hierarchical::agglomerative "agglomerative" runs
agglomerative clustering on the columns of a matrix with a linkage from the `LinkageMethod` enum
(Ward, average,
complete, single, WPGMA, median, centroid) and a metric from the `DistanceMetric` enum (Euclidean,
correlation, Manhattan, and two locality-sensitive hashing (LSH) variants of the correlation
distance for very wide data), both chosen at compile time through the template arguments and a
`DistancePolicy`. Single linkage runs the SLINK algorithm \cite sibson1973slink, the other linkages
the nearest-neighbour chain \cite benzecri1982nnchain and \cite murtagh2012overview with the
Lance-Williams update \cite lance1967linkage. The result is the chronological vector of merge steps,
each a `MergeStep` with the two cluster ids, the merge distance and the size of the new cluster. The
linkage matrix of the R and Python packages is built from it. `DendrogramUtils::cut_tree()` cuts the
merge steps into a given number of clusters. The engine clusters columns, so observations are
clustered through the transposed matrix.

\snippet{trimleft} ml_methods.cpp hac
\output{ml_methods.hac.txt}

To cluster variables instead of observations, cluster the design itself with the correlation
distance, which is what the grouped and the dependency-aware selectors do.

\snippet{trimleft} ml_methods.cpp hac-vars
\output{ml_methods.hac-vars.txt}

`use_mmap = true` keeps the intermediate distances in memory-mapped files, and a memory-mapped data
matrix is clustered in place.

## Examples

The [TRexSelector_Examples](https://github.com/trexselector/TRexSelector_Examples/tree/main)
repository holds runnable demos and Monte Carlo studies of the scalers, ridge and elastic-net
cross-validation, PCA, SVD and clustering, with their results and a README per demo, in
[R](https://github.com/trexselector/TRexSelector_Examples/tree/main/R/ml_methods),
[Python](https://github.com/trexselector/TRexSelector_Examples/tree/main/Python/ml_methods) and
[C++](https://github.com/trexselector/TRexSelector_Examples/tree/main/cpp/ml_methods).
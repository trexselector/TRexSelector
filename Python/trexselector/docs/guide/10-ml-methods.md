---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
---

# The ML methods

```{code-cell} python
:tags: [remove-cell]
# The pages build on at most two threads.
import trexselector.utils
trexselector.utils.set_num_threads(min(2, trexselector.utils.get_max_threads()))
```

## The building blocks

The machine learning (ML) methods are the numerical building blocks the selectors use: scalers,
ridge and elastic-net regression with cross-validation (the $\lambda_2$ tuners), the singular value
decomposition (SVD) and principal component analysis (PCA) behind the sparse PCA selector, and the
hierarchical clustering behind the dependency-aware and the grouped selectors. They are exposed as
classes and functions of their own under `trexselector.ml_methods`, built on Eigen
{cite}`guennebaud2010eigen`. The data are drawn with a fixed seed, as a Fortran-ordered array, which
the in-place methods share without a copy.

```{code-cell} python
import numpy as np

rng = np.random.default_rng(8)
n = 120
p = 30
X = np.asfortranarray(rng.standard_normal((n, p)))
beta = np.zeros(p)
beta[:5] = [2.0, -1.5, 0.0, 0.0, 1.0]
y = X @ beta + rng.standard_normal(n)
```

## Scalers

`ZScoreScaler` centers and divides by the sample standard deviation, `LpNormScaler` centers and
divides by the L1 or L2 norm. Both work in place on a float64 array (or a memory-mapped one) and can
undo their transformation.

```{code-cell} python
from trexselector.ml_methods import ZScoreScaler, LpNormScaler, NormType

z = ZScoreScaler()
z.fit(X)
print("centers:", np.round(z.get_centers()[:3], 3))
print("scales: ", np.round(z.get_scales()[:3], 3))
X_before = X.copy(order="F")
z.transform_inplace(X)
print("column means after:", np.round(X.mean(axis=0)[:3], 12))
print("column sd after:   ", np.round(X.std(axis=0, ddof=1)[:3], 12))
z.inverse_transform_inplace(X)
print("largest difference to the original:", np.abs(X - X_before).max())
```

`fit_transform_inplace()` does both in one call, `get_dropped_indices()` names constant columns, and
`save()` / `load()` persist the statistics.

```{code-cell} python
l2 = LpNormScaler(NormType.L2, center=True)
l2.fit_transform_inplace(X)
print(np.sqrt((X ** 2).sum(axis=0))[:3])
l2.inverse_transform_inplace(X)
```

## Ridge regression

`RidgeSolver` solves $\min_{\boldsymbol\beta} \|\boldsymbol y - \boldsymbol X \boldsymbol\beta\|_2^2
+ \lambda \|\boldsymbol\beta\|_2^2$ for one penalty {cite}`hoerl2000ridge`, by the primal or the
dual Cholesky path depending on the shape of `X`. `RidgeCV` chooses the penalty by $K$-fold
cross-validation over a geometric grid and reports the minimum and the one-standard-error rule.

```{code-cell} python
from trexselector.ml_methods import RidgeSolver, RidgeCV

print(np.round(RidgeSolver.solve(X, y, 1.0)[:5], 3))
cv = RidgeCV()
cv.fit(X, y, n_folds=5, n_lambda=50)
print("lambda_min:", round(cv.cv_min(), 4), " lambda_1se:", round(cv.cv_1se(), 4))
```

```{code-cell} python
:tags: [hide-input]
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(cv.lambdas(), cv.cv_mse(), label="CV error, 5 folds")
ax.axvline(cv.cv_min(), color="C1", linestyle="-", label="minimum")
ax.axvline(cv.cv_1se(), color="C1", linestyle="--", label="one standard error")
ax.legend()
ax.set_xscale("log")
ax.set_xlabel("Lambda")
ax.set_ylabel("CV error")
ax.set_title("Cross-validation error of the ridge penalty with the minimum and the one-standard-error choice")
plt.show()
```

`RidgeCV(low_memory=True)` runs the dual (kernel) form on column blocks of `X` read in place, with
memory $O(n^2)$ instead of $O(np)$, the route for $p \gg n$.

## The elastic net

`ElasticNet` fits the elastic-net path {cite}`zou2005elasticnet` by coordinate descent over a
glmnet-style grid {cite}`friedman2010glmnet`. `alpha = 1` is the lasso, `alpha = 0` the ridge. The
coefficients are reported on the original scale.

```{code-cell} python
from trexselector.ml_methods import ElasticNet, ElasticNetCV

en = ElasticNet()
en.fit(X, y, alpha=0.5, n_lambda=50)
B = en.coef()                     # the path can stop before n_lambda points
print("coefficient path:", B.shape)
print("grid points:     ", len(en.lambdas()))
print("converged:       ", en.converged())
print("last point, first five coefficients:", np.round(B[:5, -1], 3))
print("predictions of three rows at the last point:",
      np.round(en.predict(np.asfortranarray(X[:3]))[:, -1], 3))
```

`ElasticNetCV` selects the penalty by cross-validation on a shared grid. The strong rules
{cite}`tibshirani2012strong` speed up the sweeps.

```{code-cell} python
encv = ElasticNetCV()
encv.fit(X, y, alpha=0.5, n_folds=5, n_lambda=50)
print("lambda_min:", round(encv.cv_min(), 4), " lambda_1se:", round(encv.cv_1se(), 4))
```

Both classes accept `low_memory=True` and then run on `X` itself, a memory-mapped matrix included,
without an allocation of the size of `X`.

## SVD and PCA

`SVDSolver().compute(X, M)` returns the top-$M$ truncated SVD, by a direct, a Gram or a randomized
path depending on the shape of `X`.

```{code-cell} python
from trexselector.ml_methods import SVDSolver, PCA

svd_res = SVDSolver().compute(X, 3)
print("U:", svd_res.U.shape, " V:", svd_res.V.shape)
print("S:", np.round(svd_res.S, 3))
```

`PCA` fits the top-$M$ components. `fit()` centers and normalizes `X` in place, the same rule as the
T-solvers, and `restore()` undoes it.

```{code-cell} python
X_before = X.copy(order="F")
pca = PCA(center=True, normalize=True)
res = pca.fit(X, 3)
print("scores:  ", res.Z.shape)
print("loadings:", pca.get_loadings().shape)
print("explained variance:", np.round(pca.get_explained_variance(), 3))
print("two rows transformed:")
print(np.round(pca.transform(np.asfortranarray(X_before[:2])), 3))
pca.restore()
print("largest difference to the original:", np.abs(X - X_before).max())
```

## Hierarchical clustering

`agglomerative_cluster()` runs agglomerative clustering on the rows of a matrix with a linkage from
`LinkageMethod` (Ward, average, complete, single, WPGMA, median, centroid) and a metric from
`DistanceMetric` (Euclidean, correlation, Manhattan, and two locality-sensitive hashing (LSH)
variants of the correlation distance for very wide data). Single linkage runs the SLINK algorithm
{cite}`sibson1973slink`, the other linkages the nearest-neighbour chain
{cite}`benzecri1982nnchain,murtagh2012overview` with the Lance-Williams update
{cite}`lance1967linkage`. The result is a linkage matrix in the layout of SciPy, and `cut_tree()`
cuts it into a given number of clusters.

```{code-cell} python
from trexselector.ml_methods.clustering import (agglomerative_cluster, cut_tree, LinkageMethod,
                                                DistanceMetric)

rng9 = np.random.default_rng(9)
centers = np.repeat([-3.0, 0.0, 3.0], 20)
Z = np.asfortranarray(centers[:, None] + rng9.standard_normal((60, 4)))
linkage = agglomerative_cluster(Z, LinkageMethod.Ward, DistanceMetric.Euclidean)
print("linkage:", linkage.shape)
labels = np.array(cut_tree(linkage, num_orig_objs=60, num_clusters=3))
truth = np.repeat([0, 1, 2], 20)
table = np.array([[np.sum((labels == a) & (truth == b)) for b in range(3)] for a in np.unique(labels)])
print("found clusters (rows) against the true groups (columns):")
print(table)
```

To cluster variables instead of observations, cluster the transposed design with the correlation
distance, which is what the grouped and the dependency-aware selectors do.

```{code-cell} python
link_vars = agglomerative_cluster(np.asfortranarray(X.T), LinkageMethod.Single,
                                  DistanceMetric.Correlation)
cut_tree(link_vars, num_orig_objs=p, num_clusters=5)[:6]
```

`use_mmap=True` keeps the intermediate distances in memory-mapped files, and a memory-mapped `data`
matrix is clustered in place.

## Examples

The [TRexSelector_Examples](https://github.com/trexselector/TRexSelector_Examples/tree/main)
repository holds runnable demos and Monte Carlo studies of the scalers, ridge and elastic-net
cross-validation, PCA, SVD and clustering, with their results and a README per demo, in
[R](https://github.com/trexselector/TRexSelector_Examples/tree/main/R/ml_methods),
[Python](https://github.com/trexselector/TRexSelector_Examples/tree/main/Python/ml_methods) and
[C++](https://github.com/trexselector/TRexSelector_Examples/tree/main/cpp/ml_methods).
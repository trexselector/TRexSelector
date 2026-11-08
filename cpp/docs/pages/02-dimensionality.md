# Data dimensionality and FDR control {#dimensionality}

[TOC]

## The variable selection problem

A linear model relates a response \f$\boldsymbol{y} \in \mathbb{R}^n\f$ to a design
\f$\boldsymbol{X} \in \mathbb{R}^{n \times p}\f$ with \f$n\f$ observations and \f$p\f$ variables,

\f[
\boldsymbol{y} = \boldsymbol{X} \boldsymbol{\beta} + \boldsymbol{\varepsilon},
\f]

and in the settings of this guide only a few coefficients of \f$\boldsymbol{\beta}\f$ differ from
zero. The variables with a non-zero coefficient form the active set \f$\mathcal{A}\f$, the others are
null variables. Variable selection returns an estimate \f$\widehat{\mathcal{A}}\f$ of the active set.
A null variable in \f$\widehat{\mathcal{A}}\f$ is a false discovery, and an active variable missing
from it is a missed discovery. Prediction is a different goal. A model can predict well with many
null variables on board, while a selection with false discoveries sends a laboratory
\cite ioannidis2005false \cite macarthur2012falsepositives \cite uffelmann2021gwas
\cite tandyconnor2018dtc or an investment team \cite harvey2016crosssection after effects that are not
there.

## FDR and TPR

Let \f$1 \vee a = \max\{1, a\}\f$. The false discovery rate (FDR) and the true positive rate (TPR)
of a selection procedure are

\f[
\operatorname{FDR} = \mathbb{E}\left[\frac{|\widehat{\mathcal{A}} \setminus \mathcal{A}|}
{1 \vee |\widehat{\mathcal{A}}|}\right],
\qquad
\operatorname{TPR} = \mathbb{E}\left[\frac{|\mathcal{A} \cap \widehat{\mathcal{A}}|}
{1 \vee |\mathcal{A}|}\right].
\f]

The quantities inside the expectations are the false discovery proportion (FDP) and the true
positive proportion (TPP) of one selection, which the metrics of the library compute when the true
support is known, as in a simulation. Benjamini and Hochberg introduced the FDR as the error rate for
multiple hypothesis testing \cite benjamini1995fdr, and Benjamini reviewed its origins and its first
fifteen years in \cite benjamini2010discovering. The older family-wise error rate asks for no false
discovery at all and loses power as the number of hypotheses grows \cite holm1979sequential, while a
target FDR of 0.1 says that among the selected variables at most ten percent are expected to be false,
whether ten or ten thousand variables were screened. A good selection procedure keeps the FDR below the
target and the TPR as high as possible.

## Three regimes

What a procedure can do depends on the number of variables relative to the number of observations.

### Low dimensionality, n at least p

With at least as many observations as variables the least squares estimate exists, and every
coefficient has a standard error, a test statistic and a p-value. Multiple testing is classical
here. The Benjamini–Hochberg procedure on the \f$p\f$ p-values controls the FDR under the assumptions
of the model and the dependence conditions of the procedure \cite benjamini1995fdr
\cite benjamini2001dependency, and the T-Rex selector works here as
well.

### High dimensionality, p above n

With more variables than observations the least squares estimate does not exist, because
\f$\boldsymbol{X}^\top \boldsymbol{X}\f$ is singular, and the classical p-values go with it. Sparse
path methods like the lasso \cite tibshirani1996lasso and least angle regression \cite efron2004lars,
and greedy pursuits like orthogonal matching pursuit, build a model one variable at a time and work
for any \f$p\f$, but their stopping rules are tuned for prediction error, so they select without
controlling the FDR. The T-Rex selector wraps these solvers in its random experiments and supplies the
FDR control they lack, in finite samples and without p-values \cite machkour2025trex, and for strongly
dependent variables the dependency-aware variant T-Rex+DA restores the exchangeability its guarantee
rests on \cite machkour2025datrex.

### Highest dimensionality, p far above n

In genomics a design has hundreds of thousands or millions of variables and up to a hundred thousand
observations, and the method is the same while the computation is not. Each random experiment
enlarges the design by \f$L \cdot p\f$ dummy columns, \f$K\f$ experiments run, and the design itself
may not fit in memory. Screen-T-Rex runs the random experiments once as a fast first pass
\cite machkour2023screen, memory-mapped designs keep the data on disk \cite scheidt2023bigtrex, the
seeded dummy strategies store no dummy matrices, the greedy solvers are cheaper per step
\cite scheidt2026greedy, and the thread budgets of the library put the cores to work.

## The regimes in one map

\image html dimensionality_map.png "Figure 1: The three regimes in the plane of p and n." width=85%

The diagonal \f$p = n\f$ separates low from high dimensionality, and the dashed line
\f$p = 1000\,n\f$ marks where the computation starts to dominate, a rule of thumb rather than a
theorem. The markers show typical settings. The diagonal, not the size of \f$p\f$ alone, decides the
regime. A biobank with twenty thousand
participants and a hundred variables is a low-dimensional problem, while a laboratory study with
sixty samples and two thousand genes is a high-dimensional one. Large samples with many variables sit
near the diagonal and belong to either side, which is why the map replaces the four quadrants of
small and large samples with few and many variables.

## The T-Rex selector across the regimes

The selector is defined for any \f$n\f$ and \f$p\f$. One data-generating process with five active
variables of coefficient one at \f$n = 100\f$ and three values of \f$p\f$ walks through the regimes,
ten trials per regime with ten random experiments each, every trial with its own data and dummies.
Indices are 0-based.

\snippet{trimleft} dimensionality.cpp includes

\snippet{trimleft} dimensionality.cpp regimes
\output{dimensionality.regimes.txt}

The estimated FDR stays below the target in every regime, while the TPR falls as \f$p\f$ grows
relative to \f$n\f$. Dimensionality costs power, not the control of the false discoveries. In the
low-dimensional case the classical route exists next to the selector, and in the two others it does
not.

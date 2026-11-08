# The T-Rex selector {#trex_selector}

[TOC]

## FDR and TPR

Let \f$\widehat{\mathcal{A}} \subseteq \{1, \ldots, p\}\f$ be the set of selected variables and
\f$\mathcal{A}\f$ the set of true active variables. With \f$1 \vee a = \max\{1, a\}\f$, the false
discovery rate (FDR) and the true positive rate (TPR) are

\f[
\operatorname{FDR} = \mathbb{E}\left[\frac{|\widehat{\mathcal{A}} \setminus \mathcal{A}|}
{1 \vee |\widehat{\mathcal{A}}|}\right],
\qquad
\operatorname{TPR} = \mathbb{E}\left[\frac{|\mathcal{A} \cap \widehat{\mathcal{A}}|}
{1 \vee |\mathcal{A}|}\right].
\f]

The quantities inside the expectations are the false discovery proportion (FDP) and the true
positive proportion (TPP) of one selection. `compute_fdp()` and `compute_tpp()` of the namespace
`trex::utils::eval::rates` compute them from index sets, and their overloads from coefficient
vectors. The T-Rex selector controls the FDR at the target `tFDR` and maximizes the number of
selected variables under that constraint \cite machkour2025trex.

## The framework

\image html T-Rex_framework.png "Figure 1: The T-Rex framework." width=85%

The selector runs \f$K\f$ random experiments. In experiment \f$k\f$ the design \f$\boldsymbol{X}\f$
(\f$n \times p\f$) is enlarged by a dummy matrix with \f$L \cdot p\f$ columns drawn independently of
the data, a forward selection solver walks the solution path of the enlarged design, and the path is
terminated as soon as \f$T\f$ dummies have entered the model. The dummies are a negative control
group. A null variable and a dummy are exchangeable along the path, so the number of included dummies
estimates the number of included null variables.

\image html predictor_matrix_with_dummies.png "Figure 2: The enlarged design (predictors and dummies)." width=65%

The candidate sets of the \f$K\f$ experiments are fused by a vote. A variable is selected at voting
level \f$v\f$ when it appears in more than \f$v \cdot K\f$ of the candidate sets. The calibration
searches the voting level and the stopping time \f$T\f$ for the largest selection whose FDP estimate
stays below the target. The number of dummies is calibrated too. The L-loop starts at \f$L = 1\f$ and
raises \f$L\f$ until the FDP estimate at the reference voting level falls below the target or the cap
`max_dummy_multiplier` is reached.

## Running the selector

The headers and names of this page:

\snippet{trimleft} trex_selector.cpp includes

The data are drawn with a fixed seed, \f$n = 100\f$ observations, \f$p = 150\f$ variables and five
active variables. Eigen's column-major layout is the layout the core reads without a copy (see the
section on the design below). Indices are 0-based.

\snippet{trimleft} trex_selector.cpp data

\snippet{trimleft} trex_selector.cpp run
\output{trex_selector.run.txt}

### Seeds

`seed` seeds the dummy draws. With a fixed seed a run is reproducible on one platform. The dummies
come from the C++ standard library's random engines, whose Normal distribution differs between
standard library implementations, so a seeded run on macOS, Linux and Windows draws different
dummies and can select differently in the borderline cases. For a Monte Carlo study over many trials
pass `seed = -1` (the default): every trial then draws its own dummies from the hardware entropy
source, while a fixed seed would give every trial the same per-experiment seeds and correlate the
trials.

### The number of random experiments

`K` is the number of random experiments. The T-Rex paper uses \f$K = 20\f$, the default. A larger
\f$K\f$ makes the vote finer, since the voting grid has steps of \f$1/K\f$, at a proportional cost.

\snippet{trimleft} trex_selector.cpp voting-grid
\output{trex_selector.voting-grid.txt}

### Early stopping of the T-loop

The T-loop raises the stopping time \f$T\f$ one dummy at a time. `use_max_T_stop = true` (default)
bounds it, and `tloop_stagnation_stop` stops the loop when the FDP estimate has not improved for
`tloop_max_stagnant_steps` steps. The stagnation stop is a `std::optional<bool>` that defaults to
`std::nullopt`, which resolves to off for the path solvers (T-LARS, T-LASSO, the elastic nets, the CD
family) and to on for the greedy solvers.

## The result fields

\snippet{trimleft} trex_selector.cpp result-fields
\output{trex_selector.result-fields.txt}

`getFDPHatMat()` is the FDP estimate per stopping time (rows, up to `T_stop`) and voting level
(columns), `getRMat()` the number of selected variables per stopping time and voting level, and
`getPhiMat()` the relative occurrence of each variable over the \f$K\f$ experiments after each
stopping time. The selection sits where the estimate is below the target and the count is largest.
The rows are 0-based, so the row of the last stopping time is `T_stop - 1`.

\snippet{trimleft} trex_selector.cpp matrices
\output{trex_selector.matrices.txt}

\snippet{trimleft} trex_selector.cpp phi
\output{trex_selector.phi.txt}

\image html trex_selector_phi.png "Figure 3: Relative occurrences after the last stopping time. The true support stands out." width=80%

Relative occurrences after the last stopping time. The true support has occurrence one. The null
variables that reach the voting level are the false positives of the selection, and nearly all of
the others stay far below it.

`getPhiPrime()` is the calibrated occurrence that enters the FDP estimate (at most 1, possibly
negative), `getXMeans()`, `getXL2Norms()` and `getYMean()` are the centering and scaling constants the
selector applied, and `getNumDummies()` is the number of dummies of the last run.

## The L-loop strategies

`lloop_strategy` of TRexControlParameter says how the dummies of the L-loop are produced, with the
LLoopStrategy enum.

| Strategy | Dummies | Storage |
| --- | --- | --- |
| `STANDARD` (default) | fresh matrices at every L step | \f$K\f$ matrices |
| `HCONCAT` | the same \f$K\f$ matrices grow by \f$p\f$ columns per L step | \f$K\f$ matrices |
| `SEEDED` | the draws of `STANDARD`, re-derived from the seed at every step | none |
| `PERMUTATION` | permuted columns of the design | \f$K\f$ matrices |
| `PERMUTATION_SEEDED` | the permutation re-derived from the seed | none |
| `SKIPL` | no L-loop, the cap `max_dummy_multiplier` is used directly | \f$K\f$ matrices |

`SKIPL` skips the calibration of \f$L\f$ and runs with the cap at once. Under a large fixed \f$L\f$
without the stagnation stop its FDP estimate is anti-conservative. The expected null inflow per T step
shrinks with \f$L\f$, so the T-loop runs deep. There the same few null variables keep entering across
the random experiments, because \f$X\f$ and \f$y\f$ are the same in all of them, and they collect more
occurrence above the reference voting level than the estimate books as false. The total inflow
itself stays about right. In the Examples demo 04 (signal-to-noise ratio 5, \f$T\f$ about 145) the
strategy with twenty times \f$p\f$ dummies under T-LARS reached an FDR of 0.133 at a target of 0.1,
and 0.026 with the stagnation stop on. The selector warns from `max_dummy_multiplier = 10` on when
the stop is off, and the calibrated strategies are the safe choice.

The seeded strategies store nothing and suit designs where \f$K \cdot L \cdot p\f$ columns do not fit
in memory. `SEEDED` selects exactly as `STANDARD` for the same seed, and `PERMUTATION_SEEDED` exactly
as `PERMUTATION`. Fresh and grown dummies give the same FDR. A paired study of 500 trials in nine
cells (signal-to-noise ratio 0.5, 1 and 2 times target FDR 0.05, 0.1 and 0.2, \f$n = 300\f$,
\f$p = 1{,}000\f$, ten active variables, \f$K = 20\f$, T-LARS) found the FDR equal within Monte Carlo
error in every cell and the grown form ahead in TPR by 0.6 percentage points on average, because it
calibrates a larger \f$L\f$ and runs longer T-loops. The table lists the cells by signal-to-noise ratio
(SNR) and target FDR (tFDR).

| SNR | tFDR | FDR fresh | FDR grown | TPR fresh | TPR grown |
| --- | --- | --- | --- | --- | --- |
| 0.5 | 0.05 | 0.011 | 0.010 | 0.113 | 0.115 |
| 0.5 | 0.10 | 0.041 | 0.041 | 0.253 | 0.266 |
| 0.5 | 0.20 | 0.103 | 0.103 | 0.395 | 0.398 |
| 1.0 | 0.05 | 0.016 | 0.017 | 0.598 | 0.608 |
| 1.0 | 0.10 | 0.038 | 0.038 | 0.738 | 0.738 |
| 1.0 | 0.20 | 0.096 | 0.096 | 0.845 | 0.851 |
| 2.0 | 0.05 | 0.014 | 0.015 | 0.941 | 0.952 |
| 2.0 | 0.10 | 0.036 | 0.035 | 0.983 | 0.982 |
| 2.0 | 0.20 | 0.137 | 0.129 | 0.998 | 0.997 |

\snippet{trimleft} trex_selector.cpp lloop
\output{trex_selector.lloop.txt}

## Dummy distributions

The dummies are standard Normal by default. The `Distribution` struct of the namespace
`trex::utils::datageneration::dummygen` offers the uniform, Rademacher, Student t, Laplace, Gumbel,
triangle, uniform sphere, Mammen \cite mammen1993bootstrap, constrained sparse Rademacher and
logistic distributions through its named constructors, `Distribution::Normal()`,
`Distribution::StudentT(df)` and so on. The FDR guarantee rests on the exchangeability of dummies and
null variables along the path, which the standard Normal provides for the centered and scaled working
design.

\snippet{trimleft} trex_selector.cpp dummies
\output{trex_selector.dummies.txt}

## Scaling modes

The selector works on a centered and scaled design. `scaling_mode = ScalingMode::L2` (default)
centers `X` and `y` and divides each column of `X` by its L2 norm, `ScalingMode::ZSCORE` divides by
the sample standard deviation instead, and the dummies follow the same convention.
`ScalingMode::NONE` takes `X` and `y` as supplied, for a design that was pre-processed by the caller
or for a dictionary whose columns must not be centered (steering vectors, Fourier atoms). A design
normalized by hand the way the selector does it selects under `ScalingMode::NONE` exactly as the raw
design under `ScalingMode::L2`.

\snippet{trimleft} trex_selector.cpp scaling
\output{trex_selector.scaling.txt}

Non-finite entries and constant columns are rejected in every mode with an exception that names the
column.

## The design is used in place

The selector centers and scales `X` in place at construction and writes the original values back
when `select()` returns, so no copy of the design is made. Between construction and `select()` the
matrix holds the working design. `restore()` gives it back early when a constructed selector will not
run. The constructor takes an `Eigen::Map<Eigen::MatrixXd>` of the caller's column-major memory, so
a matrix of another layout or type is mapped or converted by the caller first. The page on working
without data copies has the full contract.

\snippet{trimleft} trex_selector.cpp inplace
\output{trex_selector.inplace.txt}

## FDR control in a small simulation

A Monte Carlo study over three target levels. Each trial draws its own data and its own dummies.

\snippet{trimleft} trex_selector.cpp mc
\output{trex_selector.mc.txt}

\image html trex_selector_mc.png "Figure 4: Estimated FDR against the target (left) and the TPR (right), 20 trials per level." width=90%

The estimated FDR stays below the target at every level, 20 trials per level. The simulations of
the T-Rex paper \cite machkour2025trex show the same over a far wider range of settings.

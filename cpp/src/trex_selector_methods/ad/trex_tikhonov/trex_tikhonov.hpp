// ===================================================================================
// trex_tikhonov.hpp
// ===================================================================================
#ifndef TREX_SELECTOR_METHODS_TREX_TIKHONOV_HPP
#define TREX_SELECTOR_METHODS_TREX_TIKHONOV_HPP
// ===================================================================================
/**
 * @file trex_tikhonov.hpp
 *
 * @brief T-Rex selector with a general Tikhonov penalty (T-Rex+Tikhonov).
 *
 * @details The method is described at the class TRexTikhonovSelector.
 *  Implementations are in trex_tikhonov.cpp.
 */
// ===================================================================================

// std includes
#include <cstddef>
#include <string>
#include <vector>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Base T-Rex
#include <trex_selector_methods/ad/trex_core/trex.hpp>

// Hierarchical clustering types (LinkageMethod for the dummy clustering)
#include <ml_methods/clustering/hierarchical/agglomerative/agglomerative_types.hpp>

// lambda_2 tuner health record
#include <ml_methods/model_selection/lambda2_cv_diagnostics.hpp>

// The pick of the power / pilot rules (selector-free)
#include <trex_selector_methods/ad/trex_tikhonov/tik_lambda2_pick.hpp>

// ===================================================================================

// Embedded into namespace trex::trex_selector_methods::ad::trex_tikhonov
namespace trex::trex_selector_methods::ad::trex_tikhonov {

// Local namespace aliases
namespace tc = trex::trex_selector_methods::ad::trex_core;
namespace disp = trex::trex_selector_methods::ad::utils::solver_dispatch;
namespace er = trex::trex_selector_methods::ad::utils::experiment_runner;
namespace hac = trex::ml_methods::clustering::hierarchical::agglomerative;
namespace ms = trex::ml_methods::model_selection;


// ===================================================================================
// Enums
// ===================================================================================

/**
 * @brief Rule that chooses the penalty weight `lambda_2` when no value is
 *        supplied (`lambda_2 < 0`). The pilot rule reads it off a few random
 *        experiments per value of `pilot_grid`. The four cross-validation
 *        rules tune on the normalized data before the random experiments.
 *        Every rule reports its pick on the working scale of the run (see
 *        the class documentation, "lambda_2").
 */
enum class TikLambda2Method {
    /** @brief IEN-geometry coordinate-descent CV (lambda1 profiled), 1SE. */
    CV_1SE_IEN_CCD,
    /** @brief IEN-geometry coordinate-descent CV (lambda1 profiled), min. */
    CV_MIN_IEN_CCD,
    /** @brief Generalized-Tikhonov SVD ridge CV (lambda1 = 0 backbone), 1SE. */
    CV_1SE_TIK_SVD,
    /** @brief Generalized-Tikhonov SVD ridge CV (lambda1 = 0 backbone), min. */
    CV_MIN_TIK_SVD,
    /** @brief Pilot rule (the default): over `pilot_grid`, the smallest
     *  lambda_2 whose pilot selection count is within
     *  `pilot_knee_tolerance` of the largest. */
    PILOT_KNEE
};

/** @brief The count the pilot rule picks on, per grid point. Both are
 *  always recorded in `Lambda2PilotResult`. */
enum class PilotStatistic {
    /** @brief The number of originals ahead of the T_pilot-th dummy in more
     *  than half of the pilot experiments — the T-Rex's own count at voting
     *  level 1/2, without the FDP gate. Zero under the global null. */
    MAJORITY,
    /** @brief The T-Rex selection rule applied to the pilot data (Phi',
     *  FDP-hat on a grid of step 1/K_pilot, the largest selection subject
     *  to FDP-hat <= tFDR over t <= T_pilot). Coarser at a small K_pilot. */
    CALIBRATED,
    /** @brief The MAJORITY count restricted to the groups of `prior_groups`
     *  in which at least `pilot_group_fraction` of the members are ahead:
     *  a group that enters as a whole contributes all of its members
     *  ahead, a group that enters as a minority contributes nothing. The
     *  false positives a group prior admits under crosstalk are scattered
     *  members of null groups correlated with an active one (a few per
     *  touched group); they inflate MAJORITY as lambda_2 grows and pull a
     *  count-based pick upward in exactly the trials that are false-
     *  positive-prone, whereas an active group crosses the fraction as a
     *  whole. Needs `prior_groups` (refused otherwise). Falls to zero when
     *  no group reaches the fraction (very low SNR): the pick is then the
     *  grid's first point. */
    GROUP_MAJORITY
};

/** @brief What the pilot of the last select() / runPilotOnly() did. */
struct Lambda2PilotResult {

    /** @brief A pilot executed (false: lambda_2 came from the caller or a
     *  CV tuner). */
    bool ran = false;

    /** @brief The lambda_2 values the pilot ran at, ascending and
     *  de-duplicated: `pilot_grid` times `grid_scale`. */
    std::vector<double> grid;

    /** @brief The factor between `pilot_grid` and `grid`: n - 1 under
     *  ScalingMode::ZSCORE, 1 otherwise. */
    double grid_scale = 1.0;

    /** @brief The count the pick was read off, per grid point: the count of
     *  `pilot_statistic` (MAJORITY, CALIBRATED or GROUP_MAJORITY). */
    std::vector<std::size_t> counts;

    /** @brief The majority count per grid point and pilot stopping time
     *  t = 1 .. T_pilot (`curve_majority[g][t - 1]`). */
    std::vector<std::vector<std::size_t>> curve_majority;

    /** @brief The mean number of originals ahead of the t-th dummy, per
     *  grid point and t (the column sums of phi_T_mat; diagnostic). */
    std::vector<std::vector<double>> curve_soft;

    /** @brief The calibrated pilot selection count per grid point. */
    std::vector<std::size_t> calibrated;

    /** @brief The group-majority count per grid point: the members ahead of
     *  the T_pilot-th dummy in the majority of experiments, summed over the
     *  groups in which at least `pilot_group_fraction` of the members are
     *  ahead. Zeros without `prior_groups`. */
    std::vector<std::size_t> group_majority;

    /** @brief The number of groups that reach `pilot_group_fraction`, per
     *  grid point. Zeros without `prior_groups`. */
    std::vector<std::size_t> groups_qualifying;

    /** @brief Wall-clock seconds per grid point. */
    std::vector<double> probe_seconds;

    /** @brief The pick diagnostics on `counts`. */
    Lambda2PickInfo pick;

    /** @brief Index of the pick into `grid`. */
    std::size_t index = 0;

    /** @brief The picked lambda_2, `grid[index]`, on the working scale. */
    double lambda2 = 0.0;

    /** @brief The pick sits on the first grid point. */
    bool at_lower_end = false;

    /** @brief The pick sits on the last grid point: widen the grid before
     *  trusting it. */
    bool at_upper_end = false;

    /** @brief Pilot experiments per grid point that ran: `pilot_K`, capped
     *  at K. */
    std::size_t K_pilot = 0;

    /** @brief Stopping time of the pilot paths that was used: `pilot_T`,
     *  capped at min(L_0, ceil(n / 2)). */
    std::size_t T_pilot = 0;

    /** @brief Number of pilot dummies, `pilot_L_multiplier` times p. */
    std::size_t L_0 = 0;

    /** @brief Seed tag under which the pilot dummies were drawn. */
    std::size_t dummy_tag = 0;

    /** @brief Wall-clock seconds of the whole pilot. */
    double seconds = 0.0;
};


/**
 * @brief The clusters of the variables in the dummy generation of a run.
 *        The fields have the names of the GVS result.
 */
struct DummyClusters {

    /** @brief Cluster of each variable (0-based, length p). Empty when the
     *  dummies are not drawn per cluster or no run was made. */
    Eigen::VectorXi groups_vec;

    /** @brief Number of clusters M (0 when there are none). */
    std::size_t max_clusters = 0;

    /** @brief Origin of the clusters: "prior" when `prior_groups` was
     *  supplied, else the linkage of the clustering, "single", "complete",
     *  "average" or "wpgma". Empty when there are no clusters. */
    std::string hc_method_used;
};


/**
 * @brief Structural form of the Tikhonov penalty (solver routing).
 */
enum class TikPenaltyStructure {
    /** @brief The supplied matrix `tikhonov_K` is the penalty. The solver
     *  works on the matrix, at a cost per update that grows with the
     *  non-zero entries of a column of K. */
    SPARSE_GENERAL,

    /** @brief K = blockdiag(I_{p_m} - 1_m 1_m^T / p_m) over `prior_groups`:
     *  the members of a group share a coefficient. The groups are required.
     *  `tikhonov_K` may be left empty, the selector builds the matrix for
     *  the cross-validation rules. With TCDTIKNET the solver runs on the
     *  O(1) group-centering hooks, the same penalty at a constant cost. */
    GROUP_CENTERING
};


// ===================================================================================
// Control Parameters
// ===================================================================================

/**
 * @brief Tikhonov-specific control parameters.
 */
struct TRexTikhonovControlParameter {

    /** @brief Structural form of the penalty (see TikPenaltyStructure).
     *  Default SPARSE_GENERAL (matrix-passing). GROUP_CENTERING requires
     *  non-empty `prior_groups` and routes the solver through the O(1)
     *  group-centering hooks. */
    TikPenaltyStructure penalty_structure = TikPenaltyStructure::SPARSE_GENERAL;

    /** @brief General Tikhonov matrix K = Gamma^T Gamma (p x p, sparse,
     *  PSD). REQUIRED for SPARSE_GENERAL: an empty matrix throws at
     *  construction. The constructor checks its content once (finite
     *  entries, symmetry, positive semi-definiteness), and the solvers of
     *  the selector take it as checked. Build from a regularization operator
     *  Gamma (m x p) via `TRexTikhonovSelector::gammaToK()`. Prefer
     *  diag(K) = 1 so the real columns carry the same ridge as the dummies
     *  (kappa_dummy = 1).
     *  OPTIONAL for GROUP_CENTERING: when empty, the selector builds
     *  K_GC = blockdiag(I - J/p_m) from `prior_groups` internally (the CV
     *  tuners consume it). A supplied matrix is compared with K_GC, and
     *  another matrix is refused at construction. The solver runs on the
     *  O(1) hooks either way. */
    Eigen::SparseMatrix<double> tikhonov_K;

    /** @brief User-supplied lambda_2.
     *  - `< 0` (default `-1.0`): not supplied, `lambda2_method` chooses it.
     *  - `== 0`: no penalty, the lasso path of the solver.
     *  - `> 0` : used as it is.
     *  The value must be finite. A supplied value refers to the scale of the
     *  standardized columns. Under ScalingMode::ZSCORE the columns have the
     *  squared norm n - 1, so the same shrinkage needs the value times n - 1.
     *  The rules report their pick on this scale. */
    double lambda_2 = -1.0;

    /** @brief Selection rule for auto-computing `lambda_2` when
     *  `lambda_2 < 0`. Default: PILOT_KNEE, the knee of the pilot profile;
     *  the CV tuners (well-posed in any p/n regime) are the alternatives. */
    TikLambda2Method lambda2_method = TikLambda2Method::PILOT_KNEE;

    /** @brief Number of folds for the lambda_2 CV, 2 to n. Default: 10. */
    int cv_n_folds = 10;

    /** @brief lambda_2 grid size for the CV_*_TIK_SVD analytic scan.
     *  Default: 100. */
    Eigen::Index cv_n_lambda = 100;

    /** @brief lambda_2 grid size for the CV_*_IEN_CCD 2D tuner. The scan is
     *  warm-started across the lambda_2 grid with glmnet-style path
     *  truncation decided by a full-data reference scan (cv.glmnet
     *  grid-once semantics), so cost grows sublinearly in this value.
     *  Default: 100 — resolution parity with the tuner's internal lambda_1
     *  path and glmnet's nlambda. */
    Eigen::Index cv_n_lambda2 = 100;

    /** @brief Length of the internal lambda_1 path of the CV_*_IEN_CCD 2D
     *  tuner. lambda_1 is a NUISANCE parameter there — profiled out, needed
     *  only to keep the contrast space well-posed — so it does not need
     *  glmnet's nlambda = 100 resolution. Cost in this value is U-shaped
     *  (coarser steps degrade the pathwise warm starts faster than the point
     *  count falls); 40 is the smallest length that still reproduces the
     *  100-point grid's lambda_2 picks exactly. Default: 40. */
    Eigen::Index cv_n_lambda1 = 40;

    /** @brief Support cap of the CV_*_IEN_CCD 2D tuner as a fraction of the
     *  smallest CV fold's training size n_train_min: the full-data
     *  reference scan stops each lambda_1 path at the first point with more
     *  than floor(cv_dfmax_ratio * n_train_min) nonzero coefficients
     *  (glmnet `dfmax`). `<= 0` disables the cap. Larger values allow larger
     *  supports at a higher tuner cost (past nnz ~ n_train the coordinate
     *  descent slows sharply). Must be finite. Default: 0.8. */
    double cv_dfmax_ratio = 0.8;

    /** @brief Seed for the deterministic CV fold-permutation RNG.
     *  - `-1` (default): derived from the T-Rex `seed` parameter.
     *  - `>= 0`: explicit fold seed; overrides `seed` for CV only. */
    int cv_seed = -1;

    // ----- Pilot rule (lambda2_method = PILOT_KNEE) -----

    /** @brief lambda_2 candidates, stated for columns of unit L2 norm.
     *  Under ScalingMode::ZSCORE the pilot runs at these values times n - 1,
     *  as the CV tuners return their pick, so one grid means the same
     *  shrinkage under every scaling mode. The pick is reported on the
     *  working scale and can be supplied as `lambda_2`. Sorted and
     *  de-duplicated internally. Keep 0 on it so that "no penalty" is an
     *  admissible answer: a prior that misfits the truth loses power
     *  monotonically and the rule then switches it off. Default:
     *  `default_lambda2_grid()`, 0, 1, 3, 10, 30, 100, 300, 1000, 3000. */
    std::vector<double> pilot_grid = default_lambda2_grid();

    /** @brief Pilot experiments per grid point (>= 2; capped at
     *  `trex_ctrl.K`). The majority count has resolution 1/pilot_K. */
    std::size_t pilot_K = 10;

    /** @brief Stopping time of the pilot paths (>= 1; capped at
     *  min(L_0, ceil(n / 2)) as the T-loop caps T). The pick reads the
     *  count at this t; the counts at every t <= pilot_T are recorded. */
    std::size_t pilot_T = 5;

    /** @brief Pilot dummies L_0 = pilot_L_multiplier * p (in
     *  [1, trex_ctrl.max_dummy_multiplier]). */
    std::size_t pilot_L_multiplier = 1;

    /** @brief Selections below the largest pilot count that still count as
     *  "on the plateau" (>= 0); 0 reads the largest count. */
    double pilot_knee_tolerance = 1.0;

    /** @brief The count the pick is read off. */
    PilotStatistic pilot_statistic = PilotStatistic::MAJORITY;

    /** @brief GROUP_MAJORITY: the share of a group's members that must be
     *  ahead for the group to count, in (0, 1]; 1/2 = at least half. */
    double pilot_group_fraction = 0.5;

    /** @brief false (default): the pilot draws its dummies from a reserved
     *  seed stream, independent of the main run's — conditional on the
     *  pick, the main run is a standard T-Rex whose dummies are independent
     *  of the pick, and its selection equals a fixed-lambda_2 run's at that
     *  value with the same seed. true: the pilot reuses the main strategy's
     *  first dummy block (HCONCAT / STANDARD / SKIPL / SEEDED; the
     *  PERMUTATION family is refused). */
    bool pilot_shared_dummies = false;

    /** @brief Dummy-coupling mode of the sparse-K penalty.
     *  Default: true (FOLDED) — dummy copy t joins its originating
     *  variable (t - p) mod p in the K-coupling, the exact generalization
     *  of GVS's layered group convention; dummy/null exchangeability then
     *  holds for arbitrary K (with K = sum_m 1_m 1_m^T / p_m the solver
     *  reproduces GVS-TCDIENET bit-exactly). Set to false
     *  (INDEPENDENT_RIDGE: decoupled unit-ridge dummy block) ONLY for
     *  validation / equivalence runs, and with solver_type TCDTIKNET, the
     *  one solver that has this mode — it preserves the K = I ≡ TCDENET
     *  collapse but deforms exchangeability whenever K has off-diagonals
     *  (conservative for positive off-diagonals such as group-mean K,
     *  anti-conservative for negative ones such as Laplacians).
     *  @deprecated Will be removed with TCDTIKNET and the coordinate descent
     *  family in the next minor release. */
    bool fold_dummy_coupling = true;

    /** @brief Cluster-aware dummy generation, as in TRexGVSSelector.
     *  Default: true. When true, dummies are drawn as per-cluster colored
     *  MVN layers on the clustering below via the DummyGenerator cluster
     *  mode — all strategies and memory mapping supported; only the
     *  PERMUTATION variants throw. Set to false ONLY for validation /
     *  equivalence runs against the canonical i.i.d. dummy path (e.g. the
     *  K = I ≡ TCDENET anchor test): with an informed K, i.i.d. dummies
     *  deform the dummy/null exchangeability the FDR calibration relies on. */
    bool cluster_dummies = true;

    /** @brief Groups of the variables: a 0-based group ID per variable
     *  (length p, IDs without gaps). They define the clusters of the
     *  dummies (with `cluster_dummies`), the penalty under GROUP_CENTERING
     *  and the count GROUP_MAJORITY of the pilot rule. Empty (default): the
     *  clusters of the dummies are found by clustering the correlations
     *  with `corr_max` and `hc_linkage`. */
    std::vector<Eigen::Index> prior_groups;

    /** @brief Correlation threshold of the automatic clustering (only when
     *  groups are not provided). The dendrogram is cut at the height
     *  1 - corr_max. With Single linkage (default) no two variables of
     *  different clusters correlate above corr_max. With Complete linkage
     *  every two variables of one cluster correlate at least at corr_max,
     *  and variables of different clusters can exceed it. Default 0.5, in
     *  [0, 1]. Ignored when `prior_groups` is set or
     *  `cluster_dummies == false`. */
    double corr_max = 0.5;

    /** @brief Linkage of the hierarchical clustering (default: Single).
     *  Supported: Single, Complete, Average (UPGMA), WPGMA (McQuitty).
     *  Ward is rejected, because its update is valid for squared Euclidean
     *  distances only and the correlation distance 1 - |cor| is not
     *  Euclidean. The centroid and median linkages are rejected as well,
     *  because they are not reducible. Ignored when `prior_groups` is set
     *  or `cluster_dummies == false`. */
    hac::LinkageMethod hc_linkage = hac::LinkageMethod::Single;

    /** @brief Base T-Rex algorithmic control parameters (nested).
     *  `solver_type` must be one of the general-K solvers: TIKNET (native
     *  LARS homotopy, the default, FOLDED coupling only) or TCDTIKNET (CCD).
     *  The CCD solver is not finally developed: once lambda_2 times K is
     *  large it can stop at its sweep budget (cd_max_sweeps) before it has
     *  converged. It then throws instead of continuing on an unconverged
     *  path. TCDTIKNET is deprecated and will be removed in the next minor
     *  release.
     *  With penalty_structure = GROUP_CENTERING, the TCDTIKNET dispatch
     *  routes to TCDIENET's O(1) group-centering hooks, so the solver's
     *  warnings (sweep budget, certification fallback) carry the tag
     *  [TCDIENET]. TIKNET consumes the internally built K_GC on its sparse
     *  path (the O(1) fast-path hooks are CCD-only). The default here
     *  overrides the base class's own default (TLARS) accordingly. All other
     *  base knobs (K, L-loop strategy, warm starts, memory mapping, CD
     *  hyperparameters, ...) apply unchanged. */
    tc::TRexControlParameter trex_ctrl = [] {
        tc::TRexControlParameter t;
        t.solver_type = disp::SolverTypeForTRex::TIKNET;
        return t;
    }();
};


// ===================================================================================
// TRexTikhonovSelector
// ===================================================================================

/**
 * @class TRexTikhonovSelector
 *
 * @brief T-Rex selector with a general Tikhonov penalty (T-Rex+Tikhonov):
 *        FDR-controlled variable selection with the penalty
 *        (lambda_2 / 2) * beta^T K beta, whose positive semi-definite matrix
 *        K carries the structure among the variables.
 *
 * @details
 *  The penalty. K = Gamma^T Gamma (p x p) for a regularization operator
 *  Gamma, see gammaToK(). Special cases:
 *    - K = I: the ridge part of the elastic net.
 *    - K = sum_m 1_m 1_m^T / p_m (group means): the penalty of the informed
 *      elastic net of TRexGVSSelector.
 *    - K = blockdiag(I - 1_m 1_m^T / p_m) (group centering,
 *      TikPenaltyStructure::GROUP_CENTERING): the members of a group share
 *      a coefficient.
 *    - a banded matrix or a graph Laplacian: neighbouring coefficients are
 *      smooth.
 *
 *  The solver. TIKNET, the homotopy solver, by default. TCDTIKNET, the
 *  coordinate-descent solver, is deprecated.
 *
 *  The dummies. By default they are drawn per cluster of correlated
 *  variables, as in TRexGVSSelector (`cluster_dummies`): the clusters are
 *  `prior_groups`, or come from a hierarchical clustering of the
 *  correlations (`corr_max`, `hc_linkage`). getDummyClusters() reports
 *  them. With an informed K, i.i.d. dummies would be less competitive than
 *  correlated null variables, which deforms the exchangeability the FDR
 *  calibration relies on. In the penalty every dummy takes the place of the
 *  variable it was drawn for (FOLDED coupling, `fold_dummy_coupling`), so
 *  dummies and null variables stay exchangeable for every K.
 *
 *  L-loop strategies. With cluster dummies STANDARD, HCONCAT, SKIPL and
 *  SEEDED are supported, and SEEDED equals STANDARD. PERMUTATION and
 *  PERMUTATION_SEEDED are refused, because a row permutation destroys the
 *  covariance of the dummies within a cluster. With
 *  `cluster_dummies = false` the two run as well.
 *
 *  lambda_2. A supplied value is used as it is. Otherwise `lambda2_method`
 *  chooses it before the L-loop, on the normalized data:
 *    - PILOT_KNEE (default), the pilot rule. For every value of
 *      `pilot_grid`, `pilot_K` random experiments run with L_0 = p dummies
 *      to the stopping time `pilot_T`, on one set of dummies, and the
 *      variables ahead of the dummies in more than half of the experiments
 *      are counted (`pilot_statistic`). The pick is the knee of this
 *      profile, the smallest lambda_2 whose count is within
 *      `pilot_knee_tolerance` of the largest. It costs about
 *      |grid| x pilot_K paths, a fraction of one selection per grid value.
 *      The pilot dummies come from a seed stream of their own, so the main
 *      run is a T-Rex run at a fixed lambda_2 whose dummies do not depend on
 *      the pick. The pick itself depends on y and lies outside the FDR
 *      guarantee of a fixed lambda_2. In simulations the pick cost nothing
 *      under the global null and at most 0.02 of true positive rate
 *      against the best fixed value under the alternatives. The record is
 *      in getLambda2PilotResult(), and tikhonov_lambda2_power_rule() runs
 *      the same rule on full selections.
 *    - CV_*_IEN_CCD: cross-validation with the lambda_1 path profiled out
 *      per lambda_2 candidate, well-posed for every p and n.
 *    - CV_*_TIK_SVD: analytic ridge cross-validation (lambda_1 = 0), exact
 *      for n > p and refused when the unpenalized null space interpolates
 *      the training folds.
 *  The scale. A supplied `lambda_2` refers to the scale of the standardized
 *  columns: under ScalingMode::ZSCORE the columns have the squared norm
 *  n - 1, so the same shrinkage needs the value times n - 1. `pilot_grid`
 *  is stated for columns of unit norm, and the pilot rule and the
 *  cross-validation rules apply the factor n - 1 themselves. Every rule
 *  reports its pick on the scale of a supplied value (getLambda2Used()).
 *
 *  A selector runs once: select() and runPilotOnly() exclude each other.
 *
 *  Architecture. Unlike TRexGVSSelector, which has its own loop over the
 *  random experiments, this class reuses the machinery of the base class
 *  (L-loop, T-loop, warm starts, the ExperimentRunner): buildRunnerConfig()
 *  hands the penalty matrix and the resolved lambda_2 to the solver
 *  dispatch, and the cluster model enters through the DummyGenerator
 *  (setClusterModel() in onSelectBegin()), which serves every strategy and
 *  the memory-mapped writers.
 *
 *  References:
 *  -----------
 *   A. N. Tikhonov, "Solution of Incorrectly Formulated Problems and the Regularization Method," Soviet Math.
 *     Dokl., vol. 4, pp. 1035-1038, 1963 (English translation of Dokl. Akad. Nauk SSSR 151 (1963), no. 3,
 *     501-504).
 *   Y. Park, L. Reichel, G. Rodriguez, and X. Yu, "Parameter Determination for Tikhonov Regularization
 *     Problems in General Form," J. Comput. Appl. Math., vol. 343, pp. 12-25, 2018,
 *     doi: 10.1016/j.cam.2018.04.049.
 *   A. J. Smola and R. Kondor, "Kernels and Regularization on Graphs," in Learning Theory and Kernel Machines
 *     (COLT/Kernel 2003), Lecture Notes in Artificial Intelligence vol. 2777, 2003, pp. 144-158.
 *   F. R. K. Chung, Spectral Graph Theory. CBMS Regional Conference Series in Mathematics, no. 92. American
 *     Mathematical Society, 1997.
 *   A. E. Hoerl and R. W. Kennard, "Ridge Regression: Biased Estimation for Nonorthogonal Problems,"
 *     Technometrics, vol. 42, no. 1, pp. 80-86, 2000 (Reprint of the 1970 article in the Special 40th
 *     Anniversary Issue).
 *   J. Machkour, M. Muma, and D. P. Palomar, "The Terminating-Random Experiments Selector: Fast
 *     High-Dimensional Variable Selection with False Discovery Rate Control," Signal Process., vol. 231,
 *     Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.
 */
class TRexTikhonovSelector : public tc::TRexSelector {

public:

    /**
     * @brief Construct a TRexTikhonovSelector.
     *
     * @param X    Design matrix (n x p), Eigen::Map, not copied. It is
     *             standardized in place at construction and restored when
     *             select() returns (the contract of the base class).
     * @param y    Response vector (n), Eigen::Map, copied internally.
     * @param tFDR Target FDR level (default: 0.1), in [0, 1].
     * @param trex_tik_ctrl Tikhonov control parameters. Under SPARSE_GENERAL
     *             they carry the p x p matrix `tikhonov_K`, under
     *             GROUP_CENTERING the groups.
     * @param seed Random seed. A negative value (default: -1) gives a run
     *             that is not reproducible.
     * @param verbose Enable verbose output (default: true).
     *
     * @throws std::invalid_argument on an empty or mis-dimensioned
     *         `tikhonov_K`, a `tikhonov_K` that is not finite, not symmetric
     *         or not positive semi-definite, a supplied `tikhonov_K` under
     *         GROUP_CENTERING that is not the centering matrix, a `lambda_2`
     *         that is not finite, a `solver_type` other than TCDTIKNET /
     *         TIKNET, malformed `prior_groups` / `corr_max`, an unsupported
     *         `hc_linkage`, wrong sizes of the cross-validation or of the
     *         pilot rule, or — with `cluster_dummies` on — a
     *         PERMUTATION-family L-loop strategy.
     */
    TRexTikhonovSelector(Eigen::Map<Eigen::MatrixXd>& X,
                         Eigen::Map<Eigen::VectorXd>& y,
                         double tFDR = 0.1,
                         TRexTikhonovControlParameter trex_tik_ctrl =
                             TRexTikhonovControlParameter(),
                         int seed = -1,
                         bool verbose = true);

    /** @brief Build K = Gamma^T Gamma from a regularization operator
     *  Gamma (m x p), pruned and compressed. */
    static Eigen::SparseMatrix<double> gammaToK(
        const Eigen::SparseMatrix<double>& gamma);

    /** @brief Build the within-group centering matrix
     *  K_GC = blockdiag(I_{p_m} - 1_m 1_m^T / p_m) from a 0-based
     *  contiguous group assignment (length p; singleton groups contribute
     *  nothing). The canonical GROUP_CENTERING penalty representation
     *  (used by the CV tuners; the solver runs on the O(1) hooks). */
    static Eigen::SparseMatrix<double> buildGroupCenteringK(
        const std::vector<Eigen::Index>& prior_groups);

    /** @brief The lambda_2 the last select() / runPilotOnly() used: the
     *  supplied value, or the pick of the pilot rule or of a
     *  cross-validation rule, on the scale of a supplied value. -1 before
     *  a run. */
    double getLambda2Used() const noexcept { return lambda2_; }

    /** @brief Health record of the lambda_2 resolution of the last
     *  `select()` call (`ran == false`: the value came from the caller). */
    const ms::lambda2_cv_diagnostics& getLambda2CVDiagnostics() const noexcept {
        return lambda2_cv_diag_;
    }

    /** @brief What the pilot rule of the last `select()` /
     *  `runPilotOnly()` did (`ran == false`: no pilot executed). */
    const Lambda2PilotResult& getLambda2PilotResult() const noexcept {
        return lambda2_pilot_;
    }

    /** @brief The clusters the dummies of the last `select()` /
     *  `runPilotOnly()` were drawn for: `prior_groups`, or the clusters the
     *  selector found by clustering the correlations. Without a run and with
     *  `cluster_dummies = false` the record is empty. */
    const DummyClusters& getDummyClusters() const noexcept {
        return dummy_clusters_;
    }

    /**
     * @brief Run ONLY the lambda_2 pilot (cluster model + pilot experiments
     *        + pick) and release X, without the main selection.
     *
     * @details For callers that want the pick alone (a study that scores it
     *   against selections it already has). A selector runs once: after
     *   `runPilotOnly()` a call of `select()` or of `runPilotOnly()` throws,
     *   and so does `runPilotOnly()` after `select()`.
     *   Requires `lambda_2 < 0` and the pilot rule as `lambda2_method`.
     *
     * @throws std::invalid_argument when the control does not request the
     *         pilot rule.
     * @throws std::runtime_error when the selector has run already or gave
     *         X back before.
     */
    const Lambda2PilotResult& runPilotOnly();

protected:

    /** @brief Builds the dummy cluster model (prior groups or
     *  correlation-HAC on the normalized X), installs it into the base
     *  DummyGenerator (cluster mode), and resolves lambda_2 (CV on the
     *  normalized X / centered y) before the L-loop runs. The base dummy
     *  hook and ExperimentRunner then run unchanged — the generator draws
     *  cluster layers for every strategy and the memory-mapped writers. */
    void onSelectBegin() override;

    /** @brief Base config plus the Tikhonov injection: sets
     *  `tikhonov_K` and overwrites `solver_params.lambda2` with the
     *  resolved penalty weight. */
    er::ExperimentRunnerConfig buildRunnerConfig(
        std::size_t num_dummies,
        std::size_t T_stop,
        bool use_warm_start,
        er::ExperimentStrategy strategy,
        std::size_t seed_factor = 0,
        std::size_t existing_cols_on_disk = 0) const override;

private:

    /** @brief lambda_2 resolution (user value or sparse-K CV tuners). */
    double computeLambda2() const;

    /** @brief The control requests the pilot rule (lambda_2 < 0 and
     *  lambda2_method = PILOT_KNEE). */
    bool pilotRequested() const noexcept;

    /** @brief The lambda_2 pilot: for every grid point, `pilot_K` random
     *  experiments with L_0 dummies to T_pilot on ONE stored dummy set,
     *  the counts of `PilotStatistic`, then the pick. Leaves the main run's
     *  dummy streams untouched (reserved tag, or the strategy's own first
     *  block when `pilot_shared_dummies`). Returns the picked lambda_2 and
     *  fills `lambda2_pilot_`. */
    double runLambda2Pilot();

    /** @brief The seed tag under which the main strategy draws its first
     *  dummy block (the `pilot_shared_dummies` option). */
    std::size_t pilotDummyTag() const;

    /** @brief The reserved seed tag of an independent pilot: the main run
     *  only uses tags 0 (prefix-stable streams; in cluster mode layer li
     *  draws from stream li + 1) and 1 .. max_dummy_multiplier. */
    static constexpr std::size_t kPilotDummyTag = std::size_t{1} << 30U;

    /** @brief Build the dummy cluster model: clusters from `prior_groups`
     *  or correlation-HAC on the normalized X_, then the per-cluster
     *  Cholesky factors (shared machinery, trex_cluster_dummies.hpp /
     *  trex_cluster_hac.hpp); installed into the base DummyGenerator via
     *  setClusterModel(). */
    void setupClusterModel();

    /** @brief Tikhonov control parameters (owns the K matrix the runner
     *  configs point into). */
    TRexTikhonovControlParameter trex_tik_ctrl_;

    /** @brief Resolved CV fold-permutation seed. */
    unsigned int resolved_cv_seed_{0};

    /** @brief Resolved penalty weight for the current select() run. */
    double lambda2_{-1.0};

    /** @brief Health record of the last `computeLambda2()` call. Mutable
     *  because `computeLambda2()` is const — it reports what it did
     *  alongside the value it returns. */
    mutable ms::lambda2_cv_diagnostics lambda2_cv_diag_;

    /** @brief What the last pilot did (`ran == false` when none ran). */
    Lambda2PilotResult lambda2_pilot_;

    /** @brief GROUP_CENTERING: group ids as VectorXi for the solver
     *  dispatch (owned; runner configs point into it). */
    Eigen::VectorXi groups_vec_i_{};

    // ----- Cluster-aware dummy state (cluster_dummies == true) -----

    /** @brief Per-cluster column indices into X (0-based). Length M.
     *  Kept for diagnostics; the Cholesky factors live inside the
     *  DummyGenerator cluster model. */
    std::vector<std::vector<Eigen::Index>> clusters_list_;

    /** @brief The clusters as reported by getDummyClusters(). */
    DummyClusters dummy_clusters_;
};

// ===================================================================================
} /* End of namespace trex::trex_selector_methods::ad::trex_tikhonov */
// ===================================================================================
#endif /* End of TREX_SELECTOR_METHODS_TREX_TIKHONOV_HPP */

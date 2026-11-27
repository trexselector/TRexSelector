// =============================================================================
// trex_spca.hpp
// ==============================================================================
#ifndef TREX_SELECTOR_METHODS_TREX_SPCA_HPP
#define TREX_SELECTOR_METHODS_TREX_SPCA_HPP
// ==============================================================================
/**
 * @file trex_spca.hpp
 *
 * @brief T-Rex sparse PCA selector (T-Rex+SPCA).
 *
 * @details
 *  The two methods of the paper: the T-Rex PCA (SPCAMode::ActiveSet) and the
 *  T-Rex Thresholded PCA (SPCAMode::Thresholded). A T-Rex+GVS run with the
 *  elastic net selects the variables of each component.
 *
 *  Implementations are in trex_spca.cpp.
 */
// ==============================================================================

// std includes
#include <vector>

// Eigen includes
#include <Eigen/Dense>

// Centering (stored as member for exception-safe restoration)
#include <ml_methods/scaler_methods/lp_norm_scaler.hpp>

// TRex core (TRexControl)
#include <trex_selector_methods/ad/trex_core/trex.hpp>

// TRex GVS (TRexGVSSelector, TRexGVSControl, GVSType)
#include <trex_selector_methods/ad/trex_gvs/trex_gvs.hpp>

// ==============================================================================

namespace trex::trex_selector_methods::ad::trex_spca {

// Namespace aliases
namespace tc         = trex::trex_selector_methods::ad::trex_core;
namespace tg         = trex::trex_selector_methods::ad::trex_gvs;
namespace std_scaler = trex::ml_methods::scaler_methods;

// ==============================================================================
// Enums
// ==============================================================================

/**
 * @brief Assembly of the sparse loading vector from the selected set A_m.
 */
enum class SPCAMode {
    /** @brief Ridge regression of the component z_m on the selected columns
     *  X_{A_m}, normalized to unit length. The T-Rex PCA of the paper. */
    ActiveSet,

    /** @brief The |A_m| largest ordinary PCA loadings in absolute value,
     *  normalized to unit length. The T-Rex Thresholded PCA of the paper. */
    Thresholded
};


// ==============================================================================
// Result Structure
// ==============================================================================

/**
 * @brief Result of the T-Rex sparse PCA selector.
 *
 * @details A component without a selected variable has a zero column in V
 *  and Z and the explained variance 0.
 */
struct TRexSPCAResult {
    /** @brief Scores of the sparse components (n x M): column m is the
     *  centered X times column m of V. */
    Eigen::MatrixXd Z;

    /** @brief Sparse loadings (p x M), each column of unit L2 norm.
     *
     *  In ActiveSet mode the non-zero entries of column m are the variables
     *  of active_sets[m]. In Thresholded mode they are the |active_sets[m]|
     *  largest ordinary PCA loadings, which can be other variables. */
    Eigen::MatrixXd V;

    /** @brief Selected variables per component (0-based), the result of the
     *  T-Rex+GVS runs. The FDR control holds for these sets. In Thresholded
     *  mode they can differ from the non-zero entries of V (see V). */
    std::vector<std::vector<Eigen::Index>> active_sets;

    /** @brief Adjusted explained variance per component (length M):
     *  r_mm^2 / (n - 1), with r_mm the m-th diagonal element of R in the QR
     *  decomposition Z = QR. The adjustment accounts for sparse components
     *  that are not orthogonal. */
    Eigen::VectorXd adjusted_ev;

    /** @brief Cumulative adjusted explained variance as a fraction of the
     *  total variance of the centered X (length M), between 0 and 1. It is
     *  not the PEV of Definition 1 of the paper, whose denominator (signal
     *  EV plus mixed EV) needs the true support. */
    Eigen::VectorXd cumulative_ev;
};


// ==============================================================================
// Control Parameters
// ==============================================================================

/**
 * @brief Control parameters of the T-Rex sparse PCA selector.
 */
struct TRexSPCAControl {

    // ====================================
    // Loading Assembly
    // ====================================

    /** @brief Assembly of the sparse loading vectors (default: ActiveSet). */
    SPCAMode mode = SPCAMode::ActiveSet;

    /** @brief Ridge penalty of the loading assembly in ActiveSet mode
     *  (default: 1e-6), in the regression
     *  `(X_A^T X_A + lambda2_ridge_loadings * I)^{-1} X_A^T z_m`.
     *  Not used in Thresholded mode.
     */
    double lambda2_ridge_loadings = 1e-6;


    // ====================================
    // T-Rex+GVS Runs
    // ====================================

    /** @brief Control of the T-Rex+GVS run of each component, with the base
     *  control in gvs_ctrl.trex_ctrl.
     *
     *  - gvs_ctrl.gvs_type must be EN (the default). The constructor refuses
     *    IEN, because the method is defined on the elastic net.
     *  - gvs_ctrl.en_solver and gvs_ctrl.trex_ctrl.solver_type have no
     *    effect: the selector replaces them by `en_solver` of this control.
     *  - gvs_ctrl.lambda_2: a negative value (default: -1) lets the rule
     *    gvs_ctrl.lambda2_method choose it, 0 is the T-LASSO without the
     *    ridge part, and a positive value is used as it is. A supplied value
     *    refers to the scale of the standardized columns. Under
     *    ScalingMode::ZSCORE the columns have the squared norm n - 1, so the
     *    same shrinkage needs the value times n - 1.
     *  - gvs_ctrl.lambda2_method: CV_1SE_CCD (default), CV_MIN_CCD,
     *    CV_1SE_SVD or CV_MIN_SVD. The constructor refuses the rules of the
     *    informed elastic net (CV_*_IEN_CCD and CV_*_TIK_SVD).
     *  - gvs_ctrl.trex_ctrl.lloop_strategy: see the class documentation.
     */
    tg::TRexGVSControl gvs_ctrl;

    /** @brief Elastic net solver of the T-Rex+GVS runs (default: TENET_AUG).
     *
     *  - TENET_AUG: the elastic net as a LASSO on row-augmented data, the
     *    form of the paper.
     *  - TENET: the elastic net on the Gram matrix, with the ridge part in
     *    the Cholesky update of X_A^T X_A + lambda_2 * I.
     *
     *  The two solve the same elastic net problem.
     */
    tg::ENSolverType en_solver = tg::ENSolverType::TENET_AUG;

    // The base control is gvs_ctrl.trex_ctrl. A second copy at this level
    // could differ from it.
};

/** @brief Old name of TRexSPCAControl.
 *  @deprecated Use TRexSPCAControl. */
using TRexSPCAControlParameter [[deprecated("Use TRexSPCAControl")]] = TRexSPCAControl;


// ==============================================================================
// TRexSPCASelector
// ==============================================================================

/**
 * @class TRexSPCASelector
 *
 * @brief T-Rex sparse PCA selector (T-Rex+SPCA): sparse principal components
 *        whose non-zero loadings are selected with FDR control. Each ordinary
 *        principal component is the response of a T-Rex+GVS run with the
 *        elastic net.
 *
 * @details
 *  Steps:
 *    1. Center X and compute its first M ordinary principal components
 *       z_1, ..., z_M (a PCA of the covariance, the columns are not scaled).
 *    2. For each component m, a T-Rex+GVS run with the elastic net and the
 *       response z_m selects the set A_m at the target FDR.
 *    3. Assemble the sparse loading vector v_m from A_m (see SPCAMode), and
 *       the scores of the sparse component, the centered X times v_m.
 *    4. Adjust the explained variance for components that are not
 *       orthogonal, by the QR decomposition of the scores.
 *
 *  A component without a selected variable has a zero column in V and Z and
 *  the explained variance 0.
 *
 *  L-loop strategies: STANDARD, HCONCAT, SKIPL and SEEDED of the T-Rex+GVS
 *  runs are supported. PERMUTATION and PERMUTATION_SEEDED are refused at
 *  construction, as by the T-Rex+GVS selector.
 *
 *  Data contract:
 *    - X is centered in place during select() and restored on return, also
 *      when select() ends with an exception. The centering happens under
 *      every scaling mode.
 *    - X is not owned. The Eigen::Map must stay valid for the lifetime of
 *      the object.
 *
 *  References:
 *  -----------
 *   J. Machkour, A. Breloy, M. Muma, D. P. Palomar, and F. Pascal, "Sparse PCA with False Discovery Rate
 *     Controlled Variable Selection," in Proc. IEEE 49th Int. Conf. Acoust., Speech, Signal Process.
 *     (ICASSP), 2024, pp. 9716-9720.
 *   J. Machkour, M. Muma, and D. P. Palomar, "The Terminating-Random Experiments Selector: Fast
 *     High-Dimensional Variable Selection with False Discovery Rate Control," Signal Process., vol. 231,
 *     Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.
 */
class TRexSPCASelector {
public:

    // ============================================================
    // Constructor / Destructor
    // ============================================================

    /**
     * @brief Construct the T-Rex sparse PCA selector.
     *
     * @param X          Design matrix (n x p), Eigen::Map, not copied. It is
     *                   centered in place during select() and restored on return.
     * @param M          Number of sparse principal components, 1 to min(n, p).
     * @param tFDR       Target FDR level of each component, in [0, 1].
     * @param spca_ctrl  Control parameters: the loading assembly, the elastic
     *                   net solver and the control of the T-Rex+GVS runs.
     * @param seed       Random seed. A negative value (default: -1) gives a run
     *                   that is not reproducible. Each component runs with its
     *                   own seed derived from it.
     * @param verbose    Print one line per component (default: false). The
     *                   T-Rex+GVS runs stay quiet.
     *
     * @throws std::invalid_argument if M is not in 1 to min(n, p), if tFDR is
     *         not in [0, 1], if gvs_ctrl.gvs_type is IEN, if
     *         gvs_ctrl.lambda2_method is a rule of the informed elastic net,
     *         or if the L-loop strategy is PERMUTATION or PERMUTATION_SEEDED.
     */
    TRexSPCASelector(Eigen::Map<Eigen::MatrixXd>& X,
             Eigen::Index M,
             double tFDR,
             TRexSPCAControl spca_ctrl = TRexSPCAControl(),
             int seed = -1,
             bool verbose = false);

    /** @brief Destructor. Restores X if select() was interrupted by an exception. */
    ~TRexSPCASelector();

    /** @brief Deleted copy constructor. */
    TRexSPCASelector(const TRexSPCASelector&) = delete;

    /** @brief Deleted copy assignment operator. */
    TRexSPCASelector& operator=(const TRexSPCASelector&) = delete;

    /** @brief Deleted move constructor. */
    TRexSPCASelector(TRexSPCASelector&&) = delete;

    /** @brief Deleted move assignment operator. */
    TRexSPCASelector& operator=(TRexSPCASelector&&) = delete;


    // ============================================================
    // Main Method
    // ============================================================

    /**
     * @brief Run the T-Rex sparse PCA selection.
     *
     * @return Const reference to the stored TRexSPCAResult (sparse PCs,
     *         loadings, active sets, explained variance); the bindings wrap
     *         it once instead of copying Z and V twice. A caller
     *         that writes `auto res = spca.select();` still gets its own copy.
     */
    const TRexSPCAResult& select();

    /**
     * @brief Retrieve the result from the last select() call.
     *
     * @return Const reference to TRexSPCAResult.
     */
    const TRexSPCAResult& getResult() const noexcept { return result_; }


private:

    // ============================================================
    // Data Members
    // ============================================================

    /** @brief Non-owning pointer to design matrix X (n x p). */
    Eigen::Map<Eigen::MatrixXd>* X_;

    /** @brief Number of sparse PCs to extract. */
    const Eigen::Index M_;

    /** @brief Target false discovery rate. */
    const double tFDR_;

    /** @brief Algorithmic control parameters. */
    TRexSPCAControl spca_ctrl_;

    /** @brief Random seed. A negative value gives a run that is not reproducible. */
    const int seed_;

    /** @brief Verbose output flag. */
    const bool verbose_;

    /** @brief Result from the last select() call. */
    TRexSPCAResult result_;

    /** @brief Centering scaler; stored as member so the destructor can call
     *  inverse_transform_inplace if select() was interrupted by an exception. */
    std_scaler::LpNormScaler scaler_;

    /** @brief True while X is in centered state (between transform and inverse_transform
     *  inside select()). Guards against leaving X permanently mutated on exceptions. */
    bool X_is_centered_{false};


    // ============================================================
    // Private Helpers
    // ============================================================

    /**
     * @brief Assemble loadings via Ridge regression on the active set.
     *
     * @param X_ref      Centered design matrix (n x p).
     * @param z_m        Ordinary PC response vector (n x 1).
     * @param active_set Selected variable indices.
     * @param v_out      Output loading vector (p x 1), written in-place.
     */
    void assembleActiveSetLoadings(const Eigen::Ref<const Eigen::MatrixXd>& X_ref,
                                   const Eigen::Ref<const Eigen::VectorXd>& z_m,
                                   const std::vector<Eigen::Index>& active_set,
                                   Eigen::Ref<Eigen::VectorXd> v_out) const;

    /**
     * @brief Assemble thresholded loadings: top |A| ordinary loadings, normalized.
     *
     * @param ord_v       Ordinary PCA loading vector (p x 1).
     * @param active_size Cardinality of the active set |A|.
     * @param v_out       Output loading vector (p x 1), written in-place.
     */
    void assembleThresholdedLoadings(const Eigen::Ref<const Eigen::VectorXd>& ord_v,
                                     Eigen::Index active_size,
                                     Eigen::Ref<Eigen::VectorXd> v_out) const;

    /**
     * @brief Compute adjusted explained variance via QR decomposition of Z.
     *
     * @param result  TRexSPCAResult to update (Z must already be assembled).
     * @param X_ref   Centered design matrix (n x p) for total variance computation.
     */
    void adjustExplainedVariance(TRexSPCAResult& result,
                                 const Eigen::Ref<const Eigen::MatrixXd>& X_ref) const;
};

/** @brief Old name of TRexSPCASelector.
 *  @deprecated Use TRexSPCASelector. */
using TRexSPCA [[deprecated("Use TRexSPCASelector")]] = TRexSPCASelector;

// ==============================================================================
} /* End of namespace trex::trex_selector_methods::ad::trex_spca */
// ==============================================================================
#endif /* End of TREX_SELECTOR_METHODS_TREX_SPCA_HPP */

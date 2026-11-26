// =====================================================================================
// trex_core_bindings.hpp - Core Python bindings for TRexSelector classes
// =====================================================================================
#ifndef TREX_PYTHON_TREX_CORE_BINDINGS_HPP
#define TREX_PYTHON_TREX_CORE_BINDINGS_HPP
// =====================================================================================
/**
 * @file trex_core_bindings.hpp
 *
 * @brief Pybind11 bindings for the core TRexSelector algorithms and configurations.
 */
// =====================================================================================

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

// Core TRexSelector includes
#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <trex_selector_methods/ad/trex_utils/trex_solver_dispatch.hpp>
#include <utils/datageneration/utils_dummygen.hpp>

// binding helpers
#include "deprecation.hpp"
#include "design_views.hpp"

// std includes
#include <numbers>

// =====================================================================================

namespace py = pybind11;
using namespace trex::trex_selector_methods::ad::trex_core;

// =====================================================================================

/**
 * @brief Design of a selector from an array that is taken as it is.
 *
 * @param X            Design array, typically not writable.
 * @param trex_control Control parameters that hold the scaling mode.
 *
 * @return Reference to the array's buffer.
 *
 * @throws py::value_error under a scaling mode that writes into a design
 *         that is not writable.
 */
inline Eigen::Ref<Eigen::MatrixXd> selector_design_ref(
    const DesignArray& X,
    const TRexControlParameter& trex_control
) {
    namespace dn = trex::trex_selector_methods::ad::utils::data_normalizer;
    return reader_design_map(X, dn::writesDesign(trex_control.scaling_mode), "X", kSelectorDesignRule);
}

// =====================================================================================
// Core Class Mapping Wrapper
// =====================================================================================

/**
 * @brief Python wrapper class for TRexSelector leveraging zero-copy references.
 *
 * @details Retains references to the underlying numpy arrays via Eigen::Map without performing
 *          unnecessary deep copies of the data, analogous to the solver wrappers.
 */
class PyTRexSelector {
protected:
    PyTRexSelector() = default;
public:
    virtual ~PyTRexSelector() = default;
protected:
    /** @brief Map for the original design matrix (X) */
    std::unique_ptr<Eigen::Map<Eigen::MatrixXd>> X_map_;

    /** @brief Map for the response vector (y) */
    std::unique_ptr<Eigen::Map<Eigen::VectorXd>> y_map_;

    /** @brief Pointer to the underlying core orchestrator instance */
    std::unique_ptr<TRexSelector> selector_;

public:
    /**
     * @brief Zero-copy constructor for the TRexSelector wrapper.
     *
     * @param X Reference to the design matrix (n × p).
     * @param y Reference to the response vector (n × 1), copied by the selector.
     * @param tFDR Target False Discovery Rate in [0, 1].
     * @param trex_control Control parameters of the algorithm.
     * @param seed Random seed (< 0 for a run that is not reproducible).
     * @param verbose Enable verbose output.
     */
    PyTRexSelector(
        Eigen::Ref<Eigen::MatrixXd> X,
        Eigen::Ref<Eigen::VectorXd> y,
        double tFDR,
        const TRexControlParameter& trex_control,
        int seed,
        bool verbose
    ) {
        X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X.data(),
                                                                  X.rows(),
                                                                      X.cols());
        y_map_ = std::make_unique<Eigen::Map<Eigen::VectorXd>>(y.data(),
                                                                  y.size());

        selector_ = std::make_unique<TRexSelector>(
            *X_map_, *y_map_, tFDR, trex_control, seed, verbose
        );
    }

    /**
     * @brief Zero-copy constructor on a design that is not writable.
     *
     * Accepted under ScalingMode::NONE, where the selector writes nothing
     * into the design, and refused under L2 and ZSCORE.
     *
     * @param X Design array, taken as it is (n × p).
     */
    PyTRexSelector(
        const DesignArray& X,
        Eigen::Ref<Eigen::VectorXd> y,
        double tFDR,
        const TRexControlParameter& trex_control,
        int seed,
        bool verbose
    ) : PyTRexSelector(selector_design_ref(X, trex_control), y, tFDR, trex_control, seed, verbose) {}

    /**
     * @brief Execute the full selection algorithm.
     *
     * @return SelectionResult Containing all solution structures, selected variables,
     *         and thresholds.
     */
    virtual TRexSelector::SelectionResult select() {
        return selector_->select();
    }

    /** @brief Give X back before select() (idempotent; no-op after select()). */
    void restore() { selector_->restore(); }

    /** @brief Number of dummies used by the last select() (the L-loop's calibration). */
    std::size_t getNumDummies() const noexcept { return selector_->getNumDummies(); }

    /** @brief Get target FDR level. */
    double getTFDR() const { return selector_->getTFDR(); }

    /** @brief Get number of random experiments K. */
    std::size_t getK() const { return selector_->getK(); }

    /** @brief Get the cap of the dummy multiplier (max_dummy_multiplier). */
    std::size_t getMaxNumDummies() const { return selector_->getMaxNumDummies(); }

    /** @brief Check whether the bound on the stopping time is on (use_max_T_stop). */
    bool getMaxTStop() const { return selector_->getMaxTStop(); }

    /** @brief Check if the stagnation stop is enabled. */
    bool getStagnationCheck() const { return selector_->getStagnationCheck(); }

    /** @brief Get means of the design matrix columns. */
    const Eigen::VectorXd& getXMeans() const { return selector_->getXMeans(); }

    /** @brief Get the scale divisors of the design matrix columns. */
    const Eigen::VectorXd& getXL2Norms() const { return selector_->getXL2Norms(); }

    /** @brief Get mean of the response vector. */
    double getYMean() const { return selector_->getYMean(); }

    /** @brief Get selected variables as a binary vector. */
    const Eigen::VectorXi& getSelectedVar() const { return selector_->getSelectedVar(); }

    /** @brief Get indices of selected variables. */
    const std::vector<std::size_t>& getSelectedIndices() const { return selector_->getSelectedIndices(); }

    /** @brief Get number of selected variables. */
    std::size_t getNumSelected() const { return selector_->getNumSelected(); }

    /** @brief Get the stopping time T_stop at which the T-loop stopped. */
    std::size_t getTStop() const { return selector_->getTStop(); }

    /** @brief Get the voting grid. */
    const Eigen::VectorXd& getVotingGrid() const { return selector_->getVotingGrid(); }

    /** @brief Get voting threshold. */
    double getVotingThreshold() const { return selector_->getVotingThreshold(); }

    /** @brief Get FDP hat matrix. */
    const Eigen::MatrixXd& getFDPHatMat() const { return selector_->getFDPHatMat(); }

    /** @brief Get Phi matrix. */
    const Eigen::MatrixXd& getPhiMat() const { return selector_->getPhiMat(); }

    /** @brief Get the numbers of selected variables (SelectionResult::R_mat). */
    const Eigen::MatrixXd& getRMat() const { return selector_->getRMat(); }

    /** @brief Get Phi prime vector. */
    const Eigen::VectorXd& getPhiPrime() const { return selector_->getPhiPrime(); }

    /** @brief Get dummy multiplier L. */
    std::size_t getDummyMultiplierL() const { return selector_->getDummyMultiplierL(); }
};


/**
 * @brief Binds the TRexSelector core configurations and main wrapper object.
 *
 * @param m_tsm The submodule to bind against.
 */
inline void bind_trex_core(py::module& m_tsm) {
    namespace disp = trex::trex_selector_methods::ad::utils::solver_dispatch;
    namespace dg = trex::utils::datageneration::dummygen;
    namespace dn = trex::trex_selector_methods::ad::utils::data_normalizer;

    // Bind ScalingMode enum
    py::enum_<dn::ScalingMode>(m_tsm, "ScalingMode",
        "Working-design convention of the T-Rex selectors. L2 and ZSCORE center X and y and "
        "scale the columns of X (dummies follow); NONE leaves X and y exactly as supplied "
        "(user-managed centering and scaling, e.g. a pre-normalized design or an intercept-free "
        "dictionary), with the dummies kept on the unit-L2 standard.")
        .value("L2",     dn::ScalingMode::L2,     "Divide each centered column by its L2 norm (unit L2 norm).")
        .value("ZSCORE", dn::ScalingMode::ZSCORE, "Divide each centered column by its sample SD (n-1); glmnet "
            "standardize=TRUE convention.")
        .value("NONE",   dn::ScalingMode::NONE,   "Leave X and y untouched (no centering, no scaling); the design "
            "should be unit-L2-commensurate since the dummies keep the unit-L2 standard. A design pre-normalized with "
            "normalize_design() selects under NONE exactly what the raw design selects under L2.")
        .export_values();

    // The library's own normalization as a user-callable helper, so a caller
    // can reproduce the selector's L2 / ZSCORE preprocessing (then run under
    // NONE), or inspect the parameters. Copies; the inputs are untouched.
    m_tsm.def("normalize_design",
        [](Eigen::MatrixXd X, Eigen::VectorXd y, dn::ScalingMode scaling_mode, double eps,
           std::optional<int> n_threads) {
            dn::NormalizationParams params;
            dn::centerY(y, params, scaling_mode);
            dn::centerAndL2NormalizeX(X, params, eps, /*verbose=*/false, scaling_mode,
                                      n_threads.value_or(trex::utils::openmp::get_num_threads()));
            py::dict out;
            out["X"]        = X;
            out["y"]        = y;
            out["X_means"]  = params.X_means;
            out["X_scales"] = params.X_l2norms;
            out["y_mean"]   = params.y_mean;
            return out;
        },
        py::arg("X"), py::arg("y"),
        py::arg("scaling_mode") = dn::ScalingMode::L2,
        py::arg("eps") = 1e-6,
        py::arg("n_threads") = py::none(),
        "Apply the selectors' data normalization to copies of X and y and return a dict with "
        "keys 'X', 'y' (the normalized copies), 'X_means', 'X_scales' (the per-column divisors: "
        "L2 norms for L2, sample SDs for ZSCORE, ones for NONE) and 'y_mean'. Under "
        "ScalingMode.NONE the copies equal the inputs and the parameters are the identity. "
        "eps is the dead-column threshold on the working L2 norm (the selectors use the solver "
        "tolerance, default 1e-6). n_threads is the thread budget of the column passes (None: "
        "the library's default budget, utils.get_num_threads(), initially 1). Raises ValueError "
        "on a non-finite entry or a dead column (constant under L2 / ZSCORE, all-zero under "
        "NONE), naming the column (0-based).");

    // Bind LLoopStrategy enum
    py::enum_<LLoopStrategy>(m_tsm, "LLoopStrategy", "L-loop strategy controlling how dummy variables are constructed "
        "in each experiment.")
        .value("SKIPL", LLoopStrategy::SKIPL, "Stored, one-shot: fixed dummies = max_dummy_multiplier * p.")
        .value("STANDARD", LLoopStrategy::STANDARD, "Stored; fresh dummies at each L-loop iteration (default).")
        .value("HCONCAT", LLoopStrategy::HCONCAT, "Stored; the dummy matrices grow by p columns per L step (faster, "
            "same FDR in simulations, not the same selection as STANDARD).")
        .value("PERMUTATION", LLoopStrategy::PERMUTATION, "Stored base dummy matrix + deterministic row permutations "
            "per experiment.")
        .value("PERMUTATION_SEEDED", LLoopStrategy::PERMUTATION_SEEDED, "Seed-derived base + row permutations per "
            "experiment; nothing stored.")
        .value("SEEDED", LLoopStrategy::SEEDED, "The storage-free twin of STANDARD: fresh dummies at every L step, "
            "re-derived from the seed per experiment and bit-identical to STANDARD for the same seed; nothing stored.")
        .export_values();

    // Bind SolverTypeForTRex enum
    py::enum_<disp::SolverTypeForTRex>(m_tsm, "SolverTypeForTRex", "Solver algorithm variants available for T-Rex.")
        .value("TLARS",      disp::SolverTypeForTRex::TLARS,      "Terminating LARS.")
        .value("TLASSO",     disp::SolverTypeForTRex::TLASSO,     "Terminating LASSO.")
        .value("TENET",      disp::SolverTypeForTRex::TENET,      "Terminating Elastic Net.")
        .value("TSTEPWISE",  disp::SolverTypeForTRex::TSTEPWISE,  "Terminating Stepwise.")
        .value("TSTAGEWISE", disp::SolverTypeForTRex::TSTAGEWISE, "Terminating Stagewise.")
        .value("TOMP",       disp::SolverTypeForTRex::TOMP,       "Terminating OMP.")
        .value("TGP",        disp::SolverTypeForTRex::TGP,        "Terminating Gradient Pursuit.")
        .value("TACGP",      disp::SolverTypeForTRex::TACGP,      "Terminating Approximate Conjugate Gradient Pursuit.")
        .value("TMP",        disp::SolverTypeForTRex::TMP,        "Terminating Matching Pursuit.")
        .value("TNCGMP",     disp::SolverTypeForTRex::TNCGMP,     "Terminating Norm-Corrected Generalized Matching "
            "Pursuit.")
        .value("TOOLS",      disp::SolverTypeForTRex::TOOLS,      "Terminating Orthogonal Least Squares.")
        .value("TAFS",       disp::SolverTypeForTRex::TAFS,       "Terminating AFS.")
        .value("TENET_AUG",  disp::SolverTypeForTRex::TENET_AUG,  "Terminating Elastic Net via augmented LASSO; "
            "accepted wherever TENET runs (LASSO inner path standard, LARS via tenet_aug_use_lars).")
        .value("TIENET",     disp::SolverTypeForTRex::TIENET,     "Terminating Informed Elastic Net, native pathwise "
            "(GVS-only, gvs_type=IEN; IEN default).")
        .value("TIENET_AUG", disp::SolverTypeForTRex::TIENET_AUG, "Terminating Informed Elastic Net via augmented "
            "LASSO "
            "(GVS-only, gvs_type=IEN).")
        .value("TIKNET",     disp::SolverTypeForTRex::TIKNET,     "Terminating Tikhonov Net, native pathwise "
            "(sparse-K, "
            "FOLDED; TRexTikhonovSelector only).")
        .export_values();

    // Bind SolverHyperparameters struct
    py::class_<disp::SolverHyperparameters>(m_tsm, "SolverHyperparameters", "Hyperparameters for the selected T-Rex "
        "solver.")
        .def(py::init<>())
        .def_readwrite("lambda2",       &disp::SolverHyperparameters::lambda2,       "L2 penalty for ENET-type "
            "solvers.")
        .def_readwrite("rho_afs",       &disp::SolverHyperparameters::rho_afs,       "Shrinkage step size for AFS "
            "solver in (0, 1].")
        .def_readwrite("ncgmp_variant", &disp::SolverHyperparameters::ncgmp_variant, "NCGMP variant: 0 = LineSearch, 1 "
            "= FullyCorrective.")
        .def_readwrite("tol",           &disp::SolverHyperparameters::tol,           "Numerical tolerance of the "
            "selector (default 1e-6): the threshold of its own comparisons and the tolerance of the solvers in its "
            "random experiments. A stand-alone solver defaults to machine epsilon (SolverOptions.tol). The "
            "two tolerances serve different purposes and differ on purpose.")
        .def_readwrite("exch_tie_alpha", &disp::SolverHyperparameters::exch_tie_alpha, "Exchangeable-tie band width "
            "of the greedy solvers (TOMP, TOOLS, TGP, TACGP, TMP, TNCGMP and TAFS) in pairwise ranking-noise sd "
            "units. -1 (default): the selector decides. TRexDASelector runs with 0.25, every other selector without "
            "the band. 0 = off, also under T-Rex+DA. The path solvers ignore the option.")
        .def_readwrite("exch_tie_floor", &disp::SolverHyperparameters::exch_tie_floor, "Minimum |correlation| for "
            "exchangeable-tie candidates in (0, 1); ignored unless exch_tie_alpha > 0.")
        .def_readwrite("tenet_aug_use_lars", &disp::SolverHyperparameters::tenet_aug_use_lars, "TENET_AUG: False "
            "(default) = standard LASSO inner path; True = optional pure-LARS inner path that never drops variables.");

    // Bind DummyDistribution (wraps C++ dummygen::Distribution)
    py::class_<dg::Distribution> dummy_dist(m_tsm, "DummyDistribution",
        "Specification of the dummy variable distribution injected by T-Rex.\n"
        "\n"
        "References\n"
        "----------\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n"
        "E. Mammen, \"Bootstrap and Wild Bootstrap for High Dimensional Linear Models,\" Ann.\n"
        "  Statist., vol. 21, no. 1, pp. 255-285, 1993, doi: 10.1214/aos/1176349025.\n");

    py::enum_<dg::Distribution::Type>(dummy_dist, "Type",
        "Dummy distribution type enumeration.")
        .value("Normal",                      dg::Distribution::Type::Normal)
        .value("Uniform",                     dg::Distribution::Type::Uniform)
        .value("Rademacher",                  dg::Distribution::Type::Rademacher)
        .value("StudentT",                    dg::Distribution::Type::StudentT)
        .value("Laplace",                     dg::Distribution::Type::Laplace)
        .value("Gumbel",                      dg::Distribution::Type::Gumbel)
        .value("Triangle",                    dg::Distribution::Type::Triangle)
        .value("UniformSphere",               dg::Distribution::Type::UniformSphere)
        .value("Mammen",                      dg::Distribution::Type::Mammen)
        .value("ConstrainedSparseRademacher", dg::Distribution::Type::ConstrainedSparseRademacher)
        .value("Logistic",                    dg::Distribution::Type::Logistic)
        .export_values();

    dummy_dist
        .def(py::init<>(), "Default constructor: Normal distribution.")
        .def_static("normal",     &dg::Distribution::Normal,     "Create Normal distribution.")
        .def_static("rademacher", &dg::Distribution::Rademacher, "Create Rademacher distribution.")
        .def_static("mammen",     &dg::Distribution::Mammen,
            py::arg("p1") = (std::sqrt(5.0) + 1.0) / (2.0 * std::sqrt(5.0)),
            py::arg("p2") = (std::sqrt(5.0) - 1.0) / (2.0 * std::sqrt(5.0)),
            "Create Mammen two-point distribution.")
        .def_static("uniform", &dg::Distribution::Uniform,
            py::arg("a") = std::sqrt(3.0),
            "Create Uniform U(-a, a) distribution.")
        .def_static("student_t", &dg::Distribution::StudentT,
            py::arg("df") = 5.0,
            "Create Student's t-distribution.")
        .def_static("laplace", &dg::Distribution::Laplace,
            py::arg("location") = 0.0, py::arg("scale") = 1.0 / std::sqrt(2.0),
            "Create Laplace distribution.")
        .def_static("gumbel", &dg::Distribution::Gumbel,
            py::arg("location") = 0.0, py::arg("scale") = 1.0,
            "Create Gumbel distribution.")
        .def_static("triangle", &dg::Distribution::Triangle,
            py::arg("a") = -std::sqrt(6.0), py::arg("b") = 0.0, py::arg("c") = std::sqrt(6.0),
            "Create Triangle distribution.")
        .def_static("uniform_sphere", &dg::Distribution::UniformSphere,
            py::arg("dim") = 3,
            "Create Uniform distribution on the unit d-sphere.")
        .def_static("constrained_sparse_rademacher", &dg::Distribution::ConstrainedSparseRademacher,
            py::arg("s") = 0.1,
            "Create Constrained Sparse Rademacher distribution: an entry is +1 or -1 with probability s / 2 each "
            "and 0 otherwise, s in (0, 1] (default 0.1), and s = 1 gives the Rademacher distribution.")
        .def_static("logistic", &dg::Distribution::Logistic,
            py::arg("location") = 0.0, py::arg("scale") = std::sqrt(3.0) / std::numbers::pi,
            "Create Logistic distribution.")
        .def("get_type", &dg::Distribution::get_type, "Get the distribution type.")
        .def_readwrite("uniform_a",                       &dg::Distribution::uniform_a,
            "Half-range of the uniform distribution U(-a, a) (default sqrt(3)).")
        .def_readwrite("student_t_df",                    &dg::Distribution::student_t_df,
            "Degrees of freedom of Student's t-distribution (default 5).")
        .def_readwrite("laplace_location",                &dg::Distribution::laplace_location,
            "Location of the Laplace distribution (default 0).")
        .def_readwrite("laplace_scale",                   &dg::Distribution::laplace_scale,
            "Scale of the Laplace distribution (default 1 / sqrt(2)).")
        .def_readwrite("gumbel_location",                 &dg::Distribution::gumbel_location,
            "Location of the Gumbel distribution (default 0).")
        .def_readwrite("gumbel_scale",                    &dg::Distribution::gumbel_scale,
            "Scale of the Gumbel distribution (default 1).")
        .def_readwrite("triangle_a",                      &dg::Distribution::triangle_a,
            "Lower bound of the triangle distribution (default -sqrt(6)).")
        .def_readwrite("triangle_b",                      &dg::Distribution::triangle_b,
            "Mode of the triangle distribution (default 0).")
        .def_readwrite("triangle_c",                      &dg::Distribution::triangle_c,
            "Upper bound of the triangle distribution (default sqrt(6)).")
        .def_readwrite("sphere_dim",                      &dg::Distribution::sphere_dim,
            "Dimension of the sphere of the uniform sphere distribution (default 3).")
        .def_readwrite("mammen_param_p1",                 &dg::Distribution::mammen_param_p1,
            "Probability of the point -(sqrt(5) - 1) / 2 of the Mammen distribution.")
        .def_readwrite("mammen_param_p2",                 &dg::Distribution::mammen_param_p2,
            "Probability of the point (sqrt(5) + 1) / 2 of the Mammen distribution.")
        .def_readwrite("constrained_sparse_rademacher_s", &dg::Distribution::constrained_sparse_rademacher_s,
            "Sparsity parameter s of the constrained sparse Rademacher distribution in (0, 1] (default 0.1).")
        .def_readwrite("logistic_location",               &dg::Distribution::logistic_location,
            "Location of the logistic distribution (default 0).")
        .def_readwrite("logistic_scale",                  &dg::Distribution::logistic_scale,
            "Scale of the logistic distribution (default sqrt(3) / pi).");

    // Bind TRexControlParameter
    py::class_<TRexControlParameter>(m_tsm, "TRexControlParameter", "Control parameters for the core T-Rex selection "
        "algorithm.")
        .def(py::init<>())
        .def_readwrite("K", &TRexControlParameter::K,
            "Number of random experiments (default 20, at least 2).")
        .def_readwrite("max_dummy_multiplier", &TRexControlParameter::max_dummy_multiplier, "Cap of the dummy "
            "multiplier L: the L-loop adds dummies up to this multiple of p (default 10).")
        .def_readwrite("use_max_T_stop", &TRexControlParameter::use_max_T_stop, "If true (default), the T-loop ends at "
            "ceil(n / 2) at the latest. If false, it can run up to the number of dummies.")
        .def_readwrite("opt_threshold", &TRexControlParameter::opt_threshold,
            "Reference voting level of the L-loop (default 0.75, in [0, 1]). The L-loop raises L until "
            "the FDP estimate at the largest voting level that does not exceed opt_threshold is at or "
            "below tFDR.")
        .def_readwrite("lloop_strategy", &TRexControlParameter::lloop_strategy,
            "L-loop strategy (default STANDARD). STANDARD draws fresh dummies at every L step, as the 1.x "
            "versions do. HCONCAT keeps the dummies of the earlier steps and appends p columns per step, "
            "which is faster. The two gave the same FDR in simulations and can select differently. See "
            "LLoopStrategy for the others.")
        .def_readwrite("tloop_stagnation_stop", &TRexControlParameter::tloop_stagnation_stop,
                       "T-loop early stopping when the selection stagnates. None (default) resolves by solver "
                       "family: off for the path solvers (TLARS, TLASSO, TENET, TENET_AUG, TIENET, TIENET_AUG "
                       "and TIKNET), which terminate on their own, as in the 1.x "
                       "versions and the T-Rex paper, and on for the greedy solvers (TSTEPWISE, TSTAGEWISE, TOMP, "
                       "TGP, TACGP, TMP, TNCGMP, TOOLS, TAFS), which otherwise run into a noise trap. Set True or "
                       "False to override.")
        .def_readwrite("tloop_max_stagnant_steps", &TRexControlParameter::tloop_max_stagnant_steps,
            "Number of T steps in a row without a growing selection that ends the T-loop under the "
            "stagnation stop (default 5, at least 1). The selection is counted at the highest voting "
            "level whose FDP estimate does not exceed tFDR.")
        .def_readwrite("parallel_rnd_experiments", &TRexControlParameter::parallel_rnd_experiments,
                       "Run the K random experiments in parallel (default False: in sequence). Inside a "
                       "parallel run every solver runs serially; the lean L-loop strategies (SEEDED, "
                       "PERMUTATION, PERMUTATION_SEEDED) always run their experiments in sequence.")
        .def_readwrite("use_memory_mapping",  &TRexControlParameter::use_memory_mapping,
            "Keep the dummy matrices in memory-mapped files instead of in memory (default False).")
        .def_readwrite("max_outer_threads",   &TRexControlParameter::max_outer_threads,
                       "Threads of the parallel K-loop (default 1; a whole number of at least 1), capped at "
                       "K. Ignored while parallel_rnd_experiments is False.")
        .def_readwrite("max_inner_threads",   &TRexControlParameter::max_inner_threads,
                       "Thread budget of the solvers' own passes over the design (default 1, serial; a whole "
                       "number of at least 1). Live when the experiments run in sequence; inside a parallel "
                       "run the solvers stay serial. A pass splits only above about a million entries.")
        .def_readwrite("solver_type",         &TRexControlParameter::solver_type,
            "Solver of the random experiments (default TLARS).")
        .def_readwrite("solver_params",       &TRexControlParameter::solver_params,       "Hyperparameters for the "
            "selected solver.")
        .def_readwrite("scaling_mode",        &TRexControlParameter::scaling_mode,        "Working-design convention: "
            "L2 (default) or ZSCORE center X and y and scale X; NONE leaves X and y as supplied (user-managed; dummies "
            "keep the unit-L2 standard). Non-finite entries and dead columns are rejected in every mode (ValueError). "
            "L2 and ZSCORE standardize X in place and restore it, so X must be writable. NONE writes nothing into X, "
            "so an array that is not writable serves (the view of a ReadOnly MemoryMappedMatrix).")
        .def_readwrite("dummy_distribution",  &TRexControlParameter::dummy_distribution,  "Distribution used to "
            "generate dummy variables.");

    // Bind TRexSelector::SelectionResult
    py::class_<TRexSelector::SelectionResult>(m_tsm, "SelectionResult", "Return value of TRexSelector.select(): "
        "selected variables and diagnostic matrices.")
        .def_readonly("selected_var", &TRexSelector::SelectionResult::selected_var, "Selection indicator vector "
            "(length "
            "p; 1 = selected, 0 = not selected).")
        .def_readonly("v_thresh", &TRexSelector::SelectionResult::v_thresh,
            "Voting threshold applied to the relative occurrences to obtain the selection.")
        .def_readonly("r_mat", &TRexSelector::SelectionResult::R_mat,
            "Number of selected variables per stopping time (rows) and voting level (columns), (T_stop - "
            "1) x len(voting_grid). The stopping time T_stop is not part of it. One row for T_stop = 1.")
        .def_readonly("T_stop", &TRexSelector::SelectionResult::T_stop,
            "Stopping time T_stop: the number of included dummies at which the T-loop stopped, so the "
            "random experiments ran until T_stop dummies were included. For T_stop > 1 it is not the "
            "stopping time of the selection: the selection is taken from the stopping times 1 to T_stop - "
            "1, mostly from T_stop - 1 and sometimes from an earlier one (see select()).")
        .def_readonly("num_dummies", &TRexSelector::SelectionResult::num_dummies,
            "Total number of dummy variables used (L * p).")
        .def_readonly("dummy_factor_L", &TRexSelector::SelectionResult::dummy_factor_L,
            "Dummy multiplier L determined by the L-loop.")
        .def_readonly("fdp_hat_mat", &TRexSelector::SelectionResult::FDP_hat_mat,
            "Estimated FDP per stopping time (rows) and voting level (columns), T_stop x "
            "len(voting_grid).")
        .def_readonly("phi_mat", &TRexSelector::SelectionResult::Phi_mat,
            "Relative occurrences of the variables over the K random experiments per stopping time "
            "(rows), T_stop x p.")
        .def_readonly("phi_prime", &TRexSelector::SelectionResult::Phi_prime,
            "Calibrated relative occurrences that enter the FDP estimate (length p), at most 1 and "
            "possibly negative.")
        .def_readonly("voting_grid", &TRexSelector::SelectionResult::voting_grid,
            "Voting levels 0.5, 0.5 + 1/K, ... below 1, closed by 1 - eps.");

    // Bind PyTRexSelector
    py::class_<PyTRexSelector>(m_tsm, "TRexSelector",
        "T-Rex Selector: FDR-controlled variable selection via dummy comparison.\n"
        "\n"
        "References\n"
        "----------\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n")
        .def(py::init<Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::VectorXd>, double, const TRexControlParameter&,
            int, bool>(),
             py::arg("X"), py::arg("y"), py::arg("tFDR") = 0.1,
             py::arg("trex_control") = TRexControlParameter(),
             py::arg("seed") = -1, py::arg("verbose") = true,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             "Construct the selector. X is used without a copy, y is copied by the selector.")
        .def(py::init<const DesignArray&, Eigen::Ref<Eigen::VectorXd>, double, const TRexControlParameter&,
            int, bool>(),
             py::arg("X").noconvert(), py::arg("y"), py::arg("tFDR") = 0.1,
             py::arg("trex_control") = TRexControlParameter(),
             py::arg("seed") = -1, py::arg("verbose") = true,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             kReaderSelectorDoc.c_str())
        .def("getNumDummies", &PyTRexSelector::getNumDummies,
             "Number of dummies used by the last select() (the L-loop's calibration).")
        .def("restore", &PyTRexSelector::restore,
             "Give X back before select(): undo the in-place normalization applied at "
             "construction and release the buffer claim (idempotent; a no-op after "
             "select(), which restores X itself). A restored selector cannot run: "
             "select() then raises, construct a new one. Without it, an unselected "
             "selector gives X back only when it is garbage-collected.")
        .def("select", &PyTRexSelector::select, py::call_guard<py::gil_scoped_release>(),
             "Run the T-Rex algorithm and return a SelectionResult. A selector runs once: select() gives X "
             "back when it ends, so a second call raises RuntimeError instead of running on a design that "
             "is no longer standardized. Construct a new selector to run again. A select() that ends with "
             "an error gives X back at once (original scale) and retires the selector as well.\n\nThe "
             "T-loop raises the stopping time until the FDP estimate at the highest voting level exceeds "
             "tFDR, or until its bound or the stagnation stop ends it, and T_stop is where it stopped. "
             "The selection is the cell with the largest number of "
             "selected variables among the stopping times 1 to T_stop - 1 (the stopping time 1 for T_stop "
             "= 1) and the voting levels whose FDP estimate does not exceed tFDR. Among equal counts the "
             "highest voting level wins, then the latest stopping time. The stopping time of the selection "
             "is therefore below T_stop, mostly by one and sometimes by more.")
        .def("getTFDR", &PyTRexSelector::getTFDR, "Return the target FDR.")
        .def("getK", &PyTRexSelector::getK, "Return K: the number of random experiments.")
        .def("getMaxNumDummies", &PyTRexSelector::getMaxNumDummies, "Return the cap of the dummy multiplier "
            "(max_dummy_multiplier), a multiplier and not a number of dummies.")
        .def("getMaxTStop", &PyTRexSelector::getMaxTStop, "Return whether use_max_T_stop is active.")
        .def("getStagnationCheck", &PyTRexSelector::getStagnationCheck, "Return whether T-loop stagnation stopping is "
            "active.")
        .def("getXMeans", &PyTRexSelector::getXMeans,
            "Return the column means of X that the selector subtracted (zeros under NONE).")
        .def("getXL2Norms", &PyTRexSelector::getXL2Norms,
            "Return the scale divisors of the columns of X: the L2 norms of the centered columns under "
            "L2, the sample standard deviations under ZSCORE, ones under NONE.")
        .def("getYMean", &PyTRexSelector::getYMean,
            "Return the mean of y that the selector subtracted (0 under NONE).")
        .def("getSelectedVar", &PyTRexSelector::getSelectedVar, "Return the binary selection indicator vector (1 = "
            "selected).")
        .def("getSelectedIndices", &PyTRexSelector::getSelectedIndices,
            "Return the 0-based indices of the selected variables.")
        .def("getNumSelected", &PyTRexSelector::getNumSelected,
            "Return the number of selected variables.")
        .def("getTStop", &PyTRexSelector::getTStop,
            "Return T_stop, the stopping time at which the T-loop stopped. The selection is taken from an "
            "earlier one, see SelectionResult.T_stop.")
        .def("getVotingGrid", &PyTRexSelector::getVotingGrid,
            "Return the voting levels 0.5, 0.5 + 1/K, ... below 1, closed by 1 - eps.")
        .def("getVotingThreshold", &PyTRexSelector::getVotingThreshold,
            "Return the voting threshold of the selection.")
        .def("getFDPHatMat", &PyTRexSelector::getFDPHatMat,
            "Return the estimated FDP per stopping time (rows) and voting level (columns).")
        .def("getPhiMat", &PyTRexSelector::getPhiMat,
            "Return the relative occurrences of the variables per stopping time (rows).")
        .def("getRMat", &PyTRexSelector::getRMat,
            "Return the number of selected variables per stopping time (rows) and voting level (columns), "
            "see SelectionResult.r_mat.")
        .def("getPhiPrime", &PyTRexSelector::getPhiPrime,
            "Return the calibrated relative occurrences that enter the FDP estimate, at most 1 and "
            "possibly negative.")
        .def("getDummyMultiplierL", &PyTRexSelector::getDummyMultiplierL, "Return L: the dummy multiplier used.");
}
// =====================================================================================
#endif /* End of TREX_PYTHON_TREX_CORE_BINDINGS_HPP */

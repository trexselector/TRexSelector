// =====================================================================================
// tsolver_bindings.hpp - Python bindings for T-Solvers
// =====================================================================================
#ifndef TREX_PYTHON_TSOLVER_BINDINGS_HPP
#define TREX_PYTHON_TSOLVER_BINDINGS_HPP
// =====================================================================================
/**
 * @file tsolver_bindings.hpp
 *
 * @brief Pybind11 bindings for T-Solvers used in the T-Rex Selector package.
 */
// =====================================================================================

 // pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

// T-Solvers includes
// LARS-based solvers
#include <tsolvers/ad/linear_model/lars_based/tlars_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tlasso_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tstepwise_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tstagewise_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tenet_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tenet_aug_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tienet_aug_solver.hpp>

// OMP-based solvers
#include <tsolvers/ad/linear_model/omp_based/tomp_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tgp_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tacgp_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tmp_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tncgmp_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tools_solver.hpp>

// AFS-based solvers
#include <tsolvers/ad/linear_model/afs_based/tafs_solver.hpp>
// Native IEN and the sparse-K solvers
#include <tsolvers/ad/linear_model/lars_based/tienet_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tiknet_solver.hpp>
// CD-based solvers

// binding helpers
#include "deprecation.hpp"
#include "design_views.hpp"

// std includes
#include <utility>

// =====================================================================================

namespace py = pybind11;

/**
 * @brief Python wrapper class for T-Solvers.
 *
 * @details Retains references to the underlying numpy arrays via Eigen::Map without performing
 *          unnecessary deep copies of the data. The wrapper operates zero-copy directly on the
 *          Python buffers.
 *
 * @tparam Solver The specific underlying T-Solver class from the C++ library.
 */
template <typename Solver>
class PySolverWrapper {
protected:

    /** @brief Map for the design matrix (X) */
    std::unique_ptr<Eigen::Map<Eigen::MatrixXd>> X_map_;

    /** @brief Map for the dummy variable matrix (D) */
    std::unique_ptr<Eigen::Map<Eigen::MatrixXd>> D_map_;

    /** @brief Map for the response vector (y) */
    std::unique_ptr<Eigen::Map<Eigen::VectorXd>> y_map_;

    /** @brief Pointer to the underlying solver instance */
    std::unique_ptr<Solver> solver_;

public:
    /**
     * @brief Construct a new Py Solver Wrapper object utilizing zero-copy mapping.
     *
     * @tparam Args Variadic template parameters for the specific solver's constructor.
     *
     * @param X Reference to the design matrix (n × p).
     * @param D Reference to the dummy variable matrix (n × num_dummies).
     * @param y Reference to the response vector (n × 1).
     */
    template <typename... Args>
    PySolverWrapper(Eigen::Ref<Eigen::MatrixXd> X,
                    Eigen::Ref<Eigen::MatrixXd> D,
                    Eigen::Ref<Eigen::VectorXd> y,
                    Args&&... args)
    {
        X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X.data(),
                                                               X.rows(),
                                                               X.cols());
        D_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(D.data(),
                                                               D.rows(),
                                                               D.cols());
        y_map_ = std::make_unique<Eigen::Map<Eigen::VectorXd>>(y.data(), y.size());

        solver_ = std::make_unique<Solver>(*X_map_, *D_map_, *y_map_, std::forward<Args>(args)...);
    }

    /**
     * @brief Execute the solver algorithm.
     *
     * @details Computes the solution path by iteratively including variables into the active set.
     *          If `T_stop` is set to 0 (or -1 in Python, mapping to max unsigned integer), the
     *          solver calculates the full solution path until all possible variables are included
     *          or the constraints prevent further additions.
     *
     * @param T_stop The maximum number of dummy variables in the active set.
     *               Pass 0 or -1 to calculate the full solution path.
     * @param early_stop Whether to break execution early immediately upon selecting `T_stop`
     *                   dummy variables. Default is true.
     */
    void executeStep(std::size_t T_stop = 0, bool early_stop = true) {
        solver_->executeStep(T_stop, early_stop);
    }

    /**
     * @brief Get the regression coefficients (Beta) for a specific step.
     *
     * @param step Specific path step to retrieve Beta for. Default is -1, returning the
     *             coefficients from the final computed step.
     *
     * @return Eigen::VectorXd containing the coefficients.
     */
    Eigen::VectorXd getBeta(int step = -1) const { return solver_->getBeta(step); }

    /**
     * @brief Get the intercept for a specific step.
     *
     * @param step Specific path step to retrieve the intercept for. Default is -1, returning
     *             the intercept from the final computed step.
     *
     * @return double representing the intercept value.
     */
    double getIntercept(int step = -1) const { return solver_->getIntercept(step); }

    /**
     * @brief Retrieve the full path history of the Beta coefficients.
     *
     * @details The solver stores the path sparsely; this densifies it for the
     * numpy exchange — O((p + num_dummies) * steps) memory at call time.
     *
     * @return Eigen::MatrixXd where each column represents the Beta coefficient state at a
     *         specific step in the sequential path computation.
     */
    Eigen::MatrixXd getBetaPath() const { return solver_->getBetaPathSparse().dense(); }

    /**
     * @brief Retrieve the string identifier for the underlying T-Solver class used.
     *
     * @return std::string representation of the solver engine.
     */
    std::string solverTypeToString() const { return solver_->solverTypeToString(); }

    // =================================================================================
    // API Functionality Expansions (Metrics, Path Info, Diagnostics, and State Mgmt)
    // =================================================================================

    /** @brief Get the estimated Mallow's Cp values at each step. */
    std::vector<double> getCp() const { return solver_->getCp(); }

    /** @brief Get the estimated Mallow's Cp values at each step with provided sigma squared. */
    std::vector<double> getCpWithSigma(double sigma_hat_sq) const {
        return static_cast<const trex::tsolvers::ad::TSolverBase*>(solver_.get())->getCp(sigma_hat_sq);
    }

    /** @brief Get the residual sum of squares at each step. */
    const std::vector<double>& getRSS() const { return solver_->getRSS(); }

    /** @brief Get the coefficient of determination (R^2) at each step. */
    const std::vector<double>& getR2() const { return solver_->getR2(); }

    /** @brief Get the degrees of freedom at each step. */
    const std::vector<std::size_t>& getDoF() const { return solver_->getDoF(); }

    /** @brief Get the residual vector at the current step. */
    Eigen::VectorXd getResiduals() const { return solver_->getResiduals(); }

    /** @brief Get the integer codes for actions taken at each step (e.g., param addition/drop). */
    const std::vector<std::vector<int>>& getActions() const { return solver_->getActions(); }

    /** @brief Get indices of variables that are active at the current step. */
    const std::vector<std::size_t>& getActives() const { return solver_->getActives(); }

    /** @brief Get indices of variables that are inactive at the current step. */
    const std::vector<std::size_t>& getInactives() const { return solver_->getInactives(); }

    /** @brief Get indices of variables that were dropped during the execution. */
    const std::vector<std::size_t>& getDroppedIndices() const { return solver_->getDroppedIndices(); }

    /** @brief Get indices of active predictor (non-dummy) variables at the current step. */
    std::vector<std::size_t> getActivePredictorIndices() const {
        return solver_->getActivePredictorIndices();
    }

    /** @brief Get indices of active dummy variables at the current step. */
    std::vector<std::size_t> getActiveDummyIndices() const {
        return solver_->getActiveDummyIndices();
    }

    /** @brief Calculate how many variables are currently active. */
    std::size_t getNumActives() const { return solver_->getNumActives(); }

    /** @brief Find where dummy variables start within the merged predictor index space. */
    std::size_t getDummyStartIndex() const { return solver_->getDummyStartIndex(); }

    /** @brief Quantify the number of path steps successfully taken thus far. */
    std::size_t getNumSteps() const { return solver_->getNumSteps(); }

    /** @brief Heap bytes the solver holds (TSolverBase::memoryBytes(); the
     *  borrowed X and D and the views of this wrapper are not counted). */
    std::size_t getMemoryBytes() const { return solver_->memoryBytes(); }

    /** @brief Check whether the solver has mapped access to data predictors. */
    bool isConnected() const { return solver_->isConnected(); }

    /** @brief Set the numerical tolerance used for internal numeric comparisons. */
    void setTolerance(double tol) { solver_->setTolerance(tol); }

    /** @brief Set the random seed for reproducible tie-breaking among dummies. */
    void setTieSeed(uint32_t seed) { solver_->setTieSeed(seed); }

    /** @brief Set the thread budget of the solver's own regions (1 serial, at least 1). */
    void setNumThreads(int n_threads) { solver_->setNumThreads(n_threads); }

    /** @brief The resolved thread budget of the solver's own regions (>= 1). */
    int getNumThreads() const { return solver_->getNumThreads(); }

    /** @brief Configure exchangeability-calibrated stochastic tie-breaking (greedy solvers). */
    void setExchangeableTie(double alpha, double floor) {
        solver_->setExchangeableTie(alpha, floor);
    }

    /** @brief Get all recorded warnings (dropped columns, collinear rejections, ...). */
    const std::vector<std::string>& getWarnings() const { return solver_->getWarnings(); }

    /** @brief Clear the recorded warnings. */
    void clearWarnings() { solver_->clearWarnings(); }

    // =================================================================================
    // Solver-Specific Getters (only valid for the corresponding solver type)
    // =================================================================================

    /** @brief Get regularization parameter sequence (the solvers with a lambda path). */
    const std::vector<double>& getLambda() const { return solver_->getLambda(); }

    /** @brief Get the current Kahan summation refresh interval (TLARS-specific). */
    std::size_t getKahanRefreshInterval() const { return solver_->getKahanRefreshInterval(); }

    /** @brief Set the Kahan summation refresh interval (TLARS-specific). */
    void setKahanRefreshInterval(std::size_t interval) { solver_->setKahanRefreshInterval(interval); }

    /** @brief Get total number of variable removals (TLASSO/TSTAGEWISE/TENET-specific). */
    std::size_t getNumRemovals() const { return solver_->getNumRemovals(); }

    /** @brief Get cycling ratio of removals to additions (TLASSO/TSTAGEWISE/TENET-specific). */
    double getCyclingRatio() const { return solver_->getCyclingRatio(); }

    /** @brief Get the Elastic Net mixing parameter alpha (TENET-specific). */
    double getAlpha() const { return solver_->getAlpha(); }

    /** @brief Get the secondary L2/group-ridge penalty lambda2 (TIENETAug-specific). */
    double getLambda2() const { return solver_->getLambda2(); }

    /** @brief Get the 0-based group id per original variable (TIENETAug-specific). */
    Eigen::VectorXi getGroups() const { return solver_->getGroups(); }

    /** @brief Get the number of disjoint groups (TIENETAug-specific). */
    std::size_t getNumGroups() const { return solver_->getNumGroups(); }

    // ---- sparse-K solver ----------------------------------------------------
    /** @brief The Tikhonov matrix K as a dense p x p array. */
    Eigen::MatrixXd getTikhonovK() const { return Eigen::MatrixXd(solver_->getTikhonovK()); }

    // =================================================================================
    // Structural Bindings (State reconnects and Cereal Serialization)
    // =================================================================================

    /**
     * @brief Reconnect the underlying native solver engine to a new batch of underlying Maps.
     */
    void reconnect(Eigen::Ref<Eigen::MatrixXd> X, Eigen::Ref<Eigen::MatrixXd> D) {
        X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X.data(),
                                                              X.rows(),
                                                               X.cols());

        D_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(D.data(),
                                                               D.rows(),
                                                               D.cols());

        solver_->reconnect(*X_map_, *D_map_);
    }

    /**
     * @brief Restore solver state parameters based on fresh target views.
     */
    void restore(Eigen::Ref<Eigen::MatrixXd> X,
                 Eigen::Ref<Eigen::MatrixXd> D,
                 Eigen::Ref<Eigen::VectorXd> y) {
        X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X.data(),
                                                               X.rows(),
                                                               X.cols());

        D_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(D.data(),
                                                               D.rows(),
                                                               D.cols());

        y_map_ = std::make_unique<Eigen::Map<Eigen::VectorXd>>(y.data(), y.size());

        solver_->restore(*X_map_, *D_map_, *y_map_);
    }

    /**
     * @brief Serialize calculation state directly to disk utilizing Cereal.
     */
    void save(const std::string& filepath) const {
        solver_->save(filepath);
    }

    /**
     * @brief Hydrate a solver explicitly from a pre-calculated Cereal output file on disk.
     */
    void load(const std::string& filepath) {
        solver_ = std::make_unique<Solver>(Solver::load(filepath, *X_map_, *D_map_));
    }
};


/**
 * @brief The in-place writer/reader contract, appended to every T-solver's
 *        constructor docstring (the C++ text lives on
 *        TSolverBase::restore()).
 */
inline const std::string kInPlaceContractDoc =
    "\n\nX and D are Fortran-ordered float64 arrays, which the solver shares without "
    "a copy. An array of another layout or type is refused with TypeError, and "
    "numpy.asfortranarray() converts it once."
    "\n\nIn-place contract: with normalize=True or intercept=True the solver is a "
    "WRITER. It centers/scales the columns of X and D IN PLACE through the views, "
    "keeps the statistics itself and holds a claim on both buffers until restore() "
    "or destruction, both of which write X and D back (restore() is idempotent and "
    "must be given the same arrays); constructing a second writer (T-solver or "
    "T-Rex selector) on a claimed buffer raises RuntimeError, because a second "
    "in-place transformation moves the data by rounding. With normalize=False and "
    "intercept=False the solver is a READER: it runs as constructed and never "
    "touches the buffers. A pipeline with several solvers on one buffer should "
    "pre-scale once in front (LpNormScaler) and construct the solvers as readers; a "
    "solver resumed with load() is a reader that expects the buffers in the saved "
    "state. y is never written (the solver centers an internal copy). A READER "
    "also accepts X and D that are not writable, for example the to_numpy() view "
    "of a MemoryMappedMatrix opened with AccessMode.ReadOnly (Fortran-contiguous "
    "float64, zero-copy); a WRITER raises ValueError for them.";

/** @brief Docstring of the solver constructors that take X and D as they are. */
inline const std::string kReaderSolverDoc =
    "Initializes the solver on X and D that are not writable, for example the to_numpy() view "
    "of a MemoryMappedMatrix opened with AccessMode.ReadOnly. Both must be Fortran-contiguous "
    "float64 arrays and are accessed zero-copy. " + kSolverDesignRule +
    " Otherwise ValueError is raised.";

/** @brief X or D of a solver from an array that is taken as it is. */
inline Eigen::Ref<Eigen::MatrixXd> solver_design_ref(
    const DesignArray& M,
    const trex::tsolvers::ad::TSolverOptions& options,
    const char* name
) {
    return reader_design_map(M, options.writesDesign(), name, kSolverDesignRule);
}

/** @brief X and D of a solver, checked in this order. Two calls inside one
 *  argument list run in an unspecified order, and the message of a refused
 *  array names the one that was checked first. */
inline std::pair<Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXd>> solver_design_refs(
    const DesignArray& X,
    const DesignArray& D,
    const trex::tsolvers::ad::TSolverOptions& options
) {
    Eigen::Ref<Eigen::MatrixXd> X_ref = solver_design_ref(X, options, "X");
    Eigen::Ref<Eigen::MatrixXd> D_ref = solver_design_ref(D, options, "D");
    return {X_ref, D_ref};
}


/**
 * @brief Attach the method set shared by every T-Solver to an existing class binding.
 *
 * @details Binds all TSolverBase-backed methods (path execution, coefficient and
 *          diagnostics getters, state management, and Cereal serialization). The
 *          constructor is intentionally NOT bound here — each solver attaches the
 *          init matching its own C++ constructor signature.
 *
 * @tparam Solver The specific underlying T-Solver class being wrapped.
 *
 * @param cls The Pybind11 class the methods are attached to.
 */
template<typename Solver>
py::class_<PySolverWrapper<Solver>>& add_standard_tsolver_methods(
    py::class_<PySolverWrapper<Solver>>& cls) {
    cls.def("executeStep", &PySolverWrapper<Solver>::executeStep,
             py::arg("T_stop") = 0,
             py::arg("early_stop") = true,
             "Executes the underlying algorithmic solve mechanics up to constraint `T_stop`."
             " Pass 0 or -1 to calculate the full sequential path.")
        .def("getBeta", &PySolverWrapper<Solver>::getBeta, py::arg("step") = -1,
             "Retrieves the active regression coefficients (Beta) for a specific step.")
        .def("getIntercept", &PySolverWrapper<Solver>::getIntercept, py::arg("step") = -1,
             "Retrieves the calculated model intercept for a specific step.")
        .def("getBetaPath", &PySolverWrapper<Solver>::getBetaPath,
             "Retrieves a historical path history block containing all calculated Beta coefficients"
             " across execution steps.")
        .def("solverTypeToString", &PySolverWrapper<Solver>::solverTypeToString,
             "Returns a string descriptor mapping for the underlying wrapped T-Solver.")
        .def("getCp", &PySolverWrapper<Solver>::getCp,
             "Get the estimated Mallow's Cp values at each step.")
        .def("getCpWithSigma", &PySolverWrapper<Solver>::getCpWithSigma, py::arg("sigma_hat_sq"),
             "Get the estimated Mallow's Cp values at each step with provided sigma squared.")
        .def("getRSS", &PySolverWrapper<Solver>::getRSS,
             "Get the residual sum of squares at each step.")
        .def("getR2", &PySolverWrapper<Solver>::getR2,
             "Get the coefficient of determination (R^2) at each step.")
        .def("getDoF", &PySolverWrapper<Solver>::getDoF,
             "Get the degrees of freedom at each step.")
        .def("getResiduals", &PySolverWrapper<Solver>::getResiduals,
             "Get the residual vector at the current step.")
        .def("getActions", &PySolverWrapper<Solver>::getActions,
             "Get the integer codes for actions taken at each step. Entries "
             "follow the R lars convention: 1-based signed indices, +(j+1) "
             "for the addition and -(j+1) for the removal of 0-based "
             "variable j (decode with abs(a) - 1).")
        .def("getActives", &PySolverWrapper<Solver>::getActives,
             "Get indices of variables that are active at the current step.")
        .def("getInactives", &PySolverWrapper<Solver>::getInactives,
             "Get indices of variables that are inactive at the current step.")
        .def("getDroppedIndices", &PySolverWrapper<Solver>::getDroppedIndices,
             "Get indices of variables that were dropped during the execution.")
        .def("getActivePredictorIndices", &PySolverWrapper<Solver>::getActivePredictorIndices,
             "Get indices of active predictor (non-dummy) variables at the current step.")
        .def("getActiveDummyIndices", &PySolverWrapper<Solver>::getActiveDummyIndices,
             "Get indices of active dummy variables at the current step.")
        .def("getNumActives", &PySolverWrapper<Solver>::getNumActives,
             "Calculate how many variables are currently active.")
        .def("getDummyStartIndex", &PySolverWrapper<Solver>::getDummyStartIndex,
             "Find where dummy variables start within the merged predictor index space.")
        .def("getNumSteps", &PySolverWrapper<Solver>::getNumSteps,
             "Quantify the number of path steps successfully taken thus far.")
        .def("getMemoryBytes", &PySolverWrapper<Solver>::getMemoryBytes,
             "Heap bytes the solver holds: its copy of y, the residual, correlations, "
             "Cholesky factor, coefficient path, diagnostics and reserved capacity. X and D "
             "are borrowed and not counted, nor are the temporaries inside a step, so the "
             "peak of a step is higher.")
        .def("isConnected", &PySolverWrapper<Solver>::isConnected,
             "Check whether the solver has mapped access to data predictors.")
        .def("setTolerance", &PySolverWrapper<Solver>::setTolerance, py::arg("tol"),
             "Set the numerical tolerance used for internal numeric comparisons.")
        .def("setTieSeed", &PySolverWrapper<Solver>::setTieSeed, py::arg("seed"),
             "Set the random seed for reproducible tie-breaking among dummy candidates.")
        .def("setNumThreads", &PySolverWrapper<Solver>::setNumThreads, py::arg("n_threads"),
             "Set the thread budget of the solver's own passes over the design: 1 serial, a whole "
             "number of at least 1 (other values raise). A pass splits only above about a million "
             "entries and never inside a parallel region of the caller.")
        .def("getNumThreads", &PySolverWrapper<Solver>::getNumThreads,
             "The resolved thread budget of the solver's own passes (at least 1).")
        .def("setExchangeableTie", &PySolverWrapper<Solver>::setExchangeableTie,
             py::arg("alpha"), py::arg("floor") = 0.5,
             "Configure exchangeability-calibrated stochastic tie-breaking for greedy "
             "candidate selection (TOMP/TAFS); alpha <= 0 turns it off.")
        .def("getWarnings", &PySolverWrapper<Solver>::getWarnings,
             "Get all recorded warnings (dropped columns, collinear rejections, ...).")
        .def("clearWarnings", &PySolverWrapper<Solver>::clearWarnings,
             "Clear the recorded warnings.")
        .def("reconnect", &PySolverWrapper<Solver>::reconnect, py::arg("X"), py::arg("D"),
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             "Reconnect the underlying native solver engine to a new batch of underlying Matrices.")
        .def("restore", &PySolverWrapper<Solver>::restore, py::arg("X"), py::arg("D"),
             py::arg("y"),
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>(),
             "Undo the solver's in-place centering/scaling of X and D (a writer, see the "
             "constructor) and release its claim on the buffers, so that a later writer "
             "may take them; y is not modified. No-op for a reader.")
        .def("save", &PySolverWrapper<Solver>::save, py::arg("filepath"),
             "Serialize calculation state directly to disk utilizing Cereal.")
        .def("load", &PySolverWrapper<Solver>::load, py::arg("filepath"),
             "Hydrate a solver explicitly from a pre-calculated Cereal output file on disk.");
    return cls;
}


/**
 * @brief Standardized Pybind11 registration for the T-Solvers.
 *
 * @details Binds the canonical `(X, D, y, options)` constructor (every
 *          input except the data lives in TSolverOptions) plus the
 *          shared method set. The IEN solvers, which also take a group
 *          vector, bind their own init and call add_standard_tsolver_methods().
 *
 * @tparam Solver The specific underlying T-Solver class being wrapped.
 *
 * @param m The Pybind11 sub-module where this solver will be exposed.
 * @param name The class name assigned in the Python wrapper space.
 * @param class_doc The class docstring: a one-line brief of the solver and
 *        its References section.
 * @param doc One sentence naming the TSolverOptions fields this solver reads.
 */
template<typename Solver>
py::class_<PySolverWrapper<Solver>> bind_standard_tsolver(py::module& m, const std::string& name,
                                                          const char* class_doc,
                                                          const std::string& doc = {}) {
    using trex::tsolvers::ad::TSolverOptions;
    py::class_<PySolverWrapper<Solver>> cls(m, name.c_str(), class_doc);
    const std::string init_doc =
        "Initializes the solver via zero-copy views of X, D and y; every other input "
        "comes from `options` (SolverOptions). " + doc + kInPlaceContractDoc;
    cls.def(py::init<Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXd>,
                      Eigen::Ref<Eigen::VectorXd>, const TSolverOptions&>(),
             py::arg("X"),
             py::arg("D"),
             py::arg("y"),
             py::arg("options") = TSolverOptions{},
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>(),
             init_doc.c_str());
    cls.def(py::init([](const DesignArray& X, const DesignArray& D, Eigen::Ref<Eigen::VectorXd> y,
                        const TSolverOptions& options) {
                 auto [X_ref, D_ref] = solver_design_refs(X, D, options);
                 return std::make_unique<PySolverWrapper<Solver>>(X_ref, D_ref, y, options);
             }),
             py::arg("X").noconvert(), py::arg("D").noconvert(), py::arg("y"),
             py::arg("options") = TSolverOptions{},
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>(),
             kReaderSolverDoc.c_str());
    add_standard_tsolver_methods<Solver>(cls);
    return cls;
}


/**
 * @brief Registration for a solver whose constructor also takes the group
 *        vector (X, D, y, groups, options): TIENET and TIENETAug.
 *        `class_doc` is the class docstring, as in bind_standard_tsolver().
 */
template<typename Solver>
py::class_<PySolverWrapper<Solver>> bind_groups_tsolver(py::module& m, const std::string& name,
                                                        const char* class_doc,
                                                        const std::string& doc = {}) {
    using trex::tsolvers::ad::TSolverOptions;
    py::class_<PySolverWrapper<Solver>> cls(m, name.c_str(), class_doc);
    const std::string init_doc =
        "Initializes the solver via zero-copy views of X, D and y with a 0-based group id per "
        "original variable (`groups`, length p); every other input comes from `options` "
        "(SolverOptions). " + doc + kInPlaceContractDoc;
    cls.def(py::init<Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXd>,
                      Eigen::Ref<Eigen::VectorXd>, const Eigen::VectorXi&, const TSolverOptions&>(),
             py::arg("X"), py::arg("D"), py::arg("y"), py::arg("groups"),
             py::arg("options") = TSolverOptions{},
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>(),
             init_doc.c_str());
    cls.def(py::init([](const DesignArray& X, const DesignArray& D, Eigen::Ref<Eigen::VectorXd> y,
                        const Eigen::VectorXi& groups, const TSolverOptions& options) {
                 auto [X_ref, D_ref] = solver_design_refs(X, D, options);
                 return std::make_unique<PySolverWrapper<Solver>>(X_ref, D_ref, y, groups, options);
             }),
             py::arg("X").noconvert(), py::arg("D").noconvert(), py::arg("y"), py::arg("groups"),
             py::arg("options") = TSolverOptions{},
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>(),
             kReaderSolverDoc.c_str());
    add_standard_tsolver_methods<Solver>(cls);
    return cls;
}


/**
 * @brief Registration for a solver whose constructor also takes the sparse
 *        Tikhonov matrix (X, D, y, tikhonov_K, options): TIKNET.
 *        K is passed as a dense p x p array and converted with
 *        sparseView() (no SciPy dependency). `class_doc` is the class
 *        docstring, as in bind_standard_tsolver().
 */
template<typename Solver>
py::class_<PySolverWrapper<Solver>> bind_tikhonov_tsolver(py::module& m, const std::string& name,
                                                          const char* class_doc,
                                                          const std::string& doc = {}) {
    using trex::tsolvers::ad::TSolverOptions;
    py::class_<PySolverWrapper<Solver>> cls(m, name.c_str(), class_doc);
    const std::string init_doc =
        "Initializes the solver via zero-copy views of X, D and y with the p x p PSD Tikhonov "
        "matrix K (a dense array; converted to sparse internally); every other input comes "
        "from `options` (SolverOptions). " + doc + kInPlaceContractDoc;
    cls.def(py::init([](Eigen::Ref<Eigen::MatrixXd> X, Eigen::Ref<Eigen::MatrixXd> D,
                        Eigen::Ref<Eigen::VectorXd> y, const Eigen::MatrixXd& tikhonov_K,
                        const TSolverOptions& options) {
                 const Eigen::SparseMatrix<double> K = tikhonov_K.sparseView();
                 return std::make_unique<PySolverWrapper<Solver>>(X, D, y, K, options);
             }),
             py::arg("X"), py::arg("D"), py::arg("y"), py::arg("tikhonov_K"),
             py::arg("options") = TSolverOptions{},
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>(),
             init_doc.c_str());
    cls.def(py::init([](const DesignArray& X, const DesignArray& D, Eigen::Ref<Eigen::VectorXd> y,
                        const Eigen::MatrixXd& tikhonov_K, const TSolverOptions& options) {
                 const Eigen::SparseMatrix<double> K = tikhonov_K.sparseView();
                 auto [X_ref, D_ref] = solver_design_refs(X, D, options);
                 return std::make_unique<PySolverWrapper<Solver>>(X_ref, D_ref, y, K, options);
             }),
             py::arg("X").noconvert(), py::arg("D").noconvert(), py::arg("y"), py::arg("tikhonov_K"),
             py::arg("options") = TSolverOptions{},
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>(),
             kReaderSolverDoc.c_str());
    add_standard_tsolver_methods<Solver>(cls);
    return cls;
}


/**
 * @brief Bind trex::tsolvers::ad::TSolverOptions with a keyword constructor
 *        (`TSolverOptions(lambda2=0.5, normalize=False)`) and read/write fields.
 */
inline void bind_solver_options(py::module& m) {
    using trex::tsolvers::ad::TSolverOptions;
    using trex::tsolvers::ad::NCGMPVariant;
    using ScalingMode = trex::tsolvers::ad::ScalingMode;

    py::class_<TSolverOptions> cls(m, "SolverOptions",
        "Every T-solver input except the data X, D and y. Each solver reads the fields it "
        "needs and ignores the rest; the defaults are the solver-level defaults. Construct "
        "with keywords (`SolverOptions(lambda2=0.5, normalize=False)`) or set attributes.");
    cls.def(py::init<>(), "All defaults.")
       .def(py::init([](py::kwargs kwargs) {
                TSolverOptions o;
                py::object self = py::cast(&o, py::return_value_policy::reference);
                for (auto item : kwargs) {
                    const std::string key = py::str(item.first);
                    if (!py::hasattr(self, key.c_str())) {
                        throw py::type_error("SolverOptions: unknown option '" + key + "'");
                    }
                    py::setattr(self, key.c_str(), item.second);
                }
                return o;
            }), "Keyword constructor: every keyword must name a SolverOptions field.")
       .def_readwrite("normalize", &TSolverOptions::normalize,
            "Scale each column of X (and D) to unit L2 norm (or unit SD, see scaling_mode), "
            "IN PLACE through the views (the solver becomes a writer and claims the buffers "
            "until restore() or destruction). False when the caller pre-scaled. Default True.")
       .def_readwrite("intercept", &TSolverOptions::intercept,
            "Center X and D in place (a writer, as for normalize) and the solver's internal "
            "copy of y (an intercept is fitted). False when the caller pre-centered. Default True.")
       .def_readwrite("precentered", &TSolverOptions::precentered,
            "The caller centered X, D and y (intercept=False on centered data, as in the T-Rex "
            "selectors). The centered columns span at most n - 1 dimensions, so the solver counts "
            "the intercept as with intercept=True: the active set is capped at n - 1 and the "
            "degrees of freedom and Cp include the intercept. The data are not touched. Redundant "
            "when intercept is True. Default False.")
       .def_readwrite("scaling_mode", &TSolverOptions::scaling_mode,
            "Column scaling applied when normalize is True: ScalingMode.L2 (default) or ZSCORE.")
       .def_readwrite("verbose", &TSolverOptions::verbose, "Print status and diagnostics. Default False.")
       .def_readwrite("tol", &TSolverOptions::tol,
            "Numerical tolerance of the solver's internal comparisons. Default: machine epsilon. A\n"
            "selector sets its own tolerance on the solvers it runs (1e-6 by default).")
       .def_readwrite("tie_seed", &TSolverOptions::tie_seed,
            "Seed of the tie-break engine. Default 0: a fixed seed, so a solver run is reproducible "
            "also where exact ties occur; < 0 seeds from std::random_device (nondeterministic).")
       .def_readwrite("n_threads", &TSolverOptions::n_threads,
            "Thread budget of the solver's own passes over the design (correlation refreshes, "
            "step-length search, standardization). Default 1 (serial); a whole number of at least "
            "1. A pass splits only above about a million entries and never inside a parallel region "
            "of the caller. The T-Rex selectors set it from max_inner_threads.")
       .def_readwrite("exch_tie_alpha", &TSolverOptions::exch_tie_alpha,
            "Exchangeable-tie band width of the greedy solvers (TOMP, TOOLS, TGP, TACGP, TMP, TNCGMP and TAFS); "
            "0 (default) = off. The path solvers ignore it.")
       .def_readwrite("exch_tie_floor", &TSolverOptions::exch_tie_floor,
            "Minimum |correlation| for exchangeable-tie candidates in [0, 1). Default 0.5.")
       .def_readwrite("lambda2", &TSolverOptions::lambda2,
            "L2 (ridge / group-ridge / Tikhonov) penalty weight, >= 0; read by the elastic-net, "
            "IEN and Tikhonov solvers. Default 0.1, as inside the selectors; 0 means no quadratic penalty.")
       .def_readwrite("rho_afs", &TSolverOptions::rho_afs,
            "TAFS shrinkage step in (0, 1]. Default 0.3, as inside the selectors; 1 is forward stepwise.")
       .def_readwrite("ncgmp_variant", &TSolverOptions::ncgmp_variant,
            "TNCGMP variant (NCGMPVariant). Default LineSearch.")
       .def_readwrite("tenet_aug_use_lars", &TSolverOptions::tenet_aug_use_lars,
            "TENETAug / TIENETAug: pure-LARS inner path that never drops variables (True) instead of the LASSO inner "
                "path (False, default).")
       .def_readwrite("tikhonov_validated", &TSolverOptions::tikhonov_validated,
            "TIKNET: True says that the Tikhonov matrix was checked by the caller (finite "
            "entries, symmetric, positive semi-definite), and the solver skips its own check. For a caller "
            "that builds many solvers on one matrix, as the T-Rex selector does. Default False.")
       .def("__repr__", [](const TSolverOptions& o) {
            return "SolverOptions(normalize=" + std::string(o.normalize ? "True" : "False") +
                   ", intercept=" + (o.intercept ? "True" : "False") +
                   ", precentered=" + (o.precentered ? "True" : "False") +
                   ", lambda2=" + std::to_string(o.lambda2) +
                   ", rho_afs=" + std::to_string(o.rho_afs) + ", ...)";
       });
}


/**
 * @brief Binds all T-Solvers to the given Pybind11 module.
 *
 * @param m The Pybind11 module where the solvers will be exposed.
 */
inline void bind_tsolvers_module(py::module& m) {

    using ScalingMode = trex::tsolvers::ad::solver_utils::preprocessing::ScalingMode;

    // Column-scaling convention used by all solver constructors. This is the
    // solver-side enum (solver_utils::preprocessing::ScalingMode); the
    // selector-level trex_selector_methods.ScalingMode mirrors it.
    py::enum_<ScalingMode>(m, "ScalingMode",
                           "Column scaling convention applied by the solver after centering.")
        .value("L2", ScalingMode::L2,
               "Divide each centered column by its L2 norm (unit L2 norm).")
        .value("ZSCORE", ScalingMode::ZSCORE,
               "Divide each centered column by its sample SD (n-1).")
        .export_values();

    // Sub-modules first: the option enums live where their solvers live and
    // TSolverOptions (bound next) refers to them.
    py::module m_lars = m.def_submodule("lars_based", "LARS-based T-Solvers");
    py::module m_omp  = m.def_submodule("omp_based",  "OMP-based T-Solvers");
    py::module m_afs  = m.def_submodule("afs_based",  "AFS-based T-Solvers");

    /**
     * @brief Variants defining the internal optimization mode utilized by the TNCGMP solver.
     */
    py::enum_<trex::tsolvers::ad::NCGMPVariant>(m_omp, "NCGMPVariant")
        .value("LineSearch",
               trex::tsolvers::ad::NCGMPVariant::LineSearch,
               "Computes updates using standard backtracking / exact line search.")
        .value("FullyCorrective",
              trex::tsolvers::ad::NCGMPVariant::FullyCorrective,
              "Computes fully corrective projections for all active set members at each step.")
        .export_values();

    // Every input except X, D, y.
    bind_solver_options(m);

    // -------------------------------------------------------------------------
    // LARS-based sub-module
    // -------------------------------------------------------------------------
    using namespace trex::tsolvers::ad::linear_model::lars_based;

    bind_standard_tsolver<TSolverLARS>(m_lars, "TLARS_Solver",
        "Terminating least angle regression (T-LARS) solver.\n"
        "\n"
        "References\n"
        "----------\n"
        "B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani, \"Least Angle Regression,\" Ann.\n"
        "  Statist., vol. 32, no. 2, pp. 407-499, 2004.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n",
        "Plain T-LARS reads no family hyperparameter.")
        .def("getLambda", &PySolverWrapper<TSolverLARS>::getLambda,
             "Get the regularization parameter sequence (max correlation per step).")
        .def("getKahanRefreshInterval", &PySolverWrapper<TSolverLARS>::getKahanRefreshInterval,
             "Get the current Kahan summation refresh interval.")
        .def("setKahanRefreshInterval", &PySolverWrapper<TSolverLARS>::setKahanRefreshInterval,
             py::arg("interval"),
             "Set the Kahan summation refresh interval (call before executeStep).");

    bind_standard_tsolver<TSolverLASSO>(m_lars, "TLASSO_Solver",
        "Terminating LASSO (T-LASSO) solver: the LARS-LASSO path with variable removals.\n"
        "\n"
        "References\n"
        "----------\n"
        "R. Tibshirani, \"Regression Shrinkage and Selection via the Lasso,\" J. R. Stat. Soc. Ser.\n"
        "  B Stat. Methodol., vol. 58, no. 1, pp. 267-288, 1996.\n"
        "B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani, \"Least Angle Regression,\" Ann.\n"
        "  Statist., vol. 32, no. 2, pp. 407-499, 2004.\n"
        "M. R. Osborne, B. Presnell, and B. A. Turlach, \"On the LASSO and Its Dual,\" J. Comput.\n"
        "  Graph. Statist., vol. 9, no. 2, pp. 319-337, 2000.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n")
        .def("getLambda", &PySolverWrapper<TSolverLASSO>::getLambda,
             "Get the regularization parameter sequence (max correlation per step).")
        .def("getNumRemovals", &PySolverWrapper<TSolverLASSO>::getNumRemovals,
             "Get the total number of variable removals during solution path.")
        .def("getCyclingRatio", &PySolverWrapper<TSolverLASSO>::getCyclingRatio,
             "Get the cycling ratio (removals / additions).");

    bind_standard_tsolver<TSolverStepwise>(m_lars, "TSTEPWISE_Solver",
        "Terminating forward stepwise regression (T-Stepwise) solver.\n"
        "\n"
        "References\n"
        "----------\n"
        "B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani, \"Least Angle Regression,\" Ann.\n"
        "  Statist., vol. 32, no. 2, pp. 407-499, 2004.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n")
        .def("getLambda", &PySolverWrapper<TSolverStepwise>::getLambda,
             "Get the regularization parameter sequence (max correlation per step).");

    bind_standard_tsolver<TSolverStagewise>(m_lars, "TSTAGEWISE_Solver",
        "Terminating incremental forward stagewise (T-Stagewise) solver via the LARS-stagewise path.\n"
        "\n"
        "References\n"
        "----------\n"
        "B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani, \"Least Angle Regression,\" Ann.\n"
        "  Statist., vol. 32, no. 2, pp. 407-499, 2004.\n"
        "T. Hastie, J. Taylor, R. Tibshirani, and G. Walther, \"Forward Stagewise Regression and\n"
        "  the Monotone Lasso,\" Electron. J. Stat., vol. 1, pp. 1-29, 2007,\n"
        "  doi: 10.1214/07-EJS004.\n"
        "R. J. Tibshirani, \"A General Framework for Fast Stagewise Algorithms,\" J. Mach. Learn.\n"
        "  Res., vol. 16, no. 78, pp. 2543-2588, 2015.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n")
        .def("getLambda", &PySolverWrapper<TSolverStagewise>::getLambda,
             "Get the regularization parameter sequence (max correlation per step).")
        .def("getNumRemovals", &PySolverWrapper<TSolverStagewise>::getNumRemovals,
             "Get the total number of coefficient removals.")
        .def("getCyclingRatio", &PySolverWrapper<TSolverStagewise>::getCyclingRatio,
             "Get the cycling ratio (removals / additions).");

    bind_standard_tsolver<TSolverENET>(m_lars, "TENET_Solver",
        "Terminating elastic net (T-ENET) solver via the native T-LARS-EN path (L1 + L2 penalty).\n"
        "\n"
        "References\n"
        "----------\n"
        "H. Zou and T. Hastie, \"Regularization and Variable Selection via the Elastic Net,\" J. R.\n"
        "  Stat. Soc. Ser. B Stat. Methodol., vol. 67, no. 2, pp. 301-320, 2005.\n"
        "B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani, \"Least Angle Regression,\" Ann.\n"
        "  Statist., vol. 32, no. 2, pp. 407-499, 2004.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n",
        "Reads options.lambda2, the secondary L2 penalty (default 0: the lasso path).")
        .def("getLambda", &PySolverWrapper<TSolverENET>::getLambda,
             "Get the regularization parameter sequence: an initial entry followed by the decreasing L1 "
             "penalty per step, so its length is the number of steps plus one.")
        .def("getAlpha", &PySolverWrapper<TSolverENET>::getAlpha,
             "Get the Elastic Net mixing parameter alpha implied by lambda2.")
        .def("getNumRemovals", &PySolverWrapper<TSolverENET>::getNumRemovals,
             "Get the total number of variable removals during solution path.")
        .def("getCyclingRatio", &PySolverWrapper<TSolverENET>::getCyclingRatio,
             "Get the cycling ratio (removals / additions).");

    bind_standard_tsolver<TSolverENETAug>(m_lars, "TENETAug_Solver",
        "Terminating elastic net (T-ENETAug) solver via row augmentation around an inner T-LASSO.\n"
        "\n"
        "References\n"
        "----------\n"
        "H. Zou and T. Hastie, \"Regularization and Variable Selection via the Elastic Net,\" J. R.\n"
        "  Stat. Soc. Ser. B Stat. Methodol., vol. 67, no. 2, pp. 301-320, 2005.\n"
        "B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani, \"Least Angle Regression,\" Ann.\n"
        "  Statist., vol. 32, no. 2, pp. 407-499, 2004.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"False Discovery Rate Control for Grouped\n"
        "  Variable Selection in High-Dimensional Linear Models Using the T-Knock Filter,\" in\n"
        "  Proc. 30th Eur. Signal Process. Conf. (EUSIPCO), 2022, pp. 892-896.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n",
        "Reads options.lambda2 and options.tenet_aug_use_lars (pure-LARS inner path).");

    // TIENETAug: Informed Elastic Net via augmented LASSO (group-ridge penalty);
    // the group vector is data, not an option.
    py::class_<PySolverWrapper<TSolverIENETAug>> cls_tienet_aug(m_lars, "TIENETAug_Solver",
        "Terminating informed elastic net (T-IENETAug) solver via row augmentation around an inner\n"
        "T-LASSO.\n"
        "\n"
        "References\n"
        "----------\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Informed Elastic Net for Fast Grouped\n"
        "  Variable Selection and FDR Control in Genomics Research,\" in Proc. IEEE 9th Int.\n"
        "  Workshop Comput. Adv. Multi-Sensor Adapt. Process. (CAMSAP), 2023, pp. 466-470.\n"
        "H. Zou and T. Hastie, \"Regularization and Variable Selection via the Elastic Net,\" J. R.\n"
        "  Stat. Soc. Ser. B Stat. Methodol., vol. 67, no. 2, pp. 301-320, 2005.\n"
        "B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani, \"Least Angle Regression,\" Ann.\n"
        "  Statist., vol. 32, no. 2, pp. 407-499, 2004.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n");
    cls_tienet_aug.def(py::init<Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXd>,
                                Eigen::Ref<Eigen::VectorXd>, const Eigen::VectorXi&,
                                const trex::tsolvers::ad::TSolverOptions&>(),
             py::arg("X"),
             py::arg("D"),
             py::arg("y"),
             py::arg("groups"),
             py::arg("options") = trex::tsolvers::ad::TSolverOptions{},
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>(),
             ("Initializes the augmented-LASSO Informed Elastic Net solver with a 0-based group"
              " id per original variable (`groups`, length p); reads options.lambda2 (the"
              " group-ridge penalty) and options.tenet_aug_use_lars." + kInPlaceContractDoc).c_str());
    cls_tienet_aug.def(py::init([](const DesignArray& X, const DesignArray& D, Eigen::Ref<Eigen::VectorXd> y,
                                   const Eigen::VectorXi& groups,
                                   const trex::tsolvers::ad::TSolverOptions& options) {
                 auto [X_ref, D_ref] = solver_design_refs(X, D, options);
                 return std::make_unique<PySolverWrapper<TSolverIENETAug>>(X_ref, D_ref, y, groups, options);
             }),
             py::arg("X").noconvert(), py::arg("D").noconvert(), py::arg("y"), py::arg("groups"),
             py::arg("options") = trex::tsolvers::ad::TSolverOptions{},
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>(),
             kReaderSolverDoc.c_str());
    add_standard_tsolver_methods<TSolverIENETAug>(cls_tienet_aug);
    cls_tienet_aug
        .def("getLambda2", &PySolverWrapper<TSolverIENETAug>::getLambda2,
             "Get the group-ridge penalty lambda2 used.")
        .def("getGroups", &PySolverWrapper<TSolverIENETAug>::getGroups,
             "Get the 0-based group id per original variable (length p).")
        .def("getNumGroups", &PySolverWrapper<TSolverIENETAug>::getNumGroups,
             "Get the number of disjoint groups.");

    // TIENET: native pathwise LARS-IEN (no augmented copies); groups are data.
    bind_groups_tsolver<TSolverIENET>(m_lars, "TIENET_Solver",
        "Terminating informed elastic net (T-IENET) solver via the native pathwise LARS-IEN algorithm.\n"
        "\n"
        "References\n"
        "----------\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Informed Elastic Net for Fast Grouped\n"
        "  Variable Selection and FDR Control in Genomics Research,\" in Proc. IEEE 9th Int.\n"
        "  Workshop Comput. Adv. Multi-Sensor Adapt. Process. (CAMSAP), 2023, pp. 466-470.\n"
        "H. Zou and T. Hastie, \"Regularization and Variable Selection via the Elastic Net,\" J. R.\n"
        "  Stat. Soc. Ser. B Stat. Methodol., vol. 67, no. 2, pp. 301-320, 2005.\n"
        "B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani, \"Least Angle Regression,\" Ann.\n"
        "  Statist., vol. 32, no. 2, pp. 407-499, 2004.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n",
        "Reads options.lambda2, the group-ridge penalty.")
        .def("getLambda", &PySolverWrapper<TSolverIENET>::getLambda,
             "Get the regularization parameter sequence (max correlation per step).")
        .def("getLambda2", &PySolverWrapper<TSolverIENET>::getLambda2,
             "Get the group-ridge penalty lambda2 used.")
        .def("getGroups", &PySolverWrapper<TSolverIENET>::getGroups,
             "Get the 0-based group id per original variable (length p).")
        .def("getNumGroups", &PySolverWrapper<TSolverIENET>::getNumGroups,
             "Get the number of disjoint groups.");

    // TIKNET: native LARS homotopy under lambda2 * beta' K beta (FOLDED dummies).
    bind_tikhonov_tsolver<TSolverTIKNET>(m_lars, "TIKNET_Solver",
        "Terminating Tikhonov net (T-TIKNET) solver via the native pathwise LARS-TIK algorithm.\n"
        "\n"
        "References\n"
        "----------\n"
        "A. N. Tikhonov, \"Solution of Incorrectly Formulated Problems and the Regularization\n"
        "  Method,\" Soviet Math. Dokl., vol. 4, pp. 1035-1038, 1963 (English translation of Dokl.\n"
        "  Akad. Nauk SSSR 151 (1963), no. 3, 501-504).\n"
        "Y. Park, L. Reichel, G. Rodriguez, and X. Yu, \"Parameter Determination for Tikhonov\n"
        "  Regularization Problems in General Form,\" J. Comput. Appl. Math., vol. 343, pp. 12-25,\n"
        "  2018, doi: 10.1016/j.cam.2018.04.049.\n"
        "A. E. Hoerl and R. W. Kennard, \"Ridge Regression: Biased Estimation for Nonorthogonal\n"
        "  Problems,\" Technometrics, vol. 42, no. 1, pp. 80-86, 2000 (Reprint of the 1970 article\n"
        "  in the Special 40th Anniversary Issue).\n"
        "Y. Miche, M. van Heeswijk, P. Bas, O. Simula, and A. Lendasse, \"TROP-ELM: A\n"
        "  Double-Regularized ELM Using LARS and Tikhonov Regularization,\" Neurocomputing,\n"
        "  vol. 74, pp. 2413-2421, 2011, doi: 10.1016/j.neucom.2010.12.042.\n"
        "B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani, \"Least Angle Regression,\" Ann.\n"
        "  Statist., vol. 32, no. 2, pp. 407-499, 2004.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n",
        "Reads options.lambda2, the quadratic penalty weight; FOLDED dummy coupling only "
        "(D.cols() a multiple of X.cols() when lambda2 > 0).")
        .def("getLambda", &PySolverWrapper<TSolverTIKNET>::getLambda,
             "Get the regularization parameter sequence (max correlation per step).")
        .def("getLambda2", &PySolverWrapper<TSolverTIKNET>::getLambda2,
             "Get the quadratic penalty weight lambda2 used.")
        .def("getTikhonovK", &PySolverWrapper<TSolverTIKNET>::getTikhonovK,
             "The Tikhonov matrix K as a dense p x p array.");

    // -------------------------------------------------------------------------
    // OMP-based sub-module
    // -------------------------------------------------------------------------
    using namespace trex::tsolvers::ad::linear_model::omp_based;

    bind_standard_tsolver<TSolverOMP>(m_omp, "TOMP_Solver",
        "Terminating orthogonal matching pursuit (T-OMP) solver.\n"
        "\n"
        "References\n"
        "----------\n"
        "Y. C. Pati, R. Rezaiifar, and P. S. Krishnaprasad, \"Orthogonal Matching Pursuit:\n"
        "  Recursive Function Approximation with Applications to Wavelet Decomposition,\" in Proc.\n"
        "  27th Asilomar Conf. Signals, Syst., Comput., 1993, pp. 40-44.\n"
        "J. A. Tropp, \"Greed is Good: Algorithmic Results for Sparse Approximation,\" IEEE Trans.\n"
        "  Inf. Theory, vol. 50, no. 10, pp. 2231-2242, 2004, doi: 10.1109/TIT.2004.834793.\n"
        "T. Blumensath and M. E. Davies, \"On the Difference Between Orthogonal Matching Pursuit\n"
        "  and Orthogonal Least Squares,\" unpublished note, draft of March 29, 2007.\n"
        "F. Scheidt and M. Muma, \"Terminating Greedy Pursuits for High-Dimensional False\n"
        "  Discovery Rate Control,\" in Proc. Eur. Signal Process. Conf. (EUSIPCO), 2026,\n"
        "  pp. 2511-2515.\n",
        "Reads options.exch_tie_alpha / exch_tie_floor, the exchangeable-tie band.");
    bind_standard_tsolver<TSolverGP>(m_omp, "TGP_Solver",
        "Terminating gradient pursuit (T-GP) solver.\n"
        "\n"
        "References\n"
        "----------\n"
        "T. Blumensath and M. E. Davies, \"Gradient Pursuits,\" IEEE Trans. Signal Process.,\n"
        "  vol. 56, no. 6, pp. 2370-2382, 2008.\n"
        "F. Scheidt and M. Muma, \"Terminating Greedy Pursuits for High-Dimensional False\n"
        "  Discovery Rate Control,\" in Proc. Eur. Signal Process. Conf. (EUSIPCO), 2026,\n"
        "  pp. 2511-2515.\n",
        "Reads options.exch_tie_alpha / exch_tie_floor, the exchangeable-tie band.");
    bind_standard_tsolver<TSolverACGP>(m_omp, "TACGP_Solver",
        "Terminating approximate conjugate gradient pursuit (T-ACGP) solver.\n"
        "\n"
        "References\n"
        "----------\n"
        "T. Blumensath and M. E. Davies, \"Gradient Pursuits,\" IEEE Trans. Signal Process.,\n"
        "  vol. 56, no. 6, pp. 2370-2382, 2008.\n"
        "F. Scheidt and M. Muma, \"Terminating Greedy Pursuits for High-Dimensional False\n"
        "  Discovery Rate Control,\" in Proc. Eur. Signal Process. Conf. (EUSIPCO), 2026,\n"
        "  pp. 2511-2515.\n",
        "Reads options.exch_tie_alpha / exch_tie_floor, the exchangeable-tie band.");
    bind_standard_tsolver<TSolverMP>(m_omp, "TMP_Solver",
        "Terminating matching pursuit (T-MP) solver.\n"
        "\n"
        "References\n"
        "----------\n"
        "S. G. Mallat and Z. Zhang, \"Matching Pursuits with Time-Frequency Dictionaries,\" IEEE\n"
        "  Trans. Signal Process., vol. 41, no. 12, pp. 3397-3415, 1993, doi: 10.1109/78.258082.\n"
        "F. Scheidt and M. Muma, \"Terminating Greedy Pursuits for High-Dimensional False\n"
        "  Discovery Rate Control,\" in Proc. Eur. Signal Process. Conf. (EUSIPCO), 2026,\n"
        "  pp. 2511-2515.\n",
        "Reads options.exch_tie_alpha / exch_tie_floor, the exchangeable-tie band.");
    bind_standard_tsolver<TSolverOOLS>(m_omp, "TOOLS_Solver",
        "Terminating orthogonal least squares (T-OLS) solver.\n"
        "\n"
        "References\n"
        "----------\n"
        "S. Chen, S. A. Billings, and W. Luo, \"Orthogonal Least Squares Methods and Their\n"
        "  Application to Non-Linear System Identification,\" Int. J. Control, vol. 50, no. 5,\n"
        "  pp. 1873-1896, 1989, doi: 10.1080/00207178908953472.\n"
        "T. Blumensath and M. E. Davies, \"On the Difference Between Orthogonal Matching Pursuit\n"
        "  and Orthogonal Least Squares,\" unpublished note, draft of March 29, 2007.\n"
        "F. Scheidt and M. Muma, \"Terminating Greedy Pursuits for High-Dimensional False\n"
        "  Discovery Rate Control,\" in Proc. Eur. Signal Process. Conf. (EUSIPCO), 2026,\n"
        "  pp. 2511-2515.\n",
        "Reads options.exch_tie_alpha / exch_tie_floor, the exchangeable-tie band.");
    bind_standard_tsolver<TSolverNCGMP>(m_omp, "TNCGMP_Solver",
        "Terminating norm-corrective generalized matching pursuit (T-NCGMP) solver.\n"
        "\n"
        "References\n"
        "----------\n"
        "F. Locatello, R. Khanna, M. Tschannen, and M. Jaggi, \"A Unified Optimization View on\n"
        "  Generalized Matching Pursuit and Frank-Wolfe,\" in Proc. 20th Int. Conf. Artif. Intell.\n"
        "  Statist. (AISTATS), Proc. Mach. Learn. Res. vol. 54, 2017, pp. 860-868.\n"
        "F. Scheidt and M. Muma, \"Terminating Greedy Pursuits for High-Dimensional False\n"
        "  Discovery Rate Control,\" in Proc. Eur. Signal Process. Conf. (EUSIPCO), 2026,\n"
        "  pp. 2511-2515.\n",
        "Reads options.ncgmp_variant (LineSearch or FullyCorrective) and options.exch_tie_alpha / "
        "exch_tie_floor, the exchangeable-tie band.");

    // -------------------------------------------------------------------------
    // AFS-based sub-module
    // -------------------------------------------------------------------------
    using namespace trex::tsolvers::ad::linear_model::afs_based;

    bind_standard_tsolver<TSolverAFS>(m_afs, "TAFS_Solver",
        "Terminating adaptive forward stepwise (T-AFS) solver.\n"
        "\n"
        "References\n"
        "----------\n"
        "I. Zhang and R. Tibshirani, \"Adaptive Forward Stepwise: A Method for High Sparsity\n"
        "  Regression,\" J. Mach. Learn. Res., vol. 27, no. 35, pp. 1-24, 2026.\n"
        "F. Scheidt and M. Muma, \"Terminating Greedy Pursuits for High-Dimensional False\n"
        "  Discovery Rate Control,\" in Proc. Eur. Signal Process. Conf. (EUSIPCO), 2026,\n"
        "  pp. 2511-2515.\n",
        "Reads options.rho_afs, the shrinkage step in (0, 1] (1: forward stepwise), and the"
        " exchangeable-tie fields.");
}
// =====================================================================================
#endif /* End of TREX_PYTHON_TSOLVER_BINDINGS_HPP */
// =====================================================================================

// =====================================================================================
// trex_da_bindings.hpp - Pybind11 bindings for TRexDASelector
// =====================================================================================
#ifndef TREX_PYTHON_TREX_DA_BINDINGS_HPP
#define TREX_PYTHON_TREX_DA_BINDINGS_HPP
// =====================================================================================
/**
 * @file trex_da_bindings.hpp
 * @brief Pybind11 bindings for TRexDASelector.
 *
 */
// =====================================================================================

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

// trex includes
#include <trex_selector_methods/ad/trex_da/trex_da.hpp>
#include "trex_core_bindings.hpp"

// =====================================================================================

namespace py = pybind11;
using namespace trex::trex_selector_methods::ad::trex_da;
using namespace trex::trex_selector_methods::ad::trex_core;

// =====================================================================================

/**
 * @brief Python wrapper class for TRexDASelector leveraging zero-copy references.
 */
class PyTRexDASelector : public PyTRexSelector {
public:
    /**
     * @brief Zero-copy constructor for the PyTRexDASelector wrapper.
     *
     * @param X Reference to the design matrix (n × p).
     * @param y Reference to the response vector (n × 1).
     * @param tFDR Target False Discovery Rate threshold.
     * @param da_control Configuration parameters of the dependency-aware correction.
     * @param trex_control General configuration parameters for the TRex algorithm.
     * @param seed Random seed for reproducible computations.
     * @param verbose If true, prompts standard out messaging during select().
     */
    PyTRexDASelector(
        Eigen::Ref<Eigen::MatrixXd> X,
        Eigen::Ref<Eigen::VectorXd> y,
        double tFDR,
        const TRexDAControlParameter& da_control,
        const TRexControlParameter& trex_control,
        int seed,
        bool verbose
    ) {
        warn_if_cd_solver_type(trex_control.solver_type);

        this->X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X.data(),
                                                                  X.rows(),
                                                                  X.cols());
        this->y_map_ = std::make_unique<Eigen::Map<Eigen::VectorXd>>(y.data(),
                                                                  y.size());

        // The refactored C++ constructor takes a single DA control struct that
        // nests the base algorithmic parameters as `trex_ctrl`. Merge the
        // separately-supplied `trex_control` into that nested field.
        TRexDAControlParameter da_ctrl = da_control;
        da_ctrl.trex_ctrl = trex_control;

        this->selector_ = std::make_unique<TRexDASelector>(
            *(this->X_map_), *(this->y_map_), tFDR, da_ctrl, seed, verbose
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
    PyTRexDASelector(
        const DesignArray& X,
        Eigen::Ref<Eigen::VectorXd> y,
        double tFDR,
        const TRexDAControlParameter& da_control,
        const TRexControlParameter& trex_control,
        int seed,
        bool verbose
    ) : PyTRexDASelector(selector_design_ref(X, trex_control), y, tFDR, da_control, trex_control, seed,
                         verbose) {}

    /**
     * @brief Get the complete DA selection result.
     */
    const TRexDASelector::DASelectionResult& getDAResult() const {
        return static_cast<TRexDASelector*>(this->selector_.get())->getDAResult();
    }
};


/**
 * @brief Binds the TRexDASelector and related configurations.
 *
 * @param m The Python module to bind against.
 */
inline void bind_trex_da(py::module& m) {
    py::enum_<DAMethod>(m, "DAMethod", "Dependency-aware (DA) methods.")
        .value("AR1", DAMethod::AR1, "Made for a first-order autoregressive dependency across the columns. A "
            "variable is compared with a window of neighbouring variables.")
        .value("EQUI", DAMethod::EQUI, "Made for variables that all share one correlation. A variable is compared "
            "with all other variables.")
        .value("BT", DAMethod::BT, "Binary tree, made for disjoint groups of highly correlated variables. The "
            "groups come from a hierarchical clustering.")
        .value("NN", DAMethod::NN, "Nearest neighbours, made for overlapping groups of highly correlated "
            "variables. The groups come from a sweep of correlation thresholds.")
        .value("PRIOR_GROUPS", DAMethod::PRIOR_GROUPS, "Groups supplied by the user, clustered within each group.")
        .export_values();

    py::enum_<BTSelectionMode>(m, "BTSelectionMode", "Cell-selection policy for the BT (binary-tree) DA method.")
        .value("FeasibleOnly", BTSelectionMode::FeasibleOnly, "Only FDP-feasible cells may be selected (default; "
            "controls realized FDR).")
        .value("RFaithful",    BTSelectionMode::RFaithful,    "Reproduces the R reference; may inflate FDR at high "
            "within-group correlation.")
        .export_values();

    py::class_<TRexDAControlParameter>(m, "TRexDAControlParameter", "Control parameters for TRexDASelector.")
        .def(py::init<>())
        .def_readwrite("method", &TRexDAControlParameter::method, "Dependency-aware method (default: BT). Ignored\n"
            "when prior_groups is set.")
        .def_readwrite("cor_coef", &TRexDAControlParameter::cor_coef, "Correlation coefficient of the AR1 and EQUI "
            "methods; setting the parameter to -2 estimates it from the data. Both methods use the magnitude of a "
            "supplied coefficient. Any other value outside [-1, 1] gives a warning and is replaced by the estimate.")
        .def_readwrite("rho_thr_DA", &TRexDAControlParameter::rho_thr_DA, "Correlation threshold of the AR1 and "
            "EQUI methods (default 0.02). AR1: the window around a variable holds the lags whose correlation "
            "abs(cor_coef)**lag reaches the threshold, range [0, 1). EQUI: the correction is skipped when "
            "abs(cor_coef) does not exceed the threshold, range [0, 1].")
        .def_readwrite("hc_linkage", &TRexDAControlParameter::hc_linkage,
            "Linkage of the hierarchical clustering: Single (default), Complete, Average or WPGMA. Ward is\n"
            "rejected, because the correlation distance 1 - |cor| is not Euclidean.")
        .def_readwrite("hc_grid_length", &TRexDAControlParameter::hc_grid_length, "Grid length for hierarchical "
            "clustering.")
        .def_readwrite("bt_selection_mode", &TRexDAControlParameter::bt_selection_mode, "Cell-selection policy for the "
            "BT method (FeasibleOnly or RFaithful).")
        .def_readwrite("prior_groups", &TRexDAControlParameter::prior_groups, "Prior grouping constraints (label "
            "vectors of length p); sub-clustering runs within their finest common refinement.")
        .def_readwrite("trex_ctrl", &TRexDAControlParameter::trex_ctrl, "Nested base T-Rex algorithmic control "
            "parameters (overridden by the separate trex_control argument to TRexDASelector).");

    py::class_<TRexDASelector::DASelectionResult, TRexSelector::SelectionResult>(m, "DASelectionResult", "Extended "
        "SelectionResult holding specific DA outputs.")
        .def_readonly("rho_thresh", &TRexDASelector::DASelectionResult::rho_thresh, "The threshold applied for rho.")
        .def_readonly("rho_grid", &TRexDASelector::DASelectionResult::rho_grid, "Grid values used for rho.")
        .def_readonly("method",           &TRexDASelector::DASelectionResult::method,           "The resulting DA "
            "method applied.")
        .def_readonly("cor_coef",          &TRexDASelector::DASelectionResult::cor_coef,          "Auto-estimated or "
            "supplied correlation coefficient.")
        .def_readonly("fdp_hat_array_bt",  &TRexDASelector::DASelectionResult::FDP_hat_array_BT,  "Per-rho-level FDP "
            "estimate matrices (BT method only).")
        .def_readonly("phi_array_bt",      &TRexDASelector::DASelectionResult::Phi_array_BT,      "Per-rho-level "
            "deflated Phi matrices (BT method only).")
        .def_readonly("r_array_bt",        &TRexDASelector::DASelectionResult::R_array_BT,        "Per-rho-level R "
            "matrices (BT method only).");

    py::class_<PyTRexDASelector, PyTRexSelector>(m, "TRexDASelector",
        "Dependency-aware T-Rex selector (T-Rex+DA): FDR-controlled variable selection for dependent\n"
        "variables, with a correlation-structure model (AR1, EQUI, BT, NN or PRIOR_GROUPS) that\n"
        "accounts for the dependencies among the variables.\n"
        "\n"
        "Notes\n"
        "-----\n"
        "The penalty of the last step. The dependency-aware relative occurrences of variable j at\n"
        "step t of the T-loop are\n"
        "\n"
        "    Phi^DA_t(j) = Psi_T(j) * Phi_t(j),   t = 1, ..., T,\n"
        "    Psi_T(j)    = 1 / (2 - min_j' |Phi_T(j) - Phi_T(j')|),\n"
        "\n"
        "where T is the last step of the evaluated path and j' runs over the variables that j is\n"
        "compared with (its group, its window or all other variables). A variable without a group\n"
        "gets Psi_T(j) = 1/2.\n"
        "\n"
        "This is a deliberate deviation from the papers below. Read literally, they recompute the\n"
        "penalty at every step, Psi_t(j) from Phi_t. The selector takes the penalty of the last\n"
        "step for all steps, with all five methods. The reason is the FDP estimator. It sums over\n"
        "the increments Phi^DA_t(j) - Phi^DA_(t-1)(j) and counts expected false discoveries only\n"
        "while no increment is negative. One penalty per variable keeps the increments\n"
        "non-negative. A penalty per step does not, because the penalty of a variable changes from\n"
        "step to step, and the FDP estimate then becomes too small.\n"
        "\n"
        "A Monte Carlo comparison of the two constructions (1,000 runs per setting on block, AR(1)\n"
        "and independent designs) supports the choice. With the method NN the per-step penalty\n"
        "exceeded a target FDR of 0.1, with 0.15 on an independent design and 0.13 on a block\n"
        "design, where the last-step penalty gave 0.00 and 0.03. No setting had the last-step\n"
        "penalty above the target and the per-step penalty below it. The per-step penalty selected\n"
        "more variables only where its estimate was invalid, and on AR(1) designs the last-step\n"
        "penalty had the higher true positive rate.\n"
        "\n"
        "The methods AR1 and EQUI used the per-step penalty up to trexselector 1.0.0, so their\n"
        "selections can differ slightly from earlier versions. BT, NN and PRIOR_GROUPS have always\n"
        "used the penalty of the last step.\n"
        "\n"
        "References\n"
        "----------\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"High-Dimensional False Discovery Rate Control\n"
        "  for Dependent Variables,\" Signal Process., vol. 234, Art. no. 109990, 2025,\n"
        "  doi: 10.1016/j.sigpro.2025.109990.\n"
        "J. Machkour, D. P. Palomar, and M. Muma, \"FDR-Controlled Portfolio Optimization for\n"
        "  Sparse Financial Index Tracking,\" in Proc. IEEE Int. Conf. Acoust., Speech, Signal\n"
        "  Process. (ICASSP), 2025, pp. 1-5, doi: 10.1109/ICASSP49660.2025.10890629.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"FDR-Controlled Sparse Index Tracking with\n"
        "  Autoregressive Stock Dependency Models,\" in Proc. Eur. Signal Process. Conf.\n"
        "  (EUSIPCO), 2024, pp. 2662-2666.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n")
        .def(py::init<Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::VectorXd>, double, const TRexDAControlParameter&,
            const TRexControlParameter&, int, bool>(),
             py::arg("X"), py::arg("y"), py::arg("tFDR") = 0.1,
             py::arg("da_control") = TRexDAControlParameter(),
             py::arg("trex_control") = TRexControlParameter(),
             py::arg("seed") = -1, py::arg("verbose") = true,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             "Construct the DA selector. X and y are accessed zero-copy via Eigen::Map.")
        .def(py::init<const DesignArray&, Eigen::Ref<Eigen::VectorXd>, double, const TRexDAControlParameter&,
            const TRexControlParameter&, int, bool>(),
             py::arg("X").noconvert(), py::arg("y"), py::arg("tFDR") = 0.1,
             py::arg("da_control") = TRexDAControlParameter(),
             py::arg("trex_control") = TRexControlParameter(),
             py::arg("seed") = -1, py::arg("verbose") = true,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             kReaderSelectorDoc.c_str())
        .def("getDAResult", &PyTRexDASelector::getDAResult, "Fetch the detailed DA selection result.");
}
// =====================================================================================
#endif /* End of TREX_PYTHON_TREX_DA_BINDINGS_HPP */

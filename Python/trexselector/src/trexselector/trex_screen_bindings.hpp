// =====================================================================================
// trex_screen_bindings.hpp - Pybind11 bindings for ScreenTRexSelector
// =====================================================================================
#ifndef TREX_PYTHON_TREX_SCREEN_BINDINGS_HPP
#define TREX_PYTHON_TREX_SCREEN_BINDINGS_HPP
// =====================================================================================
/**
 * @file trex_screen_bindings.hpp
 *
 * @brief Pybind11 bindings for ScreenTRexSelector.
 *
 */
// =====================================================================================

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

// trex includes
#include <trex_selector_methods/ad/trex_screening/trex_screening.hpp>
#include "trex_core_bindings.hpp"

// =====================================================================================
namespace py = pybind11;
using namespace trex::trex_selector_methods::ad::trex_screening;
using namespace trex::trex_selector_methods::ad::trex_core;

// =====================================================================================

/**
 * @brief Python wrapper class for ScreenTRexSelector leveraging zero-copy references.
 */
class PyScreenTRexSelector : public PyTRexSelector {
public:
    /**
     * @brief Zero-copy constructor for the PyScreenTRexSelector wrapper.
     *
     * @param X Reference to the design matrix (n × p).
     * @param y Reference to the response vector (n × 1).
     * @param screen_control Configuration parameters specific to the Screening TRex method.
     * @param trex_control General configuration parameters for the TRex algorithm.
     * @param seed Random seed for reproducible computations.
     * @param verbose If true, prompts standard out messaging during select().
     */
    PyScreenTRexSelector(
        Eigen::Ref<Eigen::MatrixXd> X,
        Eigen::Ref<Eigen::VectorXd> y,
        const ScreenTRexControlParameter& screen_control,
        const TRexControlParameter& trex_control,
        int seed,
        bool verbose
    ) {
        this->X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X.data(),
                                                                  X.rows(),
                                                                  X.cols());
        this->y_map_ = std::make_unique<Eigen::Map<Eigen::VectorXd>>(y.data(),
                                                                  y.size());

        // The refactored C++ constructor takes a single screening control
        // struct that nests the base algorithmic parameters as `trex_ctrl`.
        // Merge the separately-supplied `trex_control` into that nested field.
        ScreenTRexControlParameter screen_ctrl = screen_control;
        screen_ctrl.trex_ctrl = trex_control;

        this->selector_ = std::make_unique<ScreenTRexSelector>(
            *(this->X_map_), *(this->y_map_), screen_ctrl, seed, verbose
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
    PyScreenTRexSelector(
        const DesignArray& X,
        Eigen::Ref<Eigen::VectorXd> y,
        const ScreenTRexControlParameter& screen_control,
        const TRexControlParameter& trex_control,
        int seed,
        bool verbose
    ) : PyScreenTRexSelector(selector_design_ref(X, trex_control), y, screen_control, trex_control, seed,
                             verbose) {}

    /**
     * @brief Get the complete Screening selection result.
     */
    const ScreenTRexSelectionResult& getScreenResult() const {
        return static_cast<ScreenTRexSelector*>(this->selector_.get())->getScreenResult();
    }
};


/**
 * @brief Binds the ScreenTRexSelector and related configurations.
 *
 * @param m The Python module to bind against.
 */
inline void bind_trex_screen(py::module& m) {
    py::enum_<ScreenTRexMethod>(m, "ScreenTRexMethod", "Variants of the screening: without or with a "
        "dependency-aware filter on the relative occurrences.")
        .value("TREX", ScreenTRexMethod::TREX, "Screening without a dependency-aware filter.")
        .value("TREX_DA_AR1", ScreenTRexMethod::TREX_DA_AR1, "Screening with the dependency-aware filter for a "
            "first-order autoregressive dependency: a variable is compared with the neighbors inside the AR(1) "
            "window.")
        .value("TREX_DA_EQUI", ScreenTRexMethod::TREX_DA_EQUI, "Screening with the dependency-aware filter for "
            "equi-correlated variables: a variable is compared with all other variables.")
        .value("TREX_DA_BLOCK_EQUI", ScreenTRexMethod::TREX_DA_BLOCK_EQUI, "Screening with the dependency-aware "
            "filter for block equi-correlated variables: the same comparison inside contiguous blocks.")
        .value("TREX_DA_NN", ScreenTRexMethod::TREX_DA_NN, "Screening with the dependency-aware filter of the "
            "nearest neighbors: a variable is compared with the neighbors inside the AR(1) window whose sample "
            "correlation reaches rho_thr_DA, and a variable without such a neighbor keeps its relative "
            "occurrence. With a threshold at or above the correlation of the data no neighbor is left, and the "
            "variant selects as TREX.")
        .export_values();

    py::class_<ScreenTRexControlParameter>(m, "ScreenTRexControlParameter", "Control parameters for "
        "TRexScreeningSelector.")
        .def(py::init<>())
        .def_readwrite("use_bootstrap_CI", &ScreenTRexControlParameter::use_bootstrap_CI, "Whether the "
            "confidence-based screening (bootstrap) makes the selection (default False). The FDR estimate stays one "
            "over the number of selected variables. With a dependency-aware variant it mostly selects nothing: the "
            "filter lowers the ordinary count that the confidence level has to reach, and in simulations no level "
            "reached it in 90 % to 100 % of the runs.")
        .def_readwrite("R_boot", &ScreenTRexControlParameter::R_boot, "Number of bootstrap resamples, at least 2 "
            "(default 1000).")
        .def_readwrite("ci_grid_step", &ScreenTRexControlParameter::ci_grid_step, "Step of the grid over the "
            "confidence level, in (0, 1] (default 0.001).")
        .def_readwrite("trex_method", &ScreenTRexControlParameter::trex_method, "Variant of the screening (default "
            "TREX). See ScreenTRexMethod.")
        .def_readwrite("rho_thr_DA", &ScreenTRexControlParameter::rho_thr_DA, "Correlation threshold of the "
            "dependency-aware variants (default 0.02). TREX_DA_AR1: the window around a variable reaches "
            "ceil(log(rho_thr_DA) / log(abs(cor_coef))) columns to each side, up to the first lag whose correlation "
            "abs(cor_coef)**lag has fallen to the threshold, range [0, 1). TREX_DA_NN: the same window, of which "
            "the variables stay whose sample correlation with the variable reaches the threshold, range [0, 1). "
            "A threshold at or above the correlation of the data leaves no neighbor, the variant then selects as "
            "TREX and the filter protects no longer, so the threshold stays below the correlation that the filter "
            "is meant to guard against. At the default the variant selected as TREX_DA_AR1 in simulations. "
            "TREX_DA_EQUI and TREX_DA_BLOCK_EQUI: the filter is skipped when abs(cor_coef) does not exceed the "
            "threshold, range [0, 1].")
        .def_readwrite("cor_coef",   &ScreenTRexControlParameter::cor_coef,   "Correlation coefficient of the "
            "dependency-aware variants; setting the parameter to -2 estimates it from the data. The variants use the "
            "magnitude of a supplied coefficient. Any other value outside [-1, 1] gives a warning and is replaced "
            "by the estimate.")
        .def_readwrite("n_blocks",   &ScreenTRexControlParameter::n_blocks,   "Number of blocks of "
            "TREX_DA_BLOCK_EQUI, between 1 and the number of variables (default 5).")
        .def_readwrite("trex_ctrl",  &ScreenTRexControlParameter::trex_ctrl,  "Nested base T-Rex algorithmic control "
            "parameters (overridden by the separate trex_control argument to TRexScreeningSelector).");

    py::class_<ScreenTRexSelectionResult, TRexSelector::SelectionResult>(m, "ScreenTRexSelectionResult",
        "SelectionResult "
        ""
        "extended with screening-specific fields.")
        .def_readonly("estimated_FDR", &ScreenTRexSelectionResult::estimated_FDR, "Estimated false discovery rate, "
            "one over the number of selected variables. An empty selection reports 1. It rests on the expected "
            "number of false discoveries, about one, and is no bound on the FDP of a run (see the notes of "
            "TRexScreeningSelector).")
        .def_readonly("used_bootstrap", &ScreenTRexSelectionResult::used_bootstrap, "Whether the confidence-based "
            "screening (bootstrap) made the selection.")
        .def_readonly("confidence_level", &ScreenTRexSelectionResult::confidence_level, "Confidence level that the "
            "grid search of the bootstrap picked. 0 without the bootstrap.")
        .def_readonly("beta_mat", &ScreenTRexSelectionResult::beta_mat, "Coefficients of the variables, averaged "
            "over the K experiments (one per variable).")
        .def_readonly("dummy_betas", &ScreenTRexSelectionResult::dummy_betas, "Non-zero coefficients of the dummies, "
            "collected over the K experiments. They are the sample of the bootstrap.")
        .def_readonly("estimated_correlation", &ScreenTRexSelectionResult::estimated_correlation, "Correlation "
            "coefficient of the dependency-aware variants: the supplied value, or the estimate from the data. 0 "
            "for TREX.");

    py::class_<PyScreenTRexSelector, PyTRexSelector>(m, "TRexScreeningSelector",
        "T-Rex screening selector (T-Rex+Screening): fast FDR-controlled variable screening in one\n"
        "pass of the random experiments, until one dummy is included (T = 1), with as many dummies\n"
        "as variables (L = p).\n"
        "\n"
        "Notes\n"
        "-----\n"
        "Attributes of the base class after a screening run. The screening has no target FDR and\n"
        "no calibration, so the attributes of TRexSelector hold fixed values: T_stop is 1, L is 1,\n"
        "num_dummies is p, v_thresh is 0.5 and tFDR is 0. r_mat is a 1 x 1 matrix with the number\n"
        "of selected variables, phi_prime equals the relative occurrences, and fdp_hat_mat is a\n"
        "1 x 1 matrix with 0, because no FDP is estimated per voting level. The FDR estimate of the\n"
        "run is estimated_FDR of the result.\n"
        "\n"
        "A response without signal. The screening stops at the first dummy, so about one null\n"
        "variable is expected to enter before it. On a response of pure noise it selected at least\n"
        "one variable in a quarter to a half of the runs in simulations, with fewer than one\n"
        "selected variable on average. A small selection is therefore no evidence of a signal. The\n"
        "estimate, one over the number of selected variables, rests on that expected number of\n"
        "false discoveries and is no bound on the FDP of a run: with two or more variables selected\n"
        "on noise it is 0.5 or less, while every one of them is false. TRexBiobankScreeningSelector\n"
        "accepts a screening result only when the estimate lies in its window and runs the T-Rex\n"
        "selector otherwise, which selected nothing on noise in about 99 % of the runs.\n"
        "\n"
        "L-loop strategies. The screening takes the strategies STANDARD (default), SEEDED,\n"
        "PERMUTATION and PERMUTATION_SEEDED. SEEDED gives the results of STANDARD and keeps one\n"
        "dummy matrix in memory at a time. PERMUTATION_SEEDED gives the results of PERMUTATION\n"
        "without storing the base dummy matrix. HCONCAT and SKIPL are refused, because they\n"
        "describe the growth of the dummy matrix over the L-loop, and the screening is one pass\n"
        "at L = p.\n"
        "\n"
        "References\n"
        "----------\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"False Discovery Rate Control for Fast\n"
        "  Screening of Large-Scale Genomics Biobanks,\" in Proc. 22nd IEEE Statist. Signal\n"
        "  Process. Workshop (SSP), 2023, pp. 666-670.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n")
        .def(py::init<Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::VectorXd>, const ScreenTRexControlParameter&,
            const TRexControlParameter&, int, bool>(),
             py::arg("X"), py::arg("y"),
             py::arg("screen_control") = ScreenTRexControlParameter(),
             py::arg("trex_control") = TRexControlParameter(),
             py::arg("seed") = -1, py::arg("verbose") = true,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             "Construct the screening selector. X and y are accessed zero-copy via Eigen::Map. A negative\n"
             "seed (default: -1) gives a run that is not reproducible.")
        .def(py::init<const DesignArray&, Eigen::Ref<Eigen::VectorXd>, const ScreenTRexControlParameter&,
            const TRexControlParameter&, int, bool>(),
             py::arg("X").noconvert(), py::arg("y"),
             py::arg("screen_control") = ScreenTRexControlParameter(),
             py::arg("trex_control") = TRexControlParameter(),
             py::arg("seed") = -1, py::arg("verbose") = true,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             kReaderSelectorDoc.c_str())
        .def("getScreenResult", &PyScreenTRexSelector::getScreenResult, "Return the ScreenTRexSelectionResult after "
            "calling select().");
}

#endif // TREX_PYTHON_TREX_SCREEN_BINDINGS_HPP

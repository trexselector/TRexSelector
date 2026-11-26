// =====================================================================================
// trex_spca_bindings.hpp - Pybind11 bindings for TRexSPCA
// =====================================================================================
#ifndef TREX_PYTHON_TREX_SPCA_BINDINGS_HPP
#define TREX_PYTHON_TREX_SPCA_BINDINGS_HPP
// =====================================================================================
/**
 * @file trex_spca_bindings.hpp
 *
 * @brief Pybind11 bindings for the T-Rex sparse PCA selector (T-Rex+SPCA).
 *
 * @details Exposes TRexSPCAResult, TRexSPCAControlParameter, SPCAMode, and the
 *          PyTRexSPCASelector instance wrapper to Python.
 */
// =====================================================================================

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

// trex_spca include
#include <trex_selector_methods/ad/trex_spca/trex_spca.hpp>
#include "deprecation.hpp"
#include "design_views.hpp"

// =====================================================================================

namespace py = pybind11;

using namespace trex::trex_selector_methods::ad::trex_spca;
namespace tc = trex::trex_selector_methods::ad::trex_core;
namespace tg = trex::trex_selector_methods::ad::trex_gvs;

// =====================================================================================

/**
 * @brief Zero-copy Python wrapper for TRexSPCA.
 *
 * @details Holds an Eigen::Map pointing into the numpy array's buffer (no copy),
 *          and a TRexSPCA instance. Mirrors the pattern of PyTRexGVSSelector.
 */
class PyTRexSPCASelector {
public:
    /**
     * @brief Construct a PyTRexSPCASelector.
     *
     * @param X        Design matrix (n × p) — accessed zero-copy via Eigen::Map.
     * @param M        Number of sparse PCs to extract.
     * @param tFDR     Target false discovery rate in (0, 1).
     * @param ctrl     Algorithmic control parameters.
     * @param seed     Random seed (-1 for non-deterministic).
     * @param verbose  Enable verbose output.
     */
    PyTRexSPCASelector(
        Eigen::Ref<Eigen::MatrixXd> X,
        Eigen::Index M,
        double tFDR,
        const TRexSPCAControlParameter& ctrl,
        int seed,
        bool verbose
    ) {
        X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X.data(), X.rows(), X.cols());
        selector_ = std::make_unique<TRexSPCA>(*X_map_, M, tFDR, ctrl, seed, verbose);
    }

    /**
     * @brief Constructor on a design that is not writable: always refused,
     *        because the selector centers the design in every run.
     *
     * @param X Design array, taken as it is (n × p).
     */
    PyTRexSPCASelector(
        const DesignArray& X,
        Eigen::Index M,
        double tFDR,
        const TRexSPCAControlParameter& ctrl,
        int seed,
        bool verbose
    ) : PyTRexSPCASelector(
            Eigen::Ref<Eigen::MatrixXd>(reader_design_map(
                X, true, "X",
                "TRexSPCASelector centers the design in place and restores it on return, which "
                "needs a writable array.")),
            M, tFDR, ctrl, seed, verbose) {}

    /** @brief Run the T-Rex sparse PCA selection. */
    void select() {
        result_ = selector_->select();
        has_result_ = true;
    }

    /** @brief Return the result from the last select() call. */
    const TRexSPCAResult& getResult() const {
        if (!has_result_)
            throw std::runtime_error("Call select() before accessing results.");
        return result_;
    }

private:
    std::unique_ptr<Eigen::Map<Eigen::MatrixXd>> X_map_;
    std::unique_ptr<TRexSPCA>                    selector_;
    TRexSPCAResult                               result_;
    bool                                         has_result_{false};
};

// =====================================================================================

/**
 * @brief Bind TRexSPCA types and PyTRexSPCASelector to a Python module.
 *
 * @param m The Python module (trex_selector_methods submodule).
 */
inline void bind_trex_spca(py::module& m) {

    // =========================================================================
    // SPCAMode enum
    // =========================================================================

    py::enum_<SPCAMode>(m, "SPCAMode",
        "Assembly of the sparse loading vector from the selected set.")
        .value("ActiveSet",   SPCAMode::ActiveSet,
               "Ridge regression of the component on the selected variables, normalized to unit "
               "length. The T-Rex PCA of the paper.")
        .value("Thresholded", SPCAMode::Thresholded,
               "The largest ordinary PCA loadings in absolute value, as many as variables were "
               "selected, normalized to unit length. The T-Rex Thresholded PCA of the paper.")
        .export_values();


    // =========================================================================
    // TRexSPCAControlParameter
    // =========================================================================

    py::class_<TRexSPCAControlParameter>(m, "TRexSPCAControlParameter",
        "Control parameters of the T-Rex sparse PCA selector.")
        .def(py::init<>())
        .def_readwrite("mode",                  &TRexSPCAControlParameter::mode,
                       "Assembly of the sparse loading vectors (default: ActiveSet).")
        .def_readwrite("lambda2_ridge_loadings", &TRexSPCAControlParameter::lambda2_ridge_loadings,
                       "Ridge penalty of the loading assembly in ActiveSet mode (default: 1e-6). "
                       "Not used in Thresholded mode.")
        .def_readwrite("gvs_ctrl",              &TRexSPCAControlParameter::gvs_ctrl,
                       "Control of the T-Rex+GVS run of each component (TRexGVSControlParameter), with "
                       "the base control in gvs_ctrl.trex_ctrl.\n"
                       "\n"
                       "gvs_ctrl.gvs_type must be GVSType.EN (the default). The selector raises "
                       "ValueError for IEN, because the method is defined on the elastic net.\n"
                       "\n"
                       "gvs_ctrl.en_solver and gvs_ctrl.trex_ctrl.solver_type have no effect: the "
                       "selector replaces them by en_solver of this control.\n"
                       "\n"
                       "gvs_ctrl.lambda_2: a negative value (default: -1) lets the rule "
                       "gvs_ctrl.lambda2_method choose it, 0 is the T-LASSO without the ridge part, and "
                       "a positive value is used as it is. A supplied value refers to the scale of the "
                       "standardized columns. Under ScalingMode.ZSCORE the columns have the squared "
                       "norm n - 1, so the same shrinkage needs the value times n - 1.\n"
                       "\n"
                       "gvs_ctrl.lambda2_method: CV_1SE_CCD (default), CV_MIN_CCD, CV_1SE_SVD or "
                       "CV_MIN_SVD. The selector raises ValueError for the rules of the informed "
                       "elastic net (CV_*_IEN_CCD and CV_*_TIK_SVD).")
        .def_readwrite("en_solver",             &TRexSPCAControlParameter::en_solver,
                       "Elastic net solver of the T-Rex+GVS runs (default: TENET_AUG). TENET_AUG is the "
                       "elastic net as a LASSO on row-augmented data, the form of the paper. TENET is "
                       "the elastic net on the Gram matrix. The two solve the same elastic net problem.");


    // =========================================================================
    // TRexSPCAResult
    // =========================================================================

    py::class_<TRexSPCAResult>(m, "TRexSPCAResult",
        "Result of the T-Rex sparse PCA selector. A component without a selected variable has a "
        "zero column in V and Z and the explained variance 0.")
        .def(py::init<>())
        .def_readwrite("Z",            &TRexSPCAResult::Z,
                       "Scores of the sparse components (n x M): column m is the centered X times "
                       "column m of V.")
        .def_readwrite("V",            &TRexSPCAResult::V,
                       "Sparse loadings (p x M), each column of unit L2 norm. In ActiveSet mode the "
                       "non-zero entries of column m are the variables of active_sets[m]. In "
                       "Thresholded mode they are the largest ordinary PCA loadings, as many as "
                       "active_sets[m] has entries, which can be other variables.")
        .def_readwrite("active_sets",  &TRexSPCAResult::active_sets,
                       "Selected variables per component (list of M index lists, 0-based), the result "
                       "of the T-Rex+GVS runs. The FDR control holds for these sets. In Thresholded "
                       "mode they can differ from the non-zero entries of V.")
        .def_readwrite("adjusted_ev",  &TRexSPCAResult::adjusted_ev,
                       "Adjusted explained variance per component (length M): r_mm^2 / (n - 1), with "
                       "r_mm the m-th diagonal element of R in the QR decomposition Z = QR. The "
                       "adjustment accounts for sparse components that are not orthogonal.")
        .def_readwrite("cumulative_ev",&TRexSPCAResult::cumulative_ev,
                       "Cumulative adjusted explained variance as a fraction of the total variance of "
                       "the centered X (length M), between 0 and 1. It is not the PEV of Definition 1 "
                       "of the paper, whose denominator (signal EV plus mixed EV) needs the true "
                       "support.");


    // =========================================================================
    // PyTRexSPCASelector (instance wrapper)
    // =========================================================================

    py::class_<PyTRexSPCASelector>(m, "TRexSPCASelector",
        "T-Rex sparse PCA selector (T-Rex+SPCA): sparse principal components whose non-zero\n"
        "loadings are selected with FDR control. Each ordinary principal component is the response\n"
        "of a T-Rex+GVS run with the elastic net.\n"
        "\n"
        "Notes\n"
        "-----\n"
        "Steps. (1) X is centered and its first M ordinary principal components are computed (a\n"
        "PCA of the covariance, the columns are not scaled). (2) For each component, a T-Rex+GVS\n"
        "run with the elastic net and the component as response selects its variables at the\n"
        "target FDR. (3) The sparse loading vector is assembled from the selected set (see\n"
        "SPCAMode), and the scores of the sparse component are the centered X times the loading\n"
        "vector. (4) The explained variance is adjusted for components that are not orthogonal,\n"
        "by the QR decomposition of the scores.\n"
        "\n"
        "A component without a selected variable has a zero column in V and Z and the explained\n"
        "variance 0.\n"
        "\n"
        "L-loop strategies. STANDARD, HCONCAT, SKIPL and SEEDED of the T-Rex+GVS runs are\n"
        "supported. PERMUTATION and PERMUTATION_SEEDED are refused at construction, as by the\n"
        "T-Rex+GVS selector.\n"
        "\n"
        "References\n"
        "----------\n"
        "J. Machkour, A. Breloy, M. Muma, D. P. Palomar, and F. Pascal, \"Sparse PCA with False\n"
        "  Discovery Rate Controlled Variable Selection,\" in Proc. IEEE 49th Int. Conf. Acoust.,\n"
        "  Speech, Signal Process. (ICASSP), 2024, pp. 9716-9720.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n")
        .def(py::init<Eigen::Ref<Eigen::MatrixXd>, Eigen::Index, double,
                      const TRexSPCAControlParameter&, int, bool>(),
             py::arg("X"),
             py::arg("M"),
             py::arg("tFDR"),
             py::arg("ctrl")    = TRexSPCAControlParameter(),
             py::arg("seed")    = -1,
             py::arg("verbose") = false,
             py::keep_alive<1, 2>(),
             "Construct the selector. X (n x p) is used without a copy and is centered in place "
             "during select() and restored on return. M is the number of sparse principal "
             "components, 1 to min(n, p). tFDR is the target FDR level of each component, in "
             "[0, 1]. A negative seed (default: -1) gives a run that is not reproducible, and "
             "each component runs with its own seed derived from it. verbose prints one line per "
             "component (default: False), and the T-Rex+GVS runs stay quiet.")
        .def(py::init<const DesignArray&, Eigen::Index, double,
                      const TRexSPCAControlParameter&, int, bool>(),
             py::arg("X").noconvert(),
             py::arg("M"),
             py::arg("tFDR"),
             py::arg("ctrl")    = TRexSPCAControlParameter(),
             py::arg("seed")    = -1,
             py::arg("verbose") = false,
             py::keep_alive<1, 2>(),
             "A design that is not writable is refused with ValueError: the selector centers the "
             "design in place and restores it on return, whatever the scaling mode.")
        .def("select",    &PyTRexSPCASelector::select,
             "Run the T-Rex sparse PCA selection.")
        .def("getResult", &PyTRexSPCASelector::getResult,
             py::return_value_policy::reference_internal,
             "The TRexSPCAResult of the last select().");
}

// =====================================================================================
#endif /* End of TREX_PYTHON_TREX_SPCA_BINDINGS_HPP */

// =====================================================================================
// trex_screening_biobanks_bindings.hpp - Pybind11 bindings for BiobankScreenTRex
// =====================================================================================
#ifndef TREX_PYTHON_TREX_SCREENING_BIOBANKS_BINDINGS_HPP
#define TREX_PYTHON_TREX_SCREENING_BIOBANKS_BINDINGS_HPP
// =====================================================================================
/**
 * @file trex_screening_biobanks_bindings.hpp
 * @brief Pybind11 bindings for BiobankScreenTRex.
 *
 */
// =====================================================================================

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

// trex includes
#include <trex_selector_methods/ad/trex_screening/trex_biobank_screening.hpp>
#include "trex_core_bindings.hpp"

// =====================================================================================

namespace py = pybind11;
using namespace trex::trex_selector_methods::ad::trex_screening;
using namespace trex::trex_selector_methods::ad::trex_biobank_screening;

// =====================================================================================

/**
 * @brief Response matrix Y from an array that is taken as it is.
 *
 * The selector only reads Y, so an array that is not writable serves. An
 * array that is not two-dimensional is left to the next overload.
 *
 * @param Y Response array, typically not writable.
 *
 * @return Reference to the array's buffer.
 */
inline Eigen::Ref<Eigen::MatrixXd> response_ref(const DesignArray& Y) {
    if (Y.ndim() != 2) { throw py::reference_cast_error(); }
    return reader_design_map(Y, false, "Y", "");
}

// =====================================================================================

/**
 * @brief Python wrapper class for BiobankScreenTRex leveraging zero-copy references.
 */
class PyBiobankScreenTRex {
private:
    /** @brief Map for the design matrix (X) */
    std::unique_ptr<Eigen::Map<Eigen::MatrixXd>> X_map_;

    /** @brief Map for the response vector (y) when single phenotype mapping is used. */
    std::unique_ptr<Eigen::Map<Eigen::VectorXd>> y_map_1d_;

    /** @brief Map for the response matrix (Y) when multiple phenotype mapping is used. */
    std::unique_ptr<Eigen::Map<Eigen::MatrixXd>> Y_map_2d_;

    /** @brief Pointer to the underlying BiobankScreenTRex selector instance. */
    std::unique_ptr<BiobankScreenTRex> selector_;

public:

    /**
     * @brief Zero-copy constructor for single phenotype mapping.
     *
     * @param X Reference to the design matrix (n × p).
     * @param y Reference to the response vector (n × 1).
     * @param bio_ctrl Configuration parameters specific to the BiobankScreenTRex method.
     * @param seed Random seed for reproducible computations.
     * @param verbose If true, prompts standard out messaging during select().
     */
    PyBiobankScreenTRex(
        Eigen::Ref<Eigen::MatrixXd> X,
        Eigen::Ref<Eigen::VectorXd> y,
        const BiobankScreenTRexControl& bio_ctrl,
        int seed,
        bool verbose
    ) {
        warn_if_cd_solver_type(bio_ctrl.trex_screen_ctrl.trex_ctrl.solver_type);

        X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X.data(),
                                                                  X.rows(),
                                                                  X.cols());

        y_map_1d_ = std::make_unique<Eigen::Map<Eigen::VectorXd>>(y.data(),
                                                                     y.size());

        selector_ = std::make_unique<BiobankScreenTRex>(
            *X_map_, *y_map_1d_, bio_ctrl, seed, verbose
        );
    }


    /**
     * @brief Zero-copy constructor for multiple phenotype mapping.
     */
    PyBiobankScreenTRex(
        Eigen::Ref<Eigen::MatrixXd> X,
        Eigen::Ref<Eigen::MatrixXd> Y,
        const BiobankScreenTRexControl& bio_ctrl,
        int seed,
        bool verbose
    ) {
        warn_if_cd_solver_type(bio_ctrl.trex_screen_ctrl.trex_ctrl.solver_type);

        X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X.data(),
                                                                  X.rows(),
                                                                  X.cols());

        Y_map_2d_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(Y.data(),
                                                                     Y.rows(),
                                                                     Y.cols());

        selector_ = std::make_unique<BiobankScreenTRex>(
            *X_map_, *Y_map_2d_, bio_ctrl, seed, verbose
        );
    }

    /**
     * @brief Zero-copy constructors on a design that is not writable.
     *
     * Accepted under ScalingMode::NONE, where the inner selectors write
     * nothing into the design, and refused under L2 and ZSCORE.
     *
     * @param X Design array, taken as it is (n × p).
     */
    PyBiobankScreenTRex(
        const DesignArray& X,
        Eigen::Ref<Eigen::VectorXd> y,
        const BiobankScreenTRexControl& bio_ctrl,
        int seed,
        bool verbose
    ) : PyBiobankScreenTRex(selector_design_ref(X, bio_ctrl.trex_screen_ctrl.trex_ctrl), y, bio_ctrl, seed,
                            verbose) {}

    PyBiobankScreenTRex(
        const DesignArray& X,
        Eigen::Ref<Eigen::MatrixXd> Y,
        const BiobankScreenTRexControl& bio_ctrl,
        int seed,
        bool verbose
    ) : PyBiobankScreenTRex(selector_design_ref(X, bio_ctrl.trex_screen_ctrl.trex_ctrl), Y, bio_ctrl, seed,
                            verbose) {}

    /**
     * @brief Zero-copy constructors on a response matrix that is not
     *        writable. The selector only reads Y.
     *
     * @param Y Response array, taken as it is (n × q).
     */
    PyBiobankScreenTRex(
        Eigen::Ref<Eigen::MatrixXd> X,
        const DesignArray& Y,
        const BiobankScreenTRexControl& bio_ctrl,
        int seed,
        bool verbose
    ) : PyBiobankScreenTRex(X, response_ref(Y), bio_ctrl, seed, verbose) {}

    PyBiobankScreenTRex(
        const DesignArray& X,
        const DesignArray& Y,
        const BiobankScreenTRexControl& bio_ctrl,
        int seed,
        bool verbose
    ) : PyBiobankScreenTRex(selector_design_ref(X, bio_ctrl.trex_screen_ctrl.trex_ctrl), response_ref(Y),
                            bio_ctrl, seed, verbose) {}


    /**
     * @brief Execute the Biobank Screening mapping for the loaded phenotype.
     */
    BiobankScreenTRexResult screenPhenotype() {
        return selector_->screenPhenotype();
    }


    /**
     * @brief Execute multi-phenotype bounds sequence mapping.
     */
    std::vector<BiobankScreenTRexResult> screenPhenotypes() {
        return selector_->screenPhenotypes();
    }
};


/**
 * @brief Binds the BiobankScreenTRex and related configurations.
 *
 * @param m The Python module to bind against.
 */
inline void bind_trex_screening_biobanks(py::module& m) {
    py::class_<BiobankScreenTRexControl>(m, "BiobankScreenTRexControl", "Control parameters for "
        "TRexBiobankScreeningSelector.")
        .def(py::init<>())
        .def_readwrite("target_FDR_trex", &BiobankScreenTRexControl::target_FDR_trex, "Target FDR level of the "
            "fallback T-Rex selector (default 0.1).")
        .def_readwrite("lower_bound_FDR", &BiobankScreenTRexControl::lower_bound_FDR, "Lower bound of the accepted "
            "range of the FDR estimate of the screening (default 0.05). An estimate below it suggests too little "
            "power.")
        .def_readwrite("upper_bound_FDR", &BiobankScreenTRexControl::upper_bound_FDR, "Upper bound of the accepted "
            "range of the FDR estimate of the screening (default 0.15). An estimate above it suggests too many "
            "false discoveries.")
        .def_readwrite("trex_screen_ctrl", &BiobankScreenTRexControl::trex_screen_ctrl,
                       "Nested T-Rex+Screening control shared by both the screening and "
                       "fallback-T-Rex paths. Set the screening method via "
                       "trex_screen_ctrl.trex_method and the base algorithm (K, "
                       "solver_type, ...) via trex_screen_ctrl.trex_ctrl. Its "
                       "lloop_strategy defaults to STANDARD; accepted values are "
                       "STANDARD, SEEDED, PERMUTATION, and PERMUTATION_SEEDED "
                       "(validated at construction). Its use_bootstrap_CI has no effect "
                       "here: the biobank screening computes the ordinary and the "
                       "bootstrap variant for every phenotype.");

    py::class_<BiobankScreenTRexResult>(m, "BiobankScreenTRexResult", "Result for a single phenotype returned by "
        "screenPhenotype() or screenPhenotypes().")
        .def_readonly("phenotype_index", &BiobankScreenTRexResult::phenotype_index, "0-based column index of the "
            "screened phenotype.")
        .def_readonly("selected_indices", &BiobankScreenTRexResult::selected_indices, "0-based indices of selected "
            "predictors.")
        .def_readonly("estimated_FDR", &BiobankScreenTRexResult::estimated_FDR, "FDR level of the final selection: "
            "the estimate of the accepted screening variant, or, for a fallback result, the target level "
            "target_FDR_trex, at which the T-Rex selector controls the FDR. An estimate can equal the target, so "
            "used_fallback_trex tells the two cases apart.")
        .def_readonly("method_used", &BiobankScreenTRexResult::method_used, "Method that gave the final selection: "
            "'Screen-TRex (ordinary)', 'Screen-TRex (bootstrap-CI)' or 'T-Rex (fallback)'.")
        .def_readonly("estimated_FDR_screen_ordinary", &BiobankScreenTRexResult::estimated_FDR_screen_ordinary,
            "Estimated FDR of the ordinary screening, stored for every phenotype, a fallback result included.")
        .def_readonly("estimated_FDR_screen_bootstrap", &BiobankScreenTRexResult::estimated_FDR_screen_bootstrap,
            "Estimated FDR of the bootstrap screening, stored for every phenotype, a fallback result included.")
        .def_readonly("selected_indices_screen_ordinary", &BiobankScreenTRexResult::selected_indices_screen_ordinary,
            "Selected indices of the ordinary screening (0-based), stored for every phenotype, also when another "
            "method gave the final selection.")
        .def_readonly("selected_indices_screen_bootstrap", &BiobankScreenTRexResult::selected_indices_screen_bootstrap,
            "Selected indices of the bootstrap screening (0-based), stored for every phenotype, also when another "
            "method gave the final selection.")
        .def_readonly("used_fallback_trex", &BiobankScreenTRexResult::used_fallback_trex, "True if the fallback T-Rex "
            "selector gave the final selection. estimated_FDR then holds the target level and not an estimate.");

    py::class_<PyBiobankScreenTRex>(m, "TRexBiobankScreeningSelector",
        "Biobank screening with the T-Rex selector: per phenotype the fast screening\n"
        "(T-Rex+Screening), and the T-Rex selector as fallback when no FDR estimate of the\n"
        "screening lies in the accepted range.\n"
        "\n"
        "Notes\n"
        "-----\n"
        "Decision rule per phenotype. The screening gives two FDR estimates, the one of its\n"
        "ordinary variant and the one of its bootstrap variant. A variant is accepted when its\n"
        "estimate lies in [lower_bound_FDR, upper_bound_FDR]. Of two accepted variants the one\n"
        "with the larger estimate is taken, and the bootstrap variant on a tie. With no accepted\n"
        "variant the T-Rex selector runs at target_FDR_trex.\n"
        "\n"
        "Both estimates come from one set of K experiments. The screening runs once, and the\n"
        "bootstrap variant is derived from its coefficients.\n"
        "\n"
        "L-loop strategies. STANDARD (default), SEEDED, PERMUTATION and PERMUTATION_SEEDED, as for\n"
        "TRexScreeningSelector. HCONCAT and SKIPL are refused.\n"
        "\n"
        "References\n"
        "----------\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"False Discovery Rate Control for Fast\n"
        "  Screening of Large-Scale Genomics Biobanks,\" in Proc. 22nd IEEE Statist. Signal\n"
        "  Process. Workshop (SSP), 2023, pp. 666-670.\n"
        "F. Scheidt, J. Machkour, and M. Muma, \"Solving FDR-Controlled Sparse Regression Problems\n"
        "  with Five Million Variables on a Laptop,\" in Proc. IEEE 9th Int. Workshop Comput. Adv.\n"
        "  Multi-Sensor Adapt. Process. (CAMSAP), 2023, pp. 116-120.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n")
        .def(py::init<Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::VectorXd>, const BiobankScreenTRexControl&, int,
            bool>(),
             py::arg("X"), py::arg("y"),
             py::arg("bio_ctrl") = BiobankScreenTRexControl(),
             py::arg("seed") = -1, py::arg("verbose") = false,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             "Construct for a single phenotype (1-D y). X and y are accessed zero-copy. A negative seed "
             "(default: -1) gives a run that is not reproducible.")
        .def(py::init<Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXd>, const BiobankScreenTRexControl&, int,
            bool>(),
             py::arg("X"), py::arg("Y"),
             py::arg("bio_ctrl") = BiobankScreenTRexControl(),
             py::arg("seed") = -1, py::arg("verbose") = false,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             "Construct for multiple phenotypes (2-D Y). X and Y are accessed zero-copy. Phenotype i (0-based) "
             "runs with seed + i, so the phenotypes draw independent dummies and a single phenotype is reproduced "
             "by a single-phenotype run with that seed.")
        .def(py::init<const DesignArray&, Eigen::Ref<Eigen::VectorXd>, const BiobankScreenTRexControl&, int,
            bool>(),
             py::arg("X").noconvert(), py::arg("y"),
             py::arg("bio_ctrl") = BiobankScreenTRexControl(),
             py::arg("seed") = -1, py::arg("verbose") = false,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             ("Single phenotype (1-D y). " + kReaderSelectorDoc).c_str())
        .def(py::init<const DesignArray&, Eigen::Ref<Eigen::MatrixXd>, const BiobankScreenTRexControl&, int,
            bool>(),
             py::arg("X").noconvert(), py::arg("Y"),
             py::arg("bio_ctrl") = BiobankScreenTRexControl(),
             py::arg("seed") = -1, py::arg("verbose") = false,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             ("Multiple phenotypes (2-D Y). " + kReaderSelectorDoc).c_str())
        .def(py::init<Eigen::Ref<Eigen::MatrixXd>, const DesignArray&, const BiobankScreenTRexControl&, int,
            bool>(),
             py::arg("X"), py::arg("Y").noconvert(),
             py::arg("bio_ctrl") = BiobankScreenTRexControl(),
             py::arg("seed") = -1, py::arg("verbose") = false,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             "Multiple phenotypes (2-D Y) with a response matrix that is not writable, for example the "
             "to_numpy() view of a MemoryMappedMatrix opened with AccessMode.ReadOnly. Y must be a "
             "Fortran-contiguous float64 array and is accessed zero-copy. The selector only reads Y.")
        .def(py::init<const DesignArray&, const DesignArray&, const BiobankScreenTRexControl&, int, bool>(),
             py::arg("X").noconvert(), py::arg("Y").noconvert(),
             py::arg("bio_ctrl") = BiobankScreenTRexControl(),
             py::arg("seed") = -1, py::arg("verbose") = false,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             ("Multiple phenotypes (2-D Y) with X and Y that are not writable. The selector only reads Y. " +
              kReaderSelectorDoc).c_str())
        .def("screenPhenotype", &PyBiobankScreenTRex::screenPhenotype, py::call_guard<py::gil_scoped_release>(), "Run "
            "screening for the single loaded phenotype.")
        .def("screenPhenotypes", &PyBiobankScreenTRex::screenPhenotypes, py::call_guard<py::gil_scoped_release>(),
            "Run "
            ""
            "screening for all loaded phenotypes.");

// =====================================================================================
}
// =====================================================================================
#endif /* End of TREX_PYTHON_TREX_SCREENING_BIOBANKS_BINDINGS_HPP */

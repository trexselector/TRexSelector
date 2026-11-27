// =====================================================================================
// trex_gvs_bindings.hpp - Pybind11 bindings for TRexGVSSelector
// =====================================================================================
#ifndef TREX_PYTHON_TREX_GVS_BINDINGS_HPP
#define TREX_PYTHON_TREX_GVS_BINDINGS_HPP
// =====================================================================================
/**
 * @file trex_gvs_bindings.hpp
 * @brief Pybind11 bindings for TRexGVSSelector.
 *
 */
// =====================================================================================

// std includes
#include <optional>

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

// trex includes
#include <ml_methods/model_selection/lambda2_cv_diagnostics.hpp>
#include <trex_selector_methods/ad/trex_gvs/trex_gvs.hpp>
#include "trex_core_bindings.hpp"

// =====================================================================================

namespace py = pybind11;
using namespace trex::trex_selector_methods::ad::trex_gvs;
using namespace trex::trex_selector_methods::ad::trex_core;
namespace gvs_sd = trex::trex_selector_methods::ad::utils::solver_dispatch;

// =====================================================================================

/**
 * @brief The control of a GVS run. The C++ constructor takes one control that
 *        nests the base control. Python can pass a base control as a second
 *        argument: it then replaces the nested one, as it is. Nothing is
 *        changed in either, and the C++ constructor derives or refuses the
 *        solver.
 */
inline TRexGVSControl merge_gvs_control(
    const TRexGVSControl& gvs_control,
    const std::optional<TRexControl>& trex_control)
{
    TRexGVSControl gvs_ctrl = gvs_control;
    if (trex_control) { gvs_ctrl.trex_ctrl = *trex_control; }
    return gvs_ctrl;
}


/**
 * @brief Python wrapper class for TRexGVSSelector leveraging zero-copy references.
 */
class PyTRexGVSSelector : public PyTRexSelector {
public:
    /**
     * @brief Zero-copy constructor for the PyTRexGVSSelector wrapper.
     *
     * @param X Reference to the design matrix (n × p).
     * @param y Reference to the response vector (n × 1).
     * @param tFDR Target False Discovery Rate threshold.
     * @param gvs_control Configuration parameters specific to the Group Variable Selection method.
     * @param trex_control Base control that replaces the nested one of gvs_control, or none.
     * @param seed Random seed for reproducible computations.
     * @param verbose If true, prompts standard out messaging during select().
     */
    PyTRexGVSSelector(
        Eigen::Ref<Eigen::MatrixXd> X,
        Eigen::Ref<Eigen::VectorXd> y,
        double tFDR,
        const TRexGVSControl& gvs_control,
        const std::optional<TRexControl>& trex_control,
        int seed,
        bool verbose
    ) {
        this->X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X.data(),
                                                                  X.rows(),
                                                                  X.cols());
        this->y_map_ = std::make_unique<Eigen::Map<Eigen::VectorXd>>(y.data(),
                                                                  y.size());

        const TRexGVSControl gvs_ctrl = merge_gvs_control(gvs_control, trex_control);

        this->selector_ = std::make_unique<TRexGVSSelector>(
            *(this->X_map_), *(this->y_map_), tFDR, gvs_ctrl, seed, verbose
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
    PyTRexGVSSelector(
        const DesignArray& X,
        Eigen::Ref<Eigen::VectorXd> y,
        double tFDR,
        const TRexGVSControl& gvs_control,
        const std::optional<TRexControl>& trex_control,
        int seed,
        bool verbose
    ) : PyTRexGVSSelector(selector_design_ref(X, trex_control ? *trex_control : gvs_control.trex_ctrl), y, tFDR,
                          gvs_control, trex_control, seed, verbose) {}

    /**
     * @brief Get the complete GVS selection result.
     */
    const TRexGVSResult& getGVSResult() const {
        return static_cast<TRexGVSSelector*>(this->selector_.get())->getGVSResult();
    }

    /**
     * @brief The lambda_2 tuner health record of the last select().
     */
    const trex::ml_methods::model_selection::Lambda2CVDiagnostics&
    getLambda2CVDiagnostics() const {
        return static_cast<TRexGVSSelector*>(this->selector_.get())
            ->getLambda2CVDiagnostics();
    }
};


/**
 * @brief Binds the TRexGVSSelector and related configurations.
 *
 * @param m The Python module to bind against.
 */
inline void bind_trex_gvs(py::module& m) {
    py::enum_<GVSType>(m, "GVSType", "Penalty of the grouped selection of T-Rex+GVS.")
        .value("EN", GVSType::EN, "Elastic net. en_solver chooses the solver.")
        .value("IEN", GVSType::IEN, "Informed elastic net, with the groups in the penalty (solver via "
            "trex_control.solver_type: TIENET default, TIENET_AUG).")
        .export_values();

    py::enum_<ENSolverType>(m, "ENSolverType", "Elastic-Net solver variant used when gvs_type == EN.")
        .value("TENET",     ENSolverType::TENET,     "Gram-based Elastic Net (ridge absorbed via Cholesky update).")
        .value("TENET_AUG", ENSolverType::TENET_AUG, "Augmented-LASSO Elastic Net (row-augmented system).")
        .export_values();

    py::class_<trex::ml_methods::model_selection::Lambda2CVDiagnostics>(m, "Lambda2CVDiagnostics",
        "Health record of a lambda_2 tuner run: whether it ran, its cost, the truncation-edge count, and the "
            "boundary-pick guard fields for the returned pick (grid_index / grid_size, at_grid_edge, profile_at_edge, "
            "boundary_pick).")
        .def(py::init<>())
        .def_readonly("ran", &trex::ml_methods::model_selection::Lambda2CVDiagnostics::ran, "A tuner ran (False: "
            "lambda_2 came from the caller).")
        .def_readonly("seconds", &trex::ml_methods::model_selection::Lambda2CVDiagnostics::seconds, "Wall-clock "
            "seconds spent resolving lambda_2.")
        .def_readonly("n_profile_at_edge",
            &trex::ml_methods::model_selection::Lambda2CVDiagnostics::n_profile_at_edge, "lambda_2 "
            ""
            "values whose profiled lambda_1 optimum landed on a truncation edge (0 = healthy).")
        .def_readonly("support_cap", &trex::ml_methods::model_selection::Lambda2CVDiagnostics::support_cap, "Support "
            "cap (glmnet dfmax) applied, or -1 when disabled.")
        .def_readonly("dfmax_ratio", &trex::ml_methods::model_selection::Lambda2CVDiagnostics::dfmax_ratio, "IEN "
            "tuner: the support-cap ratio it was called with (cv_dfmax_ratio, as passed; <= 0 disabled the cap). 0 "
                "when "
            "no IEN-geometry tuner ran.")
        .def_readonly("low_dimensional", &trex::ml_methods::model_selection::Lambda2CVDiagnostics::low_dimensional,
            "p "
            ""
            "<= cap: the cap was inert.")
        // Boundary-pick guard: the chosen pick on its grid.
        .def_readonly("grid_index", &trex::ml_methods::model_selection::Lambda2CVDiagnostics::grid_index, "0-based "
            "index of the returned lambda_2 on the tuner's grid (-1: no tuner ran).")
        .def_readonly("grid_size", &trex::ml_methods::model_selection::Lambda2CVDiagnostics::grid_size, "Number of "
            "lambda_2 grid points searched (0: no tuner ran).")
        .def_readonly("at_grid_edge", &trex::ml_methods::model_selection::Lambda2CVDiagnostics::at_grid_edge, "The "
            "returned lambda_2 is the first or last grid point: the grid boundary, not a CV optimum.")
        .def_readonly("profile_at_edge", &trex::ml_methods::model_selection::Lambda2CVDiagnostics::profile_at_edge,
            "IEN "
            ""
            "tuner: the profiled lambda_1 optimum of the chosen lambda_2 sat on the truncation edge of the lambda_1 "
            "path.")
        .def_readonly("profile_at_support_cap",
            &trex::ml_methods::model_selection::Lambda2CVDiagnostics::profile_at_support_cap, "IEN "
            ""
            "tuner: profile_at_edge, and that edge is the support cap (raise cv_dfmax_ratio to move it).")
        .def_readonly("boundary_pick", &trex::ml_methods::model_selection::Lambda2CVDiagnostics::boundary_pick,
            "at_grid_edge "
            ""
            "or profile_at_edge: the selector warned once and returned the pick unchanged.");

    py::enum_<GVSLambda2Method>(m, "LambdaSelectionMethod", "Rule for auto-selecting the ridge parameter lambda_2 "
        "when lambda_2 < 0.")
        .value("CV_1SE_SVD", GVSLambda2Method::CV_1SE_SVD, "JacobiSVD ridge CV, 1-SE rule.")
        .value("CV_MIN_SVD", GVSLambda2Method::CV_MIN_SVD, "JacobiSVD ridge CV, minimum-CV-error rule.")
        .value("CV_1SE_CCD", GVSLambda2Method::CV_1SE_CCD, "Coordinate-descent ridge CV, glmnet-faithful 1-SE "
            "rule "
            "(default).")
        .value("CV_MIN_CCD", GVSLambda2Method::CV_MIN_CCD, "Coordinate-descent ridge CV, glmnet-faithful minimum "
            "rule.")
        .value("CV_1SE_IEN_CCD", GVSLambda2Method::CV_1SE_IEN_CCD, "IEN-geometry coordinate-descent CV (lambda1 "
            "profiled, group-consuming), 1-SE rule. IEN default.")
        .value("CV_MIN_IEN_CCD", GVSLambda2Method::CV_MIN_IEN_CCD, "IEN-geometry coordinate-descent CV (lambda1 "
            "profiled, group-consuming), minimum rule.")
        .value("CV_1SE_TIK_SVD", GVSLambda2Method::CV_1SE_TIK_SVD, "Generalized-Tikhonov SVD ridge CV (IEN "
            "geometry), 1-SE rule.")
        .value("CV_MIN_TIK_SVD", GVSLambda2Method::CV_MIN_TIK_SVD, "Generalized-Tikhonov SVD ridge CV (IEN "
            "geometry), minimum rule.")
        .export_values();

    py::class_<TRexGVSControl>(m, "TRexGVSControlParameter", "Control parameters for TRexGVSSelector.")
        .def(py::init<>())
        .def_readwrite("gvs_type",    &TRexGVSControl::gvs_type,    "Penalty of the grouped selection: EN "
            "(elastic net, default) or IEN (informed elastic net).")
        .def_readwrite("en_solver",   &TRexGVSControl::en_solver,   "EN solver variant when gvs_type == EN: "
            "TENET (default) or TENET_AUG.")
        .def_readwrite("corr_max",    &TRexGVSControl::corr_max,    "Correlation threshold of the automatic "
            "clustering (default 0.5). The dendrogram is cut at the height 1 - corr_max. With Single linkage "
            "(default) no two variables of different clusters correlate above corr_max. With Complete linkage every "
            "two variables of one cluster correlate at least at corr_max, and variables of different clusters can "
            "exceed it.")
        .def_readwrite("hc_linkage",  &TRexGVSControl::hc_linkage,
            "Linkage of the hierarchical clustering: Single (default), Complete, Average or WPGMA. Used when no\n"
            "groups are supplied. Ward is rejected, because the correlation distance 1 - |cor| is not Euclidean.")
        .def_readwrite("lambda_2", &TRexGVSControl::lambda_2, "Ridge L2 penalty (< 0 = auto-select, 0 = "
            "pure TLASSO, > 0 = fixed). A supplied value refers to the scale of the standardized columns. Under "
            "ScalingMode.ZSCORE the columns have the squared norm n - 1, so the same shrinkage needs the value "
            "times n - 1.")
        .def_readwrite("lambda2_method", &TRexGVSControl::lambda2_method, "Rule for auto-selecting lambda_2 "
            "when lambda_2 < 0.")
        .def_readwrite("cv_n_folds", &TRexGVSControl::cv_n_folds, "Number of folds for cross-validated "
            "lambda_2 selection.")
        .def_readwrite("cv_n_lambda", &TRexGVSControl::cv_n_lambda, "Number of points on the lambda_2 grid "
            "of the rules CV_*_SVD, CV_*_CCD and CV_*_TIK_SVD (default 1000). The rules CV_*_IEN_CCD use "
            "cv_n_lambda2.")
        .def_readwrite("cv_n_lambda1", &TRexGVSControl::cv_n_lambda1, "IEN-geometry CD tuners (CV_*_IEN_CCD): "
            "points on the profiled lambda_1 path per lambda_2 value (default 40).")
        .def_readwrite("cv_dfmax_ratio", &TRexGVSControl::cv_dfmax_ratio, "IEN-geometry CD tuners "
            "(CV_*_IEN_CCD): support cap as a fraction of the smallest fold's training size (default 0.8; <= 0 "
                "disables "
            "the cap; larger values allow larger supports at a higher tuner cost; must be finite).")
        .def_readwrite("cv_n_lambda2", &TRexGVSControl::cv_n_lambda2, "IEN-geometry CD tuners: points on the "
            "lambda_2 grid (default 100).")
        .def_readwrite("cv_seed", &TRexGVSControl::cv_seed, "Seed for the CV fold-permutation RNG (-1 = "
            "derive "
            "from the T-Rex seed).")
        .def_readwrite("prior_groups", &TRexGVSControl::prior_groups, "Manually specified group assignments: "
            "one entry per variable, numbered 0 to M - 1 without gaps. Empty (default) clusters the variables.")
        .def_readwrite("group_labels", &TRexGVSControl::group_labels, "Optional human-readable cluster names "
            "(size must equal number of clusters).")
        .def_readwrite("tenet_aug_use_lars", &TRexGVSControl::tenet_aug_use_lars, "For the augmented "
            "solvers, en_solver = TENET_AUG or gvs_type = IEN with the solver TIENET_AUG: use a pure-LARS inner "
            "solver that never drops variables, as trex(type = \"lar\") of the 1.x versions of the R package does "
            "(default False). Ignored for TENET and TIENET.")
        .def_readwrite("trex_ctrl", &TRexGVSControl::trex_ctrl, "Nested base T-Rex control, whose solver "
            "is TENET. It applies when TRexGVSSelector gets no trex_control argument, and that argument replaces it "
            "otherwise. On the EN track its solver_type stays TENET or names the solver of en_solver. On the IEN "
            "track it has to name TIENET or TIENET_AUG. TRexSPCASelector uses it directly.");

    py::class_<TRexGVSResult, TRexResult>(m, "GVSSelectionResult",
        "SelectionResult "
        ""
        "extended with GVS-specific fields.")
        .def_readonly("lambda2_used", &TRexGVSResult::lambda2_used, "lambda_2 of the last "
            "select(): the supplied value, or the pick of the cross-validation rule lambda2_method.")
        .def_readonly("gvs_type", &TRexGVSResult::gvs_type, "Penalty of the grouped selection "
            "that was applied, EN or IEN.")
        .def_readonly("max_clusters", &TRexGVSResult::max_clusters, "Number of groups M, "
            "supplied or found by the clustering.")
        .def_readonly("hc_method_used", &TRexGVSResult::hc_method_used, "Origin of the groups: "
            "\"prior\" when they were supplied, else the linkage of the clustering, \"single\", \"complete\", "
            "\"average\" or \"wpgma\".")
        .def_readonly("groups_vec", &TRexGVSResult::groups_vec, "Group assignment vector "
            "(0-based "
            "group index per predictor).")
        .def_readonly("group_labels", &TRexGVSResult::group_labels, "Names of the groups "
            "(length max_clusters), or an empty list when none were supplied.")
        .def_readonly("lambda2_cv", &TRexGVSResult::lambda2_cv, "Health record of the lambda_2 "
            "tuner (Lambda2CVDiagnostics).");

    py::class_<PyTRexGVSSelector, PyTRexSelector>(m, "TRexGVSSelector",
        "T-Rex selector for grouped variable selection (T-Rex+GVS): FDR-controlled variable\n"
        "selection for groups of highly correlated variables.\n"
        "\n"
        "Notes\n"
        "-----\n"
        "The solver. On the EN track en_solver of the GVS control chooses it, TENET by default, and\n"
        "the solver_type of the base control stays TENET or names the same solver. On the IEN track\n"
        "the solver_type of the base control names it, TIENET or TIENET_AUG. Every other solver is\n"
        "refused at construction, also TLARS, the default of TRexControlParameter, and the default\n"
        "TENET of the nested base control on the IEN track.\n"
        "\n"
        "L-loop strategies. The selector supports STANDARD, HCONCAT, SKIPL and SEEDED, and SEEDED\n"
        "equals STANDARD here. PERMUTATION and PERMUTATION_SEEDED are refused, because a row\n"
        "permutation destroys the covariance of the dummies within a cluster.\n"
        "\n"
        "References\n"
        "----------\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"False Discovery Rate Control for Grouped\n"
        "  Variable Selection in High-Dimensional Linear Models Using the T-Knock Filter,\" in\n"
        "  Proc. 30th Eur. Signal Process. Conf. (EUSIPCO), 2022, pp. 892-896.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Informed Elastic Net for Fast Grouped\n"
        "  Variable Selection and FDR Control in Genomics Research,\" in Proc. IEEE 9th Int.\n"
        "  Workshop Comput. Adv. Multi-Sensor Adapt. Process. (CAMSAP), 2023, pp. 466-470.\n"
        "H. Zou and T. Hastie, \"Regularization and Variable Selection via the Elastic Net,\" J. R.\n"
        "  Stat. Soc. Ser. B Stat. Methodol., vol. 67, no. 2, pp. 301-320, 2005.\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n")
        .def(py::init<Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::VectorXd>, double, const TRexGVSControl&,
            const std::optional<TRexControl>&, int, bool>(),
             py::arg("X"), py::arg("y"), py::arg("tFDR") = 0.1,
             py::arg("gvs_control") = TRexGVSControl(),
             py::arg("trex_control") = py::none(),
             py::arg("seed") = -1, py::arg("verbose") = true,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             "Construct the GVS selector. X and y are accessed zero-copy via Eigen::Map. Without trex_control the\n"
             "base control nested in gvs_control applies. A trex_control that is passed replaces it as it is, and\n"
             "a solver that does not fit the track raises ValueError. A negative seed (default: -1) gives a run\n"
             "that is not reproducible.")
        .def(py::init<const DesignArray&, Eigen::Ref<Eigen::VectorXd>, double, const TRexGVSControl&,
            const std::optional<TRexControl>&, int, bool>(),
             py::arg("X").noconvert(), py::arg("y"), py::arg("tFDR") = 0.1,
             py::arg("gvs_control") = TRexGVSControl(),
             py::arg("trex_control") = py::none(),
             py::arg("seed") = -1, py::arg("verbose") = true,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             kReaderSelectorDoc.c_str())
        .def("getGVSResult", &PyTRexGVSSelector::getGVSResult, "Return the GVSSelectionResult after calling select().")
        .def("getLambda2CVDiagnostics", &PyTRexGVSSelector::getLambda2CVDiagnostics, "Return the lambda_2 tuner health "
            "record of the last select().");
}
// =====================================================================================
#endif /* End of TREX_PYTHON_TREX_GVS_BINDINGS_HPP */

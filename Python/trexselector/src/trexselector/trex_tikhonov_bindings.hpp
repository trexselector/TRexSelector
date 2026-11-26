// =====================================================================================
// trex_tikhonov_bindings.hpp - Pybind11 bindings for TRexTikhonovSelector
// =====================================================================================
#ifndef TREX_PYTHON_TREX_TIKHONOV_BINDINGS_HPP
#define TREX_PYTHON_TREX_TIKHONOV_BINDINGS_HPP
// =====================================================================================
/**
 * @file trex_tikhonov_bindings.hpp
 * @brief Pybind11 bindings for TRexTikhonovSelector (general-Tikhonov T-Rex).
 *
 * @details Mirrors the R surface (rcpp_trex_tikhonov.cpp): the control struct
 *          with its nested base control merged from a separate `trex_control`
 *          argument, the lambda_2 getters (used value, CV diagnostics, pilot
 *          record) and the two K builders. The penalty matrix K crosses the
 *          boundary as a dense p x p array and is converted with
 *          sparseView() (no SciPy dependency), the same convention as the
 *          Tikhonov T-solvers.
 */
// =====================================================================================

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

// std includes
#include <optional>

// trex includes
#include <ml_methods/model_selection/lambda2_cv_diagnostics.hpp>
#include <trex_selector_methods/ad/trex_tikhonov/trex_tikhonov.hpp>
#include <trex_selector_methods/ad/trex_tikhonov/tik_lambda2_pick.hpp>
#include <trex_selector_methods/ad/trex_tikhonov/tik_lambda2_power_rule.hpp>
#include "trex_core_bindings.hpp"
#include "trex_gvs_bindings.hpp"

// =====================================================================================

namespace py = pybind11;
namespace tik    = trex::trex_selector_methods::ad::trex_tikhonov;
namespace tik_sd = trex::trex_selector_methods::ad::utils::solver_dispatch;
using trex::trex_selector_methods::ad::trex_core::TRexControlParameter;

// =====================================================================================

/**
 * @brief The control of a run. The C++ constructor takes one control that
 *        nests the base control, with TIKNET as its solver. Python can pass
 *        a base control as a second argument: it then replaces the nested
 *        one, as it is. Without it the nested control applies. Nothing is
 *        changed in either, and the C++ constructor validates the result.
 */
inline tik::TRexTikhonovControlParameter merge_tikhonov_control(
    const tik::TRexTikhonovControlParameter& tik_control,
    const std::optional<TRexControlParameter>& trex_control)
{
    tik::TRexTikhonovControlParameter tik_ctrl = tik_control;
    if (trex_control) { tik_ctrl.trex_ctrl = *trex_control; }
    return tik_ctrl;
}


/**
 * @brief Python wrapper class for TRexTikhonovSelector leveraging zero-copy references.
 */
class PyTRexTikhonovSelector : public PyTRexSelector {
public:
    /**
     * @brief Zero-copy constructor for the PyTRexTikhonovSelector wrapper.
     *
     * @param X Reference to the design matrix (n × p).
     * @param y Reference to the response vector (n × 1).
     * @param tFDR Target False Discovery Rate threshold.
     * @param tik_control Tikhonov-specific control (K or penalty structure, lambda_2 rule, dummies).
     * @param trex_control Base control that replaces the nested one of tik_control, or none.
     * @param seed Random seed for reproducible computations.
     * @param verbose If true, prompts standard out messaging during select().
     */
    PyTRexTikhonovSelector(
        Eigen::Ref<Eigen::MatrixXd> X,
        Eigen::Ref<Eigen::VectorXd> y,
        double tFDR,
        const tik::TRexTikhonovControlParameter& tik_control,
        const std::optional<TRexControlParameter>& trex_control,
        int seed,
        bool verbose
    ) {
        this->X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X.data(),
                                                                  X.rows(),
                                                                  X.cols());
        this->y_map_ = std::make_unique<Eigen::Map<Eigen::VectorXd>>(y.data(),
                                                                  y.size());
        const tik::TRexTikhonovControlParameter tik_ctrl =
            merge_tikhonov_control(tik_control, trex_control);
        this->selector_ = std::make_unique<tik::TRexTikhonovSelector>(
            *(this->X_map_), *(this->y_map_), tFDR, tik_ctrl, seed, verbose
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
    PyTRexTikhonovSelector(
        const DesignArray& X,
        Eigen::Ref<Eigen::VectorXd> y,
        double tFDR,
        const tik::TRexTikhonovControlParameter& tik_control,
        const std::optional<TRexControlParameter>& trex_control,
        int seed,
        bool verbose
    ) : PyTRexTikhonovSelector(
            selector_design_ref(X, trex_control ? *trex_control : tik_control.trex_ctrl), y, tFDR,
            tik_control, trex_control, seed, verbose) {}

    /** @brief The lambda_2 the last run used, on the scale of a supplied value; -1 before a run. */
    double getLambda2Used() const {
        return static_cast<tik::TRexTikhonovSelector*>(this->selector_.get())
            ->getLambda2Used();
    }

    /** @brief The lambda_2 tuner health record of the last select(). */
    const trex::ml_methods::model_selection::lambda2_cv_diagnostics&
    getLambda2CVDiagnostics() const {
        return static_cast<tik::TRexTikhonovSelector*>(this->selector_.get())
            ->getLambda2CVDiagnostics();
    }

    /** @brief The lambda_2 pilot record of the last select() (the pilot rule). */
    const tik::Lambda2PilotResult& getLambda2PilotResult() const {
        return static_cast<tik::TRexTikhonovSelector*>(this->selector_.get())
            ->getLambda2PilotResult();
    }

    /** @brief Run only the lambda_2 pilot and give X back. */
    tik::Lambda2PilotResult runPilotOnly() {
        return static_cast<tik::TRexTikhonovSelector*>(this->selector_.get())
            ->runPilotOnly();
    }

    /** @brief The clusters the dummies of the last run were drawn for. */
    const tik::DummyClusters& getDummyClusters() const {
        return static_cast<tik::TRexTikhonovSelector*>(this->selector_.get())
            ->getDummyClusters();
    }
};


/**
 * @brief Binds the TRexTikhonovSelector and related configurations.
 *
 * @param m The Python module to bind against (the AD selector submodule;
 *          bind_trex_gvs must have registered Lambda2CVDiagnostics before).
 */
inline void bind_trex_tikhonov(py::module& m) {
    py::enum_<tik::TikLambda2Method>(m, "TikLambda2Method",
        "Rule that chooses lambda_2 when no value is supplied (tik_control.lambda_2 < 0). The pilot rule reads it "
            "off a few random experiments per value of pilot_grid. The four cross-validation rules tune on the "
            "normalized data before the random experiments. Every rule reports its pick on the working scale of "
            "the run. The knee of the full selection profile is tikhonov_lambda2_power_rule().")
        .value("CV_1SE_IEN_CCD", tik::TikLambda2Method::CV_1SE_IEN_CCD, "Cross-validation with the lambda_1 path "
            "profiled out per lambda_2 candidate, 1-SE rule. Well-posed for every p and n.")
        .value("CV_MIN_IEN_CCD", tik::TikLambda2Method::CV_MIN_IEN_CCD, "Cross-validation with the lambda_1 path "
            "profiled out, minimum rule.")
        .value("CV_1SE_TIK_SVD", tik::TikLambda2Method::CV_1SE_TIK_SVD, "Analytic ridge cross-validation "
            "(lambda_1 = 0), 1-SE rule. Exact for n > p, refused when the unpenalized null space interpolates the "
            "training folds.")
        .value("CV_MIN_TIK_SVD", tik::TikLambda2Method::CV_MIN_TIK_SVD, "Analytic ridge cross-validation "
            "(lambda_1 = 0), minimum rule.")
        .value("PILOT_KNEE",     tik::TikLambda2Method::PILOT_KNEE,     "The pilot rule (default): a few random "
            "experiments per value of pilot_grid, then the smallest lambda_2 whose count is within "
            "pilot_knee_tolerance of the largest.")
        .export_values();

    py::enum_<tik::PilotStatistic>(m, "PilotStatistic",
        "The count the pilot rule picks on, per value of the grid. The record holds all three.")
        .value("MAJORITY",       tik::PilotStatistic::MAJORITY,       "The number of variables ahead of the "
            "T_pilot-th dummy in more than half of the pilot experiments, the count of the T-Rex selector at the "
            "voting level 1/2. Zero under the global null.")
        .value("CALIBRATED",     tik::PilotStatistic::CALIBRATED,     "The T-Rex selection rule applied to the "
            "pilot data: the largest selection whose estimated FDP stays within tFDR. Coarser at a small pilot_K.")
        .value("GROUP_MAJORITY", tik::PilotStatistic::GROUP_MAJORITY, "The MAJORITY count restricted to the groups "
            "of prior_groups in which at least pilot_group_fraction of the members are ahead. A group that enters "
            "as a whole counts with its members, a group that enters as a minority counts nothing. It needs "
            "prior_groups.")
        .export_values();

    py::enum_<tik::TikPenaltyStructure>(m, "TikPenaltyStructure",
        "Structural form of the Tikhonov penalty.")
        .value("SPARSE_GENERAL",  tik::TikPenaltyStructure::SPARSE_GENERAL,  "tik_control.tikhonov_K is the p x p PSD "
            "penalty matrix (a scipy.sparse matrix or a dense array).")
        .value("GROUP_CENTERING", tik::TikPenaltyStructure::GROUP_CENTERING, "K = blockdiag(I - J/p_m) over "
            "prior_groups: the members of a group share a coefficient. The groups are required, and tikhonov_K "
            "may be left out.")
        .export_values();

    py::class_<tik::Lambda2PickInfo>(m, "Lambda2PickInfo",
        "How the knee was read off a selection-count profile (the pilot rule and the power rule).")
        .def(py::init<>())
        .def_readonly("index",        &tik::Lambda2PickInfo::index,        "0-based index of the pick on the grid.")
        .def_readonly("max_selected", &tik::Lambda2PickInfo::max_selected, "Largest selection count on the grid.")
        .def_readonly("n_ties",       &tik::Lambda2PickInfo::n_ties,       "Grid values within the tolerance of the "
            "maximum (at tolerance 0: the values attaining it).");

    py::class_<tik::Lambda2Probe>(m, "Lambda2Probe",
        "One grid value's selection in the full lambda_2 profile of tikhonov_lambda2_power_rule().")
        .def(py::init<>())
        .def_readonly("lambda2",  &tik::Lambda2Probe::lambda2,  "The lambda_2 this probe ran at: the grid value "
            "times grid_scale of the result.")
        .def_readonly("selected", &tik::Lambda2Probe::selected, "The selected variables (0-based, ascending).")
        .def_readonly("T_stop",   &tik::Lambda2Probe::T_stop,   "The calibrated stopping time of the run.")
        .def_readonly("seconds",  &tik::Lambda2Probe::seconds,  "Wall-clock seconds of the run.");

    py::class_<tik::Lambda2PowerRuleResult>(m, "Lambda2PowerRuleResult",
        "The full lambda_2 profile and its knee: every grid value's selection, the pick and the boundary flags.")
        .def(py::init<>())
        .def_readonly("lambda2",      &tik::Lambda2PowerRuleResult::lambda2,      "The chosen lambda_2, on the "
            "working scale of the runs. Supplied as lambda_2, it reproduces the selection at the pick.")
        .def_readonly("grid_scale",   &tik::Lambda2PowerRuleResult::grid_scale,   "The factor between the supplied "
            "grid and the values the probes ran at: n - 1 under ScalingMode.ZSCORE, 1 otherwise.")
        .def_readonly("index",        &tik::Lambda2PowerRuleResult::index,        "0-based index of the pick into "
            "probes.")
        .def_readonly("probes",       &tik::Lambda2PowerRuleResult::probes,       "Every grid value's selection, "
            "ascending in lambda_2.")
        .def_readonly("pick",         &tik::Lambda2PowerRuleResult::pick,         "Lambda2PickInfo of the pick.")
        .def_readonly("at_lower_end", &tik::Lambda2PowerRuleResult::at_lower_end, "The pick sits on the first "
            "grid value.")
        .def_readonly("at_upper_end", &tik::Lambda2PowerRuleResult::at_upper_end, "The pick sits on the last grid "
            "value, a boundary solution: widen the grid before trusting it.")
        .def_readonly("seconds",      &tik::Lambda2PowerRuleResult::seconds,      "Wall-clock seconds of the whole "
            "profile.")
        .def("selection", [](const tik::Lambda2PowerRuleResult& r) { return r.selection(); },
             "The selection at the chosen lambda_2 (0-based, ascending).")
        .def("counts", &tik::Lambda2PowerRuleResult::counts, "The selection counts, ascending in lambda_2.");

    py::class_<tik::Lambda2PilotResult>(m, "Lambda2PilotResult",
        "What the lambda_2 pilot of the last select() or run_pilot_only() did: the grid, the counts and curves per "
            "grid value, the pick and the size of the pilot.")
        .def(py::init<>())
        .def_readonly("ran",               &tik::Lambda2PilotResult::ran,               "A pilot ran (False: lambda_2 "
            "came from the caller or a cross-validation rule).")
        .def_readonly("grid",              &tik::Lambda2PilotResult::grid,              "The lambda_2 values the "
            "pilot ran at, ascending and without duplicates: pilot_grid times grid_scale.")
        .def_readonly("grid_scale",        &tik::Lambda2PilotResult::grid_scale,        "The factor between "
            "pilot_grid and grid: n - 1 under ScalingMode.ZSCORE, 1 otherwise.")
        .def_readonly("counts",            &tik::Lambda2PilotResult::counts,            "The count the pick was read "
            "off, per grid value: the count of pilot_statistic.")
        .def_readonly("curve_majority",    &tik::Lambda2PilotResult::curve_majority,    "The majority count per grid "
            "value and pilot stopping time t = 1 .. T_pilot.")
        .def_readonly("curve_soft",        &tik::Lambda2PilotResult::curve_soft,        "The mean number of "
            "variables ahead of the t-th dummy, per grid value and t (a diagnostic).")
        .def_readonly("calibrated",        &tik::Lambda2PilotResult::calibrated,        "The calibrated pilot "
            "selection count per grid value.")
        .def_readonly("group_majority",    &tik::Lambda2PilotResult::group_majority,    "The group-majority count "
            "per grid value: the members ahead, summed over the groups that reach pilot_group_fraction. Zeros "
            "without prior_groups.")
        .def_readonly("groups_qualifying", &tik::Lambda2PilotResult::groups_qualifying, "The number of groups that "
            "reach pilot_group_fraction, per grid value. Zeros without prior_groups.")
        .def_readonly("probe_seconds",     &tik::Lambda2PilotResult::probe_seconds,     "Wall-clock seconds per grid "
            "value.")
        .def_readonly("pick",              &tik::Lambda2PilotResult::pick,              "Lambda2PickInfo of the pick.")
        .def_readonly("index",             &tik::Lambda2PilotResult::index,             "0-based index of the pick "
            "into grid.")
        .def_readonly("lambda2",           &tik::Lambda2PilotResult::lambda2,           "The picked lambda_2, "
            "grid[index], on the working scale.")
        .def_readonly("at_lower_end",      &tik::Lambda2PilotResult::at_lower_end,      "The pick sits on the first "
            "grid value.")
        .def_readonly("at_upper_end",      &tik::Lambda2PilotResult::at_upper_end,      "The pick sits on the last "
            "grid value: widen the grid before trusting it.")
        .def_readonly("K_pilot",           &tik::Lambda2PilotResult::K_pilot,           "Pilot experiments per grid "
            "value that ran: pilot_K, capped at K.")
        .def_readonly("T_pilot",           &tik::Lambda2PilotResult::T_pilot,           "Stopping time of the pilot "
            "paths that was used: pilot_T, capped at min(L_0, ceil(n / 2)).")
        .def_readonly("L_0",               &tik::Lambda2PilotResult::L_0,               "Number of pilot dummies, "
            "pilot_L_multiplier times p.")
        .def_readonly("dummy_tag",         &tik::Lambda2PilotResult::dummy_tag,         "Seed tag under which the "
            "pilot dummies were drawn.")
        .def_readonly("seconds",           &tik::Lambda2PilotResult::seconds,           "Wall-clock seconds of the "
            "whole pilot.");

    py::class_<tik::DummyClusters>(m, "DummyClusters",
        "The clusters of the variables for which the dummies of a run were drawn. The fields have the names of "
            "the GVS result.")
        .def(py::init<>())
        .def_readonly("groups_vec",     &tik::DummyClusters::groups_vec,     "Cluster of each variable (0-based, "
            "length p). Empty when the dummies are not drawn per cluster or no run was made.")
        .def_readonly("max_clusters",   &tik::DummyClusters::max_clusters,   "Number of clusters (0 when there are "
            "none).")
        .def_readonly("hc_method_used", &tik::DummyClusters::hc_method_used, "Origin of the clusters: \"prior\" "
            "when prior_groups was supplied, else the linkage of the clustering, \"single\", \"complete\", "
            "\"average\" or \"wpgma\". Empty when there are no clusters.");

    py::class_<tik::TRexTikhonovControlParameter>(m, "TRexTikhonovControlParameter",
        "Control parameters of the Tikhonov selector: the penalty (a p x p PSD matrix tikhonov_K, sparse or dense, "
            "or GROUP_CENTERING built from groups), the rule for lambda_2, and the dummies.")
        .def(py::init<>())
        .def_readwrite("penalty_structure", &tik::TRexTikhonovControlParameter::penalty_structure, "Structural form "
            "of the penalty: SPARSE_GENERAL (default, tikhonov_K is the penalty) or GROUP_CENTERING (K from "
            "prior_groups, which it requires).")
        .def_property("tikhonov_K",
            [](const tik::TRexTikhonovControlParameter& c) { return c.tikhonov_K; },
            [](tik::TRexTikhonovControlParameter& c, const Eigen::SparseMatrix<double>& K) {
                c.tikhonov_K = K;
            },
            "The p x p PSD penalty matrix K. It is held sparse: assign a scipy.sparse matrix of any format, or a "
                "dense array, which is converted. It is read back as a scipy.sparse.csc_matrix. "
                "Required under SPARSE_GENERAL, where the selector checks its content at construction (finite "
                "entries, symmetry, positive semi-definiteness). Optional under GROUP_CENTERING, where the "
                "selector builds the centering matrix from prior_groups when it is left out. A supplied matrix "
                "is compared with that matrix, and another matrix raises ValueError.")
        .def_readwrite("lambda_2",        &tik::TRexTikhonovControlParameter::lambda_2,        "Weight of the "
            "penalty. A negative value (default -1) lets lambda2_method choose it, 0 is no penalty (the lasso path "
            "of the solver), and a positive value is used as it is. It must be finite. A supplied value refers to "
            "the scale of the standardized columns: under ScalingMode.ZSCORE the columns have the squared norm "
            "n - 1, so the same shrinkage needs the value times n - 1. The rules report their pick on this scale.")
        .def_readwrite("lambda2_method",  &tik::TRexTikhonovControlParameter::lambda2_method,  "Rule that chooses "
            "lambda_2 when lambda_2 < 0 (default PILOT_KNEE, the pilot rule; the cross-validation rules are the "
            "alternatives).")
        .def_readwrite("cv_n_folds",      &tik::TRexTikhonovControlParameter::cv_n_folds,      "Folds of the "
            "cross-validation rules, 2 to n (default 10).")
        .def_readwrite("cv_n_lambda",     &tik::TRexTikhonovControlParameter::cv_n_lambda,     "CV_*_TIK_SVD: points "
            "on the lambda_2 grid (default 100).")
        .def_readwrite("cv_n_lambda1",    &tik::TRexTikhonovControlParameter::cv_n_lambda1,    "CV_*_IEN_CCD: points "
            "on the profiled lambda_1 path per lambda_2 value (default 40).")
        .def_readwrite("cv_dfmax_ratio",  &tik::TRexTikhonovControlParameter::cv_dfmax_ratio,  "CV_*_IEN_CCD: support "
            "cap as a fraction of the smallest fold's training size (default 0.8). A value <= 0 disables the cap, "
            "larger values allow larger supports at a higher cost. It must be finite.")
        .def_readwrite("cv_n_lambda2",    &tik::TRexTikhonovControlParameter::cv_n_lambda2,    "CV_*_IEN_CCD: points "
            "on the lambda_2 grid (default 100).")
        .def_readwrite("cv_seed",         &tik::TRexTikhonovControlParameter::cv_seed,         "Seed of the fold "
            "permutation of the cross-validation rules. -1 (default) derives it from the seed of the selector, a "
            "value >= 0 replaces it for the folds only.")
        .def_readwrite("pilot_grid",      &tik::TRexTikhonovControlParameter::pilot_grid,      "lambda_2 grid of the "
            "pilot rule, stated for columns of unit norm. Under ScalingMode.ZSCORE the pilot runs at these values "
            "times n - 1, as the CV tuners return their pick, so one grid means the same shrinkage under every "
            "scaling mode. The pick is reported on the working scale and can be supplied as lambda_2. Keep 0 on "
            "the grid so that no penalty is an admissible answer. Default: default_lambda2_grid().")
        .def_readwrite("pilot_K",         &tik::TRexTikhonovControlParameter::pilot_K,         "Pilot experiments per "
            "grid value (default 10; at least 2, capped at K of the base control).")
        .def_readwrite("pilot_T",         &tik::TRexTikhonovControlParameter::pilot_T,         "Stopping time of the "
            "pilot paths (default 5; at least 1, capped at min(L_0, ceil(n / 2))).")
        .def_readwrite("pilot_L_multiplier", &tik::TRexTikhonovControlParameter::pilot_L_multiplier, "Pilot "
            "dummies L_0 = pilot_L_multiplier times p (default 1; in [1, max_dummy_multiplier] of the base "
            "control).")
        .def_readwrite("pilot_knee_tolerance", &tik::TRexTikhonovControlParameter::pilot_knee_tolerance, "Pilot "
            "rule: selections below the largest count that still count as on the plateau (default 1; 0 reads the "
            "largest count).")
        .def_readwrite("pilot_statistic", &tik::TRexTikhonovControlParameter::pilot_statistic, "The count the pilot "
            "rule picks on: MAJORITY (default), CALIBRATED or GROUP_MAJORITY.")
        .def_readwrite("pilot_group_fraction", &tik::TRexTikhonovControlParameter::pilot_group_fraction,
            "GROUP_MAJORITY: the share of a group's members that must be ahead for the group to count, in (0, 1] "
            "(default 0.5).")
        .def_readwrite("pilot_shared_dummies", &tik::TRexTikhonovControlParameter::pilot_shared_dummies, "False "
            "(default): the pilot draws its dummies from a seed stream of its own, so the dummies of the main run "
            "do not depend on the pick. True: the pilot uses the first dummy block of the main strategy (not with "
            "the PERMUTATION strategies).")
        .def_readwrite("cluster_dummies", &tik::TRexTikhonovControlParameter::cluster_dummies, "Draw the dummies per "
            "cluster of correlated variables, as in TRexGVSSelector (default True). False gives i.i.d. dummies, "
            "for validation runs only: with an informed K they deform the exchangeability of dummies and null "
            "variables that the FDR calibration relies on.")
        .def_readwrite("prior_groups",    &tik::TRexTikhonovControlParameter::prior_groups,    "Groups of the "
            "variables: a 0-based group ID per variable (length p, IDs without gaps). They define the clusters of "
            "the dummies (with cluster_dummies), the penalty under GROUP_CENTERING and the count GROUP_MAJORITY of "
            "the pilot rule. Empty (default): the clusters of the dummies are found by clustering the "
            "correlations with corr_max and hc_linkage.")
        .def_readwrite("corr_max",        &tik::TRexTikhonovControlParameter::corr_max,        "Correlation threshold "
            "of the automatic clustering (default 0.5, in [0, 1]). The dendrogram is cut at the height "
            "1 - corr_max. With Single linkage (default) no two variables of different clusters correlate above "
            "corr_max. With Complete linkage every two variables of one cluster correlate at least at corr_max, "
            "and variables of different clusters can exceed it. Ignored when prior_groups is set or "
            "cluster_dummies is False.")
        .def_readwrite("hc_linkage",      &tik::TRexTikhonovControlParameter::hc_linkage,      "Linkage of the "
            "hierarchical clustering (ml_methods.clustering.LinkageMethod, default Single). Supported: Single, "
            "Complete, Average, WPGMA. Ward is rejected, because its update is valid for squared Euclidean "
            "distances only and the correlation distance 1 - abs(cor) is not Euclidean. Ignored when prior_groups "
            "is set or cluster_dummies is False.")
        .def_readwrite("trex_ctrl",       &tik::TRexTikhonovControlParameter::trex_ctrl,       "Nested base T-Rex "
            "control, whose solver is TIKNET, the one solver of the selector. It applies when TRexTikhonovSelector "
            "gets no trex_control argument, and that argument replaces it otherwise.");

    m.def("gamma_to_k",
          [](const Eigen::SparseMatrix<double>& gamma) {
              return tik::TRexTikhonovSelector::gammaToK(gamma);
          },
          py::arg("gamma"),
          "Build the penalty matrix K = Gamma^T Gamma from a (q x p) Tikhonov operator Gamma, a scipy.sparse "
          "matrix or a dense array. Returns a scipy.sparse.csc_matrix (p x p).");

    m.def("build_group_centering_k",
          [](const std::vector<Eigen::Index>& prior_groups) {
              return tik::TRexTikhonovSelector::buildGroupCenteringK(prior_groups);
          },
          py::arg("prior_groups"),
          "Build the GROUP_CENTERING penalty K = blockdiag(I - J/p_m) from a 0-based group assignment. Returns a "
          "scipy.sparse.csc_matrix (p x p).");

    m.def("default_lambda2_grid", &tik::default_lambda2_grid,
          "The default lambda_2 grid of the pilot rule and of the power rule: 0, 1, 3, 10, 30, 100, 300, 1000, "
          "3000. It holds 0, so that no penalty is an admissible answer.");

    m.def("gvs_lambda2_power_rule",
          [](Eigen::Ref<Eigen::MatrixXd> X, Eigen::Ref<Eigen::VectorXd> y, double tFDR,
             const TRexGVSControlParameter& gvs_control,
             const std::optional<TRexControlParameter>& trex_control,
             int seed, const std::vector<double>& grid, double knee_tolerance, bool verbose) {
              // Every probe normalizes X in place and hands it back on return,
              // so X and y are read through maps and never copied.
              Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
              Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());
              tik::Lambda2PowerRuleConfig cfg;
              cfg.grid           = grid;
              cfg.knee_tolerance = knee_tolerance;
              return tik::gvs_lambda2_power_rule(X_map, y_map, tFDR, merge_gvs_control(gvs_control, trex_control),
                                                 seed, cfg, verbose);
          },
          py::arg("X"), py::arg("y"), py::arg("tFDR") = 0.1,
          py::arg("gvs_control") = TRexGVSControlParameter(),
          py::arg("trex_control") = py::none(),
          py::arg("seed") = 0,
          py::arg("grid") = tik::default_lambda2_grid(),
          py::arg("knee_tolerance") = 1.0, py::arg("verbose") = false,
          "The power rule on the GVS selector: one FDR-controlled selection per value of grid for lambda_2, the "
          "ridge part of the elastic net or of the informed elastic net, with the same seed (common dummies), then "
          "the smallest lambda_2 whose selection is within knee_tolerance selections of the largest. The grid is "
          "stated for columns of unit norm, and under ScalingMode.ZSCORE the probes run at its values times n - 1. "
          "gvs_control and trex_control are taken as by TRexGVSSelector (without trex_control the base control "
          "nested in gvs_control applies), and lambda_2 of gvs_control is replaced per probe. X (F-contiguous "
          "float64) and y are accessed zero-copy: every probe normalizes X in place and hands it back. The pick "
          "depends on y and lies outside the fixed-lambda_2 FDR guarantee that each probe "
          "carries. A negative seed is refused.");

    m.def("tikhonov_lambda2_power_rule",
          [](Eigen::Ref<Eigen::MatrixXd> X, Eigen::Ref<Eigen::VectorXd> y, double tFDR,
             const tik::TRexTikhonovControlParameter& tik_control,
             const std::optional<TRexControlParameter>& trex_control,
             int seed, const std::vector<double>& grid, double knee_tolerance, bool verbose) {
              // Every probe normalizes X in place and hands it back on return,
              // so X and y are read through maps and never copied.
              Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
              Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());
              tik::Lambda2PowerRuleConfig cfg;
              cfg.grid           = grid;
              cfg.knee_tolerance = knee_tolerance;
              const tik::TRexTikhonovControlParameter tik_ctrl =
                  merge_tikhonov_control(tik_control, trex_control);
              return tik::tikhonov_lambda2_power_rule(X_map, y_map, tFDR, tik_ctrl, seed, cfg, verbose);
          },
          py::arg("X"), py::arg("y"), py::arg("tFDR") = 0.1,
          py::arg("tik_control") = tik::TRexTikhonovControlParameter(),
          py::arg("trex_control") = py::none(),
          py::arg("seed") = 0,
          py::arg("grid") = tik::default_lambda2_grid(),
          py::arg("knee_tolerance") = 1.0, py::arg("verbose") = false,
          "The knee of the full lambda_2 profile of the Tikhonov T-Rex: one FDR-controlled selection per grid value "
          "with the same seed (common dummies), then the smallest lambda_2 whose selection is within knee_tolerance "
          "selections of the largest (0 reads the largest selection, ties toward the smallest lambda_2). This is the "
          "reference rule behind TikLambda2Method.PILOT_KNEE at the cost of one selection per grid value; it returns "
          "every probe's selection. The grid is stated for columns of unit norm, and under ScalingMode.ZSCORE the "
          "probes run at its values times n - 1, as the pilot rule does. "
          "X (F-contiguous float64) and y are accessed zero-copy: every probe normalizes X "
          "in place and hands it back, so X is unchanged on return. The pick depends on y and lies outside the "
          "fixed-lambda_2 FDR guarantee that each probe carries.");

    py::class_<PyTRexTikhonovSelector, PyTRexSelector>(m, "TRexTikhonovSelector",
        "T-Rex selector with a general Tikhonov penalty (T-Rex+Tikhonov): FDR-controlled variable selection\n"
        "with the penalty (lambda_2 / 2) * beta^T K beta, whose positive semi-definite matrix K carries the\n"
        "structure among the variables.\n"
        "\n"
        "Notes\n"
        "-----\n"
        "The penalty. K = Gamma^T Gamma (p x p) for a regularization operator Gamma, see gamma_to_k().\n"
        "Special cases: K = I is the ridge part of the elastic net. The group-mean matrix is the penalty of\n"
        "the informed elastic net of TRexGVSSelector. The group-centering matrix\n"
        "(TikPenaltyStructure.GROUP_CENTERING) makes the members of a group share a coefficient. A banded\n"
        "matrix or a graph Laplacian makes neighbouring coefficients smooth.\n"
        "\n"
        "The solver. TIKNET, the homotopy solver, is the one solver of the selector and the default of the\n"
        "nested base control. Every other solver is refused.\n"
        "\n"
        "The dummies. By default they are drawn per cluster of correlated variables, as in TRexGVSSelector\n"
        "(cluster_dummies): the clusters are prior_groups, or come from a hierarchical clustering of the\n"
        "correlations (corr_max, hc_linkage). dummy_clusters reports them. In the penalty every dummy takes\n"
        "the place of the variable it was drawn for (the folded coupling), so dummies and null variables\n"
        "stay exchangeable for every K.\n"
        "\n"
        "L-loop strategies. With cluster dummies STANDARD, HCONCAT, SKIPL and SEEDED are supported, and\n"
        "SEEDED equals STANDARD. PERMUTATION and PERMUTATION_SEEDED are refused, because a row permutation\n"
        "destroys the covariance of the dummies within a cluster. With cluster_dummies = False the two run\n"
        "as well.\n"
        "\n"
        "lambda_2. A supplied value is used as it is. Otherwise lambda2_method chooses it before the L-loop.\n"
        "PILOT_KNEE (default), the pilot rule: for every value of pilot_grid, pilot_K random experiments run\n"
        "to the stopping time pilot_T on one set of dummies, and the variables ahead of the dummies in more\n"
        "than half of the experiments are counted. The pick is the smallest lambda_2 whose count is within\n"
        "pilot_knee_tolerance of the largest. The pilot dummies come from a seed stream of their own, so the\n"
        "main run is a T-Rex run at a fixed lambda_2 whose dummies do not depend on the pick. The pick\n"
        "itself depends on y and lies outside the FDR guarantee of a fixed lambda_2. In simulations the pick\n"
        "cost nothing under the global null and at most 0.02 of true positive rate against the best fixed\n"
        "value under the alternatives. The record is lambda2_pilot, and tikhonov_lambda2_power_rule() runs\n"
        "the same rule on full selections. The rules CV_*_IEN_CCD cross-validate with the lambda_1 path\n"
        "profiled out and are well-posed for every p and n. The rules CV_*_TIK_SVD are an analytic ridge\n"
        "cross-validation, exact for n > p.\n"
        "\n"
        "The scale. A supplied lambda_2 refers to the scale of the standardized columns: under\n"
        "ScalingMode.ZSCORE the columns have the squared norm n - 1, so the same shrinkage needs the value\n"
        "times n - 1. pilot_grid is stated for columns of unit norm, and the pilot rule and the\n"
        "cross-validation rules apply the factor n - 1 themselves. Every rule reports its pick on the scale\n"
        "of a supplied value (lambda2_used).\n"
        "\n"
        "A selector runs once: select() and run_pilot_only() exclude each other.\n"
        "\n"
        "References\n"
        "----------\n"
        "A. N. Tikhonov, \"Solution of Incorrectly Formulated Problems and the Regularization\n"
        "  Method,\" Soviet Math. Dokl., vol. 4, pp. 1035-1038, 1963 (English translation of Dokl.\n"
        "  Akad. Nauk SSSR 151 (1963), no. 3, 501-504).\n"
        "Y. Park, L. Reichel, G. Rodriguez, and X. Yu, \"Parameter Determination for Tikhonov\n"
        "  Regularization Problems in General Form,\" J. Comput. Appl. Math., vol. 343, pp. 12-25,\n"
        "  2018, doi: 10.1016/j.cam.2018.04.049.\n"
        "A. J. Smola and R. Kondor, \"Kernels and Regularization on Graphs,\" in Learning Theory\n"
        "  and Kernel Machines (COLT/Kernel 2003), Lecture Notes in Artificial Intelligence\n"
        "  vol. 2777, 2003, pp. 144-158.\n"
        "F. R. K. Chung, Spectral Graph Theory. CBMS Regional Conference Series in Mathematics,\n"
        "  no. 92. American Mathematical Society, 1997.\n"
        "A. E. Hoerl and R. W. Kennard, \"Ridge Regression: Biased Estimation for Nonorthogonal\n"
        "  Problems,\" Technometrics, vol. 42, no. 1, pp. 80-86, 2000 (Reprint of the 1970 article\n"
        "  in the Special 40th Anniversary Issue).\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n")
        .def(py::init<Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::VectorXd>, double,
                      const tik::TRexTikhonovControlParameter&, const std::optional<TRexControlParameter>&,
                      int, bool>(),
             py::arg("X"), py::arg("y"), py::arg("tFDR") = 0.1,
             py::arg("tik_control") = tik::TRexTikhonovControlParameter(),
             py::arg("trex_control") = py::none(),
             py::arg("seed") = -1, py::arg("verbose") = true,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             "Construct the Tikhonov selector. X and y are accessed zero-copy via Eigen::Map. Without "
             "trex_control the base control nested in tik_control applies, whose solver is TIKNET. A "
             "trex_control replaces it as it is, and its solver_type must be TIKNET: any other solver, "
             "TLARS included, raises ValueError.")
        .def(py::init<const DesignArray&, Eigen::Ref<Eigen::VectorXd>, double,
                      const tik::TRexTikhonovControlParameter&, const std::optional<TRexControlParameter>&,
                      int, bool>(),
             py::arg("X").noconvert(), py::arg("y"), py::arg("tFDR") = 0.1,
             py::arg("tik_control") = tik::TRexTikhonovControlParameter(),
             py::arg("trex_control") = py::none(),
             py::arg("seed") = -1, py::arg("verbose") = true,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             kReaderSelectorDoc.c_str())
        .def("getLambda2Used", &PyTRexTikhonovSelector::getLambda2Used, "The lambda_2 the last run used: the "
            "supplied value, or the pick of the pilot rule or of a cross-validation rule, on the scale of a "
            "supplied value. -1 before a run.")
        .def("getLambda2CVDiagnostics", &PyTRexTikhonovSelector::getLambda2CVDiagnostics, "Return the lambda_2 tuner "
            "health record of the last select().")
        .def("getLambda2PilotResult", &PyTRexTikhonovSelector::getLambda2PilotResult, "Return the lambda_2 pilot "
            "record "
            "of the last select().")
        .def("runPilotOnly", &PyTRexTikhonovSelector::runPilotOnly, py::call_guard<py::gil_scoped_release>(),
             "Run only the lambda_2 pilot (the cluster model, the pilot experiments and the pick), give X back and "
             "return the Lambda2PilotResult. It needs lambda_2 < 0 and lambda2_method = PILOT_KNEE (ValueError "
             "otherwise). A selector runs once: after it select() and runPilotOnly() raise RuntimeError, and so "
             "does runPilotOnly() after select().")
        .def("getDummyClusters", &PyTRexTikhonovSelector::getDummyClusters, "Return the DummyClusters of the last "
            "run: the clusters the dummies were drawn for.");
}
// =====================================================================================
#endif /* End of TREX_PYTHON_TREX_TIKHONOV_BINDINGS_HPP */

// =============================================================================
// rcpp_trex_spca.cpp - Rcpp bindings for TRexSPCASelector
// =============================================================================

// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "rcpp_trex_spca.h"
#include <utils/memmap/memory_mapped_matrix.hpp>

using namespace Rcpp;
using namespace trex::trex_selector_methods::ad::trex_spca;

namespace tc = trex::trex_selector_methods::ad::trex_core;
namespace tg = trex::trex_selector_methods::ad::trex_gvs;

/*
 * Note on Roxygen Documentation:
 * All functions exported via Rcpp in this file use @noRd. They are internal
 * C++ bindings wrapped by the TRexSPCASelector R6 class.
 */

// =============================================================================

/**
 * @brief Parses an R list into a TRexSPCAControl C++ struct.
 *
 * @param spca_ctrl_list R list of trex_spca_control(): mode, lambda2_ridge_loadings,
 *   gvs_type, lambda_2, lambda2_method, en_solver, tenet_aug_use_lars, K,
 *   max_dummy_multiplier. Other keys are not read.
 * @return Populated TRexSPCAControl instance.
 */
TRexSPCAControl parse_spca_control(const Rcpp::List& spca_ctrl_list) {
    TRexSPCAControl ctrl;

    if (spca_ctrl_list.containsElementNamed("mode")) {
        std::string mode_str = spca_ctrl_list["mode"];
        if (mode_str == "ActiveSet") {
            ctrl.mode = SPCAMode::ActiveSet;
        } else if (mode_str == "Thresholded") {
            ctrl.mode = SPCAMode::Thresholded;
        } else {
            Rcpp::stop("Unknown SPCAMode: '" + mode_str +
                       "'. Use 'ActiveSet' or 'Thresholded'.");
        }
    }

    if (spca_ctrl_list.containsElementNamed("lambda2_ridge_loadings")) {
        ctrl.lambda2_ridge_loadings = spca_ctrl_list["lambda2_ridge_loadings"];
    }

    // --- GVS sub-selector control ---
    if (spca_ctrl_list.containsElementNamed("gvs_type")) {
        std::string gvs_str = spca_ctrl_list["gvs_type"];
        if (gvs_str == "EN") {
            ctrl.gvs_ctrl.gvs_type = tg::GVSType::EN;
        } else if (gvs_str == "IEN") {
            ctrl.gvs_ctrl.gvs_type = tg::GVSType::IEN;
        } else {
            Rcpp::stop("Unknown GVSType: '" + gvs_str + "'. Use 'EN' or 'IEN'.");
        }
    }

    if (spca_ctrl_list.containsElementNamed("lambda_2")) {
        ctrl.gvs_ctrl.lambda_2 = spca_ctrl_list["lambda_2"];
    }

    if (spca_ctrl_list.containsElementNamed("en_solver")) {
        std::string solver = spca_ctrl_list["en_solver"];
        if (solver == "TENET") {
            ctrl.en_solver = tg::ENSolverType::TENET;
        } else if (solver == "TENET_AUG") {
            ctrl.en_solver = tg::ENSolverType::TENET_AUG;
        } else {
            Rcpp::stop("Unknown ENSolverType: '" + solver +
                       "'. Use 'TENET' or 'TENET_AUG'.");
        }
    }

    // LARS or LASSO inner path of the TENET_AUG solver. The selector passes
    // gvs_ctrl to the T-Rex+GVS run of each component and replaces only the
    // solver fields, so the flag reaches every run.
    if (spca_ctrl_list.containsElementNamed("tenet_aug_use_lars")) {
        ctrl.gvs_ctrl.tenet_aug_use_lars =
            Rcpp::as<bool>(spca_ctrl_list["tenet_aug_use_lars"]);
    }

    if (spca_ctrl_list.containsElementNamed("lambda2_method")) {
        std::string method_str = spca_ctrl_list["lambda2_method"];
        if (method_str == "CV_1SE_SVD") {
            ctrl.gvs_ctrl.lambda2_method = tg::LambdaSelectionMethod::CV_1SE_SVD;
        } else if (method_str == "CV_MIN_SVD") {
            ctrl.gvs_ctrl.lambda2_method = tg::LambdaSelectionMethod::CV_MIN_SVD;
        } else if (method_str == "CV_1SE_CCD") {
            ctrl.gvs_ctrl.lambda2_method = tg::LambdaSelectionMethod::CV_1SE_CCD;
        } else if (method_str == "CV_MIN_CCD") {
            ctrl.gvs_ctrl.lambda2_method = tg::LambdaSelectionMethod::CV_MIN_CCD;
        } else {
            Rcpp::stop("Unknown lambda2_method: '" + method_str +
                       "'. Use 'CV_1SE_SVD', 'CV_MIN_SVD', 'CV_1SE_CCD', or 'CV_MIN_CCD'.");
        }
    }

    // --- Base T-Rex sub-selector control (nested inside the GVS control) ---
    if (spca_ctrl_list.containsElementNamed("K")) {
        ctrl.gvs_ctrl.trex_ctrl.K = static_cast<std::size_t>(
            static_cast<int>(spca_ctrl_list["K"]));
    }
    if (spca_ctrl_list.containsElementNamed("max_dummy_multiplier")) {
        ctrl.gvs_ctrl.trex_ctrl.max_dummy_multiplier = static_cast<std::size_t>(
            static_cast<int>(spca_ctrl_list["max_dummy_multiplier"]));
    }

    return ctrl;
}


// ==============================================================================
// Shared body of the in-memory and the memory-mapped entry points: the
// selector centers X in place through the map during select() and restores
// it on return.
// ==============================================================================
static Rcpp::List spca_select_on_map(
    Eigen::Map<Eigen::MatrixXd>& X_map,
    int M,
    double tFDR,
    const Rcpp::List& spca_ctrl_list,
    int seed,
    bool verbose
) {
    TRexSPCAControl ctrl = parse_spca_control(spca_ctrl_list);

    TRexSPCASelector spca(X_map, static_cast<Eigen::Index>(M), tFDR, ctrl, seed, verbose);
    const TRexSPCAResult& res = spca.select();   // wrapped once below, not copied here

    // Convert active sets to R-friendly 1-based integer vectors
    Rcpp::List active_sets_r(static_cast<R_xlen_t>(M));
    for (int m = 0; m < M; ++m) {
        const auto& as_cpp = res.active_sets[static_cast<std::size_t>(m)];
        Rcpp::IntegerVector as_r(static_cast<R_xlen_t>(as_cpp.size()));
        for (std::size_t i = 0; i < as_cpp.size(); ++i) {
            as_r[static_cast<R_xlen_t>(i)] =
                static_cast<int>(as_cpp[i]) + 1; // 0-based → 1-based
        }
        active_sets_r[m] = as_r;
    }

    return Rcpp::List::create(
        Rcpp::Named("Z")            = res.Z,
        Rcpp::Named("V")            = res.V,
        Rcpp::Named("active_sets")  = active_sets_r,
        Rcpp::Named("adjusted_ev")  = res.adjusted_ev,
        Rcpp::Named("cumulative_ev")= res.cumulative_ev
    );
}

//' @title Run the T-Rex sparse PCA selection
//'
//' @description The ordinary PCA of the centered X, one T-Rex+GVS run with the
//'   elastic net per component, and the assembly of the sparse loadings.
//'
//' @param X Design matrix (n x p). Centered in place and restored on return.
//' @param M Number of sparse principal components, 1 to min(n, p).
//' @param tFDR Target FDR level of each component, between 0 and 1.
//' @param spca_ctrl_list The list of trex_spca_control().
//' @param seed Random seed. A negative value gives a run that is not reproducible.
//' @param verbose Whether to print one line per component.
//'
//' @return Named list:
//'   \itemize{
//'     \item \code{Z} Scores of the sparse components (n x M).
//'     \item \code{V} Sparse loadings (p x M).
//'     \item \code{active_sets} List of M integer vectors (1-based indices).
//'     \item \code{adjusted_ev} Adjusted explained variance per component (length M).
//'     \item \code{cumulative_ev} Cumulative adjusted explained variance as a fraction of
//'       the total variance (length M).
//'   }
//' @noRd
// [[Rcpp::export]]
Rcpp::List trex_spca_select(
    Rcpp::NumericMatrix X,
    int M,
    double tFDR,
    Rcpp::List spca_ctrl_list,
    int seed,
    bool verbose
) {
    // Zero-copy view over the R buffer (X is REALSXP by contract — the R6 layer
    // enforces double storage, so no silent integer→double copy happens here).
    // No prot slot needed: the TRexSPCASelector object lives and dies within this call,
    // so its map never outlives X.
    Eigen::Map<Eigen::MatrixXd> X_map(X.begin(), X.nrow(), X.ncol());
    return spca_select_on_map(X_map, M, tFDR, spca_ctrl_list, seed, verbose);
}

//' @title Sparse PCA selection on a MemoryMappedMatrix (zero-copy)
//'
//' @description The file is centered in place during the selection and
//'   restored on return ("readwrite" mapping required). Before, the R6 class
//'   read the whole file into RAM through as.matrix().
//' @noRd
// [[Rcpp::export]]
Rcpp::List trex_spca_select_mmap(
    Rcpp::XPtr<trex::utils::memmap::MemoryMappedMatrix<double>> X_ptr, // NOLINT(performance-unnecessary-value-param)
    int M,
    double tFDR,
    Rcpp::List spca_ctrl_list,
    int seed,
    bool verbose
) {
    // The selector centers the design in every run, whatever the scaling mode.
    if (X_ptr->mode() == trex::utils::memmap::AccessMode::ReadOnly) {
        Rcpp::stop("TRexSPCASelector centers the design in place and restores it on return, "
                   "so an mmap_matrix must be open in \"readwrite\" mode.");
    }
    auto X_map = X_ptr->getMap();
    return spca_select_on_map(X_map, M, tFDR, spca_ctrl_list, seed, verbose);
}

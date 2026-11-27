// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "TRexSelector_types.h"   // TSolverRcpp, the zero-copy holder
#include <fstream>
#include <memory>
#include <cereal/archives/portable_binary.hpp>
#include <tsolvers/ad/tsolver_base.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_lars.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_lasso.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_stepwise.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_enet.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_enet_aug.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_ienet_aug.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_stagewise.hpp>
#include <tsolvers/ad/linear_model/omp_based/tsolver_omp.hpp>
#include <tsolvers/ad/linear_model/omp_based/tsolver_gp.hpp>
#include <tsolvers/ad/linear_model/omp_based/tsolver_acgp.hpp>
#include <tsolvers/ad/linear_model/omp_based/tsolver_mp.hpp>
#include <tsolvers/ad/linear_model/omp_based/tsolver_ools.hpp>
#include <tsolvers/ad/linear_model/omp_based/tsolver_ncgmp.hpp>
#include <tsolvers/ad/linear_model/afs_based/tsolver_afs.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_ienet.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_tiknet.hpp>
#include <utils/memmap/memory_mapped_matrix.hpp>

// ==============================================================================

/*
 * Note on Roxygen Documentation:
 * All functions exported via Rcpp in this file are documented using standard Roxygen tags
 * but strictly include the `@noRd` tag. These are internal C++ bindings that are wrapped
 * by user-friendly R6 classes and R functions in the package namespace. Generating `.Rd`
 * manuals for these endpoints would clutter the public API.
 */

using namespace Rcpp;
using namespace trex::utils::memmap;
using namespace trex::tsolvers::ad;
using namespace trex::tsolvers::ad::linear_model::lars_based;
using namespace trex::tsolvers::ad::linear_model::omp_based;
using namespace trex::tsolvers::ad::linear_model::afs_based;

// ==============================================================================
// TSolverRcpp — zero-copy holder, defined once in TRexSelector_types.h
// ==============================================================================
//
// The holder keeps the Eigen maps over the R matrices (or the mmap file) and
// owns the solver; RcppExports.cpp instantiates the XPtr finalizer from the
// same definition. See the header for the in-place contract.

/**
 * @brief Eigen::Map of X or D of a solver, without a copy.
 *
 * @param m    R double matrix (in memory) or mmap_matrix (the handle of a
 *             memory-mapped file).
 * @param opts Solver options. A writer needs a memory-mapped file opened
 *             "readwrite", a reader takes one opened "readonly" as well.
 * @param name Argument name for the message.
 *
 * @return Eigen::Map of the R matrix or of the memory-mapped file.
 */
static Eigen::Map<Eigen::MatrixXd> design_map_of(SEXP m, const TSolverOptions& opts, const char* name) {
    if (Rf_inherits(m, "mmap_matrix")) {
        XPtr<MemoryMappedMatrix<double>> ptr(m);
        return ptr->getMapFor(opts.writesDesign());
    }
    if (!Rf_isMatrix(m) || TYPEOF(m) != REALSXP) {
        stop(std::string(name) + " must be a double matrix or an mmap_matrix.");
    }
    return {REAL(m), Rf_nrows(m), Rf_ncols(m)};
}

/** @brief Map an R double vector without copying. */
static Eigen::Map<Eigen::VectorXd> map_of(Rcpp::NumericVector& v) {
    return {v.begin(), v.size()};
}

// ==============================================================================
// TSolver Instance Creators
// ==============================================================================
//
// Every solver is constructed as (X, D, y, options) — plus the 1-based group
// vector of the IEN solvers (TIENETAug, TIENET) or the Tikhonov
// matrix of the sparse-K solver TIKNET — where
// `options` is the R list built by tsolver_options().
// parse_solver_options() maps it onto the C++ TSolverOptions; absent fields
// keep the C++ defaults.

/**
 * @brief Parse the tsolver_options() list into trex::tsolvers::ad::TSolverOptions.
 */
TSolverOptions parse_solver_options(const Rcpp::List& o) {
    TSolverOptions opts;
    auto has = [&](const char* name) { return o.containsElementNamed(name); };

    if (has("normalize")) opts.normalize = Rcpp::as<bool>(o["normalize"]);
    if (has("intercept")) opts.intercept = Rcpp::as<bool>(o["intercept"]);
    if (has("precentered")) opts.precentered = Rcpp::as<bool>(o["precentered"]);
    if (has("verbose"))   opts.verbose   = Rcpp::as<bool>(o["verbose"]);
    if (has("scaling_mode")) {
        const std::string m = Rcpp::as<std::string>(o["scaling_mode"]);
        if (m == "L2")          opts.scaling_mode = ScalingMode::L2;
        else if (m == "ZSCORE") opts.scaling_mode = ScalingMode::ZSCORE;
        else Rcpp::stop("tsolver_options: scaling_mode must be \"L2\" or \"ZSCORE\"; got " + m);
    }
    if (has("tol"))            opts.tol            = Rcpp::as<double>(o["tol"]);
    if (has("tie_seed"))       opts.tie_seed       = static_cast<long long>(Rcpp::as<double>(o["tie_seed"]));
    if (has("n_threads"))      opts.n_threads      = Rcpp::as<int>(o["n_threads"]);
    if (has("exch_tie_alpha")) opts.exch_tie_alpha = Rcpp::as<double>(o["exch_tie_alpha"]);
    if (has("exch_tie_floor")) opts.exch_tie_floor = Rcpp::as<double>(o["exch_tie_floor"]);
    if (has("lambda2"))        opts.lambda2        = Rcpp::as<double>(o["lambda2"]);
    if (has("rho_afs"))        opts.rho_afs        = Rcpp::as<double>(o["rho_afs"]);
    if (has("ncgmp_variant")) {
        const std::string v = Rcpp::as<std::string>(o["ncgmp_variant"]);
        if (v == "line_search")           opts.ncgmp_variant = NCGMPVariant::LineSearch;
        else if (v == "fully_corrective") opts.ncgmp_variant = NCGMPVariant::FullyCorrective;
        else Rcpp::stop("tsolver_options: ncgmp_variant must be \"line_search\" or \"fully_corrective\"; got " + v);
    }
    if (has("tenet_aug_use_lars")) opts.tenet_aug_use_lars = Rcpp::as<bool>(o["tenet_aug_use_lars"]);
    if (has("tikhonov_validated")) opts.tikhonov_validated = Rcpp::as<bool>(o["tikhonov_validated"]);
    return opts;
}

/** @brief 1-based R group ids -> 0-based contiguous ids (the solver contract). */
static Eigen::VectorXi groups_from_r(const IntegerVector& groups) {
    Eigen::VectorXi c_groups(groups.size());
    for (int i = 0; i < groups.size(); ++i) c_groups[i] = groups[i] - 1;
    return c_groups;
}

/** @brief A Tikhonov matrix from R: a Matrix::dgCMatrix (S4) or a dense
 *  matrix (converted via sparseView()); same contract as the Tikhonov
 *  selector's tikhonov_K. */
static Eigen::SparseMatrix<double> tsolver_as_sparse_matrix(SEXP m) {
    if (Rf_isS4(m)) {
        return Rcpp::as<Eigen::SparseMatrix<double>>(m);
    }
    const Eigen::MatrixXd dense = Rcpp::as<Eigen::MatrixXd>(m);
    return dense.sparseView();
}

// ==============================================================================
// Instance creators (X and D: an R matrix in memory or a memory-mapped file)
// ==============================================================================

//' @title Create TSolverLARS
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_lars_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverLARS(h->X_map_, h->D_map_, h->y_map_, opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

//' @title Create TSolverLASSO
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_tlasso_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverLASSO(h->X_map_, h->D_map_, h->y_map_, opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

//' @title Create TSolverStepwise
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_tstepwise_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverStepwise(h->X_map_, h->D_map_, h->y_map_, opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

//' @title Create TSolverStagewise
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_tstagewise_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverStagewise(h->X_map_, h->D_map_, h->y_map_, opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

//' @title Create TSolverENET
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_tenet_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverENET(h->X_map_, h->D_map_, h->y_map_, opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

//' @title Create TSolverENETAug
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_tenet_aug_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverENETAug(h->X_map_, h->D_map_, h->y_map_, opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

//' @title Create TSolverIENETAug
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param groups 1-based group ids, one per predictor
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_tienet_aug_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    IntegerVector groups,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverIENETAug(h->X_map_, h->D_map_, h->y_map_, groups_from_r(groups), opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

//' @title Create TSolverIENET
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param groups 1-based group ids, one per predictor
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_tienet_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    IntegerVector groups,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverIENET(h->X_map_, h->D_map_, h->y_map_, groups_from_r(groups), opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

//' @title Create TSolverTIKNET
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param tikhonov_K p x p PSD Tikhonov matrix (dgCMatrix or dense)
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_tiknet_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    SEXP tikhonov_K,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverTIKNET(h->X_map_, h->D_map_, h->y_map_, tsolver_as_sparse_matrix(tikhonov_K), opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

//' @title Create TSolverOMP
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_omp_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverOMP(h->X_map_, h->D_map_, h->y_map_, opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

//' @title Create TSolverGP
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_tgp_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverGP(h->X_map_, h->D_map_, h->y_map_, opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

//' @title Create TSolverACGP
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_tacgp_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverACGP(h->X_map_, h->D_map_, h->y_map_, opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

//' @title Create TSolverMP
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_tmp_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverMP(h->X_map_, h->D_map_, h->y_map_, opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

//' @title Create TSolverOOLS
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_tools_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverOOLS(h->X_map_, h->D_map_, h->y_map_, opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

//' @title Create TSolverNCGMP
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_tncgmp_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverNCGMP(h->X_map_, h->D_map_, h->y_map_, opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

//' @title Create TSolverAFS
//'
//' @param X Design matrix or mmap_matrix
//' @param D Dummy matrix or mmap_matrix
//' @param y Response vector
//' @param options List from tsolver_options()
//'
//' @return XPtr to TSolverRcpp
//' @noRd
// [[Rcpp::export]]
XPtr<TSolverRcpp> tsolver_tafs_create(
    SEXP X,
    SEXP D,
    Rcpp::NumericVector y,
    Rcpp::List options
) {
    const TSolverOptions opts = parse_solver_options(options);
    auto h = std::make_unique<TSolverRcpp>(design_map_of(X, opts, "X"), design_map_of(D, opts, "D"),
                                           map_of(y));
    h->solver_ = new TSolverAFS(h->X_map_, h->D_map_, h->y_map_, opts);
    // The protection slot pins X, D and y until the finalizer has run (a
    // writer's restore() writes back through the maps).
    return XPtr<TSolverRcpp>(h.release(), true, R_NilValue,
                             Rcpp::List::create(X, D, y));
}

// ==============================================================================
// Common TSolver API endpoints (called via XPtr<TSolverRcpp>)
// ==============================================================================

//' @title Execute Step in TSolver
//'
//' @param ptr TSolverRcpp holder pointer
//' @param T_stop Number of steps to compute
//' @param early_stop Whether to stop early based on constraints
//' @noRd
// [[Rcpp::export]]
void tsolver_execute_step(
    XPtr<TSolverRcpp> ptr, // NOLINT(performance-unnecessary-value-param)
    int T_stop,
    bool early_stop = true
) {
    // Treat -1 or any negative as full computation (max size_t)
    std::size_t c_T_stop = (T_stop < 0) ? -1 : static_cast<std::size_t>(T_stop);
    ptr->solver_->executeStep(c_T_stop, early_stop);
}


//' @title Get Beta from TSolver
//'
//' @param ptr TSolverRcpp holder pointer
//' @param step Step to retrieve beta for
//'
//' @return Beta vector
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd tsolver_get_beta(
    XPtr<TSolverRcpp> ptr, // NOLINT(performance-unnecessary-value-param)
    int step = -1
) {
    return ptr->solver_->getBeta(step);
}


//' @title Get Cp scores from TSolver
//'
//' @param ptr TSolverRcpp holder pointer
//' @return Numeric vector of Cp scores
//' @noRd
// [[Rcpp::export]]
NumericVector tsolver_get_cp(
    XPtr<TSolverRcpp> ptr // NOLINT(performance-unnecessary-value-param)
) {
    auto cp = ptr->solver_->getCp();
    return Rcpp::wrap(cp);
}


//' @title Get RSS from TSolver
//'
//' @param ptr TSolverRcpp holder pointer
//' @return Numeric vector of RSS scores
//' @noRd
// [[Rcpp::export]]
NumericVector tsolver_get_rss(
    XPtr<TSolverRcpp> ptr // NOLINT(performance-unnecessary-value-param)
) {
    return Rcpp::wrap(ptr->solver_->getRSS());
}


//' @title Get R2 from TSolver
//'
//' @param ptr TSolverRcpp holder pointer
//' @return Numeric vector of R2 scores
//' @noRd
// [[Rcpp::export]]
NumericVector tsolver_get_r2(
    XPtr<TSolverRcpp> ptr // NOLINT(performance-unnecessary-value-param)
) {
    return Rcpp::wrap(ptr->solver_->getR2());
}


//' @title Get Degrees of Freedom from TSolver
//'
//' @param ptr TSolverRcpp holder pointer
//' @return Numeric vector of DoF
//' @noRd
// [[Rcpp::export]]
NumericVector tsolver_get_dof(
    XPtr<TSolverRcpp> ptr // NOLINT(performance-unnecessary-value-param)
) {
    return Rcpp::wrap(ptr->solver_->getDoF());
}


//' @title Get Active Set from TSolver
//'
//' @param ptr TSolverRcpp holder pointer
//' @return Numeric vector of active variable indices
//' @noRd
// [[Rcpp::export]]
NumericVector tsolver_get_actives(
    XPtr<TSolverRcpp> ptr // NOLINT(performance-unnecessary-value-param)
) {
    return Rcpp::wrap(ptr->solver_->getActives());
}


//' @title Get Inactive Set from TSolver
//'
//' @param ptr TSolverRcpp holder pointer
//' @return Numeric vector of inactive variable indices
//' @noRd
// [[Rcpp::export]]
NumericVector tsolver_get_inactives(
    XPtr<TSolverRcpp> ptr // NOLINT(performance-unnecessary-value-param)
) {
    return Rcpp::wrap(ptr->solver_->getInactives());
}


//' @title Get Intercept from TSolver
//'
//' @param ptr TSolverRcpp holder pointer
//' @param step Step to retrieve the intercept for
//' @return Numeric intercept value
//' @noRd
// [[Rcpp::export]]
double tsolver_get_intercept(
    XPtr<TSolverRcpp> ptr, // NOLINT(performance-unnecessary-value-param)
    int step = -1
) {
    return ptr->solver_->getIntercept(step);
}


//' @noRd
// [[Rcpp::export]]
int tsolver_get_num_steps(
    XPtr<TSolverRcpp> ptr // NOLINT(performance-unnecessary-value-param)
) {
    return static_cast<int>(ptr->solver_->getNumSteps());
}


//' @noRd
// [[Rcpp::export]]
double tsolver_get_memory_bytes(
    XPtr<TSolverRcpp> ptr // NOLINT(performance-unnecessary-value-param)
) {
    // double: a byte count overflows R's 32-bit integers at 2 GiB
    return static_cast<double>(ptr->solver_->memoryBytes());
}


//' @noRd
// [[Rcpp::export]]
List tsolver_get_actions(
    XPtr<TSolverRcpp> ptr // NOLINT(performance-unnecessary-value-param)
) {
    const auto& actions = ptr->solver_->getActions();
    List result(actions.size());
    for (std::size_t i = 0; i < actions.size(); ++i)
        result[i] = Rcpp::wrap(actions[i]);
    return result;
}


// ==============================================================================
// Solver-Specific API Endpoints
// ==============================================================================

//' @noRd
// [[Rcpp::export]]
NumericVector tsolver_get_lambda(
    XPtr<TSolverRcpp> ptr // NOLINT(performance-unnecessary-value-param)
) {
    if (auto* p = dynamic_cast<TSolverLARS*>(ptr->solver_)) {
        return Rcpp::wrap(p->getLambda());
    }
    stop("tsolver_get_lambda: Solver type does not have a lambda path.");
}


//' @noRd
// [[Rcpp::export]]
int tsolver_get_num_removals(
    XPtr<TSolverRcpp> ptr // NOLINT(performance-unnecessary-value-param)
) {
    if (auto* p = dynamic_cast<TSolverStagewise*>(ptr->solver_)) {
        return static_cast<int>(p->getNumRemovals());
    }
    if (auto* p = dynamic_cast<TSolverENET*>(ptr->solver_)) {
        return static_cast<int>(p->getNumRemovals());
    }
    if (auto* p = dynamic_cast<TSolverLASSO*>(ptr->solver_)) {
        return static_cast<int>(p->getNumRemovals());
    }
    stop("tsolver_get_num_removals: Solver type does not track removals.");
}


//' @noRd
// [[Rcpp::export]]
double tsolver_get_cycling_ratio(
    XPtr<TSolverRcpp> ptr // NOLINT(performance-unnecessary-value-param)
) {
    if (auto* p = dynamic_cast<TSolverStagewise*>(ptr->solver_)) {
        return p->getCyclingRatio();
    }
    if (auto* p = dynamic_cast<TSolverENET*>(ptr->solver_)) {
        return p->getCyclingRatio();
    }
    if (auto* p = dynamic_cast<TSolverLASSO*>(ptr->solver_)) {
        return p->getCyclingRatio();
    }
    stop("tsolver_get_cycling_ratio: Solver type does not support getCyclingRatio().");
}


//' @noRd
// [[Rcpp::export]]
double tsolver_tienet_aug_get_lambda2(
    XPtr<TSolverRcpp> ptr // NOLINT(performance-unnecessary-value-param)
) {
    if (auto* p = dynamic_cast<TSolverIENETAug*>(ptr->solver_)) {
        return p->getLambda2();
    }
    stop("tsolver_tienet_aug_get_lambda2: Solver is not a TIENETAug_Solver.");
}


//' @noRd
// [[Rcpp::export]]
IntegerVector tsolver_tienet_aug_get_groups(
    XPtr<TSolverRcpp> ptr // NOLINT(performance-unnecessary-value-param)
) {
    if (auto* p = dynamic_cast<TSolverIENETAug*>(ptr->solver_)) {
        // Solver stores 0-based group ids; return 1-based (package convention).
        const Eigen::VectorXi& g = p->getGroups();
        IntegerVector out(g.size());
        for (Eigen::Index i = 0; i < g.size(); ++i) out[i] = g[i] + 1;
        return out;
    }
    stop("tsolver_tienet_aug_get_groups: Solver is not a TIENETAug_Solver.");
}


//' @noRd
// [[Rcpp::export]]
int tsolver_tienet_aug_get_num_groups(
    XPtr<TSolverRcpp> ptr // NOLINT(performance-unnecessary-value-param)
) {
    if (auto* p = dynamic_cast<TSolverIENETAug*>(ptr->solver_)) {
        return static_cast<int>(p->getNumGroups());
    }
    stop("tsolver_tienet_aug_get_num_groups: Solver is not a TIENETAug_Solver.");
}


//' @title Save TSolver to File
//'
//' @param ptr TSolverRcpp holder pointer
//' @param filename Output filename
//' @noRd
// [[Rcpp::export]]
void tsolver_save(
    XPtr<TSolverRcpp> ptr, // NOLINT(performance-unnecessary-value-param)
    const std::string& filename
) {
    std::ofstream os(filename, std::ios::binary);
    if (!os.is_open()) {
        stop("tsolver_save: Cannot open file for saving.");
    }
    cereal::PortableBinaryOutputArchive archive(os);

    TSolverBase* s = ptr->solver_;
    // Derived classes must be checked before their base classes.
    // LARS family
    if (auto* p = dynamic_cast<TSolverStagewise*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverENETAug*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverIENETAug*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverIENET*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverTIKNET*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverENET*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverLASSO*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverStepwise*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverLARS*>(s)) {
        archive(*p);
    // OMP family
    } else if (auto* p = dynamic_cast<TSolverACGP*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverGP*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverMP*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverOOLS*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverOMP*>(s)) {
        archive(*p);
    // Independent families
    } else if (auto* p = dynamic_cast<TSolverNCGMP*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverAFS*>(s)) {
        archive(*p);
    } else {
        stop("tsolver_save: Unsupported TSolver type.");
    }
}


//' @title Restore the caller's X and D and release the solver's claim
//'
//' @param ptr XPtr to TSolverRcpp
//'
//' @details Undoes the in-place centering/scaling a writer solver
//'   (normalize or intercept TRUE) applied to X and D at construction and
//'   releases its claim on the buffers, so that another writer may take
//'   them; y is never modified. No-op for a reader.
//' @noRd
// [[Rcpp::export]]
void tsolver_restore(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    ptr->solver_->restore(ptr->X_map_, ptr->D_map_, ptr->y_map_);
}

//' @title Load TSolver from File
//'
//' @param ptr TSolverRcpp holder pointer
//' @param filename Input filename
//' @noRd
// [[Rcpp::export]]
void tsolver_load(
    XPtr<TSolverRcpp> ptr, // NOLINT(performance-unnecessary-value-param)
    const std::string& filename
) {
    std::ifstream is(filename, std::ios::binary);
    if (!is.is_open()) {
        stop("tsolver_load: Cannot open file for loading.");
    }
    cereal::PortableBinaryInputArchive archive(is);

    TSolverBase* s = ptr->solver_;
    // Derived classes must be checked before their base classes.
    // LARS family
    if (auto* p = dynamic_cast<TSolverStagewise*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverENETAug*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverIENETAug*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverIENET*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverTIKNET*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverENET*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverLASSO*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverStepwise*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverLARS*>(s)) {
        archive(*p);
    // OMP family
    } else if (auto* p = dynamic_cast<TSolverACGP*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverGP*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverMP*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverOOLS*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverOMP*>(s)) {
        archive(*p);
    // Independent families
    } else if (auto* p = dynamic_cast<TSolverNCGMP*>(s)) {
        archive(*p);
    } else if (auto* p = dynamic_cast<TSolverAFS*>(s)) {
        archive(*p);
    } else {
        stop("tsolver_load: Unsupported TSolver type.");
    }
    // X_ and D_ are not serialised; reconnect to the holder's stable maps.
    ptr->solver_->reconnect(ptr->X_map_, ptr->D_map_);
}


// ==============================================================================
// Family getters of TIENET and of the sparse-K solver
// ==============================================================================

//' @title TIENET: group-ridge penalty used
//' @noRd
// [[Rcpp::export]]
double tsolver_tienet_get_lambda2(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    if (auto* p = dynamic_cast<TSolverIENET*>(ptr->solver_)) return p->getLambda2();
    stop("tsolver_tienet_get_lambda2: Solver is not a TIENET_Solver.");
}

//' @title TIENET: 1-based group ids, one per predictor
//' @noRd
// [[Rcpp::export]]
IntegerVector tsolver_tienet_get_groups(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    if (auto* p = dynamic_cast<TSolverIENET*>(ptr->solver_)) {
        const Eigen::VectorXi& g = p->getGroups();
        IntegerVector out(g.size());
        for (Eigen::Index i = 0; i < g.size(); ++i) out[i] = g[i] + 1;
        return out;
    }
    stop("tsolver_tienet_get_groups: Solver is not a TIENET_Solver.");
}

//' @title TIENET: number of groups
//' @noRd
// [[Rcpp::export]]
int tsolver_tienet_get_num_groups(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    if (auto* p = dynamic_cast<TSolverIENET*>(ptr->solver_)) return static_cast<int>(p->getNumGroups());
    stop("tsolver_tienet_get_num_groups: Solver is not a TIENET_Solver.");
}

//' @title TIKNET: quadratic penalty weight lambda_2
//' @noRd
// [[Rcpp::export]]
double tsolver_tik_get_lambda2(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    if (auto* p = dynamic_cast<TSolverTIKNET*>(ptr->solver_)) return p->getLambda2();
    stop("tsolver_tik_get_lambda2: Solver is not a TIKNET_Solver.");
}


// ==============================================================================
// The base methods shared with the Python binding, plus the family getters.
// ==============================================================================

//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd tsolver_get_residuals(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    return ptr->solver_->getResiduals();
}

//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd tsolver_get_beta_path(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    return ptr->solver_->getBetaPathSparse().dense();
}

//' @noRd
// [[Rcpp::export]]
std::vector<double> tsolver_get_cp_with_sigma(XPtr<TSolverRcpp> ptr,
    double sigma_hat_sq) { // NOLINT(performance-unnecessary-value-param)
    return ptr->solver_->getCp(sigma_hat_sq);
}

//' @noRd
// [[Rcpp::export]]
std::string tsolver_solver_type_to_string(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    return ptr->solver_->solverTypeToString();
}

//' @noRd
// NOLINTBEGIN(performance-unnecessary-value-param)
// [[Rcpp::export]]
std::vector<std::size_t> tsolver_get_dropped_indices(XPtr<TSolverRcpp> ptr) {
    return ptr->solver_->getDroppedIndices();
}
// NOLINTEND(performance-unnecessary-value-param)

//' @noRd
// NOLINTBEGIN(performance-unnecessary-value-param)
// [[Rcpp::export]]
std::vector<std::size_t> tsolver_get_active_predictor_indices(XPtr<TSolverRcpp> ptr) {
    return ptr->solver_->getActivePredictorIndices();
}
// NOLINTEND(performance-unnecessary-value-param)

//' @noRd
// NOLINTBEGIN(performance-unnecessary-value-param)
// [[Rcpp::export]]
std::vector<std::size_t> tsolver_get_active_dummy_indices(XPtr<TSolverRcpp> ptr) {
    return ptr->solver_->getActiveDummyIndices();
}
// NOLINTEND(performance-unnecessary-value-param)

//' @noRd
// [[Rcpp::export]]
int tsolver_get_num_actives(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    return static_cast<int>(ptr->solver_->getNumActives());
}

//' @noRd
// [[Rcpp::export]]
int tsolver_get_dummy_start_index(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    return static_cast<int>(ptr->solver_->getDummyStartIndex());
}

//' @noRd
// [[Rcpp::export]]
bool tsolver_is_connected(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    return ptr->solver_->isConnected();
}

//' @noRd
// [[Rcpp::export]]
void tsolver_set_tolerance(XPtr<TSolverRcpp> ptr, double tol) { // NOLINT(performance-unnecessary-value-param)
    ptr->solver_->setTolerance(tol);
}

//' @noRd
// [[Rcpp::export]]
void tsolver_set_tie_seed(XPtr<TSolverRcpp> ptr, double seed) { // NOLINT(performance-unnecessary-value-param)
    ptr->solver_->setTieSeed(static_cast<uint32_t>(seed));
}

//' @noRd
// [[Rcpp::export]]
void tsolver_set_num_threads(XPtr<TSolverRcpp> ptr, int n_threads) { // NOLINT(performance-unnecessary-value-param)
    ptr->solver_->setNumThreads(n_threads);
}

//' @noRd
// [[Rcpp::export]]
int tsolver_get_num_threads(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    return ptr->solver_->getNumThreads();
}

//' @noRd
// [[Rcpp::export]]
void tsolver_set_exchangeable_tie(XPtr<TSolverRcpp> ptr, double alpha,
    double floor) { // NOLINT(performance-unnecessary-value-param)
    ptr->solver_->setExchangeableTie(alpha, floor);
}

//' @noRd
// [[Rcpp::export]]
std::vector<std::string> tsolver_get_warnings(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    return ptr->solver_->getWarnings();
}

//' @noRd
// [[Rcpp::export]]
void tsolver_clear_warnings(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    ptr->solver_->clearWarnings();
}

// ---- family getters ----------------------------------------------------------

//' @noRd
// [[Rcpp::export]]
int tsolver_lars_get_kahan_refresh_interval(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    auto* s = dynamic_cast<TSolverLARS*>(ptr->solver_);
    if (s == nullptr) stop("get_kahan_refresh_interval: not a LARS-family solver.");
    return static_cast<int>(s->getKahanRefreshInterval());
}

//' @noRd
// [[Rcpp::export]]
void tsolver_lars_set_kahan_refresh_interval(XPtr<TSolverRcpp> ptr,
    int interval) { // NOLINT(performance-unnecessary-value-param)
    auto* s = dynamic_cast<TSolverLARS*>(ptr->solver_);
    if (s == nullptr) stop("set_kahan_refresh_interval: not a LARS-family solver.");
    s->setKahanRefreshInterval(static_cast<std::size_t>(interval));
}

//' @noRd
// [[Rcpp::export]]
double tsolver_tenet_get_alpha(XPtr<TSolverRcpp> ptr) { // NOLINT(performance-unnecessary-value-param)
    auto* s = dynamic_cast<TSolverENET*>(ptr->solver_);
    if (s == nullptr) stop("get_alpha: not a TENET solver.");
    return s->getAlpha();
}

//' @noRd
// NOLINTBEGIN(performance-unnecessary-value-param)
// [[Rcpp::export]]
Eigen::SparseMatrix<double> tsolver_tik_get_tikhonov_k(XPtr<TSolverRcpp> ptr) {
    if (auto* a = dynamic_cast<TSolverTIKNET*>(ptr->solver_)) return a->getTikhonovK();
    stop("get_tikhonov_k: not a TIKNET solver.");
}
// NOLINTEND(performance-unnecessary-value-param)


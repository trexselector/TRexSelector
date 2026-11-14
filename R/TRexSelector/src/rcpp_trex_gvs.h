// ==============================================================================
#ifndef RCPP_TREX_GVS_H
#define RCPP_TREX_GVS_H
// ==============================================================================

// std includes
#include <memory>

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// TRex includes
#include <trex_selector_methods/ad/trex_gvs/trex_gvs.hpp>
#include "rcpp_trex_selector.h"

// ==============================================================================

using namespace trex::trex_selector_methods::ad::trex_gvs;
using namespace trex::trex_selector_methods::ad::trex_core;

// ==============================================================================

/**
 * @brief The control of a GVS run: the base control goes into the nested
 *        field, and the solver follows the penalty. EN takes the solver of
 *        en_solver. IEN keeps an IEN-family solver of the base control and
 *        takes TIENET otherwise. The R level warns first on a mismatch.
 */
inline TRexGVSControlParameter merge_gvs_control(
    const TRexGVSControlParameter& gvs_control,
    const tc::TRexControlParameter& trex_control) {
    namespace gvs_sd = trex::trex_selector_methods::ad::utils::solver_dispatch;
    TRexGVSControlParameter gvs_ctrl = gvs_control;
    gvs_ctrl.trex_ctrl = trex_control;
    if (gvs_ctrl.gvs_type == GVSType::IEN) {
        const auto st = gvs_ctrl.trex_ctrl.solver_type;
        const bool ien_family =
            st == gvs_sd::SolverTypeForTRex::TIENET ||
            st == gvs_sd::SolverTypeForTRex::TIENET_AUG ||
            st == gvs_sd::SolverTypeForTRex::TCDIENET;
        if (!ien_family) {
            gvs_ctrl.trex_ctrl.solver_type =
                gvs_sd::SolverTypeForTRex::TIENET;
        }
    } else { // EN
        switch (gvs_ctrl.en_solver) {
            case ENSolverType::TENET_AUG:
                gvs_ctrl.trex_ctrl.solver_type =
                    gvs_sd::SolverTypeForTRex::TENET_AUG;
                break;
            case ENSolverType::TCDENET:
                gvs_ctrl.trex_ctrl.solver_type =
                    gvs_sd::SolverTypeForTRex::TCDENET;
                break;
            default:
                gvs_ctrl.trex_ctrl.solver_type =
                    gvs_sd::SolverTypeForTRex::TENET;
        }
    }
    return gvs_ctrl;
}


/**
 * @brief Rcpp wrapper class for TRexGVSSelector leveraging zero-copy references.
 *
 * Maps R objects (vectors, matrices, lists) to the corresponding C++ structures
 * for the Dependency-Aware Group Variable Selection T-Rex variant.
 */
class RTRexGVSSelector : public RTRexSelector {
public:
    RTRexGVSSelector(
        Eigen::Map<Eigen::MatrixXd> X,
        Eigen::Map<Eigen::VectorXd> y,
        double tFDR,
        const TRexGVSControlParameter& gvs_control,
        const tc::TRexControlParameter& trex_control,
        int seed,
        bool verbose
    ) {
        this->X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X);
        this->y_map_ = std::make_unique<Eigen::Map<Eigen::VectorXd>>(y);
        const TRexGVSControlParameter gvs_ctrl = merge_gvs_control(gvs_control, trex_control);
        this->selector_ = std::make_unique<TRexGVSSelector>(
            *(this->X_map_), *(this->y_map_), tFDR, gvs_ctrl, seed, verbose
        );
    }

    TRexGVSSelector* get() const { return static_cast<TRexGVSSelector*>(this->selector_.get()); }
};

#endif // RCPP_TREX_GVS_H

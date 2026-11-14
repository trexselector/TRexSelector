// ==============================================================================
#ifndef RCPP_TREX_GVS_H
#define RCPP_TREX_GVS_H
// ==============================================================================

// std includes
#include <memory>
#include <optional>

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
 * @brief The control of a GVS run. Without a base control the one nested in
 *        the GVS control applies, whose solver is TENET. A base control that
 *        is passed replaces the nested one as it is, and the C++ constructor
 *        derives or refuses the solver.
 */
inline TRexGVSControlParameter merge_gvs_control(
    const TRexGVSControlParameter& gvs_control,
    const std::optional<tc::TRexControlParameter>& trex_control) {
    TRexGVSControlParameter gvs_ctrl = gvs_control;
    if (trex_control) {
        gvs_ctrl.trex_ctrl = *trex_control;
    }
    return gvs_ctrl;
}


/**
 * @brief Rcpp wrapper class for TRexGVSSelector leveraging zero-copy references.
 *
 * Holds the maps of the design and the response and the selector that works
 * on them. The control is the one of the run, see merge_gvs_control().
 */
class RTRexGVSSelector : public RTRexSelector {
public:
    RTRexGVSSelector(
        Eigen::Map<Eigen::MatrixXd> X,
        Eigen::Map<Eigen::VectorXd> y,
        double tFDR,
        const TRexGVSControlParameter& gvs_ctrl,
        int seed,
        bool verbose
    ) {
        this->X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X);
        this->y_map_ = std::make_unique<Eigen::Map<Eigen::VectorXd>>(y);
        this->selector_ = std::make_unique<TRexGVSSelector>(
            *(this->X_map_), *(this->y_map_), tFDR, gvs_ctrl, seed, verbose
        );
    }

    TRexGVSSelector* get() const { return static_cast<TRexGVSSelector*>(this->selector_.get()); }
};

#endif // RCPP_TREX_GVS_H

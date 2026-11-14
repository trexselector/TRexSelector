// ==============================================================================
#ifndef RCPP_TREX_TIKHONOV_H
#define RCPP_TREX_TIKHONOV_H
// ==============================================================================

// std includes
#include <memory>
#include <optional>

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// TRex includes
#include <trex_selector_methods/ad/trex_tikhonov/trex_tikhonov.hpp>
#include "rcpp_trex_selector.h"

// ==============================================================================

using namespace trex::trex_selector_methods::ad::trex_tikhonov;
using namespace trex::trex_selector_methods::ad::trex_core;

// ==============================================================================

/**
 * @brief The control of a run. Without a base control the one nested in the
 *        Tikhonov control applies, whose solver is TIKNET. A base control
 *        that is passed replaces the nested one as it is, and the C++
 *        constructor refuses a solver outside TIKNET / TCDTIKNET.
 */
inline TRexTikhonovControlParameter merge_tikhonov_control(
    const TRexTikhonovControlParameter& tik_control,
    const std::optional<tc::TRexControlParameter>& trex_control) {
    TRexTikhonovControlParameter tik_ctrl = tik_control;
    if (trex_control) {
        tik_ctrl.trex_ctrl = *trex_control;
    }
    return tik_ctrl;
}


/**
 * @brief Rcpp wrapper class for TRexTikhonovSelector (general-Tikhonov T-Rex).
 *
 * Holds the maps of the design and the response and the selector that works
 * on them. The control is the one of the run, see merge_tikhonov_control().
 */
class RTRexTikhonovSelector : public RTRexSelector {
public:
    RTRexTikhonovSelector(
        Eigen::Map<Eigen::MatrixXd> X,
        Eigen::Map<Eigen::VectorXd> y,
        double tFDR,
        const TRexTikhonovControlParameter& tik_ctrl,
        int seed,
        bool verbose
    ) {
        this->X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X);
        this->y_map_ = std::make_unique<Eigen::Map<Eigen::VectorXd>>(y);
        this->selector_ = std::make_unique<TRexTikhonovSelector>(
            *(this->X_map_), *(this->y_map_), tFDR, tik_ctrl, seed, verbose
        );
    }

    TRexTikhonovSelector* get() const {
        return static_cast<TRexTikhonovSelector*>(this->selector_.get());
    }
};

#endif // RCPP_TREX_TIKHONOV_H

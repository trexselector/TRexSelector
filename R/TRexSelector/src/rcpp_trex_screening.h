#ifndef RCPP_TREX_SCREENING_H
#define RCPP_TREX_SCREENING_H

#include <memory>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <trex_selector_methods/ad/trex_screening/trex_screening.hpp>
#include "rcpp_trex_selector.h"

using namespace trex::trex_selector_methods::ad::trex_screening;
using namespace trex::trex_selector_methods::ad::trex_core;

/**
 * @brief Rcpp wrapper class for TRexScreeningSelector leveraging zero-copy references.
 * 
 * Exposes methods to parse R parameters and retrieve screening-centric results 
 * like bootstrap confidence levels and estimated FDR.
 */
class RTRexScreeningSelector : public RTRexSelector {
public:
    RTRexScreeningSelector(
        Eigen::Map<Eigen::MatrixXd> X,
        Eigen::Map<Eigen::VectorXd> y,
        const TRexScreeningControl& screen_control,
        const tc::TRexControl& trex_control,
        int seed,
        bool verbose
    ) {
        this->X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X);
        this->y_map_ = std::make_unique<Eigen::Map<Eigen::VectorXd>>(y);

        // The T-Rex+Screening control nests the base T-Rex control; merge the
        // separately-supplied trex_control into that nested field.
        TRexScreeningControl screen_ctrl = screen_control;
        screen_ctrl.trex_ctrl = trex_control;

        this->selector_ = std::make_unique<TRexScreeningSelector>(
            *(this->X_map_), *(this->y_map_), screen_ctrl, seed, verbose
        );
    }
    
    TRexScreeningSelector* get() const { return static_cast<TRexScreeningSelector*>(this->selector_.get()); }
};

#endif // RCPP_TREX_SCREENING_H

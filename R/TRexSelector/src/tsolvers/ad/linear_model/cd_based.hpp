// ====================================================================
// cd_based.hpp
// ====================================================================
#ifndef TSOLVERS_LINEAR_MODEL_CD_BASED_HPP
#define TSOLVERS_LINEAR_MODEL_CD_BASED_HPP
// ====================================================================
/**
 * @file cd_based.hpp
 *
 * @brief Umbrella header for CCD-based terminating solvers.
 *
 * @deprecated The coordinate descent family will be removed in the next
 *  minor release. Use the homotopy solvers of lars_based.hpp: TLASSO, TENET,
 *  TIENET and TIKNET.
 */
// ====================================================================

#include <tsolvers/ad/linear_model/cd_based/tcdlasso_solver.hpp>
#include <tsolvers/ad/linear_model/cd_based/tcdenet_solver.hpp>
#include <tsolvers/ad/linear_model/cd_based/tcdienet_solver.hpp>
#include <tsolvers/ad/linear_model/cd_based/tcdtiknet_solver.hpp>

// ====================================================================
#endif /* End of TSOLVERS_LINEAR_MODEL_CD_BASED_HPP */

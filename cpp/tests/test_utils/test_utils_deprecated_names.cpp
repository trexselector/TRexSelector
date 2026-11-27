// ========================================================================================
/**
 * @file test_utils_deprecated_names.cpp
 *
 * @brief The old type names are aliases of the new ones, and the old include paths of the T-solvers
 *        forward to the new files.
 */
// ========================================================================================

// google test includes
#include <gtest/gtest.h>

// std includes
#include <type_traits>

// project includes
#include <ml_methods/model_selection/enet_cv_ccd.hpp>
#include <ml_methods/model_selection/lambda2_cv_diagnostics.hpp>
#include <ml_methods/model_selection/ridge_cv_svd.hpp>
#include <ml_methods/ridge_regression/ridge.hpp>
#include <ml_methods/svd/svd.hpp>
#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <trex_selector_methods/ad/trex_da/trex_da.hpp>
#include <trex_selector_methods/ad/trex_gvs/trex_gvs.hpp>
#include <trex_selector_methods/ad/trex_screening/trex_biobank_screening.hpp>
#include <trex_selector_methods/ad/trex_screening/trex_screening.hpp>
#include <trex_selector_methods/ad/trex_spca/trex_spca.hpp>
#include <trex_selector_methods/ad/trex_tikhonov/trex_tikhonov.hpp>
#include <trex_selector_methods/ad/trex_utils/trex_solver_dispatch.hpp>
#include <tsolvers/ad/linear_model/afs_based/tsolver_afs.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_enet_aug.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_enet.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_ienet_aug.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_ienet.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_tiknet.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_lars.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_lasso.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_stagewise.hpp>
#include <tsolvers/ad/linear_model/lars_based/tsolver_stepwise.hpp>
#include <tsolvers/ad/linear_model/omp_based/tsolver_acgp.hpp>
#include <tsolvers/ad/linear_model/omp_based/tsolver_gp.hpp>
#include <tsolvers/ad/linear_model/omp_based/tsolver_mp.hpp>
#include <tsolvers/ad/linear_model/omp_based/tsolver_ncgmp.hpp>
#include <tsolvers/ad/linear_model/omp_based/tsolver_omp.hpp>
#include <tsolvers/ad/linear_model/omp_based/tsolver_ools.hpp>
#include <tsolvers/ad/tsolver_options.hpp>
#include <tsolvers/ad/tsolver_base.hpp>
#include <utils/datageneration/utils_dummygen.hpp>

// The old include paths of the T-solvers forward to the new files. Each one prints a note on purpose.
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-W#pragma-messages"
#endif
#include <tsolvers/ad/linear_model/afs_based/tafs_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tenet_aug_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tenet_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tienet_aug_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tienet_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tiknet_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tlars_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tlasso_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tstagewise_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tstepwise_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tacgp_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tgp_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tmp_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tncgmp_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tomp_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tools_solver.hpp>
#include <tsolvers/ad/solver_options.hpp>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

// ========================================================================================

// The old names warn on purpose. This file is the one place that uses them.
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#elif defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#endif

static_assert(std::is_same_v<trex::ml_methods::model_selection::enet_cv_ccd,
                             trex::ml_methods::model_selection::ElasticNetCV>);
static_assert(std::is_same_v<trex::ml_methods::model_selection::enet_gaussian,
                             trex::ml_methods::model_selection::ElasticNet>);
static_assert(std::is_same_v<trex::ml_methods::model_selection::lambda2_cv_diagnostics,
                             trex::ml_methods::model_selection::Lambda2CVDiagnostics>);
static_assert(std::is_same_v<trex::ml_methods::model_selection::ridge_cv_svd,
                             trex::ml_methods::model_selection::RidgeCV>);
static_assert(std::is_same_v<trex::ml_methods::ridge::RidgeSolver,
                             trex::ml_methods::ridge::Ridge>);
static_assert(std::is_same_v<trex::ml_methods::svd::SVDSolver,
                             trex::ml_methods::svd::SVD>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::trex_biobank_screening::BiobankScreenTRex,
                             trex::trex_selector_methods::ad::trex_biobank_screening::TRexBiobankScreeningSelector>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::trex_biobank_screening::BiobankScreenTRexControl,
                             trex::trex_selector_methods::ad::trex_biobank_screening::TRexBiobankScreeningControl>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::trex_biobank_screening::BiobankScreenTRexResult,
                             trex::trex_selector_methods::ad::trex_biobank_screening::TRexBiobankScreeningResult>);
static_assert(std::is_same_v<
    trex::trex_selector_methods::ad::trex_biobank_screening::TRexBiobankScreeningSelector::BiobankScreenTRexStatistics,
    trex::trex_selector_methods::ad::trex_biobank_screening::TRexBiobankScreeningStatistics>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::trex_core::TRexControlParameter,
                             trex::trex_selector_methods::ad::trex_core::TRexControl>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::trex_core::TRexSelector::SelectionResult,
                             trex::trex_selector_methods::ad::trex_core::TRexResult>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::trex_da::TRexDAControlParameter,
                             trex::trex_selector_methods::ad::trex_da::TRexDAControl>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::trex_da::TRexDASelector::DASelectionResult,
                             trex::trex_selector_methods::ad::trex_da::TRexDAResult>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::trex_gvs::TRexGVSControlParameter,
                             trex::trex_selector_methods::ad::trex_gvs::TRexGVSControl>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::trex_gvs::TRexGVSSelector::GVSSelectionResult,
                             trex::trex_selector_methods::ad::trex_gvs::TRexGVSResult>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::trex_screening::ScreenTRexControlParameter,
                             trex::trex_selector_methods::ad::trex_screening::TRexScreeningControl>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::trex_screening::ScreenTRexSelectionResult,
                             trex::trex_selector_methods::ad::trex_screening::TRexScreeningResult>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::trex_screening::ScreenTRexSelector,
                             trex::trex_selector_methods::ad::trex_screening::TRexScreeningSelector>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::trex_spca::TRexSPCA,
                             trex::trex_selector_methods::ad::trex_spca::TRexSPCASelector>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::trex_spca::TRexSPCAControlParameter,
                             trex::trex_selector_methods::ad::trex_spca::TRexSPCAControl>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::trex_tikhonov::TRexTikhonovControlParameter,
                             trex::trex_selector_methods::ad::trex_tikhonov::TRexTikhonovControl>);
static_assert(std::is_same_v<trex::trex_selector_methods::ad::utils::solver_dispatch::SolverHyperparameters,
                             trex::trex_selector_methods::ad::utils::solver_dispatch::TSolverParams>);
static_assert(std::is_same_v<trex::tsolvers::ad::SolverOptions,
                             trex::tsolvers::ad::TSolverOptions>);
static_assert(std::is_same_v<trex::tsolvers::ad::TSolver_Base,
                             trex::tsolvers::ad::TSolverBase>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::afs_based::TAFS_Solver,
                             trex::tsolvers::ad::linear_model::afs_based::TSolverAFS>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::lars_based::TENETAug_Solver,
                             trex::tsolvers::ad::linear_model::lars_based::TSolverENETAug>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::lars_based::TENET_Solver,
                             trex::tsolvers::ad::linear_model::lars_based::TSolverENET>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::lars_based::TIENETAug_Solver,
                             trex::tsolvers::ad::linear_model::lars_based::TSolverIENETAug>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::lars_based::TIENET_Solver,
                             trex::tsolvers::ad::linear_model::lars_based::TSolverIENET>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::lars_based::TIKNET_Solver,
                             trex::tsolvers::ad::linear_model::lars_based::TSolverTIKNET>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::lars_based::TLARS_Solver,
                             trex::tsolvers::ad::linear_model::lars_based::TSolverLARS>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::lars_based::TLASSO_Solver,
                             trex::tsolvers::ad::linear_model::lars_based::TSolverLASSO>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::lars_based::TSTAGEWISE_Solver,
                             trex::tsolvers::ad::linear_model::lars_based::TSolverStagewise>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::lars_based::TSTEPWISE_Solver,
                             trex::tsolvers::ad::linear_model::lars_based::TSolverStepwise>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::omp_based::TACGP_Solver,
                             trex::tsolvers::ad::linear_model::omp_based::TSolverACGP>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::omp_based::TGP_Solver,
                             trex::tsolvers::ad::linear_model::omp_based::TSolverGP>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::omp_based::TMP_Solver,
                             trex::tsolvers::ad::linear_model::omp_based::TSolverMP>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::omp_based::TNCGMP_Solver,
                             trex::tsolvers::ad::linear_model::omp_based::TSolverNCGMP>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::omp_based::TOMP_Solver,
                             trex::tsolvers::ad::linear_model::omp_based::TSolverOMP>);
static_assert(std::is_same_v<trex::tsolvers::ad::linear_model::omp_based::TOOLS_Solver,
                             trex::tsolvers::ad::linear_model::omp_based::TSolverOOLS>);
static_assert(std::is_same_v<trex::utils::datageneration::dummygen::Distribution,
                             trex::utils::datageneration::dummygen::DummyDistribution>);

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#elif defined(_MSC_VER)
#pragma warning(pop)
#endif

// ========================================================================================

// Embed into test namespace
namespace trex::test::utils {

/** @brief The checks of this file run at compile time. */
TEST(DeprecatedNamesTest, OldTypeNamesAreAliases) {
    SUCCEED();
}

} // namespace trex::test::utils

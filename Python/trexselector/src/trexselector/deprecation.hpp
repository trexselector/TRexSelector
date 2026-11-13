// =====================================================================================
// deprecation.hpp - Deprecation warnings of the Python bindings
// =====================================================================================
#ifndef TREX_PYTHON_DEPRECATION_HPP
#define TREX_PYTHON_DEPRECATION_HPP
// =====================================================================================
/**
 * @file deprecation.hpp
 *
 * @brief FutureWarning of the deprecated coordinate-descent solver family.
 *
 * @details
 *  The warning is raised by trexselector._deprecation.warn_deprecated(), which
 *  attributes it to the first caller outside the package. The bindings call
 *  the helpers below in their constructors, with the GIL held and outside the
 *  experiment loops: once where a solver class of the family is constructed
 *  and once where a solver type of the family reaches a selector.
 */
// =====================================================================================

// pybind11 includes
#include <pybind11/pybind11.h>

// std includes
#include <string>
#include <type_traits>

// trex includes
#include <tsolvers/ad/linear_model/cd_based/tcdlasso_solver.hpp>
#include <tsolvers/ad/linear_model/cd_based/tcdenet_solver.hpp>
#include <tsolvers/ad/linear_model/cd_based/tcdienet_solver.hpp>
#include <tsolvers/ad/linear_model/cd_based/tcdtiknet_solver.hpp>
#include <trex_selector_methods/ad/trex_gvs/trex_gvs.hpp>
#include <trex_selector_methods/ad/trex_utils/trex_solver_dispatch.hpp>

// =====================================================================================

namespace py = pybind11;

// =====================================================================================

/**
 * @brief Raise the FutureWarning that a name is deprecated.
 *
 * @param old_name Deprecated class or enumerator, as Python spells it.
 * @param new_name Its replacement.
 *
 * @throws py::error_already_set if the warning filter turns the warning into
 *         an error.
 */
inline void warn_deprecated(const std::string& old_name, const std::string& new_name) {
    py::module_::import("trexselector._deprecation").attr("warn_deprecated")(old_name, new_name);
}

/**
 * @brief Warn when a solver type of the coordinate-descent family reaches a
 *        selector.
 *
 * @param solver_type Solver type of the base control.
 */
inline void warn_if_cd_solver_type(
    trex::trex_selector_methods::ad::utils::solver_dispatch::SolverTypeForTRex solver_type
) {
    using SolverType = trex::trex_selector_methods::ad::utils::solver_dispatch::SolverTypeForTRex;
    switch (solver_type) {
        case SolverType::TCDLASSO:
            warn_deprecated("SolverTypeForTRex.TCDLASSO", "SolverTypeForTRex.TLASSO");
            break;
        case SolverType::TCDENET:
            warn_deprecated("SolverTypeForTRex.TCDENET", "SolverTypeForTRex.TENET");
            break;
        case SolverType::TCDIENET:
            warn_deprecated("SolverTypeForTRex.TCDIENET", "SolverTypeForTRex.TIENET");
            break;
        case SolverType::TCDTIKNET:
            warn_deprecated("SolverTypeForTRex.TCDTIKNET", "SolverTypeForTRex.TIKNET");
            break;
        default:
            break;
    }
}

/**
 * @brief Warn when the elastic-net solver choice of a selector is TCDENET.
 *
 * @param en_solver Elastic-net solver of the grouped or the SPCA control.
 */
inline void warn_if_cd_en_solver(trex::trex_selector_methods::ad::trex_gvs::ENSolverType en_solver) {
    using trex::trex_selector_methods::ad::trex_gvs::ENSolverType;
    if (en_solver == ENSolverType::TCDENET) {
        warn_deprecated("ENSolverType.TCDENET", "ENSolverType.TENET");
    }
}

/**
 * @brief Warn when a solver class of the coordinate-descent family is
 *        constructed.
 *
 * @tparam Solver The T-solver class behind the Python class.
 */
template <typename Solver>
void warn_if_cd_solver_class() {
    namespace cd = trex::tsolvers::ad::linear_model::cd_based;
    if constexpr (std::is_same_v<Solver, cd::TCDLASSO_Solver>) {
        warn_deprecated("TCDLASSO_Solver", "TLASSO_Solver");
    } else if constexpr (std::is_same_v<Solver, cd::TCDENET_Solver>) {
        warn_deprecated("TCDENET_Solver", "TENET_Solver");
    } else if constexpr (std::is_same_v<Solver, cd::TCDIENET_Solver>) {
        warn_deprecated("TCDIENET_Solver", "TIENET_Solver");
    } else if constexpr (std::is_same_v<Solver, cd::TCDTIKNET_Solver>) {
        warn_deprecated("TCDTIKNET_Solver", "TIKNET_Solver");
    }
}

// =====================================================================================
#endif /* End of TREX_PYTHON_DEPRECATION_HPP */

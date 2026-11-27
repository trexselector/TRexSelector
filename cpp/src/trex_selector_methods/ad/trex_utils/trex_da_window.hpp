// ===================================================================================
// trex_da_window.hpp
// ===================================================================================
#ifndef TREX_SELECTOR_METHODS_UTILS_TREX_DA_WINDOW_HPP
#define TREX_SELECTOR_METHODS_UTILS_TREX_DA_WINDOW_HPP
// ===================================================================================
/**
 * @file trex_da_window.hpp
 *
 * @brief Sliding window of the dependency-aware T-Rex under an AR(1) model.
 *
 * @details
 *  Single home of the window half-width of `TRexDASelector` (AR1 method) and
 *  `TRexScreeningSelector` (TREX_DA_AR1), so both selectors group the same
 *  neighbors for the same inputs, and of the rule that tells whether a
 *  supplied value is a correlation coefficient.
 *
 *  Under an AR(1) model two variables at lag k have the correlation |rho|^k.
 *  The window of a variable holds the neighbors whose correlation reaches the
 *  threshold rho_thr, which gives the half-width
 *
 *      kappa = ceil(log(rho_thr) / log(|rho|)),   rho_thr in [0, 1).
 */
// ===================================================================================

// std includes
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

// ===================================================================================

// Embedded into namespace trex::trex_selector_methods::ad::utils::da_window
namespace trex::trex_selector_methods::ad::utils::da_window {

/**
 * @brief Whether a supplied value can serve as a correlation coefficient.
 *
 * The dependency-aware selectors replace a supplied coefficient outside
 * [-1, 1] by their estimate and warn about it. A NaN is not usable either.
 *
 * @param cor_coef Supplied coefficient.
 *
 * @return True if cor_coef lies in [-1, 1].
 */
inline bool isCorrelationCoefficient(double cor_coef) noexcept {
    return cor_coef >= -1.0 && cor_coef <= 1.0;
}

/**
 * @brief Half-width of the AR(1) sliding window.
 *
 * @param cor_coef AR(1) coefficient. Its magnitude is used.
 * @param rho_thr  Correlation threshold in [0, 1).
 * @param p        Number of variables, the largest half-width returned.
 *
 * @return ceil(log(rho_thr) / log(|cor_coef|)), at least 1 and at most p.
 *         A zero coefficient gives 1 (the limit of the formula). A zero
 *         threshold or |cor_coef| >= 1 gives p (every lag reaches the
 *         threshold).
 *
 * @throws std::invalid_argument if rho_thr is outside [0, 1) or cor_coef is
 *         not a number.
 */
inline std::size_t ar1HalfWidth(double cor_coef, double rho_thr, std::size_t p) {
    if (!(rho_thr >= 0.0 && rho_thr < 1.0)) {
        throw std::invalid_argument(
            "AR(1) window: the threshold must be in [0, 1). Got: " +
             std::to_string(rho_thr));
    }

    const double rho_abs = std::abs(cor_coef);
    if (!(rho_abs >= 0.0)) {
        throw std::invalid_argument("AR(1) window: the coefficient is not a number.");
    }

    // Every lag reaches the threshold: full window.
    if (rho_thr == 0.0 || rho_abs >= 1.0) { return p; }

    // Limit of the formula as |cor_coef| -> 0+: the minimal window.
    if (rho_abs == 0.0) { return 1; }

    // Both logarithms are negative and finite here, so the ratio is positive.
    const double half_width = std::ceil(std::log(rho_thr) / std::log(rho_abs));
    return half_width < static_cast<double>(p) ? static_cast<std::size_t>(half_width) : p;
}

// ===================================================================================
} /* End of namespace trex::trex_selector_methods::ad::utils::da_window */
// ===================================================================================
#endif /* End of TREX_SELECTOR_METHODS_UTILS_TREX_DA_WINDOW_HPP */

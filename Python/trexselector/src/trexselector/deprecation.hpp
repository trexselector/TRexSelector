// =====================================================================================
// deprecation.hpp - Deprecation warnings of the Python bindings
// =====================================================================================
#ifndef TREX_PYTHON_DEPRECATION_HPP
#define TREX_PYTHON_DEPRECATION_HPP
// =====================================================================================
/**
 * @file deprecation.hpp
 *
 * @brief FutureWarning of a deprecated name of the package.
 *
 * @details
 *  The warning is raised by trexselector._deprecation.warn_deprecated(), which
 *  attributes it to the first caller outside the package. A binding calls
 *  the helper below with the GIL held and outside the experiment loops.
 */
// =====================================================================================

// pybind11 includes
#include <pybind11/pybind11.h>

// std includes
#include <string>

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

// =====================================================================================
#endif /* End of TREX_PYTHON_DEPRECATION_HPP */

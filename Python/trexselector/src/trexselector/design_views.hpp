// =====================================================================================
// design_views.hpp - Zero-copy views of a NumPy design that is not writable
// =====================================================================================
#ifndef TREX_PYTHON_DESIGN_VIEWS_HPP
#define TREX_PYTHON_DESIGN_VIEWS_HPP
// =====================================================================================
/**
 * @file design_views.hpp
 *
 * @brief Zero-copy views of a NumPy design that is not writable.
 *
 * @details
 *  A writable array reaches the bindings as a mutable Eigen::Ref. An array
 *  that is not writable (the view of a ReadOnly MemoryMappedMatrix, a
 *  numpy.memmap opened with mode "r", an array with the write flag cleared)
 *  is accepted where nothing is written into it and refused everywhere else.
 *  The library decides what writes: dn::writesDesign() for the selectors and
 *  SolverOptions::writesDesign() for the T-solvers.
 */
// =====================================================================================

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

// std includes
#include <string>

// =====================================================================================

namespace py = pybind11;

// =====================================================================================

/** @brief A float64, Fortran-contiguous array taken as it is. Bound with
 *  py::arg().noconvert(), so pybind11 never copies it. */
using DesignArray = py::array_t<double, py::array::f_style>;

/** @brief When a T-Rex selector writes into its design. */
inline const std::string kSelectorDesignRule =
    "The selector standardizes the design in place and restores it under ScalingMode.L2 and "
    "ScalingMode.ZSCORE, which needs a writable array. Under ScalingMode.NONE it writes nothing "
    "and accepts an array that is not writable.";

/** @brief When a T-solver writes into X and D. */
inline const std::string kSolverDesignRule =
    "A solver with normalize=True or intercept=True standardizes X and D in place, which needs "
    "writable arrays. With normalize=False and intercept=False it only reads them and accepts "
    "arrays that are not writable.";

/** @brief Docstring of the selector constructors that take an array as it is. */
inline const std::string kReaderSelectorDoc =
    "Construct on a design that is not writable, for example the to_numpy() view of a "
    "MemoryMappedMatrix opened with AccessMode.ReadOnly. X must be a Fortran-contiguous float64 "
    "array and is accessed zero-copy. " + kSelectorDesignRule + " Otherwise ValueError is raised.";

/**
 * @brief Eigen::Map of a design array for code that states whether it writes
 *        into it. The array lives in memory or is the view of a memory-mapped
 *        file.
 *
 * @param X      Design array, typically not writable.
 * @param writes Whether the caller writes into the design.
 * @param name   Name of the argument, for the message.
 * @param rule   Sentence that says when the caller writes.
 *
 * @return Eigen::Map of the array's buffer (a one-dimensional array is one
 *         column).
 *
 * @throws py::value_error if the caller writes and the array is not writable,
 *         or if the array has more than two dimensions.
 */
inline Eigen::Map<Eigen::MatrixXd> reader_design_map(
    const DesignArray& X,
    bool writes,
    const char* name,
    const std::string& rule
) {
    if (X.ndim() != 1 && X.ndim() != 2) {
        throw py::value_error(std::string(name) + " must be a two-dimensional array.");
    }
    if (writes && !X.writeable()) {
        throw py::value_error(
            std::string(name) + " is not writable (for example the view of a ReadOnly "
            "MemoryMappedMatrix). " + rule);
    }
    // The buffer of an array that is not writable is const: the caller reads only.
    return Eigen::Map<Eigen::MatrixXd>(
        const_cast<double*>(X.data()),
        static_cast<Eigen::Index>(X.shape(0)),
        X.ndim() == 2 ? static_cast<Eigen::Index>(X.shape(1)) : Eigen::Index{1});
}

// =====================================================================================
#endif /* End of TREX_PYTHON_DESIGN_VIEWS_HPP */

// =============================================================================
// memory_mapped_matrix.hpp
// =============================================================================
#ifndef UTILS_MEMORY_MAPPED_MATRIX_HPP
#define UTILS_MEMORY_MAPPED_MATRIX_HPP
// =============================================================================
/**
 * @file memory_mapped_matrix.hpp
 *
 * @brief Declaration of a memory-mapped matrix using Eigen and Boost.
 */
// =============================================================================

// std includes
#include <string>
#include <memory>

// boost includes
#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

// Eigen includes
#include <Eigen/Dense>

// ================================================================================

// Embedded into namespace trex::utils::memmap
namespace trex::utils::memmap {

// ================================================================================

/**
 * @brief Access mode for memory-mapped files
 *
 * - ReadOnly: Open the file in read-only mode
 * - ReadWrite: Open the file in read-write mode
 */
enum class AccessMode {
    READ_ONLY,
    READ_WRITE,

    // old names, kept for one release
    ReadOnly [[deprecated("Use READ_ONLY")]] = READ_ONLY,
    ReadWrite [[deprecated("Use READ_WRITE")]] = READ_WRITE
};


/**
 * @brief Memory-mapped matrix using Boost.Iostreams
 *
 * @details This class provides a memory-mapped matrix interface using Boost.Iostreams.
 *
 *  A file serves one writer at a time. Several selectors or T-solvers that
 *  standardize the file in place must not run on one file at the same time,
 *  neither through a second object in one process nor in several processes.
 *  Each would change the data the other works on, and the results come out
 *  wrong. To use one design in several processes at once, pre-normalize it,
 *  open the file ReadOnly and run every selector under ScalingMode::NONE.
 *
 * @tparam Scalar The scalar type of the matrix elements
 *         (e.g., float, double, std::complex<double>).
 *
 *  References:
 *  -----------
 *   F. Scheidt, J. Machkour, and M. Muma, "Solving FDR-Controlled Sparse Regression Problems with Five
 *     Million Variables on a Laptop," in Proc. IEEE 9th Int. Workshop Comput. Adv. Multi-Sensor Adapt.
 *     Process. (CAMSAP), 2023, pp. 116-120.
 */
template<typename Scalar>
class MemoryMappedMatrix {
public:
    using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    using MapType = Eigen::Map<MatrixType>;
    using ConstMapType = Eigen::Map<const MatrixType>;

    // =================================================
    // Constructors
    // =================================================
    /**
     * @brief Constructor to create or open a memory-mapped matrix.
     *
     * @param filename The path to the memory-mapped file.
     * @param rows Number of rows in the matrix.
     * @param cols Number of columns in the matrix.
     * @param mode Access mode (ReadOnly or ReadWrite).
     *
     */
    MemoryMappedMatrix(
        const std::string& filename,
        std::size_t rows,
        std::size_t cols,
        AccessMode mode = AccessMode::READ_WRITE
    );

    /** @brief Destructor: closes the memory-mapped file */
    ~MemoryMappedMatrix();

    // Prevent copying
    /** @brief Delete copy constructor */
    MemoryMappedMatrix(const MemoryMappedMatrix&) = delete;

    /** @brief Delete copy assignment operator */
    MemoryMappedMatrix& operator=(const MemoryMappedMatrix&) = delete;


    // =================================================
    // Public Methods
    // =================================================

    /** @brief Get Eigen::Map view (column-major) for writing.
     *  Refused for a file opened ReadOnly. */
    MapType getMap();

    /** @brief Get const Eigen::Map view (column-major) */
    ConstMapType getMap() const;

    /**
     * @brief Get an Eigen::Map view for code that only reads the matrix.
     *
     * @details The view has the writable Eigen::Map type that the selectors
     *  and solvers take, and it is handed out for a file opened ReadOnly as
     *  well. A write through it into such a file ends the process, so only
     *  code that writes nothing into the matrix may use it: a selector
     *  under ScalingMode::NONE, a T-solver with normalize and intercept off.
     */
    MapType getReaderMap();

    /**
     * @brief Get the view for code that states whether it writes the matrix.
     *
     * @param writes True if the code writes into the matrix.
     *
     * @return getMap() for a writer, getReaderMap() otherwise.
     *
     * @throws std::runtime_error for a writer on a file opened ReadOnly.
     */
    MapType getMapFor(bool writes);

    /** @brief Access mode the file was opened with */
    AccessMode mode() const { return mode_; }

    // Raw Pointers
    /** @brief Get raw pointer to data */
    Scalar* data();

    /** @brief Get const raw pointer to data */
    const Scalar* data() const;

    // Accessors
    /** @brief Get number of rows */
    std::size_t rows() const { return rows_; }

    /** @brief Get number of columns */
    std::size_t cols() const { return cols_; }

    /** @brief Get total size (rows * columns). Cannot wrap: the constructor
     *  already verified that rows * cols * sizeof(Scalar) fits std::size_t. */
    std::size_t size() const { return rows_ * cols_; }

    // Element Access (Eigen-like interface)
    /** @brief Read/write element access (ReadWrite mode only) — Eigen-like interface */
    Scalar& operator()(std::size_t row, std::size_t col);

    /** @brief Read-only element access — Eigen-like interface */
    Scalar operator()(std::size_t row, std::size_t col) const;

    /**
     * @brief Flush dirty pages to the backing file and drop their residency.
     *
     * @details Written pages of a read-write file mapping stay dirty in the
     * process footprint until the pager cleans them; a large mapping that is
     * written once and then only read keeps gigabytes resident (macOS
     * phys_footprint) and can drive the process into a jetsam kill. This
     * synchronously writes all dirty pages to the backing file and then
     * advises the kernel to evict the region. The mapping stays valid — the
     * next access refaults the data from disk.
     */
    void releaseResidency();


private:
    // =================================================
    // Private Members
    // =================================================

    /** @brief Path to the memory-mapped file */
    std::string filename_;

    /** @brief Number of rows */
    std::size_t rows_;

    /** @brief Number of columns */
    std::size_t cols_;

    /** @brief Access mode */
    AccessMode mode_;

    /** @brief Boost file mapping object */
    std::unique_ptr<boost::interprocess::file_mapping> file_mapping_;

    /** @brief Boost mapped region object */
    std::unique_ptr<boost::interprocess::mapped_region> mapped_region_;


    // =================================================
    // Private Methods
    // =================================================

    /** @brief Initialize the memory-mapped matrix */
    void initialize();

    /** @brief Apply access hints to the memory-mapped file */
    void applyAccessHints();

};

// ================================================================================

} /* End of namespace trex::utils::memmap */

#endif /* UTILS_MEMORY_MAPPED_MATRIX_HPP */


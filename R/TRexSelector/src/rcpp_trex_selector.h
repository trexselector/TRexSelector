// ==============================================================================
// rcpp_trex_selector.h - Safe bindings wrappers for TRex objects
// ==============================================================================
#ifndef RCPP_TREX_SELECTOR_H
#define RCPP_TREX_SELECTOR_H
// ==============================================================================
/**
 * @file rcpp_trex_selector.h
 *
 * @brief Bindings wrappers for TRex objects.
 */
// ==============================================================================

// std includes
#include <algorithm>
#include <cctype>
#include <memory>
#include <string>

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// TRex includes
#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <ml_methods/clustering/hierarchical/agglomerative/agglomerative_types.hpp>
#include <utils/memmap/memory_mapped_matrix.hpp>

// ==============================================================================

// Namespace usage for convenience
using namespace trex::trex_selector_methods::ad::trex_core;

// ==============================================================================

/**
 * @brief Eigen::Map of a memory-mapped design for a selector.
 *
 * Under ScalingMode::NONE a selector writes nothing into the design, so a
 * file opened ReadOnly serves. L2 and ZSCORE standardize the design in place
 * and restore it, which needs a file opened ReadWrite.
 *
 * @param X       Memory-mapped design.
 * @param control Control parameters that hold the scaling mode.
 *
 * @return Eigen::Map of the memory-mapped design.
 *
 * @throws std::runtime_error for a file opened ReadOnly under L2 or ZSCORE.
 */
inline Eigen::Map<Eigen::MatrixXd> mapped_design_for(
    trex::utils::memmap::MemoryMappedMatrix<double>& X,
    const TRexControlParameter& control
) {
    namespace dn = trex::trex_selector_methods::ad::utils::data_normalizer;
    return X.getMapFor(dn::writesDesign(control.scaling_mode));
}

// ==============================================================================

/**
 * @brief Rcpp wrapper class for TRexSelector leveraging zero-copy references.
 */
class RTRexSelector {
protected:

    /** @brief Default constructor */
    RTRexSelector() = default;

public:

    /** @brief Virtual destructor */
    virtual ~RTRexSelector() = default;

protected:
    /** @brief Unique ptr to Eigen::Map for the design matrix. */
    std::unique_ptr<Eigen::Map<Eigen::MatrixXd>> X_map_;

    /** @brief Unique ptr to Eigen::Map for the response vector. */
    std::unique_ptr<Eigen::Map<Eigen::VectorXd>> y_map_;

    /** @brief Unique ptr to the underlying TRexSelector instance. */
    std::unique_ptr<TRexSelector> selector_;

public:
    /** @brief Constructor
     *
     * @param X Design matrix (n x p).
     * @param y Response vector (n).
     * @param tFDR Target FDR level.
     * @param trex_control TRex control parameters.
     * @param seed Random seed.
     * @param verbose Verbosity flag.
     */
    RTRexSelector(
        Eigen::Map<Eigen::MatrixXd> X, // NOLINT (performance-unnecessary-value-param)
        Eigen::Map<Eigen::VectorXd> y, // NOLINT (performance-unnecessary-value-param)
        double tFDR,
        const TRexControlParameter& trex_control,
        int seed,
        bool verbose
    ) {
        X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X);
        y_map_ = std::make_unique<Eigen::Map<Eigen::VectorXd>>(y);

        selector_ = std::make_unique<TRexSelector>(
            *X_map_, *y_map_, tFDR, trex_control, seed, verbose
        );
    }

    /** @brief Get the underlying TRexSelector instance. */
    TRexSelector* get() const { return selector_.get(); }

    /** @brief Run the selection process. */
    virtual TRexSelector::SelectionResult select() { return selector_->select(); }

    /** @brief Give X back before select() (idempotent; no-op after select()). */
    void restore() { selector_->restore(); }

    // ============================================================
    // Public Getters
    // ============================================================

    /** @brief Get the target FDR level. */
    double getTFDR() const noexcept { return selector_->getTFDR(); }

    /** @brief Get the number of random experiments. */
    std::size_t getK() const noexcept { return selector_->getK(); }

    /** @brief Get the maximum dummy multiplier. */
    std::size_t getMaxNumDummies() const noexcept { return selector_->getMaxNumDummies(); }

    /** @brief Get the maximum T stop value. */
    std::size_t getMaxTStop() const noexcept { return selector_->getMaxTStop(); }

    /** @brief Check if stagnation check is enabled. */
    bool getStagnationCheck() const noexcept { return selector_->getStagnationCheck(); }

    /** @brief Get the selected variables. */
    const Eigen::VectorXi& getSelectedVar() const noexcept { return selector_->getSelectedVar(); }

    /** @brief Get the indices of the selected variables. */
    const std::vector<std::size_t>& getSelectedIndices() const noexcept {
        return selector_->getSelectedIndices(); }

    /** @brief Get the number of selected variables. */
    std::size_t getNumSelected() const noexcept { return selector_->getNumSelected(); }

    /** @brief Get the T stop value. */
    std::size_t getTStop() const noexcept { return selector_->getTStop(); }

    /** @brief Get the voting grid. */
    const Eigen::VectorXd& getVotingGrid() const noexcept { return selector_->getVotingGrid(); }

    /** @brief Get the voting threshold. */
    double getVotingThreshold() const noexcept { return selector_->getVotingThreshold(); }

    /** @brief Get means of the design matrix columns. */
    const Eigen::VectorXd& getXMeans() const noexcept { return selector_->getXMeans(); }

    /** @brief Get L2 norms of the design matrix columns. */
    const Eigen::VectorXd& getXL2Norms() const noexcept { return selector_->getXL2Norms(); }

    /** @brief Get mean of the response vector. */
    double getYMean() const noexcept { return selector_->getYMean(); }

    /** @brief Get FDP hat matrix. */
    const Eigen::MatrixXd& getFDPHatMat() const noexcept { return selector_->getFDPHatMat(); }

    /** @brief Get Phi matrix. */
    const Eigen::MatrixXd& getPhiMat() const noexcept { return selector_->getPhiMat(); }

    /** @brief Get selection frequency matrix R. */
    const Eigen::MatrixXd& getRMat() const noexcept { return selector_->getRMat(); }

    /** @brief Get Phi prime vector. */
    const Eigen::VectorXd& getPhiPrime() const noexcept { return selector_->getPhiPrime(); }

    /** @brief Get dummy multiplier L. */
    std::size_t getDummyMultiplierL() const noexcept { return selector_->getDummyMultiplierL(); }
};

// ==============================================================================

/**
 * @brief Linkage of an R control, from its name.
 *
 * @details The controls spell the names in lower case. The capitalized names
 *  of TRexSelector 2.0.1 are accepted as well.
 *
 * @param name  Linkage name: "single", "complete", "average" or "wpgma".
 * @param owner Selector named in the error message.
 *
 * @return The linkage method.
 */
inline trex::ml_methods::clustering::hierarchical::agglomerative::LinkageMethod
parse_hc_linkage(std::string name, const std::string& owner) {
    namespace hac = trex::ml_methods::clustering::hierarchical::agglomerative;
    std::transform(name.begin(), name.end(), name.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (name == "single")   { return hac::LinkageMethod::Single; }
    if (name == "complete") { return hac::LinkageMethod::Complete; }
    if (name == "average")  { return hac::LinkageMethod::Average; }
    if (name == "wpgma")    { return hac::LinkageMethod::WPGMA; }
    Rcpp::stop("Unsupported linkage method for " + owner + ": " + name +
               "\nSupported: single, complete, average, wpgma.");
}
// ==============================================================================
#endif /* End of RCPP_TREX_SELECTOR_H */

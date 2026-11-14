// ===================================================================================
// trex_biobank_screening.hpp
// ===================================================================================
#ifndef TREX_SELECTOR_METHODS_TREX_BIOBANK_SCREENING_HPP
#define TREX_SELECTOR_METHODS_TREX_BIOBANK_SCREENING_HPP
// ===================================================================================
/**
 * @file trex_biobank_screening.hpp
 *
 * @brief Genomics biobank screening implementing Algorithm 1 from the Screen-TRex paper.
 *
 * @details Coordinates Screen-TRex and T-Rex selectors for efficient large-scale
 *          phenotype screening of genomics biobanks.
 */
// ===================================================================================

// std includes
#include <cstddef>
#include <string>
#include <vector>

// Eigen includes
#include <Eigen/Dense>

// TRex includes
#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <trex_selector_methods/ad/trex_screening/trex_screening.hpp>

// ===================================================================================

namespace trex::trex_selector_methods::ad::trex_biobank_screening {

// Namespace aliases
namespace tc = trex::trex_selector_methods::ad::trex_core;
namespace ts = trex::trex_selector_methods::ad::trex_screening;

// ===================================================================================


// =========================================================
// Structures & Enums
// =========================================================

/** @brief Result of screening a single phenotype. */
struct BiobankScreenTRexResult {

    /** @brief Index of the screened phenotype in matrix Y. */
    std::size_t phenotype_index;

    /** @brief Final selected variable indices (from chosen method). */
    std::vector<std::size_t> selected_indices;

    /** @brief FDR level of the final selection.
     *
     * @details The estimate of the accepted screening variant, or, for a
     *  fallback result, the target level `target_FDR_trex`, at which the
     *  T-Rex selector controls the FDR. An estimate can equal the target, so
     *  `used_fallback_trex` tells the two cases apart.
     */
    double estimated_FDR;

    /** @brief Method that gave the final selection: "Screen-TRex (ordinary)",
     *  "Screen-TRex (bootstrap-CI)" or "T-Rex (fallback)". */
    std::string method_used;

    /** @brief Estimated FDR of the ordinary screening, stored for every
     *  phenotype, a fallback result included. */
    double estimated_FDR_screen_ordinary;

    /** @brief Estimated FDR of the bootstrap screening, stored for every
     *  phenotype, a fallback result included. */
    double estimated_FDR_screen_bootstrap;

    /** @brief Selected indices from ordinary Screen-TRex
     *
     * @details Stored for every phenotype, also when another method gave the
     *  final selection.
     */
    std::vector<std::size_t> selected_indices_screen_ordinary;

    /** @brief Selected indices from bootstrap Screen-TRex
     *
     * @details Stored for every phenotype, also when another method gave the
     *  final selection.
     */
    std::vector<std::size_t> selected_indices_screen_bootstrap;

    /** @brief True if the fallback T-Rex selector gave the final selection.
     *  `estimated_FDR` then holds the target level and not an estimate. */
    bool used_fallback_trex;
};


/**
 * @brief Control parameters for biobank screening (Algorithm 1).
 */
struct BiobankScreenTRexControl {

    /** @brief Target FDR level α for fallback T-Rex selector */
    double target_FDR_trex = 0.1;

    /**
     * @brief Lower bound α_l for acceptable estimated FDR.
     *
     * @details If estimated FDR < α_l, power might be too low.
     *          Algorithm 1 uses this to determine whether to accept Screen-TRex result.
     */
    double lower_bound_FDR = 0.05;

    /**
     * @brief Upper bound α_u for acceptable estimated FDR.
     *
     * @details If estimated FDR > α_u, false discoveries might be too many.
     *          Algorithm 1 uses this to determine whether to accept Screen-TRex result.
     */
    double upper_bound_FDR = 0.15;

    /** @brief Control parameters forwarded to the Screen-TRex selector (nests
     *  its own base trex_ctrl, shared by both the Screen-TRex path and the
     *  fallback plain T-Rex path — single source of truth for the base
     *  algorithmic control). Its use_bootstrap_CI has no effect here: the
     *  biobank screening computes the ordinary and the bootstrap variant for
     *  every phenotype. */
    ts::ScreenTRexControlParameter trex_screen_ctrl;

    /** @brief Default constructor. The strategy of the nested control is STANDARD, the default of
     *  the base control. */
    BiobankScreenTRexControl() {
        trex_screen_ctrl.trex_ctrl.lloop_strategy = tc::LLoopStrategy::STANDARD;
    }
};


/**
 * @brief Biobank screening with the T-Rex selector: per phenotype the fast screening
 *        (T-Rex+Screening), and the T-Rex selector as fallback when no FDR estimate of
 *        the screening lies in the accepted range.
 *
 * @details Implements Algorithm 1 of the screening paper, a way to screen thousands of
 *          phenotypes efficiently.
 *
 *  Decision rule per phenotype:
 *  ----------------------------
 *  The screening gives two FDR estimates, the one of its ordinary variant and the one of
 *  its bootstrap variant. A variant is accepted when its estimate lies in
 *  [lower_bound_FDR, upper_bound_FDR]. Of two accepted variants the one with the larger
 *  estimate is taken, and the bootstrap variant on a tie. With no accepted variant the
 *  T-Rex selector runs at target_FDR_trex.
 *
 *  Both estimates come from one set of K experiments. The screening runs once, and the
 *  bootstrap variant is derived from its coefficients.
 *
 *  L-loop strategies: STANDARD (default), SEEDED, PERMUTATION and PERMUTATION_SEEDED, as
 *  for ScreenTRexSelector. HCONCAT and SKIPL are refused.
 *
 *  References:
 *  -----------
 *   J. Machkour, M. Muma, and D. P. Palomar, "False Discovery Rate Control for Fast Screening of Large-Scale
 *     Genomics Biobanks," in Proc. 22nd IEEE Statist. Signal Process. Workshop (SSP), 2023, pp. 666-670.
 *   F. Scheidt, J. Machkour, and M. Muma, "Solving FDR-Controlled Sparse Regression Problems with Five
 *     Million Variables on a Laptop," in Proc. IEEE 9th Int. Workshop Comput. Adv. Multi-Sensor Adapt.
 *     Process. (CAMSAP), 2023, pp. 116-120.
 *   J. Machkour, M. Muma, and D. P. Palomar, "The Terminating-Random Experiments Selector: Fast
 *     High-Dimensional Variable Selection with False Discovery Rate Control," Signal Process., vol. 231,
 *     Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.
 */
class BiobankScreenTRex {
protected:
    // ======================================================================
    // Protected Members
    // ======================================================================

    /** @brief Pointer to design matrix (n x p) */
    Eigen::Map<Eigen::MatrixXd>* X_;

    /** @brief Phenotype vector for single phenotype mode (n x 1) */
    Eigen::VectorXd y_;

    /** @brief Pointer to phenotype matrix for multiple phenotype mode (n x q) */
    Eigen::Map<Eigen::MatrixXd>* Y_;

    /** @brief Control parameters for Biobank Screen-TRex */
    BiobankScreenTRexControl bio_ctrl_;

    /** @brief Random seed for reproducibility. Phenotype i runs with
     *  seed + i, see phenotypeSeed(). */
    int seed_;

    /** @brief Verbosity flag */
    bool verbose_;


public:
    // ======================================================================
    // Constructors
    // ======================================================================

    /**
     * @brief Constructor for single phenotype screening.
     *
     * @param X Design matrix (n x p)
     * @param y Phenotype vector (n x 1)
     * @param bio_ctrl Control parameters for Biobank Screen-TRex (moved)
     * @param seed Random seed. A negative value (default: -1) gives a run that is
     *             not reproducible.
     * @param verbose If true, print progress messages (default: false)
     */
    BiobankScreenTRex(
        Eigen::Map<Eigen::MatrixXd>& X,
        Eigen::Map<Eigen::VectorXd>& y,
        BiobankScreenTRexControl bio_ctrl,
        int seed = -1,
        bool verbose = false
    );


    /**
     * @brief Constructor for multiple phenotype screening.
     *
     * @param X Design matrix (n x p)
     * @param Y Phenotype matrix (n x q). Only read, one column at a time, so
     *          its memory may be read-only, for example a memory-mapped file
     *          opened ReadOnly (MemoryMappedMatrix::getReaderMap()).
     * @param bio_ctrl Control parameters for Biobank Screen-TRex (moved)
     * @param seed Random seed. A negative value (default: -1) gives a run that is
     *             not reproducible.
     *             Phenotype i (0-based) runs with seed + i, so the phenotypes
     *             draw independent dummies and a single phenotype is
     *             reproduced by a single-phenotype run with that seed.
     * @param verbose If true, print progress messages (default: false)
     */
    BiobankScreenTRex(
        Eigen::Map<Eigen::MatrixXd>& X,
        Eigen::Map<Eigen::MatrixXd>& Y,
        BiobankScreenTRexControl bio_ctrl,
        int seed = -1,
        bool verbose = false
    );


    /** @brief Virtual destructor */
    virtual ~BiobankScreenTRex() = default;


    // Delete copy and move constructors and assignment operators
    /** @brief Deleted copy constructor */
    BiobankScreenTRex(const BiobankScreenTRex&) = delete;


    /** @brief Deleted copy assignment operator */
    BiobankScreenTRex& operator=(const BiobankScreenTRex&) = delete;


    /** @brief Deleted move constructor */
    BiobankScreenTRex(BiobankScreenTRex&&) = delete;


    /** @brief Deleted move assignment operator */
    BiobankScreenTRex& operator=(BiobankScreenTRex&&) = delete;


    // ======================================================================
    // Main Methods
    // ======================================================================

    /**
     * @brief Screen a single phenotype using Biobank Screen-TRex (Algorithm 1).
     *
     * @details Procedure:
     *          1. Run bootstrap-CI Screen-TRex: get α̂_C and A_G^(C)(γ, 1)
     *          2. Run ordinary Screen-TRex: get α̂ and A_p(0.5, 1)
     *          3. Apply Algorithm 1 decision logic
     *          4. Return selected set (Screen-TRex) or run fallback T-REx
     *
     * @return Screening result with selected variables and metadata.
     *
     * @throws std::runtime_error if called with multiple phenotypes
     *         (use screenPhenotypes() instead).
     */
    BiobankScreenTRexResult screenPhenotype();


    /**
     * @brief Screen all phenotypes in matrix Y using Biobank Screen-TRex (Algorithm 1).
     *
     * @return Screening results for all phenotypes with selected variables and metadata.
     *
     * @throws std::runtime_error if called without Y matrix (use screenPhenotype() instead).
     */
    std::vector<BiobankScreenTRexResult> screenPhenotypes();


    // ======================================================================
    // Public Getters and Structures
    // ======================================================================

    /**
     * @brief Get summary statistics across all phenotypes.
     */
    struct BiobankScreenTRexStatistics {
        /** @brief Total number of phenotypes. */
        std::size_t num_phenotypes;

        /** @brief Number using ordinary Screen-TRex result. */
        std::size_t num_used_screen_ordinary = 0;

        /** @brief Number using bootstrap-CI Screen-TRex result. */
        std::size_t num_used_screen_bootstrap = 0;

        /** @brief Number of phenotypes using fallback T-REx result. */
        std::size_t num_used_trex = 0;

        /** @brief Average of `estimated_FDR` across all phenotypes: estimates
         *  for screening results and the target level for fallback results. */
        double avg_estimated_FDR;

        /** @brief Average number of selected variables across all phenotypes. */
        double avg_selected_count;
    };


    /**
     * @brief Get the Biobank Screen T-REx Result object.
     *
     * @param results Vector of individual phenotype screening results.
     *
     * @return BiobankScreenTRexStatistics object containing summary statistics.
     */
    BiobankScreenTRexStatistics getBiobankScreenTRexResult(
        const std::vector<BiobankScreenTRexResult>& results
    ) const;


protected:
    // ======================================================================
    // Protected Methods - Algorithm 1
    // ======================================================================

    /**
     * @brief Seed of one phenotype.
     *
     * @details Phenotype i runs its screening and its fallback selector with
     *  seed + i, so the phenotypes draw independent dummies. Phenotype 0 and
     *  the single-phenotype mode keep the seed itself. A negative seed is
     *  passed on, so every selector draws its own entropy.
     *
     * @param phenotype_index Index of the phenotype (0-based).
     *
     * @return The seed for that phenotype.
     */
    int phenotypeSeed(std::size_t phenotype_index) const;

    /** @brief Screen a single phenotype column.
     *
     * @param y Phenotype vector.
     * @param phenotype_index Index of the phenotype being screened.
     *
     * @return Result from Screening after applying Algorithm 1.
     */
    BiobankScreenTRexResult screenSinglePhenotype(
        const Eigen::VectorXd& y,
        std::size_t phenotype_index
    );

    /**
     * @brief Convert binary selection vector to index vector.
     *
     * @param mask Binary vector (0/1) indicating selected variables.
     *
     * @return Vector of indices where mask == 1.
     */
    std::vector<std::size_t> selectedVarsToIndices(const Eigen::VectorXi& mask) const;


    /**
     * @brief Apply Algorithm 1's decision logic (Step 2.2).
     *
     * @details Choses between:
     *      - Bootstrap-CI Screen-TRex result A_G^(C)(γ, 1) with estimate α̂_C
     *      - Ordinary Screen-TRex result A_p(0.5, 1) with estimate α̂
     *      - Fallback T-REx selector result
     *
     * @param sctrex_result_bootstrap Result from bootstrap-CI Screen-TRex.
     * @param sctrex_result_ordinary Result from ordinary Screen-TRex.
     * @param y Phenotype vector.
     * @param phenotype_index Index of the phenotype being screened.
     *
     * @return Final screening result after applying decision logic.
     */
    BiobankScreenTRexResult applyDecisionLogic(
        const ts::ScreenTRexSelectionResult& sctrex_result_bootstrap,
        const ts::ScreenTRexSelectionResult& sctrex_result_ordinary,
        const Eigen::VectorXd& y,
        std::size_t phenotype_index
    );


    /**
     * @brief Run fallback T-REx selector selector (Algorithm 1, Step 2.3).
     *
     * @param y Phenotype vector.
     * @param phenotype_idx Index of the phenotype being screened.
     *
     * @return Result from T-REx selector.
     */
    BiobankScreenTRexResult runFallbackTRexSelector(
        const Eigen::VectorXd& y,
        std::size_t phenotype_idx
    );


    /** @brief Print progress message if verbosity is enabled.
     *
     * @param msg Message to print
     */
    void printProgressMessage(const std::string& msg) const;

    /** @brief Validate the configured dummy strategy at construction.
     *
     * @details Mirrors ScreenTRexSelector's gate so an unsupported strategy
     *  fails immediately instead of mid-way through the per-phenotype loop
     *  (the per-phenotype selectors are only built inside screening).
     *
     * @throws std::invalid_argument on an unsupported lloop_strategy.
     */
    void validateBiobankStrategy() const;

};

// ===================================================================================
} /* End of namespace trex::trex_selector_methods::ad::trex_biobank_screening */
// ===================================================================================
#endif /* End of TREX_SELECTOR_METHODS_TREX_BIOBANK_SCREENING_HPP */

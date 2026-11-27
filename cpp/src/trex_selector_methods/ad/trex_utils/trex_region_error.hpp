// ===================================================================================
// trex_region_error.hpp
// ===================================================================================
#ifndef TREX_SELECTOR_METHODS_UTILS_REGION_ERROR_HPP
#define TREX_SELECTOR_METHODS_UTILS_REGION_ERROR_HPP
// ===================================================================================
/**
 * @file trex_region_error.hpp
 *
 * @brief Carries the first error of an OpenMP region out of it as kind and text.
 *
 * @details
 * An exception may not leave an OpenMP structured block. The first one is
 * recorded inside the region and thrown again after it. The record holds the
 * kind and the text of the exception and no std::exception_ptr. In a process
 * with two C++ runtimes (libc++ beside libstdc++) an exception thrown again
 * from a std::exception_ptr matches catch (...) alone, so its type and its text
 * are lost. A fresh exception of the recorded kind keeps both.
 */
// ===================================================================================

// std includes
#include <atomic>
#include <new>
#include <stdexcept>
#include <string>

// Embedded into namespace trex::trex_selector_methods::ad::utils::region_error
namespace trex::trex_selector_methods::ad::utils::region_error {

/**
 * @brief First error of a parallel region, kept as its kind and its text.
 *
 * @details
 * Usage: test raised() at the top of an iteration, call capture() in the
 * catch (...) of the iteration body, call rethrow() after the region.
 */
class RegionError {
public:
    /**
     * @brief Tells whether an error is recorded.
     * @return true once capture() has recorded an error
     */
    [[nodiscard]] bool raised() const noexcept {
        return raised_.load(std::memory_order_relaxed);
    }

    /**
     * @brief Records the exception in flight. The first call wins.
     * @note Call inside a catch block. Any thread of the region may call it.
     */
    void capture() noexcept {
        Kind kind = Kind::RUNTIME_ERROR;
        std::string text;
        try {
            try { throw; }
            catch (const std::invalid_argument& e) { kind = Kind::INVALID_ARGUMENT; text = e.what(); }
            catch (const std::domain_error& e)     { kind = Kind::DOMAIN_ERROR;     text = e.what(); }
            catch (const std::length_error& e)     { kind = Kind::LENGTH_ERROR;     text = e.what(); }
            catch (const std::out_of_range& e)     { kind = Kind::OUT_OF_RANGE;     text = e.what(); }
            catch (const std::logic_error& e)      { kind = Kind::LOGIC_ERROR;      text = e.what(); }
            catch (const std::range_error& e)      { kind = Kind::RANGE_ERROR;      text = e.what(); }
            catch (const std::overflow_error& e)   { kind = Kind::OVERFLOW_ERROR;   text = e.what(); }
            catch (const std::underflow_error& e)  { kind = Kind::UNDERFLOW_ERROR;  text = e.what(); }
            catch (const std::bad_alloc&)          { kind = Kind::BAD_ALLOC; }
            catch (const std::exception& e)        { kind = Kind::RUNTIME_ERROR;    text = e.what(); }
            catch (...) { kind = Kind::RUNTIME_ERROR; text = "unknown error in a parallel region"; }
        } catch (...) {
            // no memory for the copy of the text
            kind = Kind::BAD_ALLOC;
        }
        #pragma omp critical(trex_region_error)
        {
            if (!raised_.load(std::memory_order_relaxed)) {
                kind_ = kind;
                text_.swap(text);
                raised_.store(true, std::memory_order_relaxed);
            }
        }
    }

    /**
     * @brief Throws a fresh exception of the recorded kind with the recorded text.
     * @note Call after the region. Returns when no error is recorded.
     * @throws the standard exception class of the recorded error, std::runtime_error
     *         for an exception of another class
     */
    void rethrow() const {
        if (!raised()) { return; }
        switch (kind_) {
            case Kind::INVALID_ARGUMENT: throw std::invalid_argument(text_);
            case Kind::DOMAIN_ERROR:     throw std::domain_error(text_);
            case Kind::LENGTH_ERROR:     throw std::length_error(text_);
            case Kind::OUT_OF_RANGE:     throw std::out_of_range(text_);
            case Kind::LOGIC_ERROR:      throw std::logic_error(text_);
            case Kind::RANGE_ERROR:      throw std::range_error(text_);
            case Kind::OVERFLOW_ERROR:   throw std::overflow_error(text_);
            case Kind::UNDERFLOW_ERROR:  throw std::underflow_error(text_);
            case Kind::BAD_ALLOC:        throw std::bad_alloc();
            case Kind::RUNTIME_ERROR:    break;
        }
        throw std::runtime_error(text_);
    }

private:
    /// Standard exception classes that keep their identity across the region.
    enum class Kind {
        INVALID_ARGUMENT, DOMAIN_ERROR, LENGTH_ERROR, OUT_OF_RANGE, LOGIC_ERROR,
        RANGE_ERROR, OVERFLOW_ERROR, UNDERFLOW_ERROR, BAD_ALLOC, RUNTIME_ERROR
    };

    std::atomic<bool> raised_{false};
    Kind kind_ = Kind::RUNTIME_ERROR;
    std::string text_;
};

} // namespace trex::trex_selector_methods::ad::utils::region_error

#endif /* End of TREX_SELECTOR_METHODS_UTILS_REGION_ERROR_HPP */

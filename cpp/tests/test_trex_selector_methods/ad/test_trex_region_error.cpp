// ========================================================================================
/**
 * @file test_trex_region_error.cpp
 *
 * @brief An error raised inside an OpenMP region must arrive after the region
 *  with its standard exception class and its text. RegionError carries both as
 *  plain data and throws a fresh exception, so nothing depends on
 *  std::exception_ptr, which loses type and text in a process that holds two
 *  C++ runtimes.
 */
// ========================================================================================

// google test includes
#include <gtest/gtest.h>

// std includes
#include <new>
#include <stdexcept>
#include <string>

// project includes
#include <trex_selector_methods/ad/trex_utils/trex_region_error.hpp>

// ========================================================================================

namespace {

using trex::trex_selector_methods::ad::utils::region_error::RegionError;

/** @brief Throws `Thrown` in every iteration of a parallel loop and returns the record. */
template <typename Thrown>
void raiseInRegion(RegionError& error, const Thrown& thrown, int iterations = 8) {
    #pragma omp parallel for num_threads(2)
    for (int k = 0; k < iterations; ++k) {
        if (error.raised()) { continue; }
        try {
            throw thrown;
        } catch (...) {
            error.capture();
        }
    }
}

/** @brief Expects that `Thrown` crosses the region as the class `Expected` with its text. */
template <typename Thrown, typename Expected>
void expectKindAndText(const std::string& text) {
    RegionError error;
    raiseInRegion(error, Thrown(text));
    ASSERT_TRUE(error.raised());
    try {
        error.rethrow();
        FAIL() << "rethrow() returned";
    } catch (const Expected& e) {
        EXPECT_EQ(std::string(e.what()), text);
    } catch (...) {
        FAIL() << "another class than expected";
    }
}

/** @brief An exception class outside the standard hierarchy of the library. */
struct Foreign : std::exception {
    const char* what() const noexcept override { return "foreign text"; }
};

}  // namespace

// ========================================================================================

/** @brief Without an error the record stays empty and rethrow() returns. */
TEST(RegionError, NoErrorNoThrow) {
    RegionError error;
    EXPECT_FALSE(error.raised());
    EXPECT_NO_THROW(error.rethrow());
}

/** @brief Every standard class the library throws keeps its identity and its text. */
TEST(RegionError, StandardClassesKeepKindAndText) {
    expectKindAndText<std::invalid_argument, std::invalid_argument>("invalid argument text");
    expectKindAndText<std::domain_error, std::domain_error>("domain error text");
    expectKindAndText<std::length_error, std::length_error>("length error text");
    expectKindAndText<std::out_of_range, std::out_of_range>("out of range text");
    expectKindAndText<std::logic_error, std::logic_error>("logic error text");
    expectKindAndText<std::range_error, std::range_error>("range error text");
    expectKindAndText<std::overflow_error, std::overflow_error>("overflow error text");
    expectKindAndText<std::underflow_error, std::underflow_error>("underflow error text");
    expectKindAndText<std::runtime_error, std::runtime_error>("runtime error text");
}

/** @brief A derived class must not be caught as one of its bases of another branch. */
TEST(RegionError, InvalidArgumentIsNoRuntimeError) {
    RegionError error;
    raiseInRegion(error, std::invalid_argument("bad input"));
    EXPECT_THROW(error.rethrow(), std::logic_error);
    try {
        error.rethrow();
    } catch (const std::runtime_error&) {
        FAIL() << "invalid_argument arrived as runtime_error";
    } catch (const std::exception&) {
        SUCCEED();
    }
}

/** @brief Memory exhaustion arrives as std::bad_alloc. */
TEST(RegionError, BadAllocStaysBadAlloc) {
    RegionError error;
    raiseInRegion(error, std::bad_alloc());
    EXPECT_THROW(error.rethrow(), std::bad_alloc);
}

/** @brief Another exception class arrives as std::runtime_error with its text. */
TEST(RegionError, OtherClassBecomesRuntimeErrorWithText) {
    RegionError error;
    raiseInRegion(error, Foreign());
    try {
        error.rethrow();
        FAIL() << "rethrow() returned";
    } catch (const std::runtime_error& e) {
        EXPECT_EQ(std::string(e.what()), "foreign text");
    }
}

/** @brief An object that is no std::exception arrives as std::runtime_error. */
TEST(RegionError, NonStandardObjectBecomesRuntimeError) {
    RegionError error;
    raiseInRegion(error, 42);
    try {
        error.rethrow();
        FAIL() << "rethrow() returned";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("unknown error"), std::string::npos) << e.what();
    }
}

/** @brief The first recorded error wins, later ones do not replace it. */
TEST(RegionError, FirstErrorWins) {
    RegionError error;
    try { throw std::invalid_argument("first"); } catch (...) { error.capture(); }
    try { throw std::runtime_error("second"); } catch (...) { error.capture(); }
    try {
        error.rethrow();
        FAIL() << "rethrow() returned";
    } catch (const std::invalid_argument& e) {
        EXPECT_EQ(std::string(e.what()), "first");
    }
}

/** @brief rethrow() can be called again and throws the same error. */
TEST(RegionError, RethrowIsRepeatable) {
    RegionError error;
    try { throw std::out_of_range("index"); } catch (...) { error.capture(); }
    EXPECT_THROW(error.rethrow(), std::out_of_range);
    EXPECT_THROW(error.rethrow(), std::out_of_range);
}

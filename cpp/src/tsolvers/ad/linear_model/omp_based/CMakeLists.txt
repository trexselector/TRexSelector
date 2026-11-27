# src/tsolvers/ad/linear_model/omp_based/CMakeLists.txt
# ==============================================================================
# CMakeLists.txt for the OMP-based solvers for the linear regression model.
# ==============================================================================

# ==============================================================================
# 1. Append Source Files (.cpp - private)
# ==============================================================================
target_sources(trex_tsolvers PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/tsolver_omp.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tsolver_gp.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tsolver_acgp.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tsolver_mp.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tsolver_ncgmp.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tsolver_ools.cpp
)
# ==============================================================================


# ==============================================================================
# 2. Include Directories (.hpp - public)
# ==============================================================================
# Marked PUBLIC so that CMake tools and IDEs index them correctly.
target_sources(trex_tsolvers PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src          # project root /src
    FILES
        # OMP based solvers
        # -----------------------------------------------
        ${CMAKE_CURRENT_LIST_DIR}/tsolver_omp.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tsolver_gp.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tsolver_acgp.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tsolver_mp.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tsolver_ncgmp.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tsolver_ools.hpp
        # old include paths, forwarding headers for one release
        ${CMAKE_CURRENT_LIST_DIR}/tacgp_solver.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tgp_solver.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tmp_solver.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tncgmp_solver.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tomp_solver.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tools_solver.hpp
)
# ==============================================================================

# src/tsolvers/ad/linear_model/lars_based/CMakeLists.txt
# ==============================================================================
# CMakeLists.txt for the LARS-based solvers for the linear regression model.
# ==============================================================================


# ==============================================================================
# 1. Append Source Files (.cpp - private)
# ==============================================================================
target_sources(trex_tsolvers PRIVATE
    # LARS based solvers
    # -----------------------------------------------
    ${CMAKE_CURRENT_LIST_DIR}/tsolver_enet_aug.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tsolver_enet.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tsolver_ienet.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tsolver_ienet_aug.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tsolver_tiknet.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tsolver_lars.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tsolver_lasso.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tsolver_stagewise.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tsolver_stepwise.cpp
)
# ==============================================================================


# ==============================================================================
# 2. Append Header Files (.hpp - public)
# ==============================================================================
target_sources(trex_tsolvers PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src          # project root /src
    FILES
        # LARS based solvers
        # -----------------------------------------------
        ${CMAKE_CURRENT_LIST_DIR}/tsolver_enet_aug.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tsolver_enet.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tsolver_ienet.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tsolver_ienet_aug.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tsolver_tiknet.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tsolver_lars.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tsolver_lasso.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tsolver_stagewise.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tsolver_stepwise.hpp
        # old include paths, forwarding headers for one release
        ${CMAKE_CURRENT_LIST_DIR}/tenet_aug_solver.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tenet_solver.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tienet_aug_solver.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tienet_solver.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tiknet_solver.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tlars_solver.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tlasso_solver.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tstagewise_solver.hpp
        ${CMAKE_CURRENT_LIST_DIR}/tstepwise_solver.hpp
)
# ==============================================================================

# Example programs of the C++ user guide. Each program is compiled against the installed
# TRexSelector library and run at documentation build time; the sections of its output land in
# DOC_OUTPUT_DIR and the guide pages include them next to the code.
cmake_minimum_required(VERSION 3.23...4.3)

# Match the compiler of the library (Homebrew LLVM on macOS).
if(APPLE AND NOT DEFINED CMAKE_CXX_COMPILER AND EXISTS /opt/homebrew/opt/llvm/bin/clang++)
    set(CMAKE_C_COMPILER   /opt/homebrew/opt/llvm/bin/clang   CACHE FILEPATH "C compiler")
    set(CMAKE_CXX_COMPILER /opt/homebrew/opt/llvm/bin/clang++ CACHE FILEPATH "C++ compiler")
endif()

project(TRexSelectorGuideExamples LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(TRexSelector REQUIRED)

set(DOC_OUTPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/output" CACHE PATH "Folder of the captured outputs")

set(GUIDE_EXAMPLES
    getting_started
    dimensionality
    trex_selector
    trex_da
    trex_gvs
    screening
    trex_spca
    trex_tikhonov
    tsolvers
    ml_methods
    memory_mapping
    no_copies
    parallelism
)

set(RUN_COMMANDS)
foreach(example IN LISTS GUIDE_EXAMPLES)
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} PRIVATE TRexSelector::trex_selector_methods)
    target_compile_definitions(${example} PRIVATE DOC_OUTPUT_DIR="${DOC_OUTPUT_DIR}")
    list(APPEND RUN_COMMANDS COMMAND $<TARGET_FILE:${example}>)
endforeach()

# `cmake --build <dir> --target run_examples` runs every program in order.
add_custom_target(run_examples
    ${RUN_COMMANDS}
    DEPENDS ${GUIDE_EXAMPLES}
    COMMENT "Running the guide examples"
    VERBATIM
)

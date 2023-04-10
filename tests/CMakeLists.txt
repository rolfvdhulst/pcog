cmake_minimum_required(VERSION 3.10)
project(pcog-tests)

#----------------------------------------------------------------------------------------------------------------------
# general settings and options
#----------------------------------------------------------------------------------------------------------------------

include("../cmake/utils.cmake")
string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}" is_top_level)

if(is_top_level)
    enable_testing()
endif()

#----------------------------------------------------------------------------------------------------------------------
# tests dependencies
#----------------------------------------------------------------------------------------------------------------------

if(is_top_level)
    find_package(pcog REQUIRED)
endif()

find_package(GTest REQUIRED)

#----------------------------------------------------------------------------------------------------------------------
# tests sources
#----------------------------------------------------------------------------------------------------------------------

set(test_sources
        test_main.cpp
        LPSolverTest.cpp
        DenseSetTest.cpp
        NodeMapTest.cpp
        DenseGraphTest.cpp
        ColoringTest.cpp
        )
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${test_sources})

#----------------------------------------------------------------------------------------------------------------------
# tests target
#----------------------------------------------------------------------------------------------------------------------

add_executable(pcog-tests RbTreeTest.cpp)
target_sources(pcog-tests PRIVATE ${test_sources})

target_link_libraries(pcog-tests
        PRIVATE
        pcog::pcog
        GTest::GTest)
set_target_properties(pcog-tests PROPERTIES
        CXX_STANDARD 20)
if(PCOG_COVERAGE)
    target_code_coverage(pcog-tests)
endif()
if(NOT is_top_level)
    win_copy_deps_to_target_dir(pcog-tests pcog::pcog)
endif()

include(GoogleTest)
gtest_discover_tests(pcog-tests)

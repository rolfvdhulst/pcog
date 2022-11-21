cmake_minimum_required(VERSION 3.10)
project(pcog VERSION 0.1
        DESCRIPTION "Exact solution methods for difficult combinatorial optimization problems in graphs"
        LANGUAGES CXX)

include(cmake/utils.cmake)
#Check if we are included as a subdirectory of another cmake
string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}" is_top_level)

#PCOG_SHARED_LIBS option can be used to force shared/static build
option(PCOG_BUILD_TESTS "Turn on to compile the tests" OFF)
option(PCOG_INSTALL "Generate target for installing the library" OFF)
option(PCOG_BUILD_DOCS "Build documentation" OFF)

set_if_undefined(PCOG_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/mylib" CACHE STRING
        "Install path for pcog package-related CMake files")

if(DEFINED PCOG_SHARED_LIBS)
    set(BUILD_SHARED_LIBS  ${PCOG_SHARED_LIBS})
endif()

if(NOT DEFINED CMAKE_BUILD_TYPE AND NOT DEFINED CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set_if_undefined(CMAKE_CXX_VISIBILITY_PRESET hidden)
set_if_undefined(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

add_library(pcog) #initialized below
add_library(pcog::pcog ALIAS pcog)

#Dependencies:


find_package(SoPlex REQUIRED)
#Setup target

include(GenerateExportHeader)
set(export_file_name "export_shared.h")

if(NOT BUILD_SHARED_LIBS)
    set(export_file_name "export_static.h")
endif()
generate_export_header(pcog EXPORT_FILE_NAME include/pcog/${export_file_name})

set(sources
        include/pcog/export.h
        include/pcog/pcog.hpp
        src/pcog.cpp
        include/pcog/LPSolver.hpp
        src/LPSolver.cpp
        )
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${sources})


include(CMakePackageConfigHelpers)

target_sources(pcog PRIVATE ${sources})
target_compile_definitions(pcog PUBLIC "$<$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>:PCOG_STATIC_DEFINE>")

target_link_libraries(pcog
        PUBLIC ${SOPLEX_LIBRARIES}
        )
target_include_directories(pcog
        PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>"
        ${SOPLEX_INCLUDE_DIRS}
        )

set_target_properties(pcog PROPERTIES
        SOVERSION ${PROJECT_VERSION_MAJOR}
        VERSION ${PROJECT_VERSION}
        CXX_STANDARD 20
        )

if(PCOG_INSTALL AND NOT CMAKE_SKIP_INSTALL_RULES)

    configure_package_config_file(cmake/pcog-config.cmake.in pcog-config.cmake
            INSTALL_DESTINATION "${PCOG_INSTALL_CMAKEDIR}")

    write_basic_package_version_file(pcog-config-version.cmake
            COMPATIBILITY SameMajorVersion)

    install(TARGETS pcog EXPORT pcog_export
            RUNTIME COMPONENT pcog
            LIBRARY COMPONENT pcog NAMELINK_COMPONENT pcog-dev
            ARCHIVE COMPONENT pcog-dev
            INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
    install(DIRECTORY include/
            TYPE INCLUDE
            COMPONENT pcog-dev)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/include/pcog/${export_file_name}"
            COMPONENT pcog-dev
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/pcog")

    set(targets_file "pcog-shared-targets.cmake")

    if(NOT BUILD_SHARED_LIBS)
        set(targets_file "pcog-static-targets.cmake")
    endif()

    install(EXPORT pcog_export
            COMPONENT pcog-dev
            FILE "${targets_file}"
            DESTINATION "${PCOG_INSTALL_CMAKEDIR}"
            NAMESPACE pcog::)

    install(FILES
            "${CMAKE_CURRENT_BINARY_DIR}/mylib-config.cmake"
            "${CMAKE_CURRENT_BINARY_DIR}/mylib-config-version.cmake"
            COMPONENT pcog-dev
            DESTINATION "${PCOG_INSTALL_CMAKEDIR}")

    #Necessary for debug version on windows
    if(MSVC)
        set(pdb_file "")
        set(pdb_file_destination "")

        if(BUILD_SHARED_LIBS)
            set(pdb_file "$<TARGET_PDB_FILE:pcog>")
            set(pdb_file_destination "${CMAKE_INSTALL_BINDIR}")
        else()
            # TARGET_PDB_FILE does not work for pdb file generated for static library build, determining it manually
            set(pdb_file "$<TARGET_FILE_DIR:pcog>/$<TARGET_FILE_PREFIX:pcog>$<TARGET_FILE_BASE_NAME:pcog>.pdb")
            set(pdb_file_destination "${CMAKE_INSTALL_LIBDIR}")
        endif()

        install(FILES "${pdb_file}"
                COMPONENT pcog-dev
                CONFIGURATIONS Debug RelWithDebInfo
                DESTINATION "${pdb_file_destination}"
                OPTIONAL)
    endif()
endif()


if(PCOG_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
if(PCOG_BUILD_DOCS)
    find_package(Doxygen REQUIRED)
    doxygen_add_docs(docs include)
endif()
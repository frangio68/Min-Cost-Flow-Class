# --------------------------------------------------------------------------- #
#    Main CMake file for MCFClass                                             #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.8)
# ----- Project description and directories --------------------------------- #
# TODO: set proper version
project(MCFClass VERSION 1.0.0)

# ----- Requirements -------------------------------------------------------- #
# find_package() is used to find libraries and automatically configure them.
# It requires either a package configuration module (***Config.cmake) or a
# find module (Find***.cmake), usually provided by libraries.

# If we are configuring from the umbrella project, we assume that the valid
# settings are there.
# If we are building this project alone we use its settings but we check
# also for settings in its parent directory.
# TODO: Check if this behaviour is correct
if (NOT SMS++_CONFIG_ALL)
    include(CMakeCustom.txt OPTIONAL)
    include(../CMakeCustom.txt OPTIONAL)
endif ()

# This is to use the find modules we provide
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(CPLEX)

# ----- Library ------------------------------------------------------------- #
add_library(MCFClass STATIC)
target_compile_features(MCFClass PUBLIC cxx_std_17)

target_sources(MCFClass PRIVATE
               MCFSimplex/MCFSimplex.C
               RelaxIV/RelaxIV.C
               SPTree/SPTree.C
               )

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake), that must be different
# for local build and installation.
target_include_directories(MCFClass PUBLIC
                           $<INSTALL_INTERFACE:include/MCFClass>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/MCFClass>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/MCFClone>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/MCFSimplex>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/OPTUtils>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/RelaxIV>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/SPTree>
                           )

if (CPLEX_FOUND)
    target_sources(MCFClass PRIVATE MCFCplex/MCFCplex.C)
    target_include_directories(
            MCFClass PUBLIC
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/MCFCplex>)
    target_link_libraries(MCFClass PUBLIC CPLEX::Cplex)
endif ()

# Alias for testing
add_library(MCFClass::MCFClass ALIAS MCFClass)

# ----- Install instructions ------------------------------------------------ #
# For installing the library and all its CMake modules on the system.
# Not required for local builds (see below).

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/MCFClass)

# Install library
install(TARGETS MCFClass
        EXPORT MCFClass-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install headers
install(FILES
        ${PROJECT_SOURCE_DIR}/MCFClass/MCFClass.h
        ${PROJECT_SOURCE_DIR}/MCFClone/MCFClone.h
        ${PROJECT_SOURCE_DIR}/MCFCplex/MCFCplex.h
        ${PROJECT_SOURCE_DIR}/MCFSimplex/MCFSimplex.h
        ${PROJECT_SOURCE_DIR}/OPTUtils/OPTUtils.h
        ${PROJECT_SOURCE_DIR}/RelaxIV/RelaxIV.h
        ${PROJECT_SOURCE_DIR}/SPTree/SPTree.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MCFClass)

# Install target import file
install(EXPORT MCFClass-targets
        FILE MCFClassTargets.cmake
        NAMESPACE MCFClass::
        DESTINATION ${INSTALL_CONFIGDIR})

# Generate package version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate package configuration file
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/MCFClassConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfig.cmake
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR})

# Install package version and configuration files
# TODO: Check if we should install also find modules
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfigVersion.cmake
        # ${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPLEX.cmake
        DESTINATION ${INSTALL_CONFIGDIR})

# ----- Export from the build tree ------------------------------------------ #
# This is needed so other projects can find the library in a local build.
# It does not affect installation.

export(EXPORT MCFClass-targets
       FILE ${CMAKE_CURRENT_BINARY_DIR}/MCFClassTargets.cmake
       NAMESPACE MCFClass::)

# Export also the find modules for the local config file to use it
configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPLEX.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/FindCPLEX.cmake
               COPYONLY)

# FIXME: I think that exporting this is a bad idea
#configure_file(${CMAKE_CURRENT_LIST_DIR}/CMakeCustom.txt
#               ${CMAKE_CURRENT_BINARY_DIR}/CMakeCustom.txt
#               COPYONLY)

# Register package in user's package registry (~/.cmake/registry)
export(PACKAGE MCFClass)

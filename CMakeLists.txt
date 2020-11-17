# --------------------------------------------------------------------------- #
#    Main CMake file for MCFClass                                             #
#                                                                             #
#    This file allows one to build the library using CMake.                   #
#    To do so, you can use the following commands:                            #
#                                                                             #
#        $ cmake -S <source-path> -B <build-path>                             #
#        $ cmake --build <build-path>                                         #
#                                                                             #
#    The following command also installs the library in the system:           #
#                                                                             #
#        $ cmake --build <build-path> --target install                        #
#                                                                             #
#    Note: this file tries to respect the principles of Modern CMake,         #
#          keep it in mind when editing it.                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.11)
# TODO: set proper version
project(MCFClass VERSION 1.0.0 LANGUAGES CXX)

# ----- Settings ------------------------------------------------------------ #
# An option is an ON/OFF user-settable cache variable
option(MCFClass_USE_CPLEX "Whether MCFClass will use CPLEX" ON)

# ----- Requirements -------------------------------------------------------- #
# This adds the cmake directory to the find module search paths,
# allowing us to use our find modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if (MCFClass_USE_CPLEX)
    find_package(CPLEX) # Needed for MCFCplex
endif ()

# This removes from the search paths the cmake directory we added above.
list(REMOVE_AT CMAKE_MODULE_PATH -1)

# Update cache variables
# FIXME: Not sure if this is a best practice
if (NOT CPLEX_FOUND)
    set(MCFClass_USE_CPLEX OFF)
endif ()

# ----- Configuration header ------------------------------------------------ #
# This will generate a *Config.h header in the build directory.
configure_file(cmake/MCFClassConfig.h.in MCFClassConfig.h)

# ----- Library ------------------------------------------------------------- #
# With the BUILD_SHARED_LIBS variable we can specify if the library will
# be STATIC or SHARED, so no reason to do it now.
add_library(MCFClass)
target_compile_features(MCFClass PUBLIC cxx_std_17)

# When adding source files with target_sources(), PRIVATE means that the files
# should only be added to this library, whereas PUBLIC means they should be
# added to this library and to any target that links to it.
# INTERFACE can be used for sources that should not be added to this library
# but should be added to anything that links to it.
target_sources(MCFClass PRIVATE
               MCFSimplex/MCFSimplex.C
               RelaxIV/RelaxIV.C
               SPTree/SPTree.C)

# When using target_include_directories(), PUBLIC means that any targets
# that link to this target also need that include directory.
# Other options are PRIVATE (only affect the current target, not dependencies),
# and INTERFACE (only needed for dependencies).

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake).
# This means that if a target finds this library in its build directory
# will look into the BUILD_INTERFACE path for its headers, if it finds it
# installed in the system will look into the INSTALL_INTERFACE path.
target_include_directories(
        MCFClass PUBLIC
        $<INSTALL_INTERFACE:include/MCFClass>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MCFClass>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MCFClone>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MCFSimplex>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/OPTUtils>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/RelaxIV>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SPTree>)

# When linking other targets to the library with target_link_libraries(),
# PRIVATE means that the libraries will be linked only to this library,
# PUBLIC means they will be linked also to the targets that depend on this
# library, INTERFACE means they will be linked only to the targets that depend
# on this library.
if (MCFClass_USE_CPLEX)
    message(STATUS "MCFClass: CPLEX found, configuring MCFCplex")
    target_sources(MCFClass PRIVATE MCFCplex/MCFCplex.C)
    target_include_directories(
            MCFClass PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MCFCplex>)
    target_link_libraries(MCFClass PUBLIC CPLEX::Cplex)
endif ()

# This alias is defined so that executables in this same project can use
# the library with this notation.
add_library(MCFClass::MCFClass ALIAS MCFClass)

# ----- Install instructions ------------------------------------------------ #
# The following commands are used when installing the library
# and its CMake configuration files on the system.
# They are not required for local builds (see below).

include(GNUInstallDirs)

# Install the library
install(TARGETS MCFClass
        EXPORT MCFClass-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install the headers
install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/MCFClass/MCFClass.h
        ${CMAKE_CURRENT_SOURCE_DIR}/MCFClone/MCFClone.h
        ${CMAKE_CURRENT_SOURCE_DIR}/MCFSimplex/MCFSimplex.h
        ${CMAKE_CURRENT_SOURCE_DIR}/OPTUtils/OPTUtils.h
        ${CMAKE_CURRENT_SOURCE_DIR}/RelaxIV/RelaxIV.h
        ${CMAKE_CURRENT_SOURCE_DIR}/SPTree/SPTree.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MCFClass)

if (MCFClass_USE_CPLEX)
    install(FILES
            ${CMAKE_CURRENT_SOURCE_DIR}/MCFCplex/MCFCplex.h
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MCFClass)
endif ()

# Install the auto-generated configuration header (see above).
install(FILES ${PROJECT_BINARY_DIR}/MCFClassConfig.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MCFClass)

# (Generate and) install the target import file, that allows other
# CMake projects to import the target.
install(EXPORT MCFClass-targets
        FILE MCFClassTargets.cmake
        NAMESPACE MCFClass::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MCFClass)

# Generate the package version file, that allows other
# CMake projects to know the version.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate the package configuration file, that allows other
# CMake projects to find the library with find_package().
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/MCFClassConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MCFClass)

# Install the package version and configuration files.
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MCFClass)

# Shipping custom find modules should not be encouraged because third
# party developers should provide them.
# We ship FindCPLEX because CPLEX does not on its own.
install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPLEX.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MCFClass)


# ----- Export from the build tree ------------------------------------------ #
# The following commands are needed so other projects can find the library
# in a local build directory. They use CMake's User Package Registry
# (located in ~/.cmake/registry) to do so.
# See: https://cmake.org/cmake/help/latest/manual/cmake-packages.7.html
# They don't affect the installation.

# CMAKE_EXPORT_PACKAGE_REGISTRY is OFF by default.
# It may be ON if the project is included by the SMS++ umbrella project.
if (CMAKE_EXPORT_PACKAGE_REGISTRY)

    # Generate the target import file, that allows other
    # CMake projects to import the target.
    export(EXPORT MCFClass-targets
           FILE ${CMAKE_CURRENT_BINARY_DIR}/MCFClassTargets.cmake
           NAMESPACE MCFClass::)

    # Export the find modules for the local configuration file to use them.
    configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPLEX.cmake
                   ${CMAKE_CURRENT_BINARY_DIR}/FindCPLEX.cmake
                   COPYONLY)

    # Export the package in user's package registry
    export(PACKAGE MCFClass)
endif ()

# ----- Subdirectories ------------------------------------------------------ #
if (BUILD_TESTING)
    # Just consuming BUILD_TESTING to avoid warnings
endif ()

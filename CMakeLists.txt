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
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.19)
cmake_policy(VERSION 3.12)
project(MCFClass
        VERSION 1.0.0
        DESCRIPTION "C++ solvers for MCF problems"
        HOMEPAGE_URL http://groups.di.unipi.it/optimize/Software/MCF.html
        LANGUAGES C CXX)

# These variables make the code harder to read but easier to change.
set(modName ${PROJECT_NAME})
set(modNamespace "MCFClass")

# This adds the cmake directory to the module search paths,
# allowing us to use our modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ----- Settings ------------------------------------------------------------ #
# This creates an ENABLE_TESTING option (default: ON) and enables the testing
# using the ctest executable.
# See: https://cmake.org/cmake/help/latest/manual/ctest.1.html
#include(CTest)

# Sets the default build type (if none was specified).
# See: https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html
include(BuildType)

# An option is an ON/OFF user-settable cache variable
option(MCFClass_USE_CPLEX "Whether ${modName} will use CPLEX or not." ON)

option(MCFClass_READ_PATHS "Whether ${modName} will read locations for
                            dependencies or not." ON)

# ----- Requirements -------------------------------------------------------- #
if (MCFClass_READ_PATHS)
    # Read the dependencies default location from a file written in a
    # bash-style syntax since it is shared with the other hand-written Makefiles,
    # so we define a function to convert these macros written as `lhs = rhs` into
    # their equivalent cmake variants, i.e., `set(lhs rhs)`, avoiding comments.
    macro(makefile_macros_to_cmake file)
        FILE(STRINGS ${file} paths)
        foreach (path IN LISTS paths)
            if (path) # ignore empty lines, i.e., \n
                # makes empty the comment lines, i.e. those start with #...
                string(REGEX REPLACE "^[ \t]*[#].*" "" path ${path})
                if (path) # ... and then ignore them

                    string(REGEX MATCH ".*[=]" lhs "${path}") # retrieve the lhs
                    string(REGEX REPLACE "=" "" lhs ${lhs}) # remove =
                    # remove tab or white spaces from head
                    string(STRIP "${lhs}" lhs)

                    # if it's not already defined before in custom makefile-paths...
                    if (NOT DEFINED ${lhs})
                        # then read its rhs...
                        string(REGEX MATCH "[=].*" rhs "${path}") # retrieve the rhs
                        # remove comments starting with # from tails, if any
                        string(REGEX REPLACE "[#].*" "" rhs ${rhs})
                        string(REGEX REPLACE "=" "" rhs ${rhs}) # remove =
                        string(REGEX REPLACE "\"" "" rhs ${rhs}) # remove ", if any
                        # remove tab or white space from head or tail
                        string(STRIP "${rhs}" rhs)

                        # ... and if it is not a default path in Unix...
                        if (NOT (rhs STREQUAL /usr OR
                                 rhs STREQUAL /usr/local))
                            # ... then set the one defined in makefile-default-paths
                            set("${lhs}" ${rhs})
                            message(STATUS "Specify ${lhs}: " ${rhs})
                        endif ()
                    endif ()

                    unset(lhs)
                    unset(rhs)
                endif ()
            endif ()
        endforeach ()

        unset(path)
        unset(paths)
    endmacro ()

    # firstly, set the custom makefile-paths, if any...
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/extlib/makefile-paths")
        makefile_macros_to_cmake(extlib/makefile-paths)
    endif ()
    # ... then set the default other ones ignoring those already previously set,
    # according to the specific OS
    if (UNIX)
        if (APPLE)
            makefile_macros_to_cmake(extlib/makefile-default-paths-macos)
        else ()
            makefile_macros_to_cmake(extlib/makefile-default-paths-linux)
        endif ()
    elseif (WIN32)
        makefile_macros_to_cmake(extlib/makefile-default-paths-win)
    endif ()

    # Retrieve the system architecture
    if (CMAKE_SIZEOF_VOID_P EQUAL 8)
        if (UNIX)
            set(ARCH x86-64)
        elseif (WIN32)
            set(ARCH x64)
        endif ()
    else ()
        set(ARCH x86)
    endif ()
endif ()

if (MCFClass_USE_CPLEX)
    # Needed for MCFCplex
    find_package(CPLEX)
    if (NOT CPLEX_FOUND)
        set(MCFClass_USE_CPLEX OFF)
    endif ()
endif ()

# ----- Configuration header ------------------------------------------------ #
# This will generate a *Config.h header in the build directory.
configure_file(cmake/${modName}Config.h.in ${modName}Config.h)

# ----- Library ------------------------------------------------------------- #
# With the BUILD_SHARED_LIBS variable we can specify if the library will
# be STATIC or SHARED, so no reason to do it now.
add_library(${modName})
target_compile_features(${modName} PUBLIC cxx_std_20)

# When adding source files with target_sources(), PRIVATE means that the files
# should only be added to this library, whereas PUBLIC means they should be
# added to this library and to any target that links to it.
# INTERFACE can be used for sources that should not be added to this library
# but should be added to anything that links to it.
# Note: do not GLOB files here.
target_sources(${modName} PRIVATE
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
        ${modName} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MCFClass>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MCFClone>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MCFSimplex>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/OPTUtils>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/RelaxIV>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SPTree>
        $<INSTALL_INTERFACE:include/${modNamespace}>)

# When linking other targets to the library with target_link_libraries(),
# PRIVATE means that the libraries will be linked only to this library,
# PUBLIC means they will be linked also to the targets that depend on this
# library, INTERFACE means they will be linked only to the targets that depend
# on this library.
if (MCFClass_USE_CPLEX)
    message(STATUS "${modName}: CPLEX found, configuring MCFCplex")

    target_sources(${modName} PRIVATE MCFCplex/MCFCplex.C)
    target_include_directories(
            ${modName} PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MCFCplex>)
    target_link_libraries(${modName} PUBLIC CPLEX::Cplex)
endif ()

# This alias is defined so that executables in this same project can use
# the library with this notation.
add_library(${modNamespace}::${modName} ALIAS ${modName})

# ----- Subdirectories ------------------------------------------------------ #
if (BUILD_TESTING)
    # Just consuming BUILD_TESTING to avoid warnings
endif ()

# ----- Install instructions ------------------------------------------------ #
# The following commands are used when installing the library
# and its CMake configuration files on the system.
# They are not required for local builds (see below).

include(GNUInstallDirs)

# Install the library
install(TARGETS ${modName}
        EXPORT ${modName}Targets
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
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

if (MCFClass_USE_CPLEX)
    install(FILES
            ${CMAKE_CURRENT_SOURCE_DIR}/MCFCplex/MCFCplex.h
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})
endif ()

# Install the auto-generated configuration header (see above).
install(FILES ${PROJECT_BINARY_DIR}/${modName}Config.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

# (Generate and) install the target import file, that allows other
# CMake projects to import the target.
install(EXPORT ${modName}Targets
        NAMESPACE ${modNamespace}::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Generate the package version file, that allows other
# CMake projects to know the version.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate the package configuration file, that allows other
# CMake projects to find the library with find_package().
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/${modName}Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Install the package version and configuration files.
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Shipping custom find modules should not be encouraged because third
# party developers should provide for them.
file(GLOB findMods ${CMAKE_CURRENT_LIST_DIR}/cmake/Find*)
install(FILES ${findMods}
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Install the README and LICENSE files.
install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/README.md
        ${CMAKE_CURRENT_LIST_DIR}/LICENSE
        DESTINATION ${CMAKE_INSTALL_DATADIR}/${modName})

# --------------------------------------------------------------------------- #
# Remove from the search paths the cmake directory we added at the beginning.
list(REMOVE_AT CMAKE_MODULE_PATH -1)

# --------------------------------------------------------------------------- #

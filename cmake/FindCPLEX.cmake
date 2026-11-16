# --------------------------------------------------------------------------- #
#    CMake find module for CPLEX                                              #
#                                                                             #
#    This module finds CPLEX include directories and libraries.               #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(CPLEX [version] [EXACT] [REQUIRED])                     #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        CPLEX_FOUND         - True if headers are found                      #
#        CPLEX_INCLUDE_DIRS  - Include directories                            #
#        CPLEX_LIBRARIES     - Libraries to be linked                         #
#        CPLEX_DLL           - The found runtime DLL, where there is one      #
#        CPLEX_VERSION       - Version number                                 #
#                                                                             #
#    This module reads hints about search locations from variables:           #
#                                                                             #
#        CPLEX_ROOT          - Custom path to CPLEX                           #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        CPLEX::Cplex                                                         #
#                                                                             #
#    This find module is provided because CPLEX does not provide              #
#    a CMake configuration file on its own.                                   #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- ILOG directories and lib suffixes ----------------------------------- #
# The directories where the installers of the systems put the ILOG products,
# and the directories of the library in the layouts of the platforms CPLEX
# is built for, the release and the debug one where they differ (the
# libraries for the MSVC runtime): only those of the system at hand are there

set(CPLEX_ILOG_DIRS /opt/ibm/ILOG /opt/IBM/ILOG /Applications
        "C:/IBM/ILOG" "C:/Program Files/IBM/ILOG")
set(CPLEX_LIB_PATH_SUFFIXES
        lib/x86-64_linux/static_pic
        lib/x86-64_sles10_4.1/static_pic
        lib/arm64_osx/static_pic
        lib/x86-64_osx/static_pic
        lib/x86-64_darwin9_gcc4.0/static_pic
        lib/x64_windows_msvc14/stat_mda
        lib/x64_windows_vs2017/stat_mda)
set(CPLEX_LIB_PATH_SUFFIXES_DEBUG
        lib/x64_windows_msvc14/stat_mdd
        lib/x64_windows_vs2017/stat_mdd)

# ----- Find the path to CPLEX ---------------------------------------------- #
# This takes the greatest CPLEX_Studio* found in the ILOG directories

foreach (dir ${CPLEX_ILOG_DIRS})
    file(GLOB CPLEX_DIRS "${dir}/CPLEX_Studio*")
    if (NOT IS_DIRECTORY "${CPLEX_ROOT}")
        if (NOT "${CPLEX_ROOT}" STREQUAL "")
            message(STATUS "Specified CPLEX: ${CPLEX_ROOT} not found")
        endif ()
        list(SORT CPLEX_DIRS)
        list(REVERSE CPLEX_DIRS)
        if (CPLEX_DIRS)
            list(GET CPLEX_DIRS 0 CPLEX_ROOT)
            message(STATUS "Using CPLEX: ${CPLEX_ROOT}")
            break()
        else ()
            set(CPLEX_ROOT CPLEX_ROOT-NOTFOUND)
        endif ()
    else ()
        break()
    endif ()
endforeach ()

# ----- Requirements -------------------------------------------------------- #
# This sets the variable CMAKE_THREAD_LIBS_INIT, see:
# https://cmake.org/cmake/help/latest/module/FindThreads.html
find_package(Threads QUIET)

# Drop the cached paths if one of them does not exist any more, so that
# they are all searched again from the same place
set(_cplex_cached CPLEX_INCLUDE_DIR CPLEX_LIBRARY CPLEX_LIBRARY_DEBUG
        CPLEX_DLL CPLEX_VERSION)
foreach (_var ${_cplex_cached})
    foreach (_path ${${_var}})
        if (IS_ABSOLUTE "${_path}" AND NOT EXISTS "${_path}")
            set(_cplex_stale TRUE)
        endif ()
    endforeach ()
endforeach ()
if (_cplex_stale)
    message(STATUS "Cached CPLEX paths not found, searching again")
    foreach (_var ${_cplex_cached})
        unset(${_var} CACHE)
    endforeach ()
endif ()
unset(_cplex_cached)
unset(_cplex_stale)
unset(_var)
unset(_path)

# The cached version is the one of the include directory it was read from
if (NOT "${CPLEX_INCLUDE_DIR}" STREQUAL
        "${_CPLEX_VERSION_INCLUDE_DIR}")
    unset(CPLEX_VERSION CACHE)
endif ()

# Check if already in cache
if (CPLEX_INCLUDE_DIR AND CPLEX_LIBRARY AND CPLEX_LIBRARY_DEBUG AND
        CPLEX_VERSION)
    set(CPLEX_FOUND TRUE)
endif ()

if (NOT CPLEX_FOUND)

    set(CPLEX_DIR ${CPLEX_ROOT}/cplex)

    # ----- Find the CPLEX include directory -------------------------------- #
    # The libraries are taken from CPLEX_DIR, hence so is the include
    # directory
    unset(CPLEX_INCLUDE_DIR CACHE)
    find_path(CPLEX_INCLUDE_DIR
            NAMES ilcplex/cplex.h
            PATHS ${CPLEX_DIR}/include
            DOC "CPLEX include directory.")

    # ----- Parse the version ----------------------------------------------- #
    if (CPLEX_INCLUDE_DIR)
        file(STRINGS "${CPLEX_INCLUDE_DIR}/ilcplex/cpxconst.h"
                _cplex_version_lines
                REGEX "#define CPX_VERSION_(VERSION|RELEASE|MODIFICATION)")
        foreach (_part VERSION RELEASE MODIFICATION)
            string(REGEX REPLACE ".*CPX_VERSION_${_part} *\([0-9]*\).*"
                    "\\1" _cplex_${_part} "${_cplex_version_lines}")
        endforeach ()
        set(CPLEX_VERSION
                "${_cplex_VERSION}.${_cplex_RELEASE}.${_cplex_MODIFICATION}")
        # the name of the library, without or with the version (e.g.,
        # cplex2211 for 22.1.1)
        set(_cplex_names cplex
                cplex${_cplex_VERSION}${_cplex_RELEASE}${_cplex_MODIFICATION}
                cplex${_cplex_VERSION}${_cplex_RELEASE}0)
        unset(_cplex_version_lines)
        unset(_cplex_VERSION)
        unset(_cplex_RELEASE)
        unset(_cplex_MODIFICATION)
        unset(_part)
    endif ()

    # ----- Find the CPLEX library ------------------------------------------ #
    # The debug library is the release one where there is no other; where
    # the library is a DLL, the one found in lib/ is its import library, and
    # the DLL is in the bin/ of its layout
    if (_cplex_names)
        find_library(CPLEX_LIBRARY
                NAMES ${_cplex_names}
                PATHS ${CPLEX_DIR}
                PATH_SUFFIXES ${CPLEX_LIB_PATH_SUFFIXES}
                DOC "CPLEX library.")
        find_library(CPLEX_LIBRARY_DEBUG
                NAMES ${_cplex_names}
                PATHS ${CPLEX_DIR}
                PATH_SUFFIXES ${CPLEX_LIB_PATH_SUFFIXES_DEBUG}
                NO_DEFAULT_PATH
                DOC "CPLEX debug library.")
        if (NOT CPLEX_LIBRARY_DEBUG)
            set(CPLEX_LIBRARY_DEBUG ${CPLEX_LIBRARY}
                    CACHE FILEPATH "CPLEX debug library." FORCE)
        endif ()

        list(TRANSFORM _cplex_names PREPEND "${CMAKE_SHARED_LIBRARY_PREFIX}")
        list(TRANSFORM _cplex_names APPEND "${CMAKE_SHARED_LIBRARY_SUFFIX}")
        list(REMOVE_AT _cplex_names 0)
        find_file(CPLEX_DLL
                NAMES ${_cplex_names}
                PATHS ${CPLEX_DIR}/bin/x64_win64
                NO_DEFAULT_PATH
                DOC "CPLEX runtime DLL.")
        unset(_cplex_names)
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET, REQUIRED and version-related
    # options passed to find_package(). It also sets <PackageName>_FOUND if
    # REQUIRED_VARS are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/
    # FindPackageHandleStandardArgs.html
    find_package_handle_standard_args(
            CPLEX
            REQUIRED_VARS CPLEX_LIBRARY CPLEX_LIBRARY_DEBUG CPLEX_INCLUDE_DIR
            VERSION_VAR CPLEX_VERSION)
endif ()

# A found package keeps its version in the cache, with the include
# directory it was read from
if (CPLEX_FOUND)
    set(CPLEX_VERSION "${CPLEX_VERSION}" CACHE INTERNAL
            "CPLEX version.")
    set(_CPLEX_VERSION_INCLUDE_DIR "${CPLEX_INCLUDE_DIR}"
            CACHE INTERNAL "Include directory of the CPLEX version.")
endif ()

# ----- Export the target --------------------------------------------------- #
# A DLL is imported with its import library, any other library as it is
if (CPLEX_FOUND)
    set(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR})
    set(CPLEX_LIBRARIES ${CMAKE_THREAD_LIBS_INIT} ${CMAKE_DL_LIBS})

    if (NOT TARGET CPLEX::Cplex)
        if (CPLEX_DLL)
            add_library(CPLEX::Cplex SHARED IMPORTED)
            set_target_properties(
                    CPLEX::Cplex PROPERTIES
                    IMPORTED_IMPLIB "${CPLEX_LIBRARY}"
                    IMPORTED_IMPLIB_DEBUG "${CPLEX_LIBRARY_DEBUG}"
                    IMPORTED_LOCATION "${CPLEX_DLL}"
                    IMPORTED_LOCATION_DEBUG "${CPLEX_DLL}")
        else ()
            add_library(CPLEX::Cplex UNKNOWN IMPORTED)
            set_target_properties(
                    CPLEX::Cplex PROPERTIES
                    IMPORTED_LOCATION "${CPLEX_LIBRARY}"
                    IMPORTED_LOCATION_DEBUG "${CPLEX_LIBRARY_DEBUG}")
        endif ()
        set_target_properties(
                CPLEX::Cplex PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${CPLEX_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${CPLEX_LIBRARIES}")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(CPLEX_INCLUDE_DIR
        CPLEX_LIBRARY
        CPLEX_LIBRARY_DEBUG
        CPLEX_DLL
        CPLEX_VERSION)

# --------------------------------------------------------------------------- #

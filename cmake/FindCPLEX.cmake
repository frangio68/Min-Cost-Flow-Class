# --------------------------------------------------------------------------- #
#    CMake find module for CPLEX Studio                                       #
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
#        CPLEX_VERSION       - Version number                                 #
#                                                                             #
#    This module reads hints about search locations from variables:           #
#                                                                             #
#        CPLEX_STUDIO_DIR    - Custom path to CPLEX Studio                    #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        CPLEX::Cplex                                                         #
#                                                                             #
#    This find module is provided because CPLEX does not provide              #
#    a CMake configuration file on its own.                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Find ILOG directories and lib suffixes ------------------------------ #
# Based on the OS and architecture, generate:
# - a list of possible ILOG directories
# - a list of possible lib suffixes to find the library

if (APPLE)
    # macOS (usually /Applications)
    set(CPLEX_ILOG_DIRS /Applications)
    set(CPLEX_LIB_PATH_SUFFIXES
            lib/${ARCH}_darwin9_gcc4.0/static_pic
            lib/${ARCH}_osx/static_pic)
else ()
    # Other Unix-based systems (usually /opt/ibm/ILOG)
    set(CPLEX_ILOG_DIRS /opt/ibm/ILOG /opt/IBM/ILOG)
    set(CPLEX_LIB_PATH_SUFFIXES
            lib/${ARCH}_sles10_4.1/static_pic
            lib/${ARCH}_linux/static_pic)
endif ()

# ----- Find the path to CPLEX Studio --------------------------------------- #
# This takes the greatest CPLEX_Studio* found in the ILOG directories
# TODO: Sort properly, now 129 is considered > than 1210

if (NOT CPLEX_STUDIO_DIR)
    foreach (dir ${CPLEX_ILOG_DIRS})
        file(GLOB CPLEX_STUDIO_DIRS "${dir}/CPLEX_Studio*")
        list(SORT CPLEX_STUDIO_DIRS)
        list(REVERSE CPLEX_STUDIO_DIRS)
        if (CPLEX_STUDIO_DIRS)
            list(GET CPLEX_STUDIO_DIRS 0 CPLEX_STUDIO_DIR_)
            message(STATUS "Found CPLEX Studio: ${CPLEX_STUDIO_DIR_}")
            break()
        endif ()
    endforeach ()

    if (NOT CPLEX_STUDIO_DIR_)
        set(CPLEX_STUDIO_DIR_ CPLEX_STUDIO_DIR-NOTFOUND)
    endif ()
    # Set the path in the cache
    set(CPLEX_STUDIO_DIR ${CPLEX_STUDIO_DIR_})
endif ()

# ----- Requirements -------------------------------------------------------- #
# This sets the variable CMAKE_THREAD_LIBS_INIT, see:
# https://cmake.org/cmake/help/latest/module/FindThreads.html
find_package(Threads QUIET)

# Check if already in cache
if (CPLEX_INCLUDE_DIR AND CPLEX_LIBRARY AND CPLEX_LIBRARY_DEBUG)
    set(CPLEX_FOUND TRUE)
else ()

    # ----- Find the CPLEX include directory -------------------------------- #
    set(CPLEX_DIR ${CPLEX_STUDIO_DIR}/cplex)
    # Note that find_path() creates a cache entry
    find_path(CPLEX_INCLUDE_DIR
              NAMES ilcplex/cplex.h
              PATHS ${CPLEX_DIR}/include
              DOC "CPLEX include directory.")

    # ----- Find the CPLEX library ------------------------------------------ #
    # Note that find_library() creates a cache entry
    find_library(CPLEX_LIBRARY
                 NAMES cplex
                 PATHS ${CPLEX_DIR}
                 PATH_SUFFIXES ${CPLEX_LIB_PATH_SUFFIXES}
                 DOC "CPLEX library.")
    set(CPLEX_LIBRARY_DEBUG ${CPLEX_LIBRARY})

    # ----- Parse the version ----------------------------------------------- #
    if (CPLEX_INCLUDE_DIR)
        file(STRINGS
             "${CPLEX_INCLUDE_DIR}/ilcplex/cpxconst.h"
             _cplex_version_lines REGEX "#define CPX_VERSION_(VERSION|RELEASE|MODIFICATION)")

        string(REGEX REPLACE ".*CPX_VERSION_VERSION *\([0-9]*\).*" "\\1" _cplex_version_major "${_cplex_version_lines}")
        string(REGEX REPLACE ".*CPX_VERSION_RELEASE *\([0-9]*\).*" "\\1" _cplex_version_minor "${_cplex_version_lines}")
        string(REGEX REPLACE ".*CPX_VERSION_MODIFICATION *\([0-9]*\).*" "\\1" _cplex_version_patch "${_cplex_version_lines}")

        set(CPLEX_VERSION "${_cplex_version_major}.${_cplex_version_minor}.${_cplex_version_patch}")
        unset(_cplex_version_lines)
        unset(_cplex_version_major)
        unset(_cplex_version_minor)
        unset(_cplex_version_patch)
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET, REQUIRED and version-related
    # options passed to find_package(). It also sets <PackageName>_FOUND if
    # REQUIRED_VARS are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
    find_package_handle_standard_args(
            CPLEX
            REQUIRED_VARS CPLEX_LIBRARY CPLEX_LIBRARY_DEBUG CPLEX_INCLUDE_DIR
            VERSION_VAR CPLEX_VERSION)
endif ()

# ----- Export the target --------------------------------------------------- #
if (CPLEX_FOUND)
    set(CPLEX_INCLUDE_DIRS "${CPLEX_INCLUDE_DIR}")
    set(CPLEX_LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

    # See: https://cmake.org/cmake/help/latest/module/CheckLibraryExists.html
    check_library_exists(m floor "" HAVE_LIBM)
    if (HAVE_LIBM)
        set(CPLEX_LINK_LIBRARIES ${CPLEX_LINK_LIBRARIES} m)
    endif ()

    if (UNIX)
        # Required under Unix since 12.8
        set(CPLEX_LINK_LIBRARIES ${CPLEX_LINK_LIBRARIES} dl)
    endif ()

    if (NOT TARGET CPLEX::Cplex)
        add_library(CPLEX::Cplex STATIC IMPORTED)
        set_target_properties(
                CPLEX::Cplex PROPERTIES
                IMPORTED_LOCATION "${CPLEX_LIBRARY}"
                IMPORTED_LOCATION_DEBUG "${CPLEX_LIBRARY_DEBUG}"
                INTERFACE_INCLUDE_DIRECTORIES "${CPLEX_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "${CPLEX_LINK_LIBRARIES}")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(CPLEX_INCLUDE_DIR
                 CPLEX_LIBRARY
                 CPLEX_LIBRARY_DEBUG
                 CPLEX_VERSION)

# --------------------------------------------------------------------------- #

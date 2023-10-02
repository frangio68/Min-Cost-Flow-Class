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
#                              Niccolo' Iardella                              #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Find ILOG directories and lib suffixes ------------------------------ #
# Based on the OS and architecture, generate:
# - a list of possible ILOG directories
# - a list of possible lib suffixes to find the library

if (UNIX)
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
else ()
    # Windows (usually C:/Program Files/IBM/ILOG)
    set(CPLEX_ILOG_DIRS "C:/Program Files/IBM/ILOG")
    if (ARCH STREQUAL "x86")
        set(CPLEX_ILOG_DIRS
                "C:/Program Files (x86)/IBM/ILOG" ${CPLEX_ILOG_DIRS})
    endif ()

    # Amended for VS and its various toolsets
    # https://cmake.org/cmake/help/v3.11/variable/MSVC_VERSION.html
    # Can use GREATER_EQUAL instead of the mess below if cmake version >= 3.7
    #if (NOT (MSVC_VERSION LESS 1930)) # VS 17.0 (v143 toolset)
    #    set(CPLEX_LIB_PATH_SUFFIXES
    #            lib/${ARCH}_windows_msvc17/stat_mda)
    #    set(CPLEX_LIB_PATH_SUFFIXES_DEBUG
    #            lib/${ARCH}_windows_msvc17/stat_mdd)
    #elseif (NOT (MSVC_VERSION LESS 1920)) # VS 16.0 (v142 toolset)
    #    set(CPLEX_LIB_PATH_SUFFIXES
    #            lib/${ARCH}_windows_msvc16/stat_mda)
    #    set(CPLEX_LIB_PATH_SUFFIXES_DEBUG
    #            lib/${ARCH}_windows_msvc16/stat_mdd)
    #elseif (NOT (MSVC_VERSION LESS 1910)) # VS 15.0 (v141 toolset)
    #    set(CPLEX_LIB_PATH_SUFFIXES
    #            lib/${ARCH}_windows_msvc15/stat_mda)
    #    set(CPLEX_LIB_PATH_SUFFIXES_DEBUG
    #            lib/${ARCH}_windows_msvc15/stat_mdd)
    #else
    if (NOT (MSVC_VERSION LESS 1900)) # VS 14.0 (v140 toolset)
        set(CPLEX_LIB_PATH_SUFFIXES
                lib/${ARCH}_windows_msvc14/stat_mda)
        set(CPLEX_LIB_PATH_SUFFIXES_DEBUG
                lib/${ARCH}_windows_msvc14/stat_mdd)
    elseif (NOT (MSVC_VERSION LESS 1800)) # VS 12.0 (v120 toolset)
        set(CPLEX_LIB_PATH_SUFFIXES
                lib/${ARCH}_windows_msvc12/stat_mda)
        set(CPLEX_LIB_PATH_SUFFIXES_DEBUG
                lib/${ARCH}_windows_msvc12/stat_mdd)
    elseif (NOT (MSVC_VERSION LESS 1700)) # VS 11.0 (v110 toolset)
        set(CPLEX_LIB_PATH_SUFFIXES
                lib/${ARCH}_windows_msvc11/stat_mda)
        set(CPLEX_LIB_PATH_SUFFIXES_DEBUG
                lib/${ARCH}_windows_msvc11/stat_mdd)
    elseif (NOT (MSVC_VERSION LESS 1600)) # VS 10.0 (v100 toolset)
        set(CPLEX_LIB_PATH_SUFFIXES
                lib/${ARCH}_windows_msvc10/stat_mda)
        set(CPLEX_LIB_PATH_SUFFIXES_DEBUG
                lib/${ARCH}_windows_msvc10/stat_mdd)
    elseif (NOT (MSVC_VERSION LESS 1500)) # VS  9.0 (v90 toolset)
        set(CPLEX_LIB_PATH_SUFFIXES
                lib/${ARCH}_windows_msvc09/stat_mda)
        set(CPLEX_LIB_PATH_SUFFIXES_DEBUG
                lib/${ARCH}_windows_msvc09/stat_mdd)
    endif ()
endif ()

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

# Check if already in cache
if (CPLEX_INCLUDE_DIR AND CPLEX_LIBRARY AND CPLEX_LIBRARY_DEBUG)
    set(CPLEX_FOUND TRUE)
else ()

    set(CPLEX_DIR ${CPLEX_ROOT}/cplex)

    # ----- Find the CPLEX include directory -------------------------------- #
    # Note that find_path() creates a cache entry
    find_path(CPLEX_INCLUDE_DIR
              NAMES ilcplex/cplex.h
              PATHS ${CPLEX_DIR}/include
              DOC "CPLEX include directory.")

    if (UNIX)
        # ----- Find the CPLEX library -------------------------------------- #
        # Note that find_library() creates a cache entry
        find_library(CPLEX_LIBRARY
                     NAMES cplex
                     PATHS ${CPLEX_DIR}
                     PATH_SUFFIXES ${CPLEX_LIB_PATH_SUFFIXES}
                     DOC "CPLEX library.")
        set(CPLEX_LIBRARY_DEBUG ${CPLEX_LIBRARY})

    elseif (NOT CPLEX_LIBRARY)

        # ----- Macro: find_win_cplex_library ------------------------------- #
        # On Windows the version is appended to the library name which cannot be
        # handled by find_library, so here a macro to search manually.
        macro(find_win_cplex_library var path_suffixes)
            foreach (s ${path_suffixes})
                file(GLOB CPLEX_LIBRARY_CANDIDATES "${CPLEX_DIR}/${s}/cplex*.lib")
                if (CPLEX_LIBRARY_CANDIDATES)
                    list(GET CPLEX_LIBRARY_CANDIDATES 0 ${var})
                    break()
                endif ()
            endforeach ()
            if (NOT ${var})
                set(${var} NOTFOUND)
            endif ()
        endmacro ()

        # Library
        find_win_cplex_library(CPLEX_LIB "${CPLEX_LIB_PATH_SUFFIXES}")
        set(CPLEX_LIBRARY ${CPLEX_LIB})

        # Debug library
        find_win_cplex_library(CPLEX_LIB "${CPLEX_LIB_PATH_SUFFIXES_DEBUG}")
        set(CPLEX_LIBRARY_DEBUG ${CPLEX_LIB})

        # DLL
        if (CPLEX_LIBRARY MATCHES ".*/(cplex.*)\\.lib")
            file(GLOB CPLEX_DLL_ "${CPLEX_DIR}/bin/*/${CMAKE_MATCH_1}.dll")
            set(CPLEX_DLL ${CPLEX_DLL_})
        endif ()
    endif ()

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

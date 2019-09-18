cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_VERBOSE_MAKEFILE OFF)

# ----- Project description and directories ------------------------------------
project(MCFClass DESCRIPTION "C++ solvers for Minimum Cost Flow Problems")

# ----- Requirements -----------------------------------------------------------
# This is to find FindCPLEX.cmake in current folder
# TODO: Move FindCPLEX.cmake in a proper directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}")
find_package(CPLEX REQUIRED)

# ----- Includes ---------------------------------------------------------------
set(PROJECT_INCLUDE_DIR
    ${PROJECT_SOURCE_DIR}/MCFClass
    ${PROJECT_SOURCE_DIR}/MCFClone
    ${PROJECT_SOURCE_DIR}/MCFCplex
    ${PROJECT_SOURCE_DIR}/MCFSimplex
    ${PROJECT_SOURCE_DIR}/OPTUtils
    ${PROJECT_SOURCE_DIR}/RelaxIV
    ${PROJECT_SOURCE_DIR}/SPTree
    )
include_directories(${PROJECT_INCLUDE_DIR})
include_directories(${CPLEX_INCLUDE_DIR}/ilcplex)

# ----- SOURCE FILES -----------------------------------------------------------
set(${PROJECT_NAME}_SOURCES
    MCFCplex/MCFCplex.C
    MCFSimplex/MCFSimplex.C
    RelaxIV/RelaxIV.C
    SPTree/SPTree.C
    )

add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES})

# ----- Install (optional) -----------------------------------------------------
set(${PROJECT_NAME}_HEADERS
    ${PROJECT_SOURCE_DIR}/MCFClass/MCFClass.h
    ${PROJECT_SOURCE_DIR}/MCFClone/MCFClone.h
    ${PROJECT_SOURCE_DIR}/MCFCplex/MCFCplex.h
    ${PROJECT_SOURCE_DIR}/MCFSimplex/MCFSimplex.h
    ${PROJECT_SOURCE_DIR}/OPTUtils/OPTUtils.h
    ${PROJECT_SOURCE_DIR}/RelaxIV/RelaxIV.h
    ${PROJECT_SOURCE_DIR}/SPTree/SPTree.h
    )

set_target_properties(${PROJECT_NAME} PROPERTIES
                      PUBLIC_HEADER "${${PROJECT_NAME}_HEADERS}"
                      )

install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        PUBLIC_HEADER DESTINATION include/${PROJECT_NAME}
        )

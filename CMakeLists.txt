cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_VERBOSE_MAKEFILE OFF)

project(MCFClass DESCRIPTION "C++ solvers for Minimum Cost Flow Problems")

# This is to find FindCPLEX.cmake in current folder
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}")

# Set include directories
set(PROJECT_INCLUDE_DIR
    ${PROJECT_SOURCE_DIR}/MCFClass
    ${PROJECT_SOURCE_DIR}/MCFClone
    ${PROJECT_SOURCE_DIR}/MCFCplex
    ${PROJECT_SOURCE_DIR}/MCFSimplex
    ${PROJECT_SOURCE_DIR}/OPTUtils
    ${PROJECT_SOURCE_DIR}/RelaxIV
    )

# Set source files
set(MCFClass_SOURCES
    MCFCplex/MCFCplex.C
    MCFSimplex/MCFSimplex.C
    RelaxIV/RelaxIV.C
    SPTree/SPTree.C
    )

# Find CPLEX

find_package(CPLEX)
include_directories(${CPLEX_INCLUDE_DIR}/ilcplex)

include_directories("${PROJECT_INCLUDE_DIR}")
add_library(${PROJECT_NAME} STATIC ${MCFClass_SOURCES})

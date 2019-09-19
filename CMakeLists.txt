cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_VERBOSE_MAKEFILE OFF)

# ----- Project description and directories ------------------------------------
project(MCFClass DESCRIPTION "C++ solvers for Minimum Cost Flow Problems")

set(PROJECT_INCLUDE_DIR
    ${PROJECT_SOURCE_DIR}/MCFClass
    ${PROJECT_SOURCE_DIR}/MCFClone
    ${PROJECT_SOURCE_DIR}/MCFCplex
    ${PROJECT_SOURCE_DIR}/MCFSimplex
    ${PROJECT_SOURCE_DIR}/OPTUtils
    ${PROJECT_SOURCE_DIR}/RelaxIV
    ${PROJECT_SOURCE_DIR}/SPTree
    )

# This is to find FindCPLEX.cmake in cmake folder
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ----- Requirements -----------------------------------------------------------
#find_package(CPLEX REQUIRED)

# ----- Library ----------------------------------------------------------------
set(${PROJECT_NAME}_SOURCES
    MCFCplex/MCFCplex.C
    MCFSimplex/MCFSimplex.C
    RelaxIV/RelaxIV.C
    SPTree/SPTree.C
    )

add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_INCLUDE_DIR})
#target_include_directories(${PROJECT_NAME} PUBLIC ${CPLEX_INCLUDE_DIR}/ilcplex)

# FIXME: Workaround for CPLEX
if (OVERRIDE_CPLEX_INC)
    target_include_directories(${PROJECT_NAME} PUBLIC ${OVERRIDE_CPLEX_INC}/ilcplex)
else ()
    find_package(CPLEX REQUIRED)
    target_include_directories(${PROJECT_NAME} PUBLIC ${CPLEX_INCLUDE_DIR}/ilcplex)
endif ()


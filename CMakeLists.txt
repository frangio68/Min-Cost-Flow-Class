# ----- Project description and directories ------------------------------------
cmake_minimum_required(VERSION 3.1)
project(MCFClass
        VERSION 1.0.0
        DESCRIPTION "C++ solvers for Minimum Cost Flow Problems")

# This is to find FindCPLEX.cmake in cmake folder
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ----- Requirements -----------------------------------------------------------
find_package(CPLEX REQUIRED)

# ----- Library ----------------------------------------------------------------
set(MCFClass_SOURCES
    MCFCplex/MCFCplex.C
    MCFSimplex/MCFSimplex.C
    RelaxIV/RelaxIV.C
    SPTree/SPTree.C
    )

# Add target
add_library(MCFClass STATIC ${MCFClass_SOURCES})

# Alias for testing
add_library(MCFClass::MCFClass ALIAS MCFClass)

# Target properties
target_compile_features(MCFClass PUBLIC cxx_std_17)
target_include_directories(MCFClass
                           PUBLIC
                           $<INSTALL_INTERFACE:include/MCFClass>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/MCFClass>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/MCFClone>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/MCFCplex>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/MCFSimplex>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/OPTUtils>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/RelaxIV>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/SPTree>
                           )

target_link_libraries(MCFClass PUBLIC cplex-library)

# ----- Install instructions ---------------------------------------------------
include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/MCFClass)

# Install library
install(TARGETS MCFClass
        EXPORT MCFClass-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

# Install headers
set(MCFClass_HEADERS
    ${PROJECT_SOURCE_DIR}/MCFClass/MCFClass.h
    ${PROJECT_SOURCE_DIR}/MCFClone/MCFClone.h
    ${PROJECT_SOURCE_DIR}/MCFCplex/MCFCplex.h
    ${PROJECT_SOURCE_DIR}/MCFSimplex/MCFSimplex.h
    ${PROJECT_SOURCE_DIR}/OPTUtils/OPTUtils.h
    ${PROJECT_SOURCE_DIR}/RelaxIV/RelaxIV.h
    ${PROJECT_SOURCE_DIR}/SPTree/SPTree.h
    )

install(FILES ${MCFClass_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MCFClass
        )

# Export target in *Targets.cmake
install(EXPORT MCFClass-targets
        FILE MCFClassTargets.cmake
        NAMESPACE MCFClass::
        DESTINATION ${INSTALL_CONFIGDIR}
        )

# Create *ConfigVersion.cmake
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)

# Create *Config.cmake
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/MCFClassConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfig.cmake
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

# Install the Config and ConfigVersion
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/MCFClassConfigVersion.cmake
        DESTINATION ${INSTALL_CONFIGDIR}
        )

# Install also FindCPLEX.cmake
install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPLEX.cmake
        DESTINATION ${INSTALL_CONFIGDIR}
        )

# ----- Export from the build tree ---------------------------------------------
configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPLEX.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/FindCPLEX.cmake
               COPYONLY)

export(EXPORT MCFClass-targets
       FILE ${CMAKE_CURRENT_BINARY_DIR}/MCFClassTargets.cmake
       NAMESPACE MCFClass::)

# Register package in user's package registry
export(PACKAGE MCFClass)

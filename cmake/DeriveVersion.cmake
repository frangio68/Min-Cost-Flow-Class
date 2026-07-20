# --------------------------------------------------------------------------- #
# DeriveVersion.cmake                                                          #
#                                                                              #
# Derive the module version from git, so that cutting a release is just        #
# `git tag x.y.z` with NO version number written by hand anywhere: the version #
# in project(), in the generated <Module>Config.h and in the tools' --version  #
# all follow the tag.                                                          #
#                                                                              #
# Resolution order (the first that yields a numeric x[.y[.z]] wins):           #
#   1. a working git clone: `git describe --tags --abbrev=0`;                  #
#   2. a source tarball (e.g. a GitLab release, NOT a git repo): the VERSION   #
#      file, which git fills in at `git archive` time via the                  #
#      `VERSION export-subst` .gitattributes entry (its content is             #
#      "$Format:%(describe:tags,abbrev=0)$", substituted with the tag by git); #
#   3. otherwise the sentinel "0.0.0" (no version info available).             #
#                                                                              #
# The leading x[.y[.z]] is extracted from whatever git returns, so a describe  #
# with a distance suffix ("1.2.3-4-gdeadbee", e.g. an archive of a non-tag     #
# commit) still yields "1.2.3" rather than the sentinel.                       #
#                                                                              #
# Call it BEFORE project():                                                    #
#     include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/DeriveVersion.cmake)           #
#     smspp_derive_version(MOD_VERSION)                                        #
#     project(<Module> VERSION ${MOD_VERSION} ...)                            #
# --------------------------------------------------------------------------- #

function(smspp_derive_version out_var)
    set(_raw "")

    # 1) working git clone
    find_package(Git QUIET)
    if(GIT_FOUND)
        execute_process(
                COMMAND "${GIT_EXECUTABLE}" describe --tags --abbrev=0
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                OUTPUT_VARIABLE _raw
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
                RESULT_VARIABLE _result)
        if(NOT _result EQUAL 0)
            set(_raw "")
        endif()
    endif()

    # 2) tarball: the VERSION file filled in by `git archive` (export-subst);
    #    in a non-archived copy it still holds the literal "$Format:...$", which
    #    just fails the numeric match below and falls through to the sentinel
    if(NOT _raw AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION")
        file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" _raw)
        string(STRIP "${_raw}" _raw)
    endif()

    # extract the leading x[.y[.z]] (tolerating a leading 'v' and a trailing
    # "-<n>-g<sha>" describe suffix); 3) sentinel if nothing matched
    if(_raw MATCHES "^v?([0-9]+(\\.[0-9]+)*)")
        set(_version "${CMAKE_MATCH_1}")
    else()
        set(_version "0.0.0")
    endif()

    set(${out_var} "${_version}" PARENT_SCOPE)
endfunction()

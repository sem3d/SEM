# Captures the current git hash + dirty status and expands sem_git_version.F90.in.
# Invoked as a pre-build custom command (see CMakeLists.txt) so it re-runs whenever
# .git/HEAD or .git/index change -- not just at cmake configure time.
find_package(Git QUIET)

set(GIT_HASH "unknown")
set(GIT_DIRTY "unknown")

if(GIT_FOUND AND EXISTS "${SEMDIR}/.git")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short=12 HEAD
        WORKING_DIRECTORY ${SEMDIR}
        OUTPUT_VARIABLE GIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE GIT_RESULT
    )
    if(NOT GIT_RESULT EQUAL 0 OR NOT GIT_HASH)
        set(GIT_HASH "unknown")
    endif()

    execute_process(
        COMMAND ${GIT_EXECUTABLE} status --porcelain
        WORKING_DIRECTORY ${SEMDIR}
        OUTPUT_VARIABLE GIT_DIRTY_OUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(GIT_DIRTY_OUT)
        set(GIT_DIRTY "dirty")
    else()
        set(GIT_DIRTY "clean")
    endif()
endif()

configure_file(${SEMCOMMON}/sem_git_version.F90.in ${DST} @ONLY)

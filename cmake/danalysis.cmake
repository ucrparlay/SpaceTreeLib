# -------------------------------------------------------------------
#           DAnalysis: Dan's static analysis framework
#
# - Supports cppcheck, clang-tidy, and include-what-you-use
# - Creates targets for individual files and for analyzing all files
#
# -------------------------------------------------------------------

# User options to enable static analysis tools
option(PSTP_ENABLE_CPPCHECK     "Enable cppcheck static analysis"      OFF)
option(PSTP_ENABLE_CLANG_TIDY   "Enable clang-tidy static analysis"    OFF)
option(PSTP_ENABLE_IWYU         "Enable include-what-you-use analysis" OFF)

function(CONFIGURE_DANALYSIS)
  set(options OPTIONAL "")
  set(oneValueArgs PROJECT_PREFIX)
  set(multiValueArgs INCLUDES FILES IGNORE_DEFINITIONS)
  cmake_parse_arguments(CONFIGURE_DANALYSIS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Deal with multiple include paths
  list(TRANSFORM CONFIGURE_DANALYSIS_INCLUDES PREPEND "-I")
  
  # Deal with ignored definitions
  list(TRANSFORM CONFIGURE_DANALYSIS_IGNORE_DEFINITIONS PREPEND "-U")

  set(STATIC_ANALYSIS_TARGETS "")

  # ---------------------------- cppcheck -----------------------------

  if(PSTP_ENABLE_CPPCHECK)
    # Find cppcheck installation
    find_program(
      CPPCHECK_EXE
      NAMES "cppcheck"
      DOC "Path to cppcheck executable"
    )

    if(NOT CPPCHECK_EXE)
      message(FATAL_ERROR "cppcheck enabled but could not find cppcheck executable")
    else()
      message(STATUS "cppcheck:                     Enabled (${CPPCHECK_EXE})")

      # Args for cppcheck. Enable all levels of errors
      set(CPPCHECK_ARGS
        "--enable=warning,performance,style,portability"
        "--suppress=missingIncludeSystem"
        "${CONFIGURE_DANALYSIS_INCLUDES}"
        "--error-exitcode=1"
        "--inline-suppr"
        "--suppressions-list=${CMAKE_SOURCE_DIR}/.cppcheck-suppress"
        "--language=c++"
        "--force"
        "--cppcheck-build-dir=${PROJECT_SOURCE_DIR}/build"
        "--check-level=exhaustive"
        "${CONFIGURE_DANALYSIS_IGNORE_DEFINITIONS}"
      )

      # Create a target for each library header that runs cppcheck on that file
      set(CPPHECK_TARGETS "")
      foreach(LIB_FILE ${CONFIGURE_DANALYSIS_FILES})
        string(REGEX REPLACE "/" "-" TARGET_NAME ${LIB_FILE})
        set(CPPCHECK_TARGET_NAME "cppcheck-${TARGET_NAME}")
        add_custom_target(
          ${CPPCHECK_TARGET_NAME}
          COMMAND ${CPPCHECK_EXE} ${CPPCHECK_ARGS}
          ${LIB_FILE}
        )
      set(CPPCHECK_TARGETS ${CPPCHECK_TARGETS} "${CPPCHECK_TARGET_NAME}")
      endforeach()

      # Custom target to run cppcheck on all library files sequentially
      add_custom_target(pstp_cppcheck)
      add_dependencies(pstp_cppcheck ${CPPCHECK_TARGETS})

      # Custom target to run cppcheck on all library files at once
      add_custom_target(
        pstp_cppcheck_all
        COMMAND ${CPPCHECK_EXE} ${CPPCHECK_ARGS}
        ${CONFIGURE_DANALYSIS_FILES}
      )

      set(STATIC_ANALYSIS_TARGETS ${STATIC_ANALYSIS_TARGETS} pstp_cppcheck_all)

    endif()
  else()
    message(STATUS "cppcheck:                     Disabled (enable with -DENABLE_CPPCHECK=On)")
  endif()

  # --------------------------- clang-tidy ----------------------------

  if(PSTP_ENABLE_CLANG_TIDY)
    # Find clang-tidy installation
    find_program(
      CLANG_TIDY_EXE
      NAMES "clang-tidy"
      DOC "Path to clang-tidy executable"
    )

    if(NOT CLANG_TIDY_EXE)
      message(FATAL_ERROR "clang-tidy enabled but could not find clang-tidy executable")
    else()
      message(STATUS "clang-tidy:                   Enabled (${CLANG_TIDY_EXE})")

      # Args for clang-tidy.
      set(CLANG_TIDY_ARGS
        "-p" 
        "${PROJECT_SOURCE_DIR}/build"
        "--fix-notes"
        "extra-arg=-std=c++20"
        "--fix-errors"
      )
      set(CLANG_TIDY_COMPILER_ARGS
        "-p"
        "${PROJECT_SOURCE_DIR}/build"
        "-xc++"
        "-std=c++20"
        "${CONFIGURE_DANALYSIS_INCLUDES}"
      )

      # Create a target for each library header that runs clang-tidy on that file
      set(CLANG_TIDY_TARGETS "")
      foreach(LIB_FILE ${CONFIGURE_DANALYSIS_FILES})
        string(REGEX REPLACE "/" "-" TARGET_NAME ${LIB_FILE})
        set(CLANG_TIDY_TARGET_NAME "clang-tidy-${TARGET_NAME}")
        add_custom_target(
          ${CLANG_TIDY_TARGET_NAME}
          COMMAND ${CLANG_TIDY_EXE} ${CLANG_TIDY_ARGS}
          ${LIB_FILE}
          --
          ${CLANG_TIDY_COMPILER_ARGS}
        )
        set(CLANG_TIDY_TARGETS ${CLANG_TIDY_TARGETS} "${CLANG_TIDY_TARGET_NAME}")
      endforeach()

      # Custom target to run clang-tidy on all library files sequentially
      add_custom_target(pstp_clang_tidy)
      add_dependencies(pstp_clang_tidy ${CLANG_TIDY_TARGETS})

      # Custom target to run pstp_clang-tidy on all library files at once
      add_custom_target(
        pstp_clang_tidy_all
        COMMAND ${CLANG_TIDY_EXE} ${CLANG_TIDY_ARGS}
        ${CONFIGURE_DANALYSIS_FILES}
        --
        ${CLANG_TIDY_COMPILER_ARGS}
      )

      set(STATIC_ANALYSIS_TARGETS ${STATIC_ANALYSIS_TARGETS} pstp_clang_tidy_all)

    endif()
  else()
    message(STATUS "clang-tidy:                   Disabled (enable with -DENABLE_CLANG_TIDY=On)")
  endif()

  # ----------------------- include-what-you-use --------------------------

  if(PSTP_ENABLE_IWYU)
    # Find include-what-you-use executable
    find_program(
      IWYU_EXE
      NAMES "include-what-you-use" "iwyu"
      DOC "Path to include-what-you-use executable"
    )
    
    if(NOT IWYU_EXE)
      message(FATAL_ERROR "include-what-you-use enabled but could not find include-what-you-use executable")
    else()
      message(STATUS "include-what-you-use:         Enabled (${IWYU_EXE})")
     
      # Args for IWYU
      set(IWYU_ARGS 
        "${CONFIGURE_DANALYSIS_INCLUDES}"
        "-w"
        "-std=c++17"
        "-xc++"
      )

      # Create a target for each library header that runs infer on that file
      # Notice the strange exit code? We return $? - 2 since IWYU returns
      # 2 + number of errors. Makes sense, right?
      set(IWYU_TARGETS "")
      foreach(LIB_FILE ${CONFIGURE_DANALYSIS_FILES})
        string(REGEX REPLACE "/" "-" TARGET_NAME ${LIB_FILE})
        set(IWYU_TARGET_NAME "iwyu-${TARGET_NAME}")
        add_custom_target(
          ${IWYU_TARGET_NAME}
          COMMAND ${IWYU_EXE} ${IWYU_ARGS}
          ${LIB_FILE}
          || exit `expr $$? - 2`
        )
        set(IWYU_TARGETS ${IWYU_TARGETS} "${IWYU_TARGET_NAME}")
      endforeach()

      # Custom target to run include-what-you-use on all library files sequentially
      add_custom_target(pstp_iwyu)
      add_dependencies(pstp_iwyu ${IWYU_TARGETS})
      
      # include-what-you-use doesn't support compiling multiple files at once
      add_custom_target(pstp_iwyu_all)
      add_dependencies(pstp_iwyu_all pstp_iwyu)

      set(STATIC_ANALYSIS_TARGETS ${STATIC_ANALYSIS_TARGETS} pstp_iwyu_all)

    endif()
  else()
    message(STATUS "include-what-you-use:         Disabled (enable with -DENABLE_IWYU=On)")
  endif()

  # ---------------------------- all linters ------------------------------

  # Target to run all enabled linters
  if(PSTP_ENABLE_CPPCHECK OR PSTP_ENABLE_CLANG_TIDY OR PSTP_ENABLE_IWYU)
    add_custom_target(pstp_analysis_all)
    add_dependencies(pstp_analysis_all ${STATIC_ANALYSIS_TARGETS})
  endif()
endfunction()

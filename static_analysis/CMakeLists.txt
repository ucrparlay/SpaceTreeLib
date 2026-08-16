include(danalysis)

# Globbed rather than listed: the hand-written list had gone stale, omitting
# every file under p_tree_impl/. dependence/cpam is a deferred fork and
# space_filling_curve is vendored, so both stay out.
file(GLOB_RECURSE psi_analysis_files CONFIGURE_DEPENDS
  ${PSI_INCLUDE_DIR}/psi/*.h
  ${PSI_INCLUDE_DIR}/psi/*.hpp
)
list(FILTER psi_analysis_files EXCLUDE REGEX "/dependence/(cpam|space_filling_curve)/")

configure_danalysis(
  INCLUDES
    ${PSI_INCLUDE_DIR}
  FILES
    ${psi_analysis_files}
  IGNORE_DEFINITIONS
    PARLAY_CILK
    PARLAY_CILKPLUS
    PARLAY_OPENCILK
    BSD
    __APPLE__
    __MACH__
    __unix__
    __unix
    unix
    _WIN32
    __WIN32__
    WIN32
    __CYGWIN__
    __CYGWIN32__
    HW_MEMSIZE
    HW_PHYSMEM64
    HW_MEMSIZE
    HW_REALMEM
    HW_PHYSMEM
    HW_PYSMEM
    CTL_HW
    _SC_AIX_REALMEM
    _SC_PHYS_PAGES
    _SC_PAGESIZE
    _SC_PAGE_SIZE
)

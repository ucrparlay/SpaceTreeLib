cmake_minimum_required(VERSION 3.15)

project(
	PSI
	VERSION 0.1
	DESCRIPTION "spatial partition tree"
	HOMEPAGE_URL ""
	LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # works
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --gcc-toolchain=/opt/rh/gcc-toolset-14/root/usr/")
# add_compile_options(-Wall -Wextra -Wshadow -Wconversion -Wunused -Wsign-conversion -Wfloat-equal -Wpedantic -Wuninitialized -Wformat -Wformat-security -Wmissing-braces -Wparentheses -Wsequence-point -Wswitch -Wswitch-enum -Wstrict-aliasing -Wdeprecated-declarations -Wreturn-type -Wunused-result)
# add_compile_options(-Wall -Wextra -Werror)
# NOTE: warnings are scoped per-target below, not globally.
# The three example targets instantiate all three trees, so -Wall -Wextra on
# them covers the whole library. tests/ and the baselines carry -w until they
# are cleaned separately (they contribute ~109 out-of-scope warning sites).
set(PSI_WARNINGS -Wall -Wextra)

option(DEBUG "" OFF)
option(CILKPLUS "" OFF)
option(OPENCILK "" OFF)
option(SERIAL "" OFF)
option(STDALLOC "" OFF)
option(CGAL "" OFF)
option(JEMA "" OFF)
option(CCP "" OFF)
option(DISABLE_BATCH_DELETE_SIZE_OPT "" OFF)

# PSI is header-only; a consumer wants the library, not the benchmark suite.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	set(psi_is_top_level ON)
else()
	set(psi_is_top_level OFF)
endif()
option(PSI_BUILD_BENCHMARKS "Build the test and baseline drivers" ${psi_is_top_level})
option(PSI_BUILD_TESTS "Build the unit tests" ${psi_is_top_level})
option(PSI_BUILD_EXAMPLES "Build the three example programs" ${psi_is_top_level})

# ON reproduces the flags every measurement was taken with. OFF gives a build
# that runs on a machine other than the one that compiled it.
option(PSI_NATIVE_ARCH "Tune for the building machine's CPU" ON)

# The library's only error mechanism is assert(), and the only configuration
# that kept them was -DDEBUG=ON, which is unoptimised. That made a checked run
# 10-50x slower than a normal one, which is why the correctness script was
# impractical to run. This keeps the optimised flags and the checks.
option(PSI_ASSERTS "Keep assert() in an optimised build" OFF)


if(CILKPLUS)
	add_compile_options(-DPARLAY_CILKPLUS -DCILK -fcilkplus)
elseif(OPENCILK)
	add_compile_options(-DPARLAY_OPENCILK -DCILK -fopencilk)
elseif(SERIAL)
	add_compile_options(-DPARLAY_SEQUENTIAL)
else()
	add_compile_options(-pthread)
endif()

# NOTE: this ablation switch is orthogonal to the parallel backend. It used to
# sit in the elseif chain above, which silently dropped -pthread whenever it
# was enabled -- so any A/B measurement of it also differed in that flag.
if(DISABLE_BATCH_DELETE_SIZE_OPT)
	add_compile_options(-DDISABLE_BATCH_DELETE_SIZE_OPT)
endif()

if(DEBUG)
	set(psi_build_type "Debug")
	add_compile_options(-g)
else()
	set(psi_build_type "Release")
	add_compile_options(-O3)
endif()

# Do not silently discard the standard knob: the cache used to say Debug while
# the compile line said -O3 -DNDEBUG.
if(NOT CMAKE_BUILD_TYPE)
	if(psi_is_top_level OR NOT CMAKE_CONFIGURATION_TYPES)
		set(CMAKE_BUILD_TYPE "${psi_build_type}")
	endif()
	# Only PSI's own directory gets -O3 above; a consumer's targets are
	# compiled with their own flags, and an unoptimised spatial index looks
	# broken rather than slow.
	if(NOT psi_is_top_level)
		message(STATUS
		  "PSI: consuming project set no CMAKE_BUILD_TYPE. Build your "
		  "targets with optimisation (-O3, and -march=native if you are "
		  "measuring) or PSI will be an order of magnitude slower.")
	endif()
elseif(NOT CMAKE_BUILD_TYPE STREQUAL psi_build_type)
	message(WARNING
	  "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} is being kept; -DDEBUG=${DEBUG} "
	  "would have selected ${psi_build_type}. Use -DDEBUG=ON for the "
	  "assert-enabled configuration.")
endif()

if(PSI_ASSERTS)
	string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELEASE
	       "${CMAKE_CXX_FLAGS_RELEASE}")
	string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO
	       "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
	add_compile_options(-UNDEBUG)
	message(STATUS "PSI: assertions kept in an optimised build")
endif()

if(PSI_NATIVE_ARCH)
	add_compile_options(-march=native)
	# -mcx16 is an x86 option; the old guard tested the OS, so Linux on
	# aarch64 could not configure at all.
	if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|i[3-6]86)$")
		add_compile_options(-mcx16)
	endif()
endif()

if(CCP)
	message(STATUS "CCP=ON")
	add_compile_options(-DCCP)
endif()

if(STDALLOC)
	add_compile_options(-DPARLAY_USE_STD_ALLOC)
elseif(JEMA)
	include_directories("/usr/local/lib")
	set(CMAKE_BUILD_RPATH "/usr/local/lib")
endif()

foreach(psi_submodule parlaylib libmorton)
	if(NOT EXISTS "${PSI_SOURCE_DIR}/include/${psi_submodule}/CMakeLists.txt")
		message(FATAL_ERROR
		  "include/${psi_submodule} is empty -- this is a git submodule.\n"
		  "Run: git submodule update --init --recursive")
	endif()
endforeach()

# The minimums the project is actually built and tested with.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14)
	message(FATAL_ERROR "PSI needs g++ 14 or newer (found ${CMAKE_CXX_COMPILER_VERSION})")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19)
	message(FATAL_ERROR "PSI needs clang 19 or newer (found ${CMAKE_CXX_COMPILER_VERSION})")
endif()

message(STATUS "\n\n---------------------------- PARLAY -----------------------------")
add_subdirectory(include/parlaylib)
message(STATUS "---------------------------- END PARLAY -----------------------------\n\n")

message(STATUS "\n\n---------------------------- LIBMORTON -----------------------------")
set(BUILD_TESTING OFF CACHE BOOL "Disable the libmorton test")
add_subdirectory(include/libmorton)
message(STATUS "---------------------------- END LIBMORTON -----------------------------\n\n")

# NOTE: clear CMAKE_MODULE_PATH to avoid finding parlay's cmake files
list(REMOVE_ITEM CMAKE_MODULE_PATH "${PARLAY_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${PSI_SOURCE_DIR}/cmake")

# NOTE: begin of psi library
message(STATUS "\n\n---------------------------- PSI -----------------------------")

# NOTE: PSI
add_library(PSI INTERFACE)
# One include root, one spelling: every header inside include/psi includes its
# siblings as "psi/...", so include/ is the only directory a consumer needs.
set(PSI_INCLUDE_DIR "${PSI_SOURCE_DIR}/include")
target_link_libraries(PSI INTERFACE parlay libmorton)
target_include_directories(PSI INTERFACE
  $<BUILD_INTERFACE:${PSI_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(PSI INTERFACE cxx_std_20)

add_library(PSI::PSI ALIAS PSI)

# Installable so a downstream project can find_package(PSI) instead of
# vendoring the whole repository.
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS PSI EXPORT PSITargets)
install(DIRECTORY "${PSI_SOURCE_DIR}/include/psi"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
install(EXPORT PSITargets
        FILE PSITargets.cmake
        NAMESPACE PSI::
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/PSI")

configure_package_config_file(
  "${PSI_SOURCE_DIR}/cmake/PSIConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/PSIConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/PSI")
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/PSIConfigVersion.cmake"
  COMPATIBILITY SameMajorVersion)
install(FILES "${PROJECT_BINARY_DIR}/PSIConfig.cmake"
              "${PROJECT_BINARY_DIR}/PSIConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/PSI")

add_subdirectory(static_analysis)
message(STATUS "---------------------------- END PSI -----------------------------\n\n")

if(PSI_BUILD_BENCHMARKS)
	add_executable(kd_test ${PROJECT_SOURCE_DIR}/tests/kd_test.cpp)
	target_link_libraries(kd_test PUBLIC PSI)

	add_executable(p_test ${PROJECT_SOURCE_DIR}/tests/p_test.cpp)
	target_link_libraries(p_test PUBLIC PSI)

	add_executable(baselines ${PROJECT_SOURCE_DIR}/tests/baselines.cpp)
	target_link_libraries(baselines PUBLIC PSI)

	if(JEMA)
		target_link_libraries(kd_test PUBLIC libjemalloc.so)
		target_link_libraries(p_test PUBLIC libjemalloc.so)
		target_link_libraries(baselines PUBLIC libjemalloc.so)
	endif()

	add_executable(data_generator ${PROJECT_SOURCE_DIR}/tests/data_generate.cpp)
	target_link_libraries(data_generator PUBLIC PSI)


	add_executable(data_washer ${PROJECT_SOURCE_DIR}/tests/data_wash.cpp)
	target_link_libraries(data_washer PUBLIC PSI)

	find_package(Boost REQUIRED)
	add_executable(boost_rtree ${PROJECT_SOURCE_DIR}/tests/boost_rtree.cpp)
	target_link_libraries(boost_rtree PUBLIC PSI ${Boost_LIBRARIES})
	target_include_directories(boost_rtree PUBLIC ${Boost_INCLUDE_DIRS})
	target_compile_options(boost_rtree PRIVATE -DBOOST_MATH_DISABLE_FLOAT128 -DBOOST_CSTDFLOAT_NO_LIBQUADMATH_SUPPORT)
endif()

if(PSI_BUILD_TESTS)
	enable_testing()
	add_subdirectory(tests/unit)
endif()

if(PSI_BUILD_EXAMPLES)
	add_executable(kd_tree_example ${PROJECT_SOURCE_DIR}/example/kd_tree_main.cpp)
	target_link_libraries(kd_tree_example PUBLIC PSI)
	target_include_directories(kd_tree_example PUBLIC ${PROJECT_SOURCE_DIR}/example)

	add_executable(orth_tree_example ${PROJECT_SOURCE_DIR}/example/orth_tree_main.cpp)
	target_link_libraries(orth_tree_example PUBLIC PSI)
	target_include_directories(orth_tree_example PUBLIC ${PROJECT_SOURCE_DIR}/example)

	add_executable(p_tree_example ${PROJECT_SOURCE_DIR}/example/p_tree_main.cpp)
	target_link_libraries(p_tree_example PUBLIC PSI)
	target_include_directories(p_tree_example PUBLIC ${PROJECT_SOURCE_DIR}/example)
endif()

# NOTE: the warning policy is applied at the very bottom of this file, after
# every target (including the CGAL-only kd_ccp / p_ccp) has been declared.

# CGAL
if(CGAL AND PSI_BUILD_BENCHMARKS)
	message(STATUS "\n\n---------------------------- CGAL -----------------------------")
	message(STATUS "To disable using -DCGAL=OFF")
	find_package(CGAL REQUIRED)

	# TBB
	find_package(TBB)

	# copied from CGAL default cmake
	if(TBB_FOUND AND NOT TARGET CGAL::TBB_support)

		if(NOT TARGET Threads::Threads)
			find_package(Threads REQUIRED)
		endif()

		add_library(CGAL::TBB_support INTERFACE IMPORTED)
		set_target_properties(CGAL::TBB_support PROPERTIES
			INTERFACE_COMPILE_DEFINITIONS "CGAL_LINKED_WITH_TBB;NOMINMAX"
			INTERFACE_INCLUDE_DIRECTORIES "${TBB_INCLUDE_DIRS}"
			INTERFACE_LINK_LIBRARIES "TBB::tbb;TBB::tbbmalloc;Threads::Threads")
	endif()

	list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})

	# executables
	# add_executable(cgal ${PROJECT_SOURCE_DIR}/tests/cgal.cpp)
	# target_link_libraries(cgal PUBLIC PSI CGAL::CGAL CGAL::TBB_support)
	# target_include_directories(cgal PUBLIC ${PROJECT_SOURCE_DIR}/include )


	add_executable(kd_ccp ${PROJECT_SOURCE_DIR}/tests/kd_ccp.cpp)
	target_link_libraries(kd_ccp PUBLIC PSI CGAL::CGAL CGAL::TBB_support)
	target_include_directories(kd_ccp PUBLIC ${PROJECT_SOURCE_DIR}/include)

	add_executable(p_ccp ${PROJECT_SOURCE_DIR}/tests/p_ccp.cpp)
	target_link_libraries(p_ccp PUBLIC PSI CGAL::CGAL CGAL::TBB_support)
	target_include_directories(p_ccp PUBLIC ${PROJECT_SOURCE_DIR}/include)


	message(STATUS "\n\n---------------------------- END CGAL -----------------------------")
endif()

# NOTE: warning policy. Applied here, after every target exists -- kd_ccp and
# p_ccp are declared inside the if(CGAL) block above, so an earlier foreach
# would silently skip them and they would inherit no warning flag at all.
#
# The three examples instantiate KdTree, OrthTree and PTree between them, so
# they are where library warnings surface. They are NOT warning-free: ~40
# warnings remain, all from vendored code that constraint forbids touching
# (space_filling_curve/hilbert.hpp 18, dependence/cpam 15, and libstdc++
# headers inlined from cpam 7). Zero come from base_tree_impl/, kd_tree_impl/,
# orth_tree_impl/ or p_tree_impl/. Do NOT add -Werror here until hilbert.hpp
# is moved out of splitter.h and the vendored roots are marked SYSTEM.
#
# Everything else stays quiet: tests/ alone accounts for ~79 warning sites that
# are not library code. Escalating past -Wall -Wextra (see the commented-out
# list near the top of this file) costs ~104 sites concentrated in the frozen
# index arithmetic; that is a separate, separately-measured job.
foreach(psi_warn_target kd_tree_example orth_tree_example p_tree_example)
	if(TARGET ${psi_warn_target})
		target_compile_options(${psi_warn_target} PRIVATE ${PSI_WARNINGS})
	endif()
endforeach()
foreach(psi_quiet_target kd_test p_test baselines boost_rtree
                         data_washer data_generator kd_ccp p_ccp)
	if(TARGET ${psi_quiet_target})
		target_compile_options(${psi_quiet_target} PRIVATE -w)
	endif()
endforeach()

# Project-owned sources only. include/parlaylib and include/libmorton are
# submodules and must never be reformatted.
file(GLOB_RECURSE all_files CONFIGURE_DEPENDS
  ${PROJECT_SOURCE_DIR}/include/psi/*.h
  ${PROJECT_SOURCE_DIR}/include/psi/*.hpp
  ${PROJECT_SOURCE_DIR}/include/baselines/*.h
  ${PROJECT_SOURCE_DIR}/include/baselines/*.hpp
  ${PROJECT_SOURCE_DIR}/include/baselines/*.cpp
  ${PROJECT_SOURCE_DIR}/tests/*.h
  ${PROJECT_SOURCE_DIR}/tests/*.cpp
  ${PROJECT_SOURCE_DIR}/example/*.h
  ${PROJECT_SOURCE_DIR}/example/*.cpp
)


FIND_PROGRAM(CLANG_FORMAT "clang-format")

IF(CLANG_FORMAT)
	add_custom_target(
		format
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
		COMMAND ${CLANG_FORMAT}
		-style=file
		-i
		${all_files}
	)
endif()
